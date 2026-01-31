import tkinter as tk
from tkinter import ttk, messagebox
import customtkinter

import datetime
import threading
import subprocess
import os
import numpy as np  # for rise/set sampling

import queue
import time

# ZMQ for spectrum streaming (optional)
try:
    import zmq
except ImportError:
    zmq = None

# Matplotlib for live spectrum plot
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import io
import urllib.request

try:
    from PIL import Image, ImageTk
except ImportError:
    Image = None
    ImageTk = None


# Timezone handling: use stdlib zoneinfo if available, else backports
try:
    from zoneinfo import ZoneInfo
except ImportError:
    from backports.zoneinfo import ZoneInfo

import astropy.units as u
from astropy.coordinates import (
    SkyCoord,
    AltAz,
    EarthLocation,
    get_body,
)
from astropy.time import Time

from ATATools import ata_control as ac

# ======================================================
# CONFIG / CONSTANTS
# ======================================================

# Predefined observation profiles: name -> frequency (Hz)
OBSERVATION_PROFILES = {
    "Custom": None,
    "21cm Milky Way - 1420.405 MHz": 1420.405e6,
    "Tianwen-1 Mars Orbiter - 8431.0 MHz": 8431.0e6,
}

# Directory where USRP test scripts *might* live (currently unused but kept)
HOME_DIR = os.path.expanduser("~")

# Spectrum / streaming config
SPECTRUM_ZMQ_ADDRESS = "tcp://arise-sink.hcro.org:5550"
SPECTRUM_FFT_SIZE = 8192        # nominal, but we'll use whatever comes in
SPECTRUM_SAMP_RATE = 3.84e6     # Hz
SPECTRUM_SCRIPT = os.path.join(
    os.path.dirname(__file__),
    "GUI_Stream_21cm.py"
)

# Camera config
CAMERA_PAGE_URL = "http://10.3.0.30/camera/index.html?id=342&imagepath=%2Fmjpg%2Fvideo.mjpg&size=1#/video"
CAMERA_MJPEG_URL = "http://10.3.0.30/mjpg/video.mjpg"



# ======================================================
# ASTRO UTILS
# ======================================================

def compute_altitude_and_rise_set(coord, location, horizon_deg=18.0, n_steps=240):
    """
    Given an ICRS SkyCoord, compute:
      - current altitude (deg)
      - whether it's currently above `horizon_deg`
      - approximate next rise and set times (when crossing `horizon_deg`).

    We do a brute-force sampling over +/- 12 hours from now at 1h resolution
    in local time, and interpolate to find the horizon crossings. This is not
    super precise, but good enough for GUI feedback.
    """
    now_utc = datetime.datetime.utcnow().replace(tzinfo=datetime.timezone.utc)
    # Convert to astropy Time object for convenience
    t0 = Time(now_utc)

    # Current altitude
    altaz_now = coord.transform_to(AltAz(obstime=t0, location=location))
    alt_deg = altaz_now.alt.deg

    # We'll sample times in the next +/- 12 hours
    hours = np.linspace(-12, 12, n_steps)
    times = t0 + hours * u.hour

    alt_list = []
    for tt in times:
        altaz = coord.transform_to(AltAz(obstime=tt, location=location))
        alt_list.append(altaz.alt.deg)
    alt_arr = np.array(alt_list)

    # We want approximate times when alt_arr crosses horizon_deg.
    # We'll compute sign of alt - horizon, and look for zeros.
    sign = np.sign(alt_arr - horizon_deg)
    # Replace zeros with previous sign to avoid degenerate transitions
    for i in range(1, len(sign)):
        if sign[i] == 0:
            sign[i] = sign[i - 1]
    crossings = []
    for i in range(len(sign) - 1):
        if sign[i] * sign[i + 1] < 0:
            # There's a crossing between i and i+1
            # We'll do a simple linear interpolation in time.
            a1, a2 = alt_arr[i], alt_arr[i + 1]
            t1, t2 = times[i], times[i + 1]
            if a2 != a1:
                frac = (horizon_deg - a1) / (a2 - a1)
            else:
                frac = 0.5
            tcross = t1 + frac * (t2 - t1)
            crossings.append(tcross)

    # Classify the next rising and setting times relative to now.
    # alt_arr[-1], alt_arr[0], alt_deg might help determine order, but
    # we can just pick the next crossing in future for each change from
    # below->above and above->below by re-checking alt around crossing.
    alt_at_cross = []
    for t_c in crossings:
        alt_c = coord.transform_to(AltAz(obstime=t_c, location=location)).alt.deg
        alt_at_cross.append(alt_c)

    # We'll define "rise" as crossing from below horizon to above in the sampled array.
    rises = []
    sets = []
    for i, t_c in enumerate(crossings):
        # We'll look at alt slightly before and after the crossing index
        # using the sign array as a guide. However, to keep it simple,
        # we can check the index 'i' in hours: near hours[i].
        idx = i
        if idx < 0 or idx >= len(hours) - 1:
            continue
        alt_before = alt_arr[idx]
        alt_after = alt_arr[idx + 1]
        if alt_before < horizon_deg <= alt_after:
            rises.append(t_c)
        elif alt_before > horizon_deg >= alt_after:
            sets.append(t_c)

    # Now filter so we only take future events
    future_rises = [t for t in rises if t > t0]
    future_sets = [t for t in sets if t > t0]

    next_rise = future_rises[0] if future_rises else None
    next_set = future_sets[0] if future_sets else None

    is_up = alt_deg >= horizon_deg
    return alt_deg, is_up, next_rise, next_set


# ======================================================
# ATA CONTROL WRAPPERS (single place where ac.* is used)
# ======================================================

# Antenna set that can be chosen in the GUI (for now just 1a)
AVAILABLE_ANTENNAS = ['1a']

# Current antenna selection used by the control wrappers.
# This will be updated from the "Select Antennas" dropdown.
antennas = AVAILABLE_ANTENNAS.copy()

# ATA site (approx Hat Creek / ATA)
ATA_LOCATION = EarthLocation.from_geodetic(
    lon=-121.47 * u.deg,
    lat=40.82 * u.deg,
    height=986 * u.m
)

# Local timezone for logs & display
HCRO_TZ = ZoneInfo("America/Los_Angeles")


def set_global_antennas(new_ants):
    """
    Safely update the global 'antennas' list.

    Called whenever the GUI user changes antenna selection.
    """
    global antennas
    antennas = list(new_ants)


def ata_reserve_antennas():
    """
    Reserve antennas for ATA-GR pipeline by moving them
    from group 'none' to group 'atagr'.
    """
    ac.move_ant_group(antennas, 'none', 'atagr')
    return f"Antennas {antennas} moved from 'none' to 'atagr' (reserved)."



def ata_park_antennas():
    """
    Park the antennas.

    Uses ac.park_antennas(antennas).
    """
    ac.park_antennas(antennas)
    return f"Antennas {antennas} parked."


def ata_release_antennas():
    """
    Release antennas from ATA-GR pipeline.

    Uses ac.move_ant_group(antennas, 'atagr', 'none').
    """
    ac.move_ant_group(antennas, 'atagr', 'none')
    return f"Antennas {antennas} moved from 'atagr' to 'none' (released)."


def ata_track_radec(ra_str, dec_str):
    """
    Track RA/Dec given as strings, e.g. '18:18:18', '55:55:55'.

    Uses ac.track_source(antennas, radec=[ra_hours, dec_deg]).
    """
    coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg), frame="icrs")
    ra_deg = coord.ra.deg
    dec_deg = coord.dec.deg
    ra_hours = ra_deg / 15.0
    ac.track_source(antennas, radec=[ra_hours, dec_deg])
    return f"Tracking RA={ra_str}, Dec={dec_str} on antennas {antennas}."


def ata_track_altaz(alt_str, az_str):
    """
    Track Alt/Az given as strings in degrees, e.g. '45.0', '180.0'.

    Uses ac.track_source(antennas, altaz=[alt_deg, az_deg]).
    """
    alt_deg = float(alt_str)
    az_deg = float(az_str)
    ac.track_source(antennas, altaz=[alt_deg, az_deg])
    return f"Tracking Alt={alt_deg:.2f}, Az={az_deg:.2f} on antennas {antennas}."


def ata_track_galactic(l_str, b_str):
    """
    Track Galactic coords given as strings in degrees.

    l_str: Galactic longitude in deg
    b_str: Galactic latitude in deg

    Internally this is converted to ICRS RA/Dec and then ac.track_source is called.
    """
    gal = SkyCoord(l=float(l_str) * u.deg, b=float(b_str) * u.deg,
                   frame="galactic")
    icrs = gal.icrs
    ra_hours = icrs.ra.deg / 15.0
    dec_deg = icrs.dec.deg
    ac.track_source(antennas, radec=[ra_hours, dec_deg])
    return (f"Tracking Galactic l={l_str} deg, b={b_str} deg "
            f"(ICRS RA={icrs.ra.to_string(unit=u.hourangle, sep=':')}, "
            f"Dec={icrs.dec.to_string(sep=':')}) on antennas {antennas}.")


def ata_track_source_name(name_str):
    """
    Track source by name, e.g. '3C286' or 'M31'.

    We use astropy's SkyCoord.from_name to resolve the name to ICRS RA/Dec.
    """
    coord = SkyCoord.from_name(name_str)
    ra_hours = coord.ra.deg / 15.0
    dec_deg = coord.dec.deg
    ac.track_source(antennas, radec=[ra_hours, dec_deg])
    return (f"Tracking source '{name_str}' at RA={coord.ra.to_string(unit=u.hourangle, sep=':')}, "
            f"Dec={coord.dec.to_string(sep=':')} on antennas {antennas}.")


def ata_set_frequency(freq_mhz):
    """
    Set observing frequency (in MHz) for the current antennas.

    We call ac.set_freq(freq_mhz, antennas, lo='d').
    """
    lo = 'd'
    ac.set_freq(freq_mhz, antennas, lo)
    return f"Set frequency to {freq_mhz:.6f} MHz on antennas {antennas} (LO={lo})."


def ata_autotune():
    """
    Run autotune on current antennas.

    Uses ac.autotune(antennas).
    """
    ac.autotune(antennas)
    return f"Autotune complete on antennas {antennas}."


def ata_track_radec_degrees(ra_deg, dec_deg):
    """
    Start tracking RA/Dec (degrees) with the ATA.

    Uses ac.track_source(antennas, radec=[ra_hours, dec_deg]).
    """
    ra_hours = ra_deg / 15.0
    ac.track_source(antennas, radec=[ra_hours, dec_deg])
    return (f"Tracking RA={ra_deg:.3f} deg, Dec={dec_deg:.3f} deg "
            f"(RA={ra_hours:.3f} hours) on antennas {antennas}.")


def ata_get_status():
    """
    Get ASCII status from ATA and return as string.

    Uses ac.get_ascii_status().
    """
    return ac.get_ascii_status()


def _run_shell_command_capture_lines(command):
    """
    Run a shell command, capturing stdout and stderr.

    Returns:
      (stdout_lines, stderr_lines, return_code)
    where each of stdout_lines / stderr_lines is a list of strings (with newlines
    stripped).

    We use this to run USRP test commands in a way that is robust to non-zero
    exit codes and that matches the user's shell usage closely (PATH, etc.).
    """
    try:
        proc = subprocess.run(
            command,
            text=True,
            capture_output=True,
            shell=True,   # Use shell so "usrp_test.py" gets found on PATH
            check=False,  # Don't raise on non-zero exit
        )
    except Exception as e:
        return [f"Error running {command}: {e}"], [], -1

    stdout_lines = proc.stdout.splitlines()
    stderr_lines = proc.stderr.splitlines()
    return stdout_lines, stderr_lines, proc.returncode


def usrp_check_levels():
    """
    Run 'usrp_test.py' and return only the interesting lines for the GUI log.

    Expected interesting lines:
        Channel 0 level -16.19 dB
        ...
        Note: channel levels should be around -16 dB in this test

    We deliberately do NOT enforce a zero exit code; we mimic os.popen-like
    behavior and just parse whatever stdout we get. We also invoke the command
    exactly as you do in the shell (`usrp_test.py`) via the shell so $PATH is
    used.
    """
    try:
        proc = subprocess.run(
            "usrp_test.py",
            text=True,
            capture_output=True,
            check=False,   # DO NOT raise on non-zero exit
            shell=True,    # so that $PATH is used, matching the shell usage
        )
    except Exception as e:
        return [f"Error running usrp_test.py: {e}"]

    stdout_lines = proc.stdout.splitlines()
    stderr_lines = proc.stderr.splitlines()

    # If the script wrote stderr, show a short note plus those lines.
    if stderr_lines:
        return ["usrp_test.py stderr:"] + stderr_lines

    # Otherwise, return whatever stdout we got.  This matches the original
    # os.popen+read behaviour, where you saw the full output.
    if stdout_lines:
        return stdout_lines
    else:
        return [f"usrp_test.py finished with return code {proc.returncode}, but produced no output."]


def usrp_reset_clocking():
    """
    Run 'usrp_reset_clocking.py' and return a short status message.

    The GUI only needs to know that the command completed successfully, or see
    any stderr for debugging if something went wrong.
    """
    try:
        proc = subprocess.run(
            "usrp_reset_clocking.py",
            text=True,
            capture_output=True,
            check=False,
            shell=True,
        )
    except Exception as e:
        return [f"Error running usrp_reset_clocking.py: {e}"]

    stderr_lines = proc.stderr.splitlines()
    if proc.returncode == 0 and not stderr_lines:
        return ["USRP clocking reset success."]
    else:
        if stderr_lines:
            return ["usrp_reset_clocking.py stderr:"] + stderr_lines
        return [f"usrp_reset_clocking.py finished with return code {proc.returncode}."]


# ======================================================
# GUI
# ======================================================

class ATAObservationGUI:

    def __init__(self, root):
        self.root = root
        self.root.title("Allen Telescope Array GNURadio Observation GUI")
        self.root.geometry("1800x1200")

        customtkinter.set_appearance_mode("dark")
        customtkinter.set_default_color_theme("blue")

        # Coordinate mode: radec, altaz, galactic, name
        self.coord_mode = tk.StringVar(value="radec")
        # Cache for coordinate entries by mode (for the two numeric fields)
        self.coord_cache = {
            "radec": ["", ""],
            "altaz": ["", ""],
            "galactic": ["", ""],
        }
        self.last_coord_mode = "radec"

        # Selected profile – explicitly default to Custom
        self.selected_profile = tk.StringVar(value="Custom")

        # Selected antennas (as text in the dropdown)
        self.antenna_select_var = tk.StringVar(value="1a")

        # --- Spectrum streaming state ---
        self.spectrum_running = False
        self.spectrum_thread = None
        self.spectrum_stop_event = threading.Event()
        self.spectrum_queue = queue.Queue(maxsize=10)
        self.spectrum_context = None
        self.spectrum_socket = None
        self.spectrum_first_frame_received = False
        self.stream_process = None
        # --- Camera streaming state ---
        self.camera_running = False
        self.camera_thread = None
        self.camera_stop_event = threading.Event()
        self.camera_photo = None
        self.camera_label = None
        self.camera_status_label = None


        self._build_layout()

        # Ensure global antennas initialised from dropdown
        self._update_antennas_from_selection(self.antenna_select_var.get())

        # Start periodic tasks
        self._update_time_info()
        self._update_spectrum_plot()

    # ---------- Layout ----------

    def _build_layout(self):
        self.main_frame = customtkinter.CTkFrame(master=self.root)
        self.main_frame.pack(fill="both", expand=True)

        # Top info frame: location & time
        self.info_frame = customtkinter.CTkFrame(master=self.main_frame)
        self.info_frame.pack(fill="x", pady=(10, 0), padx=10)

        self._build_info_bar(self.info_frame)

        # Middle: left control panel and right tabs (Status / Spectrum / Camera)
        self.middle_frame = customtkinter.CTkFrame(master=self.main_frame)
        self.middle_frame.pack(fill="both", expand=True, padx=10, pady=10)

        self.left_frame = customtkinter.CTkFrame(master=self.middle_frame, width=400)
        self.left_frame.pack(side="left", fill="y", padx=(0, 10))

        self.right_frame = customtkinter.CTkFrame(master=self.middle_frame)
        self.right_frame.pack(side="right", fill="both", expand=True)

        # Bottom: log + progress
        self.log_frame = customtkinter.CTkFrame(master=self.main_frame)
        self.log_frame.pack(fill="x", pady=(0, 10), padx=10)

        # Build sections
        self._build_left_controls(self.left_frame)
        self._build_right_tabs(self.right_frame)
        self._build_log_area(self.log_frame)

    def _build_left_controls(self, parent):
        # ---- Antenna control ----
        antenna_label = customtkinter.CTkLabel(
            parent,
            text="Antenna Control",
            font=("Arial", 18, "bold")
        )
        antenna_label.pack(pady=(5, 5))

        button_frame = customtkinter.CTkFrame(parent)
        button_frame.pack(fill="x", pady=(0, 10))

        select_label = customtkinter.CTkLabel(
            button_frame,
            text="Select Antenna(s)"
        )
        select_label.pack(side="left", padx=5, pady=5)

        antenna_select = customtkinter.CTkComboBox(
            button_frame,
            values=AVAILABLE_ANTENNAS,
            variable=self.antenna_select_var,
            command=self.on_antennas_changed,
            width=90
        )
        antenna_select.pack(side="left", padx=5, pady=5)

        reserve_btn = customtkinter.CTkButton(
            button_frame,
            text="Reserve Antenna(s)",
            command=self.on_reserve_clicked
        )
        reserve_btn.pack(side="left", padx=5, pady=5)

        release_btn = customtkinter.CTkButton(
            button_frame,
            text="Release Antenna(s)",
            command=self.on_release_clicked
        )
        release_btn.pack(side="left", padx=5, pady=5)

        park_btn = customtkinter.CTkButton(
            button_frame,
            text="Park Antenna(s)",
            command=self.on_park_clicked
        )
        park_btn.pack(side="left", padx=5, pady=5)

        # ---- Frequency / Observation profile ----
        freq_frame = customtkinter.CTkFrame(parent)
        freq_frame.pack(fill="x", pady=(10, 10))

        freq_label = customtkinter.CTkLabel(
            freq_frame,
            text="Set Frequency / Observation Profile",
            font=("Arial", 16, "bold")
        )
        freq_label.pack(anchor="w", pady=(0, 5))

        freq_inner = customtkinter.CTkFrame(freq_frame)
        freq_inner.pack(fill="x", pady=(0, 5))

        # Left: profile dropdown
        profile_label = customtkinter.CTkLabel(
            freq_inner,
            text="Profile:"
        )
        profile_label.grid(row=0, column=0, padx=5, pady=5, sticky="w")

        profile_combo = customtkinter.CTkComboBox(
            freq_inner,
            values=list(OBSERVATION_PROFILES.keys()),
            variable=self.selected_profile,
            command=self.on_profile_changed,
            width=260,
        )
        profile_combo.grid(row=0, column=1, padx=5, pady=5, sticky="w")

        # Right: frequency entry
        freq_entry_label = customtkinter.CTkLabel(
            freq_inner,
            text="Frequency"
        )
        freq_entry_label.grid(row=1, column=0, padx=5, pady=5, sticky="w")

        self.freq_entry = customtkinter.CTkEntry(
            freq_inner,
            width=140,
            placeholder_text=""
        )
        self.freq_entry.grid(row=1, column=1, padx=5, pady=5, sticky="w")

        mhz_label = customtkinter.CTkLabel(
            freq_inner,
            text="MHz"
        )
        mhz_label.grid(row=1, column=2, padx=5, pady=5, sticky="w")

        set_freq_btn = customtkinter.CTkButton(
            freq_inner,
            text="Set Frequency",
            command=self.on_set_frequency_clicked
        )
        set_freq_btn.grid(row=1, column=3, padx=5, pady=5, sticky="w")

        # ---- Pointing / Coordinates ----
        coord_frame = customtkinter.CTkFrame(parent)
        coord_frame.pack(fill="x", pady=(10, 10))

        coord_label = customtkinter.CTkLabel(
            coord_frame,
            text="Pointing / Coordinate System",
            font=("Arial", 16, "bold")
        )
        coord_label.pack(anchor="w", pady=(0, 5))

        # Coordinate mode radio buttons
        mode_frame = customtkinter.CTkFrame(coord_frame)
        mode_frame.pack(fill="x", pady=(0, 5))

        for text, mode in [
            ("RA/Dec", "radec"),
            ("Alt/Az", "altaz"),
            ("Galactic", "galactic"),
            ("Source Name", "name"),
        ]:
            rb = customtkinter.CTkRadioButton(
                mode_frame,
                text=text,
                variable=self.coord_mode,
                value=mode,
                command=self.on_coord_mode_changed,
            )
            rb.pack(side="left", padx=5)

        # Coordinate entry fields
        coord_input_frame = customtkinter.CTkFrame(coord_frame)
        coord_input_frame.pack(fill="x", pady=(5, 5))

        # Left side: dynamic field area
        self.coord_fields_frame = customtkinter.CTkFrame(coord_input_frame)
        self.coord_fields_frame.pack(side="left", fill="x", expand=True)

        # RA/Dec entries
        self.ra_label = customtkinter.CTkLabel(self.coord_fields_frame, text="RA (H:M:S)")
        self.ra_entry = customtkinter.CTkEntry(
            self.coord_fields_frame,
            width=160,
            placeholder_text="18:18:18"
        )
        self.dec_label = customtkinter.CTkLabel(self.coord_fields_frame, text="Dec (D:M:S)")
        self.dec_entry = customtkinter.CTkEntry(
            self.coord_fields_frame,
            width=160,
            placeholder_text="55:55:55"
        )

        # Alt/Az entries
        self.alt_label = customtkinter.CTkLabel(self.coord_fields_frame, text="Alt (deg)")
        self.alt_entry = customtkinter.CTkEntry(
            self.coord_fields_frame,
            width=160,
            placeholder_text="45.0"
        )
        self.az_label = customtkinter.CTkLabel(self.coord_fields_frame, text="Az (deg)")
        self.az_entry = customtkinter.CTkEntry(
            self.coord_fields_frame,
            width=160,
            placeholder_text="180.0"
        )

        # Galactic entries
        self.l_label = customtkinter.CTkLabel(self.coord_fields_frame, text="Galactic l (deg)")
        self.l_entry = customtkinter.CTkEntry(
            self.coord_fields_frame,
            width=160,
            placeholder_text="120.0"
        )
        self.b_label = customtkinter.CTkLabel(self.coord_fields_frame, text="Galactic b (deg)")
        self.b_entry = customtkinter.CTkEntry(
            self.coord_fields_frame,
            width=160,
            placeholder_text="0.0"
        )

        # Source name entry (for 'name' mode only)
        self.name_label = customtkinter.CTkLabel(
            self.coord_fields_frame, text="Source Name"
        )
        self.name_entry = customtkinter.CTkEntry(
            self.coord_fields_frame,
            width=180,
            placeholder_text="3C286"
        )
        # name_label/name_entry will be gridded only in "name" mode

        # Right side: fixed-size "Check Source Location" button
        self.check_source_btn = customtkinter.CTkButton(
            coord_input_frame,
            text="Check Source Location",
            command=self.on_check_source_location,
            width=180
        )
        self.check_source_btn.pack(side="right", padx=10, pady=5)

        # Initially build RA/Dec fields
        self._show_coord_fields("radec")

        # Row under coordinates: track button + autotune
        track_frame = customtkinter.CTkFrame(coord_frame)
        track_frame.pack(fill="x", pady=(5, 0))

        track_btn = customtkinter.CTkButton(
            track_frame,
            text="Track Source",
            command=self.on_track_clicked
        )
        track_btn.pack(side="left", padx=5, pady=5)

        autotune_btn = customtkinter.CTkButton(
            track_frame,
            text="Autotune",
            command=self.on_autotune_clicked
        )
        autotune_btn.pack(side="left", padx=5, pady=5)

        park_btn2 = customtkinter.CTkButton(
            track_frame,
            text="Park Antenna(s)",
            command=self.on_park_clicked
        )
        park_btn2.pack(side="left", padx=5, pady=5)

        # Information label for Check Source Location
        self.source_info_label = customtkinter.CTkLabel(
            coord_frame,
            text="",
            justify="left"
        )
        self.source_info_label.pack(anchor="w", padx=5, pady=(5, 0))

        # ---- USRP Utilities ----
        usrp_frame = customtkinter.CTkFrame(parent)
        usrp_frame.pack(fill="x", pady=(10, 10))

        usrp_label = customtkinter.CTkLabel(
            usrp_frame,
            text="USRP Utilities",
            font=("Arial", 16, "bold")
        )
        usrp_label.pack(anchor="w", pady=(0, 5))

        usrp_btn_frame = customtkinter.CTkFrame(usrp_frame)
        usrp_btn_frame.pack(fill="x")

        check_usrp_btn = customtkinter.CTkButton(
            usrp_btn_frame,
            text="Check Channel Levels",
            command=self.on_usrp_check_clicked
        )
        check_usrp_btn.pack(side="left", padx=5, pady=5)

        reset_usrp_btn = customtkinter.CTkButton(
            usrp_btn_frame,
            text="Reset Channel Clocking",
            command=self.on_usrp_reset_clicked
        )
        reset_usrp_btn.pack(side="left", padx=5, pady=5)

    def _build_info_bar(self, parent):
        # Location label
        loc_label = customtkinter.CTkLabel(
            parent,
            text="Hat Creek Radio Observatory (ATA)",
            font=("Arial", 16, "bold")
        )
        loc_label.pack(side="left", padx=10)

        # Time labels
        self.local_time_label = customtkinter.CTkLabel(
            parent,
            text="Local time (Hat Creek): --",
            font=("Arial", 14)
        )
        self.local_time_label.pack(side="right", padx=10)

        self.utc_time_label = customtkinter.CTkLabel(
            parent,
            text="UTC: --",
            font=("Arial", 14)
        )
        self.utc_time_label.pack(side="right", padx=10)

    def _build_right_tabs(self, parent):
        notebook = ttk.Notebook(parent)
        notebook.pack(fill="both", expand=True)

        # Status tab
        status_frame = customtkinter.CTkFrame(notebook)
        notebook.add(status_frame, text="Status")

        self.status_text = tk.Text(
            status_frame, wrap="word", height=20, width=80
        )
        self.status_text.pack(
            fill="both", expand=True, padx=5, pady=5
        )
        # Make status read-only
        self.status_text.configure(state="disabled")

        refresh_btn = customtkinter.CTkButton(
            status_frame,
            text="Show Antenna Status",
            command=self.on_refresh_status_clicked
        )
        refresh_btn.pack(pady=5)

        # Spectrum tab
        spectrum_frame = customtkinter.CTkFrame(notebook)
        notebook.add(spectrum_frame, text="Spectrum")
        self._build_spectrum_tab(spectrum_frame)

        # Camera tab
        camera_frame = customtkinter.CTkFrame(notebook)
        notebook.add(camera_frame, text="Camera")
        self._build_camera_tab(camera_frame)

        self.right_notebook = notebook

    def _build_spectrum_tab(self, parent):
        # Top control strip
        control_frame = customtkinter.CTkFrame(parent)
        control_frame.pack(side="top", fill="x", padx=10, pady=5)

        start_btn = customtkinter.CTkButton(
            control_frame,
            text="Start Data Stream",
            command=self.start_data_stream,
        )
        start_btn.grid(row=0, column=0, sticky="w", padx=5, pady=2)

        stop_btn = customtkinter.CTkButton(
            control_frame,
            text="Stop Data Stream",
            command=self.stop_data_stream,
        )
        stop_btn.grid(row=0, column=1, sticky="w", padx=5, pady=2)

        customtkinter.CTkLabel(control_frame, text="ZMQ Address").grid(
            row=1, column=0, sticky="w", padx=5, pady=2
        )
        self.spec_address_var = tk.StringVar(value=SPECTRUM_ZMQ_ADDRESS)
        customtkinter.CTkEntry(
            control_frame,
            textvariable=self.spec_address_var,
            width=260,
        ).grid(row=1, column=1, columnspan=2, sticky="w", padx=5, pady=2)

        # Status text
        self.spectrum_status_label = customtkinter.CTkLabel(
            control_frame,
            text="Data stream not running",
            text_color="red",
        )
        self.spectrum_status_label.grid(
            row=2, column=0, columnspan=3, sticky="w", padx=5, pady=2
        )

        # Figure area
        plot_frame = customtkinter.CTkFrame(parent)
        plot_frame.pack(side="top", fill="both", expand=True, padx=10, pady=5)

        self.spec_figure = Figure(figsize=(6, 4), dpi=100)
        self.spec_ax = self.spec_figure.add_subplot(111)
        self.spec_ax.set_title("Live Spectrum")
        self.spec_ax.set_xlabel("Frequency (MHz)")
        self.spec_ax.set_ylabel("Power (dB)")
        self.spec_line, = self.spec_ax.plot([], [], linewidth=0.7)

        self.spec_canvas = FigureCanvasTkAgg(self.spec_figure, master=plot_frame)
        self.spec_canvas.draw()
        self.spec_canvas.get_tk_widget().pack(fill="both", expand=True)

    def _build_camera_tab(self, parent):
        """Build the camera tab with connect/disconnect and live view."""
        control_frame = customtkinter.CTkFrame(parent)
        control_frame.pack(side="top", fill="x", padx=10, pady=5)

        connect_btn = customtkinter.CTkButton(
            control_frame,
            text="Connect Camera",
            command=self.start_camera,
        )
        connect_btn.grid(row=0, column=0, padx=5, pady=5, sticky="w")

        disconnect_btn = customtkinter.CTkButton(
            control_frame,
            text="Disconnect Camera",
            command=self.stop_camera,
        )
        disconnect_btn.grid(row=0, column=1, padx=5, pady=5, sticky="w")

        self.camera_status_label = customtkinter.CTkLabel(
            control_frame,
            text="Camera not connected",
            text_color="red",
        )
        self.camera_status_label.grid(
            row=1, column=0, columnspan=2, sticky="w", padx=5, pady=(0, 5)
        )

        image_frame = customtkinter.CTkFrame(parent, width=620, height=500)
        image_frame.pack(side="top", fill="x", expand=False, padx=10, pady=5)
        image_frame.pack_propagate(False)

        # We use a plain Tk label here because PhotoImage works directly with it.
        self.camera_label = tk.Label(image_frame, bg="black")
        self.camera_label.pack(fill="both", expand=True)

    def _build_log_area(self, parent):
        # Progress indicator + status line
        top_frame = customtkinter.CTkFrame(parent)
        top_frame.pack(fill="x", pady=(5, 0))

        # Progress bar for long-running operations (slews, autotune, etc.)
        self.progress_bar = customtkinter.CTkProgressBar(
            top_frame,
            mode="indeterminate",
            width=200,
        )
        self.progress_bar.pack(side="left", padx=(5, 5))
        # Ensure it's visually empty at startup
        self.progress_bar.set(0.0)

        # Text label next to the bar
        self.progress_label = customtkinter.CTkLabel(
            top_frame,
            text="Idle",
            anchor="w"
        )
        self.progress_label.pack(side="left", padx=5)

        save_log_btn = customtkinter.CTkButton(
            top_frame,
            text="Save Log to File",
            command=self.on_save_log_clicked
        )
        save_log_btn.pack(side="right", padx=5)

        # Log text area
        self.log_text = tk.Text(
            parent,
            wrap="word",
            height=15,
        )
        self.log_text.pack(fill="both", expand=True, padx=5, pady=5)
        self.log_text.configure(state="disabled")

        # Tag styles for up/down messages in source location
        self.log_text.tag_config("up", foreground="green")
        self.log_text.tag_config("down", foreground="red")


    # ---------- Helpers ----------

    def _parse_antenna_selection(self, text: str):
        """
        Parse a comma-separated antenna selection string into a list.

        For now the combo only ever returns '1a', but this is written so that
        in the future if the control is changed to allow multi-selection or
        user-editable values like '1a, 1f, 5c', we get ['1a','1f','5c'].
        """
        parts = [p.strip() for p in text.split(",") if p.strip()]
        return parts if parts else AVAILABLE_ANTENNAS.copy()

    def _update_antennas_from_selection(self, choice: str):
        """
        Update global antennas based on the combo box selection.
        """
        ants = self._parse_antenna_selection(choice)
        set_global_antennas(ants)

    def log(self, message, tag=None):
        """
        Append a message to the log text area with a timestamp.
        Most recent message goes at the top.
        """
        now = datetime.datetime.now(tz=HCRO_TZ)
        timestr = now.strftime("%Y-%m-%d %H:%M:%S %Z")
        full_msg = f"[{timestr}] {message}\n"

        self.log_text.configure(state="normal")
        if tag:
            self.log_text.insert("1.0", full_msg, tag)
        else:
            self.log_text.insert("1.0", full_msg)
        self.log_text.configure(state="disabled")

    def run_with_progress(self, description, worker_func, log_in_status=False):
        """
        Run worker_func in a background thread while updating the progress
        bar/label, then log its result back on the Tk thread.

        worker_func: callable taking no args.
          It may return:
            - None
            - a single string
            - a list/tuple of strings

        If log_in_status is True, returned strings are appended to the
        antenna status text widget instead of the main log.
        """
        # Start progress bar + label immediately on the GUI thread
        try:
            if hasattr(self, "progress_bar") and self.progress_bar is not None:
                # Indeterminate spinner
                self.progress_bar.configure(mode="indeterminate")
                self.progress_bar.start()
        except Exception:
            pass

        try:
            if hasattr(self, "progress_label") and self.progress_label is not None:
                self.progress_label.configure(text=f"{description} ...")
        except Exception:
            pass

        def finish(err=None, result=None):
            # This runs back on the Tk thread via root.after
            try:
                if hasattr(self, "progress_bar") and self.progress_bar is not None:
                    self.progress_bar.stop()
                    # Reset to empty bar
                    self.progress_bar.set(0.0)
            except Exception:
                pass

            try:
                if hasattr(self, "progress_label") and self.progress_label is not None:
                    self.progress_label.configure(text="Idle")
            except Exception:
                pass

            if err is not None:
                # Log and show an error dialog
                self.log(f"{description} failed: {err}", tag="error")
                try:
                    messagebox.showerror("Error", f"{description} failed:\n{err}")
                except Exception:
                    pass
                return

            # Normalize result to a list of strings
            if result is None:
                items = []
            elif isinstance(result, (list, tuple)):
                items = [s for s in result if s is not None]
            else:
                items = [result]

            if not items:
                # Worker returned nothing; stay quiet
                return

            if log_in_status and hasattr(self, "status_text"):
                wrote_to_status = False
                try:
                    self.status_text.configure(state="normal")
                    for s in items:
                        self.status_text.insert("end", str(s) + "\n")
                    self.status_text.see("end")
                    self.status_text.configure(state="disabled")
                    wrote_to_status = True
                except Exception:
                    wrote_to_status = False

                if not wrote_to_status:
                    for s in items:
                        self.log(str(s))
            else:
                # Normal case: log results to main log
                for s in items:
                    self.log(str(s))

        def worker():
            try:
                result = worker_func()
                err = None
            except Exception as e:
                result = None
                err = e
            # Bounce back to the Tk thread
            self.root.after(
                0,
                lambda err=err, result=result: finish(err=err, result=result)
            )

        # Launch the worker in a background thread
        t = threading.Thread(target=worker, daemon=True)
        t.start()


    def _update_time_info(self):
        """
        Periodically update local time (Hat Creek) and UTC.
        """
        now_utc = datetime.datetime.now(datetime.timezone.utc)
        now_local = now_utc.astimezone(HCRO_TZ)

        self.local_time_label.configure(
            text=f"Local time (Hat Creek): {now_local.strftime('%Y-%m-%d %H:%M:%S %Z')}"
        )
        self.utc_time_label.configure(
            text=f"UTC: {now_utc.strftime('%Y-%m-%d %H:%M:%S %Z')}"
        )

        # Schedule next update
        self.root.after(1000, self._update_time_info)

    # ---------- Coordinate Fields Management ----------

    def _show_coord_fields(self, mode):
        """
        Show the appropriate coordinate input fields for the selected mode.
        Preserve previous entries by caching.
        """
        # Cache current mode entries (only numeric fields)
        if self.last_coord_mode == "radec":
            self.coord_cache["radec"] = [
                self.ra_entry.get(),
                self.dec_entry.get()
            ]
        elif self.last_coord_mode == "altaz":
            self.coord_cache["altaz"] = [
                self.alt_entry.get(),
                self.az_entry.get()
            ]
        elif self.last_coord_mode == "galactic":
            self.coord_cache["galactic"] = [
                self.l_entry.get(),
                self.b_entry.get()
            ]

        # Clear the frame
        for widget in self.coord_fields_frame.winfo_children():
            widget.grid_forget()

        if mode == "radec":
            # Restore cached values
            ra_val, dec_val = self.coord_cache["radec"]
            self.ra_entry.delete(0, tk.END)
            self.ra_entry.insert(0, ra_val)
            self.dec_entry.delete(0, tk.END)
            self.dec_entry.insert(0, dec_val)

            self.ra_label.grid(row=0, column=0, padx=5, pady=2, sticky="w")
            self.ra_entry.grid(row=0, column=1, padx=5, pady=2, sticky="w")
            self.dec_label.grid(row=1, column=0, padx=5, pady=2, sticky="w")
            self.dec_entry.grid(row=1, column=1, padx=5, pady=2, sticky="w")

        elif mode == "altaz":
            alt_val, az_val = self.coord_cache["altaz"]
            self.alt_entry.delete(0, tk.END)
            self.alt_entry.insert(0, alt_val)
            self.az_entry.delete(0, tk.END)
            self.az_entry.insert(0, az_val)

            self.alt_label.grid(row=0, column=0, padx=5, pady=2, sticky="w")
            self.alt_entry.grid(row=0, column=1, padx=5, pady=2, sticky="w")
            self.az_label.grid(row=1, column=0, padx=5, pady=2, sticky="w")
            self.az_entry.grid(row=1, column=1, padx=5, pady=2, sticky="w")

        elif mode == "galactic":
            l_val, b_val = self.coord_cache["galactic"]
            self.l_entry.delete(0, tk.END)
            self.l_entry.insert(0, l_val)
            self.b_entry.delete(0, tk.END)
            self.b_entry.insert(0, b_val)

            self.l_label.grid(row=0, column=0, padx=5, pady=2, sticky="w")
            self.l_entry.grid(row=0, column=1, padx=5, pady=2, sticky="w")
            self.b_label.grid(row=1, column=0, padx=5, pady=2, sticky="w")
            self.b_entry.grid(row=1, column=1, padx=5, pady=2, sticky="w")

        elif mode == "name":
            self.name_entry.delete(0, tk.END)

            self.name_label.grid(row=0, column=0, padx=5, pady=2, sticky="w")
            self.name_entry.grid(row=0, column=1, padx=5, pady=2, sticky="w")

        self.last_coord_mode = mode

    # ---------- Callbacks: Left Panel ----------

    def on_coord_mode_changed(self):
        mode = self.coord_mode.get()
        self._show_coord_fields(mode)

    def on_profile_changed(self, choice):
        freq_hz = OBSERVATION_PROFILES.get(choice)
        if freq_hz is None:
            # Custom
            return
        freq_mhz = freq_hz / 1e6
        self.freq_entry.delete(0, tk.END)
        self.freq_entry.insert(0, f"{freq_mhz:.6f}")
        self.log(f"Profile selected: {choice} ({freq_mhz:.6f} MHz)")

    def on_set_frequency_clicked(self):
        txt = self.freq_entry.get().strip()
        if not txt:
            messagebox.showerror("Error", "Please enter a frequency in MHz.")
            return
        try:
            freq_mhz = float(txt)
        except ValueError:
            messagebox.showerror("Error", "Invalid frequency. Please enter a number in MHz.")
            return

        def do_set_freq():
            return [ata_set_frequency(freq_mhz)]

        self.run_with_progress(
            f"Setting frequency to {freq_mhz:.6f} MHz",
            do_set_freq
        )

    def on_track_clicked(self):
        mode = self.coord_mode.get()

        if mode == "radec":
            ra = self.ra_entry.get().strip()
            dec = self.dec_entry.get().strip()
            if not ra or not dec:
                messagebox.showerror("Error", "Please enter RA and Dec.")
                return

            def do_track():
                return [ata_track_radec(ra, dec)]

        elif mode == "altaz":
            alt = self.alt_entry.get().strip()
            az = self.az_entry.get().strip()
            if not alt or not az:
                messagebox.showerror("Error", "Please enter Alt and Az.")
                return

            def do_track():
                return [ata_track_altaz(alt, az)]

        elif mode == "galactic":
            l = self.l_entry.get().strip()
            b = self.b_entry.get().strip()
            if not l or not b:
                messagebox.showerror("Error", "Please enter Galactic l and b.")
                return

            def do_track():
                return [ata_track_galactic(l, b)]

        elif mode == "name":
            name = self.name_entry.get().strip()
            if not name:
                messagebox.showerror("Error", "Please enter a source name.")
                return

            def do_track():
                return [ata_track_source_name(name)]

        else:
            messagebox.showerror("Error", "Unknown pointing mode.")
            return

        self.run_with_progress("Tracking source", do_track)

    def on_autotune_clicked(self):
        def do_autotune():
            return [ata_autotune()]

        self.run_with_progress("Running autotune", do_autotune)

    def on_check_source_location(self):
        """
        Check whether the currently entered coordinates or source name
        are above the 18-deg horizon, and when they rise/set.
        """
        mode = self.coord_mode.get()
        coord = None
        label = ""

        try:
            if mode == "radec":
                ra = self.ra_entry.get().strip()
                dec = self.dec_entry.get().strip()
                if not ra or not dec:
                    messagebox.showerror("Error", "Please enter RA and Dec.")
                    return
                coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame="icrs")
                label = f"RA={ra}, Dec={dec}"

            elif mode == "altaz":
                alt_str = self.alt_entry.get().strip()
                az_str = self.az_entry.get().strip()
                if not alt_str or not az_str:
                    messagebox.showerror("Error", "Please enter Alt and Az.")
                    return
                alt_deg = float(alt_str)
                az_deg = float(az_str)
                # Convert from AltAz at current time/location to ICRS for the check
                now_utc = datetime.datetime.now(datetime.timezone.utc)
                t = Time(now_utc)
                altaz = SkyCoord(
                    alt=alt_deg * u.deg,
                    az=az_deg * u.deg,
                    frame=AltAz(obstime=t, location=ATA_LOCATION),
                )
                coord = altaz.icrs
                label = f"Alt={alt_deg:.1f} deg, Az={az_deg:.1f} deg"

            elif mode == "galactic":
                l_str = self.l_entry.get().strip()
                b_str = self.b_entry.get().strip()
                if not l_str or not b_str:
                    messagebox.showerror("Error", "Please enter Galactic l and b.")
                    return
                l_deg = float(l_str)
                b_deg = float(b_str)
                gal = SkyCoord(l=l_deg * u.deg, b=b_deg * u.deg, frame="galactic")
                coord = gal.icrs
                label = f"l={l_deg:.1f} deg, b={b_deg:.1f} deg"

            elif mode == "name":
                name = self.name_entry.get().strip()
                if not name:
                    messagebox.showerror("Error", "Please enter a source name.")
                    return
                coord = SkyCoord.from_name(name)
                label = f"Source '{name}'"
            else:
                messagebox.showerror("Error", "Unknown coordinate mode.")
                return
        except Exception as e:
            messagebox.showerror("Error", f"Failed to interpret coordinates/source: {e}")
            return

        alt_deg, is_up, rise_time, set_time = compute_altitude_and_rise_set(
            coord, ATA_LOCATION, horizon_deg=18.0
        )

        # Format times in local time zone
        def fmt_time(t):
            if t is None:
                return "N/A"
            dt_utc = t.to_datetime(timezone=datetime.timezone.utc)
            dt_local = dt_utc.astimezone(HCRO_TZ)
            return dt_local.strftime("%Y-%m-%d %H:%M:%S %Z")

        rise_str = fmt_time(rise_time)
        set_str = fmt_time(set_time)

        if is_up:
            status_str = f"{label} is UP (alt={alt_deg:.1f} deg)"
            tag = "up"
            icon = "\u2705"  # green check
        else:
            status_str = f"{label} is NOT up (alt={alt_deg:.1f} deg)"
            tag = "down"
            icon = "\u274C"  # red X

        info = (
            f"{icon} {status_str}\n"
            f"   Next rise (>=18 deg): {rise_str}\n"
            f"   Next set (<=18 deg):  {set_str}"
        )

        self.source_info_label.configure(text=info)
        self.log(status_str, tag=tag)

    # ---------- Callbacks: Antenna ----------

    def on_antennas_changed(self, choice: str):
        """
        Callback when the 'Select Antenna(s)' dropdown changes.
        """
        self._update_antennas_from_selection(choice)
        self.log(f"Selected antennas: {antennas}")

    def on_reserve_clicked(self):
        self.run_with_progress("Reserving Antennas", ata_reserve_antennas)

    def on_park_clicked(self):
        self.run_with_progress("Parking Antennas", ata_park_antennas)

    def on_release_clicked(self):
        self.run_with_progress("Releasing Antennas", ata_release_antennas)

    def on_refresh_status_clicked(self):
        """
        Fetch and display the full ATA ASCII status (ac.get_ascii_status()).
        """
        def do_status():
            # Directly mirror: print(ac.get_ascii_status())
            full_status = ata_get_status()
            return [full_status]

        self.run_with_progress(
            "Fetching antenna status",
            do_status,
            log_in_status=True
        )

    def on_usrp_check_clicked(self):
        def do_check():
            return usrp_check_levels()

        self.run_with_progress(
            "Checking USRP Channel Levels",
            do_check
        )

    def on_usrp_reset_clicked(self):
        def do_reset():
            return usrp_reset_clocking()

        self.run_with_progress(
            "Resetting USRP Clocking / Channel Levels",
            do_reset
        )

    # ---------- Spectrum streaming ----------

    def start_data_stream(self):
        """Start the GNU Radio spectrum script and ZMQ subscriber."""
        if zmq is None:
            messagebox.showerror(
                "Missing Dependency",
                "pyzmq is not installed; cannot start data stream.",
            )
            return

        if self.spectrum_running:
            self.log("Data stream already running.")
            return

        address = self.spec_address_var.get().strip() or SPECTRUM_ZMQ_ADDRESS

        # Launch the GNU Radio script as a subprocess
        try:
            self.stream_process = subprocess.Popen(
                ["python", SPECTRUM_SCRIPT],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            self.log(f"Started spectrum script: {SPECTRUM_SCRIPT}")
        except Exception as e:
            messagebox.showerror(
                "Error", f"Failed to start spectrum script:\n{e}"
            )
            return

        # Start ZMQ subscriber in a background thread
        self.spectrum_stop_event.clear()
        self.spectrum_thread = threading.Thread(
            target=self._spectrum_worker, args=(address,), daemon=True
        )
        self.spectrum_thread.start()

        self.spectrum_running = True
        self.spectrum_status_label.configure(
            text=f"Data stream running ({address})", text_color="green"
        )

    def _spectrum_worker(self, address):
        """Worker thread: receive spectrum vectors via ZMQ and queue them."""
        ctx = None
        sock = None
        try:
            ctx = zmq.Context()
            sock = ctx.socket(zmq.SUB)
            sock.connect(address)
            sock.setsockopt(zmq.SUBSCRIBE, b"")

            self.spectrum_context = ctx
            self.spectrum_socket = sock

            while not self.spectrum_stop_event.is_set():
                try:
                    msg = sock.recv(flags=zmq.NOBLOCK)
                except zmq.Again:
                    time.sleep(0.01)
                    continue

                # Convert bytes -> float32 array
                try:
                    arr = np.frombuffer(msg, dtype=np.float32)
                except ValueError:
                    continue

                if arr.size == 0:
                    continue

                if arr.size != SPECTRUM_FFT_SIZE:
                    pass

                try:
                    if not self.spectrum_queue.full():
                        self.spectrum_queue.put_nowait(arr.copy())
                except queue.Full:
                    pass

        except Exception as e:
            self.log(f"Error in spectrum worker: {e}")
        finally:
            if sock is not None:
                sock.close(0)
            if ctx is not None:
                ctx.term()
            self.spectrum_running = False

    def stop_data_stream(self):
        """Stop ZMQ subscriber and the GNU Radio script process."""
        if not self.spectrum_running and self.stream_process is None:
            self.log("Data stream not running.")
            return

        self.spectrum_stop_event.set()
        if self.spectrum_thread is not None and self.spectrum_thread.is_alive():
            self.spectrum_thread.join(timeout=2.0)
        self.spectrum_thread = None

        if self.stream_process is not None:
            try:
                self.stream_process.terminate()
                try:
                    self.stream_process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    self.stream_process.kill()
                self.log("Spectrum script stopped.")
            except Exception as e:
                self.log(f"Error stopping spectrum script: {e}")
            self.stream_process = None

        self.spectrum_running = False
        self.spectrum_status_label.configure(
            text="Data stream not running", text_color="red"
        )
        self.log("Data stream stopped.")

    def _update_spectrum_plot(self):
        """Poll the queue and update the spectrum plot (runs in Tk thread)."""
        latest = None
        try:
            while True:
                latest = self.spectrum_queue.get_nowait()
        except queue.Empty:
            pass

        if latest is not None and latest.size > 0:
            n = latest.size
            samp_rate = SPECTRUM_SAMP_RATE

            # Center frequency from main Frequency entry (MHz)
            center_mhz = 0.0
            txt = self.freq_entry.get().strip()
            if txt:
                try:
                    center_mhz = float(txt)
                except ValueError:
                    center_mhz = 0.0
            else:
                # If empty but profile is non-custom, use that profile freq
                profile = self.selected_profile.get()
                freq_hz = OBSERVATION_PROFILES.get(profile)
                if freq_hz is not None:
                    center_mhz = freq_hz / 1e6
                else:
                    center_mhz = 0.0

            center_hz = center_mhz * 1e6
            freqs_hz = (center_hz - samp_rate / 2.0) + np.arange(n) * samp_rate / n
            freqs_mhz = freqs_hz / 1e6

            # Convert power to dB for plotting (assume linear PSD)
            psd_linear = latest.astype(np.float64)
            psd_db = 10.0 * np.log10(psd_linear + 1e-12)

            self.spec_line.set_data(freqs_mhz, psd_db)
            self.spec_ax.relim()
            self.spec_ax.autoscale_view()
            self.spec_canvas.draw_idle()

        # Schedule the next update (5 Hz)
        self.root.after(200, self._update_spectrum_plot)

    # ---------- Camera streaming ----------

    def start_camera(self):
        """Start the MJPEG camera stream and display it in the Camera tab."""
        if Image is None or ImageTk is None:
            messagebox.showerror(
                "Missing Dependency",
                "Pillow (PIL) is required for camera display.",
            )
            return

        if self.camera_running:
            self.log("Camera stream already running.")
            return

        self.camera_stop_event.clear()
        self.camera_running = True
        if self.camera_status_label is not None:
            self.camera_status_label.configure(
                text="Connecting to camera...", text_color="orange"
            )

        t = threading.Thread(target=self._camera_worker, daemon=True)
        self.camera_thread = t
        t.start()

    def _camera_worker(self):
        url = CAMERA_MJPEG_URL
        try:
            stream = urllib.request.urlopen(url, timeout=5)
        except Exception as e:
            # Capture the message now so we don't reference 'e' after the except block
            msg = f"Could not open MJPEG stream at {url}: {e}"
            self.root.after(0, lambda msg=msg: self._camera_error(msg))
            self.camera_running = False
            return


        bytes_buffer = b""
        while not self.camera_stop_event.is_set():
            try:
                chunk = stream.read(1024)
                if not chunk:
                    break
                bytes_buffer += chunk
                start = bytes_buffer.find(b"\xff\xd8")
                end = bytes_buffer.find(b"\xff\xd9", start + 2)
                if start != -1 and end != -1 and end > start:
                    jpg = bytes_buffer[start:end + 2]
                    bytes_buffer = bytes_buffer[end + 2:]
                    try:
                        img = Image.open(io.BytesIO(jpg))
                        img = img.convert("RGB")
                    except Exception:
                        continue

                    def update_image(image=img):
                        if self.camera_label is None:
                            return
                        # Resize to fit label, if we have a sensible size
                        w = self.camera_label.winfo_width()
                        h = self.camera_label.winfo_height()
                        if w > 20 and h > 20:
                            image_resized = image.resize((w, h))
                        else:
                            image_resized = image
                        photo = ImageTk.PhotoImage(image_resized)
                        self.camera_photo = photo  # keep a reference
                        self.camera_label.configure(image=photo)
                        if self.camera_status_label is not None:
                            self.camera_status_label.configure(
                                text=f"Camera streaming from {CAMERA_MJPEG_URL}",
                                text_color="green",
                            )

                    self.root.after(0, update_image)
            except Exception as e:
                self.root.after(0, lambda: self._camera_error(
                    f"Error reading camera stream: {e}"
                ))
                break

        try:
            stream.close()
        except Exception:
            pass
        self.camera_running = False
        self.root.after(0, self._camera_stopped_ui)

    def _camera_error(self, msg: str):
        self.log(msg)
        if self.camera_status_label is not None:
            self.camera_status_label.configure(
                text="Camera error", text_color="red"
            )
        messagebox.showerror("Camera Error", msg)

    def _camera_stopped_ui(self):
        if self.camera_status_label is not None:
            self.camera_status_label.configure(
                text="Camera not connected", text_color="red"
            )

    def stop_camera(self):
        if not self.camera_running:
            self.log("Camera stream is not running.")
            if self.camera_status_label is not None:
                self.camera_status_label.configure(
                    text="Camera not connected", text_color="red"
                )
            return

        self.camera_stop_event.set()
        if self.camera_thread is not None and self.camera_thread.is_alive():
            self.camera_thread.join(timeout=1.0)
        self.camera_thread = None
        self.camera_running = False
        if self.camera_status_label is not None:
            self.camera_status_label.configure(
                text="Camera stopped", text_color="red"
            )

    # ---------- Callbacks: Antenna ----------

    def on_save_log_clicked(self):
        """
        Save the current contents of the log window to a timestamped text file.
        """
        now = datetime.datetime.now(tz=HCRO_TZ)
        filename = now.strftime("ata_observation_log_%Y%m%d_%H%M%S.txt")
        filepath = os.path.join(os.getcwd(), filename)

        try:
            content = self.log_text.get("1.0", tk.END)
            with open(filepath, "w") as f:
                f.write(content)
            messagebox.showinfo(
                "Log Saved",
                f"Log saved to {filepath}"
            )
        except Exception as e:
            messagebox.showerror(
                "Error",
                f"Failed to save log: {e}"
            )

    # Main entry point

def main():
    root = tk.Tk()
    app = ATAObservationGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
