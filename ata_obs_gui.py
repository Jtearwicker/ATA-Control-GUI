import tkinter as tk
from tkinter import ttk, messagebox
import customtkinter

import datetime
import threading
import subprocess
import os
import numpy as np  # for rise/set sampling

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
# CONFIG / GLOBALS
# ======================================================

# Antenna set that can be chosen in the GUI (for now just 1a)
AVAILABLE_ANTENNAS = ['1a']

# Current antenna selection used by the control wrappers.
# This will be updated from the "Select antennas" dropdown.
antennas = AVAILABLE_ANTENNAS.copy()

# ATA site (approx Hat Creek / ATA)
ATA_LOCATION = EarthLocation.from_geodetic(
    lon=-121.47 * u.deg,
    lat=40.82 * u.deg,
    height=986 * u.m
)

# Local timezone for logs & display
LOCAL_TZ = ZoneInfo("America/Los_Angeles")

# Observation profiles: name -> center frequency in Hz (None = use manual entry)
OBSERVATION_PROFILES = {
    "Custom": None,
    "21cm Milky Way - 1420.405 MHz": 1420.405e6,
    "Tianwen-1 Mars Orbiter - 8431.0 MHz": 8431.0e6,
}

# Directory where USRP test scripts *might* live (currently unused but kept)
HOME_DIR = os.path.expanduser("~")


# ======================================================
# ASTRO UTILS
# ======================================================

def compute_altitude_and_rise_set(coord, location, horizon_deg=18.0, n_steps=240):
    """
    Given an ICRS SkyCoord, compute:
      - current altitude (deg)
      - whether it's currently above `horizon_deg`
      - approximate next rise and set times (crossings of horizon_deg)
        within the next 24 hours, using a simple sampling approach.

    Returns: (alt_now_deg, is_up, rise_time, set_time)
      alt_now_deg : float
      is_up       : bool
      rise_time   : astropy.time.Time or None
      set_time    : astropy.time.Time or None
    """
    now = Time.now()

    # Altitude now
    altaz_now = coord.transform_to(AltAz(obstime=now, location=location))
    alt_now_deg = altaz_now.alt.deg
    is_up = alt_now_deg >= horizon_deg

    # Sample altitudes over next 24 hours
    hours = np.linspace(0.0, 24.0, n_steps)
    times = now + hours * u.hour
    altaz_grid = coord.transform_to(AltAz(obstime=times, location=location))
    alt_grid = altaz_grid.alt.deg

    diff = alt_grid - horizon_deg
    rise_time = None
    set_time = None

    for i in range(len(times) - 1):
        d1 = diff[i]
        d2 = diff[i + 1]
        if d1 == 0:
            continue

        # Crossing from below to above -> rise
        if d1 < 0 and d2 > 0:
            frac = (horizon_deg - alt_grid[i]) / (alt_grid[i + 1] - alt_grid[i])
            t_cross = times[i] + frac * (times[i + 1] - times[i])
            if rise_time is None:
                rise_time = t_cross

        # Crossing from above to below -> set
        elif d1 > 0 and d2 < 0:
            frac = (horizon_deg - alt_grid[i]) / (alt_grid[i + 1] - alt_grid[i])
            t_cross = times[i] + frac * (times[i + 1] - times[i])
            if set_time is None:
                set_time = t_cross

    return alt_now_deg, is_up, rise_time, set_time


# ======================================================
# ATA / USRP CONTROL WRAPPERS
# ======================================================

def ata_reserve_antennas():
    """
    Reserve antennas into group 'atagr' if they are free.

    Uses ac.move_ant_group(antennas, 'none', 'atagr').
    Returns a descriptive status string.
    """
    # Check if already reserved in 'atagr'
    atagr_list = str(ac.list_antenna_group('atagr'))
    already_in_atagr = all(ant in atagr_list for ant in antennas)
    if already_in_atagr:
        return f"WARNING: {antennas} already reserved in group 'atagr'."

    # Check if in 'none' group (free)
    none_list = str(ac.list_antenna_group('none'))
    all_in_none = all(ant in none_list for ant in antennas)
    if not all_in_none:
        return (
            f"WARNING: Some antennas {antennas} are not in 'none' group. "
            "Cannot safely move to 'atagr'."
        )

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
    Release antennas from 'atagr' back to 'none'.

    Uses ac.move_ant_group(antennas, 'atagr', 'none').
    """
    # If already in 'none', do nothing
    none_list = str(ac.list_antenna_group('none'))
    all_in_none = all(ant in none_list for ant in antennas)
    if all_in_none:
        return f"Antennas {antennas} already in group 'none' (released)."

    # Otherwise, move from 'atagr' to 'none'
    ac.move_ant_group(antennas, 'atagr', 'none')
    return f"Antennas {antennas} moved from 'atagr' to 'none' (released)."


def ata_get_status_text():
    """
    Return a multi-line string describing current ATA status.

    Uses ac.get_ascii_status().
    """
    return ac.get_ascii_status()


def ata_set_frequency(freq_mhz):
    """
    Set observing frequency to freq_mhz (MHz) only, without changing attenuation
    or RF switch settings.
    Uses:
        ac.set_freq(freq_mhz, antennas, lo)
    """
    lo = "d"
    ac.set_freq(freq_mhz, antennas, lo)
    return f"Frequency set to {freq_mhz:.6f} MHz on antennas {antennas} with LO '{lo}'."


def ata_autotune():
    """
    Run autotune on the current antennas.

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
    return (
        f"Tracking RA = {ra_hours:.6f} h, Dec = {dec_deg:.6f} deg "
        f"on antennas {antennas}."
    )


def usrp_check_levels():
    """
    Run 'usrp_test.py' and return only the channel level lines and the note.

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
            shell=True     # use shell so PATH resolution matches your CLI
        )
    except Exception as e:
        return [f"Error running usrp_test.py: {e}"]

    stdout_lines = proc.stdout.splitlines() if proc.stdout else []
    stderr_lines = proc.stderr.splitlines() if proc.stderr else []

    # Filter out the interesting lines from stdout
    filtered = [
        line for line in stdout_lines
        if line.startswith("Channel ") or line.startswith("Note:")
    ]

    if filtered:
        return filtered

    # If no interesting lines, fall back to stderr or a generic message
    if stderr_lines:
        return ["usrp_test.py stderr:"] + stderr_lines

    return [f"usrp_test.py finished with return code {proc.returncode}, no output."]


def usrp_reset_clocking():
    """
    Run 'usrp_reset_clocking.py'.

    On success (return code 0), ignore stdout and just return:
        'USRP clocking reset success.'

    On failure, return stderr lines if available, otherwise a generic
    non-zero return-code message.
    """
    try:
        proc = subprocess.run(
            "usrp_reset_clocking.py",
            text=True,
            capture_output=True,
            check=False,
            shell=True
        )
    except Exception as e:
        return [f"Error running usrp_reset_clocking.py: {e}"]

    stderr_lines = proc.stderr.splitlines() if proc.stderr else []

    if proc.returncode == 0:
        # Success: don't spam UHD chatter, just report success.
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
        self.root.title("ATA Observation GUI")
        self.root.geometry("1400x900")

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

        self._build_layout()

        # Ensure global antennas initialised from dropdown
        self._update_antennas_from_selection(self.antenna_select_var.get())

    # ---------- Layout ----------

    def _build_layout(self):
        # Main vertical layout: top main_frame + bottom log_frame
        self.main_frame = customtkinter.CTkFrame(master=self.root)
        self.main_frame.pack(side="top", fill="both", expand=True)

        self.log_frame = customtkinter.CTkFrame(master=self.root, height=150)
        self.log_frame.pack(side="bottom", fill="x")

        # Inside main_frame: left controls + right side (status)
        self.left_frame = customtkinter.CTkFrame(
            master=self.main_frame, width=360
        )
        self.left_frame.pack(
            side="left", fill="y", padx=10, pady=10
        )

        self.right_frame = customtkinter.CTkFrame(master=self.main_frame)
        self.right_frame.pack(
            side="right", fill="both", expand=True, padx=10, pady=10
        )

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

        # Label + antenna selection dropdown
        antenna_select_label = customtkinter.CTkLabel(
            button_frame, text="Select Antenna(s)"
        )
        antenna_select_label.pack(side="left", padx=(0, 5), pady=5)

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

        # ---- USRP backend control ----
        usrp_frame = customtkinter.CTkFrame(parent)
        usrp_frame.pack(fill="x", pady=(0, 10))

        usrp_label = customtkinter.CTkLabel(
            usrp_frame,
            text="USRP Backend",
            font=("Arial", 16, "bold")
        )
        usrp_label.pack(anchor="w", pady=(5, 5))

        usrp_btn_frame = customtkinter.CTkFrame(usrp_frame)
        usrp_btn_frame.pack(fill="x", pady=(0, 5))

        usrp_test_btn = customtkinter.CTkButton(
            usrp_btn_frame,
            text="Check channel levels",
            command=self.on_usrp_check_clicked
        )
        usrp_test_btn.pack(side="left", padx=5, pady=5)

        usrp_reset_btn = customtkinter.CTkButton(
            usrp_btn_frame,
            text="Reset channel clocking",
            command=self.on_usrp_reset_clicked
        )
        usrp_reset_btn.pack(side="left", padx=5, pady=5)

        # ---- Frequency / profile ----
        freq_frame = customtkinter.CTkFrame(parent)
        freq_frame.pack(fill="x", pady=(5, 10))

        freq_label = customtkinter.CTkLabel(
            freq_frame,
            text="Observation Profile / Frequency",
            font=("Arial", 16, "bold")
        )
        freq_label.pack(anchor="w", pady=(0, 5))

        profile_combo = customtkinter.CTkComboBox(
            freq_frame,
            values=list(OBSERVATION_PROFILES.keys()),
            variable=self.selected_profile,
            command=self.on_profile_changed
        )
        profile_combo.pack(fill="x", pady=(0, 5))

        # Frequency row (no attenuation)
        freq_subframe = customtkinter.CTkFrame(freq_frame)
        freq_subframe.pack(fill="x")

        freq_name_label = customtkinter.CTkLabel(freq_subframe, text="Frequency")
        freq_name_label.pack(side="left", padx=(0, 5), pady=5)

        # BLANK by default
        self.freq_entry = customtkinter.CTkEntry(
            freq_subframe,
            placeholder_text="",
            width=100
        )
        self.freq_entry.pack(side="left", padx=(0, 5), pady=5)

        freq_unit_label = customtkinter.CTkLabel(freq_subframe, text="MHz")
        freq_unit_label.pack(side="left", padx=(0, 10), pady=5)

        apply_freq_btn = customtkinter.CTkButton(
            freq_subframe,
            text="Set frequency",
            command=self.on_set_frequency_clicked
        )
        apply_freq_btn.pack(side="left", padx=5, pady=5)

        # separate autotune button
        autotune_btn = customtkinter.CTkButton(
            freq_frame,
            text="Autotune",
            command=self.on_autotune_clicked
        )
        autotune_btn.pack(pady=(5, 0))

        # ---- Coordinate / pointing ----
        coord_frame = customtkinter.CTkFrame(parent)
        coord_frame.pack(fill="x", pady=(10, 10))

        coord_label = customtkinter.CTkLabel(
            coord_frame,
            text="Pointing / Coordinates",
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
            ("Source name", "name"),
        ]:
            rb = customtkinter.CTkRadioButton(
                mode_frame,
                text=text,
                variable=self.coord_mode,
                value=mode,
                command=self.on_coord_mode_changed
            )
            rb.pack(side="left", padx=(0, 5))

        # Coordinate input frame: left = fields, right = button
        coord_input_frame = customtkinter.CTkFrame(coord_frame)
        coord_input_frame.pack(fill="x", pady=(5, 5))

        # Left side: labels + entries
        self.coord_fields_frame = customtkinter.CTkFrame(coord_input_frame)
        self.coord_fields_frame.pack(side="left", fill="x", expand=True)

        self.coord1_label = customtkinter.CTkLabel(
            self.coord_fields_frame, text="RA (h:m:s)"
        )
        self.coord1_label.grid(
            row=0, column=0, sticky="w", padx=(0, 5), pady=2
        )

        self.coord1_entry = customtkinter.CTkEntry(
            self.coord_fields_frame,
            width=140,
            placeholder_text="12:34:56"
        )
        self.coord1_entry.grid(
            row=0, column=1, sticky="w", padx=(0, 5), pady=2
        )

        self.coord2_label = customtkinter.CTkLabel(
            self.coord_fields_frame, text="Dec (d:m:s)"
        )
        self.coord2_label.grid(
            row=1, column=0, sticky="w", padx=(0, 5), pady=2
        )

        self.coord2_entry = customtkinter.CTkEntry(
            self.coord_fields_frame,
            width=140,
            placeholder_text="+12:34:56"
        )
        self.coord2_entry.grid(
            row=1, column=1, sticky="w", padx=(0, 5), pady=2
        )

        self.name_label = customtkinter.CTkLabel(
            self.coord_fields_frame, text="Source name"
        )
        self.name_entry = customtkinter.CTkEntry(
            self.coord_fields_frame,
            width=180,
            placeholder_text="3C286"
        )
        # name_label/name_entry will be gridded only in "name" mode

        # Right side: fixed-size "Check source location" button
        self.check_source_btn = customtkinter.CTkButton(
            coord_input_frame,
            text="Check source location",
            command=self.on_check_source_location_clicked,
            width=180
        )
        self.check_source_btn.pack(side="right", padx=(10, 0), pady=2)

        # Track / Park buttons
        track_button_frame = customtkinter.CTkFrame(coord_frame)
        track_button_frame.pack(fill="x", pady=(5, 0))

        track_btn = customtkinter.CTkButton(
            track_button_frame,
            text="Track target",
            command=self.on_track_clicked
        )
        track_btn.pack(side="left", padx=5, pady=5)

        park_btn = customtkinter.CTkButton(
            track_button_frame,
            text="Park Antenna(s)",
            command=self.on_park_clicked
        )
        park_btn.pack(side="left", padx=5, pady=5)

        # Time / LST info
        self.info_label = customtkinter.CTkLabel(
            coord_frame, text="", font=("Arial", 12)
        )
        self.info_label.pack(anchor="w", pady=(5, 0))
        self._update_time_info()

        # Make initial pointing UI consistent with default mode
        self.on_coord_mode_changed()

    def _build_right_tabs(self, parent):
        notebook = ttk.Notebook(parent)
        notebook.pack(fill="both", expand=True)

        # Status tab only
        status_frame = customtkinter.CTkFrame(notebook)
        notebook.add(status_frame, text="Status")

        self.status_text = tk.Text(
            status_frame, wrap="word", height=20, width=80
        )
        self.status_text.pack(
            fill="both", expand=True, padx=5, pady=5
        )

        refresh_btn = customtkinter.CTkButton(
            status_frame,
            text="Show antenna status",
            command=self.on_refresh_status_clicked
        )
        refresh_btn.pack(pady=5)

    def _build_log_area(self, parent):
        # Progress indicator + status line
        top_frame = customtkinter.CTkFrame(parent)
        top_frame.pack(fill="x", pady=(5, 0))

        self.progress_label = customtkinter.CTkLabel(
            top_frame, text="", anchor="w"
        )
        self.progress_label.pack(side="left", padx=5, pady=5)

        self.progressbar = ttk.Progressbar(
            top_frame, mode="indeterminate", length=200
        )
        self.progressbar.pack(side="right", padx=5, pady=5)

        # Log text
        self.log_text = tk.Text(parent, wrap="word", height=8)
        self.log_text.pack(
            fill="both", expand=True, padx=5, pady=(0, 5)
        )

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

    def _update_antennas_from_selection(self, selection_text: str):
        """
        Update the global 'antennas' list from the selection string.
        """
        global antennas
        antennas = self._parse_antenna_selection(selection_text)

    def log(self, msg):
        # Local time in Pacific, with timezone label
        now_local = datetime.datetime.now(LOCAL_TZ)
        timestamp = now_local.strftime("%Y-%m-%d %H:%M:%S %Z")
        line = f"[{timestamp}] {msg}\n"

        # Insert newest log at the top
        self.log_text.insert("1.0", line)
        self.log_text.see("1.0")
        self.log_text.update_idletasks()

    def run_with_progress(self, description, func, callback=None, log_result=True):
        """
        Run `func` in a background thread while showing an
        indeterminate progress bar. `func` returns a string or a
        list/tuple of strings to log. If callback is provided, it will
        be called with the result on the main thread after logging
        (unless log_result is False, in which case only the callback
        handles the result).
        """

        def worker():
            try:
                result = func()
                error = None
            except Exception as e:
                result = None
                error = e

            def on_done():
                self.progressbar.stop()
                self.progress_label.configure(text="")
                if error is not None:
                    # Errors still get logged regardless of log_result
                    self.log(f"{description} FAILED: {error}")
                    messagebox.showerror(
                        "Error", f"{description} failed:\n{error}"
                    )
                else:
                    if log_result:
                        if isinstance(result, str):
                            self.log(result)
                        elif isinstance(result, (list, tuple)):
                            for line in result:
                                self.log(line)
                        else:
                            self.log(f"{description} completed.")
                    if callback is not None:
                        callback(result)

            self.root.after(0, on_done)

        self.progress_label.configure(text=description)
        self.progressbar.start(10)
        threading.Thread(target=worker, daemon=True).start()

    def _update_time_info(self):
        """
        Show local civil time (Pacific) and LST (Local Sidereal Time).
        """
        try:
            now_astropy = Time.now()
            lst = now_astropy.sidereal_time('apparent', longitude=ATA_LOCATION.lon)
            local_now = datetime.datetime.now(LOCAL_TZ)
            text = (
                f"Local time (America/Los_Angeles): "
                f"{local_now.strftime('%Y-%m-%d %H:%M:%S %Z')}   |   "
                f"LST (Local Sidereal Time): "
                f"{lst.to_string(unit=u.hour, sep=':', precision=0, pad=True)}"
            )
        except Exception as e:
            text = f"Time error: {e}"

        self.info_label.configure(text=text)
        # re-schedule
        self.root.after(10000, self._update_time_info)

    # ---------- Callbacks: Antenna ----------

    def on_antennas_changed(self, choice: str):
        """
        Callback when the 'Select antennas' dropdown changes.
        """
        self._update_antennas_from_selection(choice)
        self.log(f"Selected antennas: {antennas}")

    def on_reserve_clicked(self):
        self.run_with_progress("Reserving antennas", ata_reserve_antennas)

    def on_park_clicked(self):
        self.run_with_progress("Parking antennas", ata_park_antennas)

    def on_release_clicked(self):
        self.run_with_progress("Releasing antennas", ata_release_antennas)

    # ---------- Callbacks: USRP backend ----------

    def on_usrp_check_clicked(self):
        def do_check():
            return usrp_check_levels()

        self.run_with_progress(
            "Checking USRP channel levels",
            do_check
        )

    def on_usrp_reset_clicked(self):
        def do_reset():
            return usrp_reset_clocking()

        self.run_with_progress(
            "Resetting USRP clocking / channel levels",
            do_reset
        )

    # ---------- Callbacks: Frequency / profiles ----------

    def on_profile_changed(self, *_):
        profile = self.selected_profile.get()
        freq_hz = OBSERVATION_PROFILES.get(profile)
        if freq_hz is not None:
            self.freq_entry.delete(0, tk.END)
            self.freq_entry.insert(0, f"{freq_hz / 1e6:.6f}")
        else:
            # Custom: leave frequency blank
            self.freq_entry.delete(0, tk.END)

    def on_set_frequency_clicked(self):
        text = self.freq_entry.get().strip()
        if not text:
            messagebox.showwarning(
                "No frequency", "Please enter a frequency in MHz."
            )
            return
        try:
            freq_mhz = float(text)
        except ValueError:
            messagebox.showerror(
                "Invalid frequency",
                f"Could not parse '{text}' as a number."
            )
            return

        def do_set():
            return ata_set_frequency(freq_mhz)

        self.run_with_progress(
            f"Setting frequency to {freq_mhz:.6f} MHz",
            do_set
        )

    def on_autotune_clicked(self):
        self.run_with_progress("Autotuning", ata_autotune)

    # ---------- Callbacks: Coordinate modes ----------

    def on_coord_mode_changed(self):
        """
        Handle switching between RA/Dec, Alt/Az, Galactic, and Source-name modes.

        - When switching away from a coordinate mode, cache the current entries.
        - When switching back to that mode, restore the cached values.
        - RA/Dec / Alt/Az / Galactic share the two numeric entry boxes.
        - Source-name mode hides those boxes and shows only the name entry.
        """
        old_mode = self.last_coord_mode
        new_mode = self.coord_mode.get()

        # Save current entries for the old coordinate mode (if applicable)
        if old_mode in self.coord_cache:
            self.coord_cache[old_mode][0] = self.coord1_entry.get().strip()
            self.coord_cache[old_mode][1] = self.coord2_entry.get().strip()

        # Configure UI for the new mode
        if new_mode == "radec":
            # Show numeric fields
            self.coord1_label.grid(row=0, column=0, sticky="w", padx=(0, 5), pady=2)
            self.coord1_entry.grid(row=0, column=1, sticky="w", padx=(0, 5), pady=2)
            self.coord2_label.grid(row=1, column=0, sticky="w", padx=(0, 5), pady=2)
            self.coord2_entry.grid(row=1, column=1, sticky="w", padx=(0, 5), pady=2)
            # Hide source-name row
            self.name_label.grid_remove()
            self.name_entry.grid_remove()

            self.coord1_label.configure(text="RA (h:m:s)")
            self.coord2_label.configure(text="Dec (d:m:s)")
            self.coord1_entry.configure(placeholder_text="12:34:56")
            self.coord2_entry.configure(placeholder_text="+12:34:56")

        elif new_mode == "altaz":
            self.coord1_label.grid(row=0, column=0, sticky="w", padx=(0, 5), pady=2)
            self.coord1_entry.grid(row=0, column=1, sticky="w", padx=(0, 5), pady=2)
            self.coord2_label.grid(row=1, column=0, sticky="w", padx=(0, 5), pady=2)
            self.coord2_entry.grid(row=1, column=1, sticky="w", padx=(0, 5), pady=2)
            self.name_label.grid_remove()
            self.name_entry.grid_remove()

            self.coord1_label.configure(text="Alt (deg)")
            self.coord2_label.configure(text="Az (deg)")
            self.coord1_entry.configure(placeholder_text="45.0")
            self.coord2_entry.configure(placeholder_text="180.0")

        elif new_mode == "galactic":
            self.coord1_label.grid(row=0, column=0, sticky="w", padx=(0, 5), pady=2)
            self.coord1_entry.grid(row=0, column=1, sticky="w", padx=(0, 5), pady=2)
            self.coord2_label.grid(row=1, column=0, sticky="w", padx=(0, 5), pady=2)
            self.coord2_entry.grid(row=1, column=1, sticky="w", padx=(0, 5), pady=2)
            self.name_label.grid_remove()
            self.name_entry.grid_remove()

            self.coord1_label.configure(text="l (deg)")
            self.coord2_label.configure(text="b (deg)")
            self.coord1_entry.configure(placeholder_text="120.0")
            self.coord2_entry.configure(placeholder_text="-5.0")

        elif new_mode == "name":
            # Hide numeric fields entirely; show only source-name entry
            self.coord1_label.grid_remove()
            self.coord1_entry.grid_remove()
            self.coord2_label.grid_remove()
            self.coord2_entry.grid_remove()

            self.name_label.grid(row=0, column=0, sticky="w", padx=(0, 5), pady=2)
            self.name_entry.grid(row=0, column=1, sticky="w", padx=(0, 5), pady=2)

        # Restore cached coordinate entries for the new mode (if it uses them)
        if new_mode in self.coord_cache:
            c1, c2 = self.coord_cache[new_mode]
            self.coord1_entry.delete(0, tk.END)
            self.coord1_entry.insert(0, c1)
            self.coord2_entry.delete(0, tk.END)
            self.coord2_entry.insert(0, c2)
        else:
            # For 'name' mode, numeric entries are hidden; leave their contents as-is.
            pass

        self.last_coord_mode = new_mode

    def _resolve_target_for_visibility(self, mode, coord1, coord2, name):
        """
        Resolve the user inputs into an ICRS SkyCoord plus a human-readable label.
        Used by the 'Check source location' feature.

        For 'Source name':
          - If the name matches a Solar System body, use ephemerides (get_body)
            at the current time and ATA location.
          - Otherwise, use SkyCoord.from_name (Simbad/etc.).
        """
        if mode == "name":
            if not name:
                raise ValueError("Please enter a source name to look up.")

            name_lower = name.strip().lower()
            body_map = {
                "sun": "sun",
                "moon": "moon",
                "mercury": "mercury",
                "venus": "venus",
                "mars": "mars",
                "jupiter": "jupiter",
                "saturn": "saturn",
                "uranus": "uranus",
                "neptune": "neptune",
            }

            if name_lower in body_map:
                # Solar System body: use ephemeris-based position now, at ATA
                now = Time.now()
                body = get_body(body_map[name_lower], now, ATA_LOCATION)
                c = body.icrs
                label = f"Solar-system body '{name}'"
                return c, label
            else:
                # Generic name: Simbad / Vizier resolver
                c = SkyCoord.from_name(name)
                label = f"Source '{name}'"
                return c.icrs, label

        if mode == "radec":
            if not coord1 or not coord2:
                raise ValueError("RA and Dec are required.")
            c = SkyCoord(
                coord1, coord2,
                unit=(u.hourangle, u.deg),
                frame="icrs"
            )
            label = f"RA={coord1}, Dec={coord2}"
            return c.icrs, label

        if mode == "altaz":
            if not coord1 or not coord2:
                raise ValueError("Alt and Az are required.")
            alt = float(coord1)
            az = float(coord2)
            now = Time.now()
            altaz = AltAz(
                obstime=now,
                location=ATA_LOCATION,
                alt=alt * u.deg,
                az=az * u.deg
            )
            c = altaz.transform_to("icrs")
            label = f"Alt={alt} deg, Az={az} deg"
            return c.icrs, label

        if mode == "galactic":
            if not coord1 or not coord2:
                raise ValueError("l and b are required.")
            l = float(coord1)
            b = float(coord2)
            c = SkyCoord(
                l=l * u.deg, b=b * u.deg, frame="galactic"
            ).transform_to("icrs")
            label = f"l={l} deg, b={b} deg"
            return c.icrs, label

        raise ValueError(f"Unknown coordinate mode '{mode}'.")

    def on_check_source_location_clicked(self):
        """
        Check whether the current coordinates / source name are up,
        and compute approximate rise/set times (above 18°).
        """
        mode = self.coord_mode.get()
        coord1 = self.coord1_entry.get().strip()
        coord2 = self.coord2_entry.get().strip()
        name = self.name_entry.get().strip()

        try:
            coord_icrs, label = self._resolve_target_for_visibility(
                mode, coord1, coord2, name
            )
        except ValueError as e:
            messagebox.showerror("Coordinate error", str(e))
            return
        except Exception as e:
            messagebox.showerror(
                "Coordinate error",
                f"Could not interpret coordinates or resolve source:\n{e}"
            )
            return

        def do_check():
            alt_now_deg, is_up, rise_time, set_time = compute_altitude_and_rise_set(
                coord_icrs, ATA_LOCATION, horizon_deg=18.0
            )
            lines = []

            # Explicitly log the ICRS RA/Dec being used
            ra_str = coord_icrs.ra.to_string(unit=u.hour, sep=':', precision=2, pad=True)
            dec_str = coord_icrs.dec.to_string(unit=u.deg, sep=':', precision=2, alwayssign=True, pad=True)
            lines.append(f"{label}: ICRS coordinates = RA {ra_str}, Dec {dec_str}")

            lines.append(f"{label}: Elevation now = {alt_now_deg:.2f} deg")

            if is_up:
                lines.append(f"{label} is UP (above 18°).")
            else:
                lines.append(f"{label} is NOT up (below 18°).")

            # Rise time message (local time in America/Los_Angeles)
            if rise_time is not None:
                rt_local = rise_time.to_datetime(timezone=LOCAL_TZ)
                lines.append(
                    f"Next rise above 18°: "
                    f"{rt_local.strftime('%Y-%m-%d %H:%M:%S %Z')}"
                )
            else:
                if is_up:
                    lines.append(
                        "No next rise above 18° in the next 24 hours "
                        "(source already above threshold for this window)."
                    )
                else:
                    lines.append(
                        "No rise above 18° in the next 24 hours "
                        "(source stays below threshold)."
                    )

            # Set time message (local time in America/Los_Angeles)
            if set_time is not None:
                st_local = set_time.to_datetime(timezone=LOCAL_TZ)
                lines.append(
                    f"Next set below 18°: "
                    f"{st_local.strftime('%Y-%m-%d %H:%M:%S %Z')}"
                )
            else:
                if is_up:
                    lines.append(
                        "No set below 18° in the next 24 hours "
                        "(source stays above threshold)."
                    )
                else:
                    lines.append(
                        "No next set below 18° in the next 24 hours "
                        "(source already below threshold for this window)."
                    )

            return lines

        self.run_with_progress(
            f"Checking source visibility for {label}", do_check
        )

    def on_track_clicked(self):
        mode = self.coord_mode.get()
        coord1 = self.coord1_entry.get().strip()
        coord2 = self.coord2_entry.get().strip()
        name = self.name_entry.get().strip()

        # Source-name mode: use ATA catalog directly, no RA/Dec required here
        if mode == "name":
            if not name:
                messagebox.showerror(
                    "Source name error",
                    "Please enter a source name to look up in the ATA catalog."
                )
                return

            def do_track_name():
                ac.track_source(antennas, source=name)
                return f"Tracking catalog source '{name}' on antennas {antennas}."

            self.run_with_progress(
                f"Tracking catalog source '{name}'", do_track_name
            )
            return

        # All other modes use Astropy coordinate transforms
        try:
            if mode == "radec":
                if not coord1 or not coord2:
                    raise ValueError("RA and Dec are required.")
                c = SkyCoord(
                    coord1, coord2,
                    unit=(u.hourangle, u.deg),
                    frame="icrs"
                )
            elif mode == "altaz":
                if not coord1 or not coord2:
                    raise ValueError("Alt and Az are required.")
                alt = float(coord1)
                az = float(coord2)
                now = Time.now()
                altaz = AltAz(obstime=now, location=ATA_LOCATION,
                              alt=alt * u.deg, az=az * u.deg)
                c = altaz.transform_to("icrs")
            elif mode == "galactic":
                if not coord1 or not coord2:
                    raise ValueError("l and b are required.")
                l = float(coord1)
                b = float(coord2)
                c = SkyCoord(
                    l=l * u.deg, b=b * u.deg, frame="galactic"
                ).transform_to("icrs")
            else:
                raise ValueError(f"Unknown coordinate mode '{mode}'.")

        except Exception as e:
            messagebox.showerror(
                "Coordinate error",
                f"Could not interpret coordinates:\n{e}"
            )
            return

        ra_deg = c.ra.deg
        dec_deg = c.dec.deg

        def do_track():
            return ata_track_radec_degrees(ra_deg, dec_deg)

        self.run_with_progress(
            f"Tracking RA={c.ra.to_string(unit=u.hour, sep=':')}, "
            f"Dec={c.dec.to_string(unit=u.deg, sep=':')}",
            do_track
        )

    # ---------- Callbacks: Status ----------

    def on_refresh_status_clicked(self):
        def do_status():
            return ata_get_status_text()

        def update_status(text):
            if not isinstance(text, str):
                return
            # Full status in text box only (no logging to bottom console)
            self.status_text.delete("1.0", tk.END)
            self.status_text.insert("1.0", text)

        self.run_with_progress(
            "Refreshing ATA status",
            do_status,
            callback=update_status,
            log_result=False
        )


# ======================================================
# MAIN
# ======================================================

if __name__ == "__main__":
    root = tk.Tk()
    app = ATAObservationGUI(root)
    root.mainloop()
