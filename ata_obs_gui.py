import tkinter as tk
from tkinter import ttk, messagebox
import customtkinter
from tkinterweb import HtmlFrame

import datetime
import threading
from zoneinfo import ZoneInfo  # for local time

import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy.time import Time

from ATATools import ata_control as ac


# ======================================================
# CONFIG / GLOBALS
# ======================================================

# For now, single-antenna operation; change here later to generalize
antennas = ['1a']

# ATA site (approx Hat Creek / ATA)
ATA_LOCATION = EarthLocation.from_geodetic(
    lon=-121.47 * u.deg,
    lat=40.82 * u.deg,
    height=986 * u.m
)

# Local timezone for logs & display
LOCAL_TZ = ZoneInfo("America/Los_Angeles")

# Simple profiles: name -> center frequency in Hz (None = use manual entry)
OBSERVATION_PROFILES = {
    "21cm Milky Way (1420.405 MHz)": 1420.405e6,
    "Mars @ 8431 MHz": 8431e6,
    "Custom": None,
}

# Original camera page + direct MJPEG stream
CAMERA_PAGE_URL = (
    "http://10.3.0.30/camera/index.html"
    "?id=342&imagepath=%2Fmjpg%2Fvideo.mjpg&size=1#/video"
)
CAMERA_STREAM_URL = "http://10.3.0.30/mjpg/video.mjpg"


# ======================================================
# ATA CONTROL WRAPPERS
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


def ata_set_frequency_hz(freq_hz):
    """
    Set observing frequency to freq_hz (float, Hz) and autotune.

    Uses ac.set_freq(freq_mhz, antennas, 'd') and ac.autotune(antennas).
    """
    freq_mhz = freq_hz / 1e6
    # RFCB LO 'd' as in the hydrogen-line notebook
    ac.set_freq(freq_mhz, antennas, 'd')
    ac.autotune(antennas)
    return (
        f"Frequency set to {freq_mhz:.6f} MHz on antennas {antennas} using LO 'd', "
        "and autotuned."
    )


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


def ata_stop_tracking():
    """
    There is no explicit 'stop tracking' in the notebook examples.

    For now, this is a no-op with a log message to keep behavior explicit.
    If you later define a proper stop command, wire it here.
    """
    return "No explicit stop-tracking command defined; use 'Park' to stop motion."


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

        # Selected profile
        self.selected_profile = tk.StringVar(
            value=list(OBSERVATION_PROFILES.keys())[0]
        )

        # Camera URLs
        self.camera_url = CAMERA_PAGE_URL
        self.camera_stream_url = CAMERA_STREAM_URL

        self._build_layout()

    # ---------- Layout ----------

    def _build_layout(self):
        # Main vertical layout: top main_frame + bottom log_frame
        self.main_frame = customtkinter.CTkFrame(master=self.root)
        self.main_frame.pack(side="top", fill="both", expand=True)

        self.log_frame = customtkinter.CTkFrame(master=self.root, height=150)
        self.log_frame.pack(side="bottom", fill="x")

        # Inside main_frame: left controls + right side (status/camera)
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
            text=f"Antenna Control {antennas}",
            font=("Arial", 18, "bold")
        )
        antenna_label.pack(pady=(5, 5))

        button_frame = customtkinter.CTkFrame(parent)
        button_frame.pack(fill="x", pady=(0, 10))

        reserve_btn = customtkinter.CTkButton(
            button_frame,
            text="Reserve",
            command=self.on_reserve_clicked
        )
        reserve_btn.pack(side="left", padx=5, pady=5)

        park_btn = customtkinter.CTkButton(
            button_frame,
            text="Park",
            command=self.on_park_clicked
        )
        park_btn.pack(side="left", padx=5, pady=5)

        release_btn = customtkinter.CTkButton(
            button_frame,
            text="Release",
            command=self.on_release_clicked
        )
        release_btn.pack(side="left", padx=5, pady=5)

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

        freq_subframe = customtkinter.CTkFrame(freq_frame)
        freq_subframe.pack(fill="x")

        self.freq_entry = customtkinter.CTkEntry(
            freq_subframe,
            placeholder_text="Frequency [MHz]",
            width=140
        )
        self.freq_entry.pack(side="left", padx=(0, 5), pady=5)

        apply_freq_btn = customtkinter.CTkButton(
            freq_subframe,
            text="Set Frequency",
            command=self.on_set_frequency_clicked
        )
        apply_freq_btn.pack(side="left", padx=5, pady=5)

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

        coord_input_frame = customtkinter.CTkFrame(coord_frame)
        coord_input_frame.pack(fill="x", pady=(5, 5))

        self.coord1_label = customtkinter.CTkLabel(
            coord_input_frame, text="RA (h:m:s)"
        )
        self.coord1_label.grid(
            row=0, column=0, sticky="w", padx=(0, 5), pady=2
        )

        self.coord1_entry = customtkinter.CTkEntry(
            coord_input_frame,
            width=140,
            placeholder_text="e.g. 12:34:56"
        )
        self.coord1_entry.grid(
            row=0, column=1, sticky="w", padx=(0, 5), pady=2
        )

        self.coord2_label = customtkinter.CTkLabel(
            coord_input_frame, text="Dec (d:m:s)"
        )
        self.coord2_label.grid(
            row=1, column=0, sticky="w", padx=(0, 5), pady=2
        )

        self.coord2_entry = customtkinter.CTkEntry(
            coord_input_frame,
            width=140,
            placeholder_text="+12:34:56"
        )
        self.coord2_entry.grid(
            row=1, column=1, sticky="w", padx=(0, 5), pady=2
        )

        self.name_label = customtkinter.CTkLabel(
            coord_input_frame, text="Source name"
        )
        self.name_label.grid(
            row=2, column=0, sticky="w", padx=(0, 5), pady=2
        )

        self.name_entry = customtkinter.CTkEntry(
            coord_input_frame,
            width=180,
            placeholder_text="e.g. 3C286"
        )
        self.name_entry.grid(
            row=2, column=1, sticky="w", padx=(0, 5), pady=2
        )

        # Track / stop buttons
        track_button_frame = customtkinter.CTkFrame(coord_frame)
        track_button_frame.pack(fill="x", pady=(5, 0))

        track_btn = customtkinter.CTkButton(
            track_button_frame,
            text="Track target",
            command=self.on_track_clicked
        )
        track_btn.pack(side="left", padx=5, pady=5)

        stop_track_btn = customtkinter.CTkButton(
            track_button_frame,
            text="Stop tracking",
            command=self.on_stop_tracking_clicked
        )
        stop_track_btn.pack(side="left", padx=5, pady=5)

        # Time / LST info
        self.info_label = customtkinter.CTkLabel(
            coord_frame, text="", font=("Arial", 12)
        )
        self.info_label.pack(anchor="w", pady=(5, 0))
        self._update_time_info()

        self.root.after(10000, self._update_time_info)  # periodic updates

    def _build_right_tabs(self, parent):
        notebook = ttk.Notebook(parent)
        notebook.pack(fill="both", expand=True)

        # Status tab
        status_frame = customtkinter.CTkFrame(notebook)
        notebook.add(status_frame, text="Status")

        # Parsed table for antenna 1a (new)
        self.status_tree = ttk.Treeview(status_frame, show="headings", height=1)
        self.status_tree.pack(fill="x", padx=5, pady=5)

        self.status_text = tk.Text(
            status_frame, wrap="word", height=20, width=80
        )
        self.status_text.pack(
            fill="both", expand=True, padx=5, pady=5
        )

        refresh_btn = customtkinter.CTkButton(
            status_frame,
            text="Refresh status",
            command=self.on_refresh_status_clicked
        )
        refresh_btn.pack(pady=5)

        # Camera tab
        camera_frame = customtkinter.CTkFrame(notebook)
        notebook.add(camera_frame, text="Camera")

        self.camera_html_frame = HtmlFrame(camera_frame)
        self.camera_html_frame.pack(
            fill="both", expand=True, padx=5, pady=5
        )

        camera_button_frame = customtkinter.CTkFrame(camera_frame)
        camera_button_frame.pack(fill="x", pady=(0, 5))

        connect_btn = customtkinter.CTkButton(
            camera_button_frame,
            text="Connect camera",
            command=self.on_camera_connect_clicked
        )
        connect_btn.pack(side="left", padx=5, pady=5)

        disconnect_btn = customtkinter.CTkButton(
            camera_button_frame,
            text="Disconnect camera",
            command=self.on_camera_disconnect_clicked
        )
        disconnect_btn.pack(side="left", padx=5, pady=5)

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

    def log(self, msg):
        # Local time in Pacific, with timezone label
        now_local = datetime.datetime.now(LOCAL_TZ)
        timestamp = now_local.strftime("%Y-%m-%d %H:%M:%S %Z")
        line = f"[{timestamp}] {msg}\n"
        self.log_text.insert("end", line)
        self.log_text.see("end")
        self.log_text.update_idletasks()

    def run_with_progress(self, description, func, callback=None):
        """
        Run `func` in a background thread while showing an
        indeterminate progress bar. `func` returns a string or a
        list/tuple of strings to log. If callback is provided, it will
        be called with the result on the main thread after logging.
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
                    self.log(f"{description} FAILED: {error}")
                    messagebox.showerror(
                        "Error", f"{description} failed:\n{error}"
                    )
                else:
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
        Show local civil time (Pacific) and LST (sidereal).
        """
        try:
            now_astropy = Time.now()
            lst = now_astropy.sidereal_time('apparent', longitude=ATA_LOCATION.lon)
            local_now = datetime.datetime.now(LOCAL_TZ)
            text = (
                f"Local time: {local_now.strftime('%Y-%m-%d %H:%M:%S %Z')}   |   "
                f"LST (sidereal): {lst.to_string(unit=u.hour, sep=':', precision=0, pad=True)}"
            )
        except Exception as e:
            text = f"Time error: {e}"

        self.info_label.configure(text=text)
        # re-schedule
        self.root.after(10000, self._update_time_info)

    # ---------- Callbacks: Antenna ----------

    def on_reserve_clicked(self):
        self.run_with_progress("Reserving antennas", ata_reserve_antennas)

    def on_park_clicked(self):
        self.run_with_progress("Parking antennas", ata_park_antennas)

    def on_release_clicked(self):
        self.run_with_progress("Releasing antennas", ata_release_antennas)

    # ---------- Callbacks: Frequency / profiles ----------

    def on_profile_changed(self, *_):
        profile = self.selected_profile.get()
        freq_hz = OBSERVATION_PROFILES.get(profile)
        if freq_hz is not None:
            self.freq_entry.delete(0, tk.END)
            self.freq_entry.insert(0, f"{freq_hz / 1e6:.6f}")

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

        freq_hz = freq_mhz * 1e6

        def do_set_freq():
            return ata_set_frequency_hz(freq_hz)

        self.run_with_progress(
            f"Setting frequency to {freq_mhz:.6f} MHz",
            do_set_freq
        )

    # ---------- Callbacks: Coordinate modes ----------

    def on_coord_mode_changed(self):
        mode = self.coord_mode.get()
        if mode == "radec":
            self.coord1_label.configure(text="RA (h:m:s)")
            self.coord2_label.configure(text="Dec (d:m:s)")
        elif mode == "altaz":
            self.coord1_label.configure(text="Alt (deg)")
            self.coord2_label.configure(text="Az (deg)")
        elif mode == "galactic":
            self.coord1_label.configure(text="l (deg)")
            self.coord2_label.configure(text="b (deg)")
        elif mode == "name":
            # For name mode, RA/Dec fields are not required; they can be left as-is
            self.coord1_label.configure(text="RA (h:m:s)")
            self.coord2_label.configure(text="Dec (d:m:s)")

    def on_track_clicked(self):
        mode = self.coord_mode.get()
        coord1 = self.coord1_entry.get().strip()
        coord2 = self.coord2_entry.get().strip()
        name = self.name_entry.get().strip()

        # --- NEW: Source-name mode uses ATA catalog directly ---
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
            msg = ata_track_radec_degrees(ra_deg, dec_deg)
            if name:
                return f"{msg} (Source name field: {name})"
            return msg

        self.run_with_progress(
            f"Tracking RA={c.ra.to_string(unit=u.hour, sep=':')}, "
            f"Dec={c.dec.to_string(unit=u.deg, sep=':')}",
            do_track
        )

    def on_stop_tracking_clicked(self):
        self.run_with_progress("Stopping tracking", ata_stop_tracking)

    # ---------- Callbacks: Status / Camera ----------

    def on_refresh_status_clicked(self):
        def do_status():
            return ata_get_status_text()

        def update_status(text):
            if not isinstance(text, str):
                return

            # Show full status in text box
            self.status_text.delete("1.0", tk.END)
            self.status_text.insert("1.0", text)

            # Parse header + row for 1a
            header_tokens, row_tokens = self._parse_ascii_status_for_antenna(
                text, "1a"
            )
            if header_tokens and row_tokens:
                # Configure treeview
                self.status_tree.delete(*self.status_tree.get_children())
                self.status_tree["columns"] = header_tokens
                for col in header_tokens:
                    self.status_tree.heading(col, text=col)
                    self.status_tree.column(col, width=80, anchor="center")
                self.status_tree.insert("", "end", values=row_tokens)

        self.run_with_progress("Refreshing ATA status", do_status, callback=update_status)

    @staticmethod
    def _parse_ascii_status_for_antenna(ascii_status: str, antenna: str):
        """
        Simple parser for ac.get_ascii_status().

        Heuristic:
        - Find a header line whose first token starts with 'ant'.
        - Find a row whose first token is the requested antenna.
        - Split both lines on whitespace and zip them together.
        """
        lines = ascii_status.splitlines()
        header_tokens = None
        row_tokens = None

        for i, line in enumerate(lines):
            stripped = line.strip()
            if not stripped:
                continue
            tokens = stripped.split()
            # header line
            if header_tokens is None and tokens and tokens[0].lower().startswith("ant"):
                header_tokens = tokens
            # antenna row
            if tokens and tokens[0].lower() == antenna.lower():
                row_tokens = tokens
                # Grab header from previous line if needed
                if header_tokens is None and i > 0:
                    prev = lines[i - 1].strip().split()
                    if prev and prev[0].lower().startswith("ant"):
                        header_tokens = prev
                break

        if header_tokens and row_tokens:
            if len(row_tokens) > len(header_tokens):
                row_tokens = row_tokens[: len(header_tokens)]
            elif len(row_tokens) < len(header_tokens):
                row_tokens = row_tokens + [""] * (len(header_tokens) - len(row_tokens))
            return header_tokens, row_tokens

        return None, None

    def on_camera_connect_clicked(self):
        """
        Instead of loading the complex camera page (which may rely on JS),
        embed a simple <img> pointing directly at the MJPEG stream.
        """
        try:
            html = f"""
            <html>
              <body style="margin:0; padding:0; background-color:black;">
                <img src="{self.camera_stream_url}" style="width:100%; height:auto;" />
              </body>
            </html>
            """
            self.camera_html_frame.load_html(html)
            self.log("Camera connected via MJPEG stream.")
        except Exception as e:
            messagebox.showerror(
                "Camera error",
                f"Failed to load camera stream:\n{e}"
            )
            self.log(f"Camera connection failed: {e}")

    def on_camera_disconnect_clicked(self):
        try:
            self.camera_html_frame.load_html(
                "<html><body><h3 style='color:white; background:black;'>Camera disconnected.</h3></body></html>"
            )
            self.log("Camera disconnected.")
        except Exception as e:
            self.log(f"Camera disconnect failed: {e}")


# ======================================================
# MAIN
# ======================================================

if __name__ == "__main__":
    root = tk.Tk()
    app = ATAObservationGUI(root)
    root.mainloop()
