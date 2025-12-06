import datetime
import threading
import tkinter as tk
from tkinter import ttk, messagebox

# Third-party libs used at the ATA
try:
    import customtkinter as ctk
except ImportError:  # Fallback so the file can at least be imported elsewhere
    ctk = None

try:
    from ATATools import ata_control as ac
except ImportError:
    ac = None  # At the ATA this must import correctly

import pytz
import astropy.units as u
from astropy.coordinates import SkyCoord, Galactic
from astropy.time import Time

# Optional camera dependencies (used for MJPEG streaming)
try:
    import cv2
    from PIL import Image, ImageTk
except ImportError:
    cv2 = None
    Image = None
    ImageTk = None

# ---- Constants ----

# Hat Creek / ATA location (from HCRO / ATA coordinates)
ATA_LOCATION = (40.8178, -121.473, 986.0)  # lat [deg], lon [deg], height [m]

# Local civil timezone for log timestamps and display
LOCAL_TZ = pytz.timezone("America/Los_Angeles")

# Direct MJPEG stream for the site camera, derived from the HTML URL you sent
CAMERA_MJPEG_URL = "http://10.3.0.30/mjpg/video.mjpg"


class ATAObserverGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("ATA Observation Console")
        self.antennas = ["1a"]

        if ctk is not None:
            ctk.set_appearance_mode("dark")
            ctk.set_default_color_theme("blue")

        # --- Top info bar (time, UTC, LST) ---
        self.info_frame = ttk.Frame(root)
        self.info_frame.pack(side="top", fill="x", padx=5, pady=5)

        self.time_label = ttk.Label(self.info_frame, text="", font=("DejaVu Sans Mono", 10))
        self.time_label.pack(side="left", padx=5)

        # --- Main layout: left controls, right status/camera ---
        self.main_frame = ttk.Frame(root)
        self.main_frame.pack(side="top", fill="both", expand=True)

        self.left_frame = ttk.Frame(self.main_frame)
        self.left_frame.pack(side="left", fill="both", expand=True, padx=5, pady=5)

        self.right_frame = ttk.Frame(self.main_frame)
        self.right_frame.pack(side="right", fill="both", expand=True, padx=5, pady=5)

        # ---- Left: antenna, profile, pointing ----
        self._build_antenna_frame(self.left_frame)
        self._build_observing_frame(self.left_frame)
        self._build_pointing_frame(self.left_frame)

        # ---- Right: status + camera ----
        self._build_status_frame(self.right_frame)
        self._build_camera_frame(self.right_frame)

        # ---- Bottom: progress + log ----
        self.bottom_frame = ttk.Frame(root)
        self.bottom_frame.pack(side="bottom", fill="x", padx=5, pady=5)

        self.progress = ttk.Progressbar(self.bottom_frame, mode="indeterminate")
        self.progress.pack(side="top", fill="x", expand=True, padx=5, pady=2)

        self.log_text = tk.Text(
            self.bottom_frame,
            height=8,
            width=120,
            font=("DejaVu Sans Mono", 9),
        )
        self.log_text.pack(side="bottom", fill="x", expand=False, padx=5, pady=2)

        # Task state
        self.current_task = None

        # Camera state
        self.camera_capture = None
        self.camera_running = False

        # Start clock updates
        self.update_clock()

    # ------------------------------------------------------------------
    # UI building blocks
    # ------------------------------------------------------------------

    def _build_antenna_frame(self, parent):
        lf = ttk.LabelFrame(parent, text="Antenna & Backend")
        lf.pack(fill="x", padx=5, pady=5)

        ant_label = ttk.Label(lf, text=f"Antennas: {', '.join(self.antennas)}")
        ant_label.grid(row=0, column=0, columnspan=4, sticky="w", padx=5, pady=2)

        self.reserve_button = ttk.Button(lf, text="Reserve", command=self.on_reserve)
        self.reserve_button.grid(row=1, column=0, padx=5, pady=2, sticky="ew")

        self.park_button = ttk.Button(lf, text="Park", command=self.on_park)
        self.park_button.grid(row=1, column=1, padx=5, pady=2, sticky="ew")

        self.release_button = ttk.Button(lf, text="Release", command=self.on_release)
        self.release_button.grid(row=1, column=2, padx=5, pady=2, sticky="ew")

        self.autotune_button = ttk.Button(lf, text="Autotune", command=self.on_autotune)
        self.autotune_button.grid(row=1, column=3, padx=5, pady=2, sticky="ew")

        # Frequency & attenuation controls
        ttk.Label(lf, text="Freq (MHz):").grid(row=2, column=0, sticky="e", padx=5, pady=2)
        self.freq_var = tk.StringVar(value="1420.405")
        ttk.Entry(lf, textvariable=self.freq_var, width=10).grid(row=2, column=1, sticky="w", padx=5, pady=2)

        ttk.Label(lf, text="LO:").grid(row=2, column=2, sticky="e", padx=5, pady=2)
        self.lo_var = tk.StringVar(value="d")
        ttk.Entry(lf, textvariable=self.lo_var, width=5).grid(row=2, column=3, sticky="w", padx=5, pady=2)

        ttk.Label(lf, text="Atten (dB):").grid(row=3, column=0, sticky="e", padx=5, pady=2)
        self.atten_var = tk.StringVar(value="20")
        ttk.Entry(lf, textvariable=self.atten_var, width=5).grid(row=3, column=1, sticky="w", padx=5, pady=2)

        self.set_freq_button = ttk.Button(lf, text="Set Freq", command=self.on_set_freq)
        self.set_freq_button.grid(row=3, column=2, columnspan=2, padx=5, pady=4, sticky="ew")

    def _build_observing_frame(self, parent):
        lf = ttk.LabelFrame(parent, text="Observation Profile")
        lf.pack(fill="x", padx=5, pady=5)

        ttk.Label(lf, text="Profile:").grid(row=0, column=0, sticky="e", padx=5, pady=2)
        self.profile_var = tk.StringVar(value="Custom")

        profiles = ["Custom", "21cm Milky Way (1420.405 MHz)", "Mars (8431 MHz)"]
        self.profile_menu = ttk.OptionMenu(
            lf, self.profile_var, profiles[0], *profiles, command=self.on_profile_changed
        )
        self.profile_menu.grid(row=0, column=1, columnspan=3, sticky="ew", padx=5, pady=2)

    def _build_pointing_frame(self, parent):
        lf = ttk.LabelFrame(parent, text="Pointing & Tracking")
        lf.pack(fill="x", padx=5, pady=5)

        # Pointing mode radio buttons
        ttk.Label(lf, text="Mode:").grid(row=0, column=0, sticky="w", padx=5, pady=2)
        self.pointing_mode = tk.StringVar(value="radec")

        modes = [
            ("RA/Dec", "radec"),
            ("Alt/Az", "altaz"),
            ("Galactic", "gal"),
            ("Source name", "name"),
            ("Park", "park"),
        ]
        col = 1
        for text, val in modes:
            ttk.Radiobutton(lf, text=text, variable=self.pointing_mode, value=val).grid(
                row=0, column=col, sticky="w", padx=2, pady=2
            )
            col += 1

        # RA/Dec input
        ttk.Label(lf, text="RA (hh:mm:ss):").grid(row=1, column=0, sticky="e", padx=5, pady=2)
        self.ra_var = tk.StringVar()
        ttk.Entry(lf, textvariable=self.ra_var, width=15).grid(row=1, column=1, sticky="w", padx=5, pady=2)

        ttk.Label(lf, text="Dec (dd:mm:ss):").grid(row=1, column=2, sticky="e", padx=5, pady=2)
        self.dec_var = tk.StringVar()
        ttk.Entry(lf, textvariable=self.dec_var, width=15).grid(row=1, column=3, sticky="w", padx=5, pady=2)

        # Alt/Az input
        ttk.Label(lf, text="Alt (deg):").grid(row=2, column=0, sticky="e", padx=5, pady=2)
        self.alt_var = tk.StringVar()
        ttk.Entry(lf, textvariable=self.alt_var, width=10).grid(row=2, column=1, sticky="w", padx=5, pady=2)

        ttk.Label(lf, text="Az (deg):").grid(row=2, column=2, sticky="e", padx=5, pady=2)
        self.az_var = tk.StringVar()
        ttk.Entry(lf, textvariable=self.az_var, width=10).grid(row=2, column=3, sticky="w", padx=5, pady=2)

        # Galactic input
        ttk.Label(lf, text="l (deg):").grid(row=3, column=0, sticky="e", padx=5, pady=2)
        self.glon_var = tk.StringVar()
        ttk.Entry(lf, textvariable=self.glon_var, width=10).grid(row=3, column=1, sticky="w", padx=5, pady=2)

        ttk.Label(lf, text="b (deg):").grid(row=3, column=2, sticky="e", padx=5, pady=2)
        self.glat_var = tk.StringVar()
        ttk.Entry(lf, textvariable=self.glat_var, width=10).grid(row=3, column=3, sticky="w", padx=5, pady=2)

        # Source name (catalog)
        ttk.Label(lf, text="Source name (catalog):").grid(row=4, column=0, sticky="e", padx=5, pady=2)
        self.source_name_var = tk.StringVar()
        ttk.Entry(lf, textvariable=self.source_name_var, width=20).grid(
            row=4, column=1, columnspan=3, sticky="w", padx=5, pady=2
        )

        self.track_button = ttk.Button(lf, text="Track", command=self.on_track)
        self.track_button.grid(row=5, column=0, columnspan=2, padx=5, pady=4, sticky="ew")

    def _build_status_frame(self, parent):
        lf = ttk.LabelFrame(parent, text="Array Status (antenna 1a)")
        lf.pack(fill="both", expand=True, padx=5, pady=5)

        self.status_button = ttk.Button(lf, text="Refresh status", command=self.on_status)
        self.status_button.pack(side="top", anchor="w", padx=5, pady=2)

        # Parsed table for antenna 1a
        self.status_tree = ttk.Treeview(lf, columns=(), show="headings", height=5)
        self.status_tree.pack(fill="x", padx=5, pady=2)

        # Raw status (filtered) below
        self.status_raw = tk.Text(lf, height=10, font=("DejaVu Sans Mono", 8))
        self.status_raw.pack(fill="both", expand=True, padx=5, pady=2)

    def _build_camera_frame(self, parent):
        lf = ttk.LabelFrame(parent, text="Site Camera")
        lf.pack(fill="both", expand=True, padx=5, pady=5)

        self.camera_label = ttk.Label(lf, text="Camera disconnected", anchor="center")
        self.camera_label.pack(fill="both", expand=True, padx=5, pady=5)

        btn_frame = ttk.Frame(lf)
        btn_frame.pack(fill="x", padx=5, pady=2)

        self.camera_connect_btn = ttk.Button(btn_frame, text="Connect", command=self.on_camera_connect)
        self.camera_connect_btn.pack(side="left", padx=5)

        self.camera_disconnect_btn = ttk.Button(btn_frame, text="Disconnect", command=self.on_camera_disconnect)
        self.camera_disconnect_btn.pack(side="left", padx=5)

        if cv2 is None or ImageTk is None:
            self.camera_label.configure(
                text="Camera preview requires OpenCV + Pillow.\n"
                     "Install them or use the browser view."
            )
            self.camera_connect_btn.configure(state="disabled")

    # ------------------------------------------------------------------
    # Time / logging / async tasks
    # ------------------------------------------------------------------

    def update_clock(self):
        """Update local time, UTC, and LST (sidereal) in the top bar."""
        try:
            utc_now = datetime.datetime.now(datetime.timezone.utc)
            local_now = utc_now.astimezone(LOCAL_TZ)

            # Sidereal time at ATA
            lat_deg, lon_deg, height_m = ATA_LOCATION
            try:
                from astropy.coordinates import EarthLocation
                location = EarthLocation(
                    lat=lat_deg * u.deg,
                    lon=lon_deg * u.deg,
                    height=height_m * u.m,
                )
                t = Time(utc_now, location=location)
                lst = t.sidereal_time("apparent")
                lst_str = lst.to_string(unit=u.hour, sep=":", precision=0, pad=True)
            except Exception:
                lst_str = "N/A"

            self.time_label.configure(
                text=(
                    f"Local time: {local_now:%Y-%m-%d %H:%M:%S %Z}   "
                    f"UTC: {utc_now:%H:%M:%S}   "
                    f"LST (sidereal): {lst_str}"
                )
            )
        except Exception as e:
            self.time_label.configure(text=f"Time error: {e}")

        self.root.after(1000, self.update_clock)

    def log(self, msg: str):
        now_local = datetime.datetime.now(LOCAL_TZ)
        ts = now_local.strftime("%Y-%m-%d %H:%M:%S %Z")
        line = f"[{ts}] {msg}\n"
        self.log_text.insert("end", line)
        self.log_text.see("end")

    def run_long_task(self, fn, description: str, callback=None):
        """Run ATA operations in a background thread with a spinner."""
        if ac is None:
            messagebox.showerror("ATATools missing", "ATATools.ata_control could not be imported.")
            return

        if self.current_task is not None:
            messagebox.showwarning("Busy", "Another operation is still running.")
            return

        self.current_task = description
        self.progress.start()
        self.log(f"{description}...")

        def worker():
            error = None
            result = None
            try:
                result = fn()
            except Exception as e:
                error = e

            def done():
                self.progress.stop()
                self.current_task = None
                if error:
                    self.log(f"{description} FAILED: {error}")
                    messagebox.showerror("Error", f"{description} failed:\n{error}")
                else:
                    self.log(f"{description} completed.")
                if callback is not None:
                    callback(result, error)

            self.root.after(0, done)

        threading.Thread(target=worker, daemon=True).start()

    # ------------------------------------------------------------------
    # ATA control handlers
    # ------------------------------------------------------------------

    def on_reserve(self):
        def task():
            return ac.move_ant_group(self.antennas, "none", "atagr")

        self.run_long_task(task, "Reserving antennas")

    def on_release(self):
        def task():
            return ac.move_ant_group(self.antennas, "atagr", "none")

        self.run_long_task(task, "Releasing antennas")

    def on_park(self):
        def task():
            return ac.park_antennas(self.antennas)

        self.run_long_task(task, "Parking antennas")

    def on_autotune(self):
        def task():
            ac.autotune(self.antennas)

        self.run_long_task(task, "Running autotune")

    def on_set_freq(self):
        try:
            freq = float(self.freq_var.get())
        except ValueError:
            messagebox.showerror("Invalid frequency", "Please enter a valid frequency in MHz.")
            return

        lo = self.lo_var.get().strip()
        if lo not in ("a", "b", "c", "d"):
            messagebox.showerror("Invalid LO", "LO must be one of a, b, c, d.")
            return

        try:
            att = float(self.atten_var.get())
        except ValueError:
            messagebox.showerror("Invalid attenuation", "Please enter a valid attenuation in dB.")
            return

        def task():
            ac.set_freq(freq, self.antennas, lo)
            # Switch matrix + attenuators
            ac.rf_switch_thread(self.antennas)
            ac.set_atten_thread(
                [[f"{ant}x", f"{ant}y"] for ant in self.antennas],
                [[att, att] for _ in self.antennas],
            )

        self.run_long_task(task, f"Setting frequency to {freq} MHz, LO {lo}, attenuation {att} dB")

    def on_profile_changed(self, *_):
        profile = self.profile_var.get()
        if profile.startswith("21cm"):
            self.freq_var.set("1420.405")
        elif profile.startswith("Mars"):
            self.freq_var.set("8431")
        # User still clicks "Set Freq" to actually apply.

    def on_track(self):
        mode = self.pointing_mode.get()

        if mode == "park":
            self.on_park()
            return

        # Name mode: use ATA catalog directly, no RA/Dec required
        if mode == "name":
            src = self.source_name_var.get().strip()
            if not src:
                messagebox.showerror("Missing source name", "Please enter a catalog source name.")
                return

            def task():
                ac.track_source(self.antennas, source=src)

            self.run_long_task(task, f"Tracking catalog source '{src}'")
            return

        # Coordinate-based modes need astropy's EarthLocation
        try:
            from astropy.coordinates import EarthLocation
        except ImportError:
            messagebox.showerror("Astropy missing", "Astropy is required for coordinate transforms.")
            return

        lat_deg, lon_deg, height_m = ATA_LOCATION
        location = EarthLocation(
            lat=lat_deg * u.deg,
            lon=lon_deg * u.deg,
            height=height_m * u.m,
        )
        utc_now = datetime.datetime.now(datetime.timezone.utc)
        obstime = Time(utc_now, location=location)

        if mode == "radec":
            ra_str = self.ra_var.get().strip()
            dec_str = self.dec_var.get().strip()
            if not ra_str or not dec_str:
                messagebox.showerror("Missing coordinates", "Enter both RA and Dec.")
                return

            try:
                sky = SkyCoord(ra=ra_str, dec=dec_str, unit=("hourangle", "deg"))
            except Exception as e:
                messagebox.showerror("Invalid RA/Dec", f"Could not parse RA/Dec: {e}")
                return

            radec = [sky.ra.hour, sky.dec.deg]

            def task():
                ac.track_source(self.antennas, radec=radec)

            desc = (
                "Tracking RA="
                + sky.ra.to_string(unit=u.hour, sep=":")
                + " Dec="
                + sky.dec.to_string(unit=u.deg, sep=":")
            )
            self.run_long_task(task, desc)

        elif mode == "altaz":
            alt_str = self.alt_var.get().strip()
            az_str = self.az_var.get().strip()
            if not alt_str or not az_str:
                messagebox.showerror("Missing coordinates", "Enter both Alt and Az.")
                return

            try:
                alt = float(alt_str)
                az = float(az_str)
            except ValueError:
                messagebox.showerror("Invalid Alt/Az", "Alt and Az must be numbers in degrees.")
                return

            azel = [az, alt]  # ata_control expects [az, el]

            def task():
                ac.track_source(self.antennas, azel=azel)

            self.run_long_task(task, f"Tracking Az={az:.2f}°, El={alt:.2f}°")

        elif mode == "gal":
            l_str = self.glon_var.get().strip()
            b_str = self.glat_var.get().strip()
            if not l_str or not b_str:
                messagebox.showerror("Missing coordinates", "Enter both Galactic l and b.")
                return

            try:
                glon = float(l_str)
                glat = float(b_str)
            except ValueError:
                messagebox.showerror("Invalid Galactic", "l and b must be numbers in degrees.")
                return

            try:
                gal = SkyCoord(l=glon * u.deg, b=glat * u.deg, frame="galactic")
                sky = gal.transform_to("icrs")
            except Exception as e:
                messagebox.showerror("Transform error", f"Could not convert Galactic → RA/Dec: {e}")
                return

            radec = [sky.ra.hour, sky.dec.deg]

            def task():
                ac.track_source(self.antennas, radec=radec)

            desc = f"Tracking Galactic l={glon:.2f}°, b={glat:.2f}° (RA/Dec converted)"
            self.run_long_task(task, desc)

    # ------------------------------------------------------------------
    # Status handling
    # ------------------------------------------------------------------

    def on_status(self):
        def task():
            return ac.get_ascii_status()

        def callback(ascii_status, error):
            if error or not ascii_status:
                return

            lines = ascii_status.splitlines()

            # Raw view: header + 1a row, if we can find them
            filtered = []
            for line in lines:
                low = line.lower()
                if "ant" in low and ("az" in low or "el" in low):
                    filtered.append(line)
                elif " 1a " in f" {low} " or low.startswith("1a "):
                    filtered.append(line)
            if not filtered:
                filtered = lines  # fallback

            self.status_raw.configure(state="normal")
            self.status_raw.delete("1.0", "end")
            self.status_raw.insert("end", "\n".join(filtered))
            self.status_raw.configure(state="disabled")

            # Parsed table for antenna 1a
            header_tokens, row_tokens = self._parse_ascii_status_for_antenna(ascii_status, "1a")
            if header_tokens and row_tokens:
                self.status_tree.delete(*self.status_tree.get_children())
                self.status_tree["columns"] = header_tokens
                for col in header_tokens:
                    self.status_tree.heading(col, text=col)
                    self.status_tree.column(col, width=80, anchor="center")
                self.status_tree.insert("", "end", values=row_tokens)

        self.run_long_task(task, "Refreshing array status", callback=callback)

    @staticmethod
    def _parse_ascii_status_for_antenna(ascii_status: str, antenna: str):
        """
        Very generic parser for ac.get_ascii_status().

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

    # ------------------------------------------------------------------
    # Camera handling
    # ------------------------------------------------------------------

    def on_camera_connect(self):
        if cv2 is None or ImageTk is None:
            messagebox.showerror(
                "Camera dependencies missing",
                "OpenCV (cv2) and Pillow are required for the embedded camera view.",
            )
            return

        if self.camera_running:
            return

        self.log("Connecting to site camera...")
        cap = cv2.VideoCapture(CAMERA_MJPEG_URL)
        if not cap.isOpened():
            self.log("Failed to open camera stream.")
            messagebox.showerror("Camera error", "Could not open MJPEG stream.")
            return

        self.camera_capture = cap
        self.camera_running = True
        self.camera_label.configure(text="")

        self._update_camera_frame()

    def _update_camera_frame(self):
        if not self.camera_running or self.camera_capture is None:
            return

        ret, frame = self.camera_capture.read()
        if not ret:
            # Try again shortly without spamming the log
            self.root.after(500, self._update_camera_frame)
            return

        # Convert BGR -> RGB and display in the label
        frame_rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
        img = Image.fromarray(frame_rgb)

        # Resize to the current label size (or a default)
        w = self.camera_label.winfo_width() or 640
        h = self.camera_label.winfo_height() or 360
        img = img.resize((w, h))

        photo = ImageTk.PhotoImage(img)
        self.camera_label.image = photo
        self.camera_label.configure(image=photo)

        # Schedule next frame
        self.root.after(100, self._update_camera_frame)

    def on_camera_disconnect(self):
        if not self.camera_running:
            return
        self.camera_running = False
        if self.camera_capture is not None:
            self.camera_capture.release()
            self.camera_capture = None
        self.camera_label.configure(text="Camera disconnected", image="")
        self.camera_label.image = None
        self.log("Camera disconnected.")


def main():
    root = tk.Tk()
    app = ATAObserverGUI(root)
    root.geometry("1300x900")
    root.mainloop()


if __name__ == "__main__":
    main()
