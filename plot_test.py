import tkinter as tk
from tkinter import *
from tkinter import ttk
import customtkinter
import os
import subprocess
import numpy as np
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, Angle
from astropy.time import Time
from math import *
import pyautogui
from tkinterweb import HtmlFrame
import datetime
from pytz import timezone
import pytz
from ATATools import ata_control as ac
from matplotlib import image 
from matplotlib import pyplot as plt
from PIL import Image, ImageTk, ImageDraw
import time
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Create the main root
root = customtkinter.CTk()
root.title("Easy ATA GUI")
root.geometry("1200x800")
customtkinter.set_appearance_mode("Light")

# Control tabs
control_frame = customtkinter.CTkFrame(master=root, height=350, width=600)
control_frame.pack(side=TOP, pady=10)

# Galaxy pointing image frame
image_frame = customtkinter.CTkFrame(master=root, height=350, width=600)
image_frame.pack(side=TOP, pady=10)

# Create a frame for the terminal output
terminal_frame = customtkinter.CTkFrame(master=root)
terminal_frame.pack(side=TOP, pady=10)

# Create a text widget to display the terminal output
terminal_text = customtkinter.CTkTextbox(master=terminal_frame, height=400, width=1200, font=("DejaVu Sans Mono", 12))
terminal_text.pack()

# Load the image
milky_way_image = Image.open("MWimg.jpg")
image_draw = ImageDraw.Draw(milky_way_image)
photo = ImageTk.PhotoImage(milky_way_image)
image_label = tk.Label(image_frame, image=photo)
image_label.pack()

def update_image():
    global photo
    photo = ImageTk.PhotoImage(milky_way_image)
    image_label.config(image=photo)
    image_label.image = photo

def run_test_command(command):
    terminal_text.insert("0.0", "Testing USRPs...\n")
    stream = os.popen(command)
    out = stream.read()
    terminal_text.insert("0.0", out)

def run_reset_command(command):
    terminal_text.insert("0.0", "Resetting USRPs clocking...\n")
    stream = os.popen(command)
    out = stream.read()
    terminal_text.insert("0.0", out)

def run_server_command(command):
    return os.popen(command)

antennas = ['1a']
freq = "1420.406"
def activate_antenna_clicked():
    run_reset_command("/opt/ata-flowgraphs/usrp_reset_clocking.py")
    run_server_command("python /home/vgajjar/reu-2023/Hydrogen_line/server.py")
    ant_free = str(ac.list_antenna_group('none'))
    if ant_free.find('1a') == -1:
        terminal_text.insert("0.0", "WARNING: Antenna 1a has already been reserved.\n")
    else:
        ac.move_ant_group(antennas, 'none', 'atagr')
        terminal_text.insert("0.0", "Antenna 1a has been reserved.\n")
    ac.set_freq(freq, antennas, 'd')
    ac.autotune(antennas)
    terminal_text.insert("0.0", "Frequency set to 1420.406 MHz and autotuned.\n")
    terminal_text.insert("0.0", ac.get_ascii_status()[:348] + "\n")
    time.sleep(45)
    terminal_text.insert("0.0", "Calibration complete!\n")

def show_ant_status_clicked():
    terminal_text.insert("0.0", ac.get_ascii_status()[:348] + "\n")

def shut_down_antenna_clicked():
    run_server_command("pkill -f \"python /home/vgajjar/reu-2023/Hydrogen_line/server.py\"")
    terminal_text.insert("0.0", "Disconnected from server.\n")
    ac.park_antennas(antennas)
    terminal_text.insert("0.0", "Antenna 1a has been parked.\n")
    ant_free = str(ac.list_antenna_group('none'))
    if ant_free.find('1a') != -1:
        terminal_text.insert("0.0", "Antenna 1a has already been released\n")
    else:
        ac.move_ant_group(antennas, 'atagr', 'none')
        terminal_text.insert("0.0", "Antenna 1a has been released.\n")

# Target availability calculator
targets = [[0,0],[10,0],[20,0],[30,0],[40,0],[50,0],[60,0],[70,0],[80,0],
[90,0],[100,0],[110,0],[120,0],[130,0],[140,0],[150,0],[160,0],[170,0],
[180,0],[190,0],[200,0],[210,0],[220,0],[230,0],[240,0],[250,0],[260,0],
[270,0],[280,0],[290,0],[300,0],[310,0],[320,0],[340,0],[350,0]]

ATA_location = EarthLocation(lat=40.817 * u.deg, lon=-121.47 * u.deg, height=3235 * u.m)
utc = pytz.timezone('UTC')
now = utc.localize(datetime.datetime.utcnow())
la = pytz.timezone('America/Los_Angeles')
obs_time = now.astimezone(la)
alt_az = AltAz(location=ATA_location, obstime=obs_time)

def radec2alt(RADEC):
    coord = SkyCoord(RADEC[0] * u.deg, RADEC[1] * u.deg)
    aa = coord.transform_to(alt_az)
    altitude = aa.alt.value
    return altitude

def ga2equ(ga):
    l = radians(ga[0])
    b = radians(ga[1])
    pole_ra = radians(192.859508)
    pole_dec = radians(27.128336)
    posangle = radians(122.932-90.0)   
    ra = atan2( (cos(b)*cos(l-posangle)), (sin(b)*cos(pole_dec) - cos(b)*sin(pole_dec)*sin(l-posangle)) ) + pole_ra
    dec = asin( cos(b)*cos(pole_dec)*sin(l-posangle) + sin(b)*sin(pole_dec) )
    return np.array([degrees(ra), degrees(dec)])

pxlib = np.array([[2800, 500],[2233, 650],[1666, 800],[1100, 950],[833, 1443],[566, 1936],[300, 2430],[350, 2903],[400, 3376],[450, 3850],
                                    [733, 4140],[1016, 4430],[1300, 4720],[1550, 4846],[1800, 4973],[2050, 5100],[2300, 5140],[2550, 5180],[2800, 5220],[3033, 5180],
                                    [3266, 5140],[3500, 5100],[3766, 4973],[4033, 4846],[4300, 4720],[4583, 4430],[4866, 4140],[5150, 3850],[5200, 3376],[5250, 2903],
                                    [5300, 2430],[5026, 1936],[4753, 1443],[4480, 950],[3920, 800],[3360, 650]])

avail_targets = []

def list_avail_targets_clicked():
    global avail_targets
    avail_targets = []
    for i in range(0, 35):
        dd_radec = ga2equ(targets[i])
        c = SkyCoord(ra=dd_radec[0]*u.deg, dec=dd_radec[1]*u.deg)
        elevation = radec2alt(ga2equ(targets[i]))

        if elevation > 20:
            terminal_text.insert("0.0", f"Galactic longitude {targets[i][0]} has an elevation of {elevation:.2f} degrees above the horizon.\n")
            avail_targets.append(targets[i][0])

    # Example visualization of available targets
    min_pos = int(min(avail_targets) / 10)
    max_pos = int(max(avail_targets) / 10)
    px_min = pxlib[min_pos]
    px_max = pxlib[max_pos]
    vis_min_x = [2800, px_min[0]]
    vis_min_y = [3850, px_min[1]]
    vis_max_x = [2800, px_max[0]]
    vis_max_y = [3850, px_max[1]]
    image_draw.line((vis_min_x[0], vis_min_y[0], vis_min_x[1], vis_min_y[1]), fill=(255, 0, 0), width=3)
    image_draw.line((vis_max_x[0], vis_max_y[0], vis_max_x[1], vis_max_y[1]), fill=(0, 255, 0), width=3)
    update_image()

def track_source_clicked():
    vis_pos = int(np.array(avail_targets).min() / 10)
    px_vis = pxlib[vis_pos]
    vis_x = [2800, px_vis[0]]
    vis_y = [3850, px_vis[1]]
    image_draw.line((vis_x[0], vis_y[0], vis_x[1], vis_y[1]), fill=(0, 0, 255), width=3)
    update_image()

# Control frame widgets
activate_button = customtkinter.CTkButton(master=control_frame, text="Activate Antenna 1a", command=activate_antenna_clicked)
activate_button.grid(row=0, column=0, padx=5, pady=5)

status_button = customtkinter.CTkButton(master=control_frame, text="Show Antenna Status", command=show_ant_status_clicked)
status_button.grid(row=0, column=1, padx=5, pady=5)

shut_down_button = customtkinter.CTkButton(master=control_frame, text="Shut Down Antenna 1a", command=shut_down_antenna_clicked)
shut_down_button.grid(row=0, column=2, padx=5, pady=5)

list_avail_targets_button = customtkinter.CTkButton(master=control_frame, text="List Available Targets", command=list_avail_targets_clicked)
list_avail_targets_button.grid(row=1, column=0, padx=5, pady=5)

track_source_button = customtkinter.CTkButton(master=control_frame, text="Track Source", command=track_source_clicked)
track_source_button.grid(row=1, column=1, padx=5, pady=5)

# Entry for Galactic coordinates
coordinate_entry = customtkinter.CTkEntry(master=control_frame, width=200, placeholder_text="Enter Galactic Coordinate")
coordinate_entry.grid(row=1, column=2, padx=5, pady=5)

# Run the Tkinter main loop
root.mainloop()
