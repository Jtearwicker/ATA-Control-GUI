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
import numpy as np
import pyautogui
from tkinterweb import HtmlFrame
import datetime
from pytz import timezone
import pytz
from ATATools import ata_control as ac
from matplotlib import image 
from matplotlib import pyplot as plt
from PIL import Image, ImageTk
import time
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Create the main root
root = customtkinter.CTk()
root.title("Easy ATA GUI")
root.geometry("1200x800")
customtkinter.set_appearance_mode("Light")

# Create the main frame
main_frame = customtkinter.CTkFrame(master=root)
main_frame.pack(fill=BOTH, expand=True)

# Control frame on the left
control_frame = customtkinter.CTkFrame(master=main_frame, width=300)
control_frame.pack(side=LEFT, fill=Y, padx=10, pady=10)

# Image frame on the right
image_frame = customtkinter.CTkFrame(master=main_frame, width=900, height=400)
image_frame.pack(side=TOP, fill=BOTH, expand=True, padx=10, pady=10)

# Terminal frame at the bottom
terminal_frame = customtkinter.CTkFrame(master=main_frame)
terminal_frame.pack(side=BOTTOM, fill=BOTH, expand=True, padx=10, pady=10)

# Create a text widget to display the terminal output
terminal_text = customtkinter.CTkTextbox(master=terminal_frame, height=400, width=1200, font=("DejaVu Sans Mono", 12))
terminal_text.pack(fill=BOTH, expand=True)

# Load and resize the image
image_path = "MWimg.jpg"
img = Image.open(image_path)
img_tk = ImageTk.PhotoImage(img)

# Label to display the image
image_label = Label(image_frame, image=img_tk)
image_label.pack(fill=BOTH, expand=True)

def run_test_command(command):
    terminal_text.insert(END, "Testing USRPs...\n")
    stream = os.popen(command)
    out = stream.read()
    terminal_text.insert(END, out)

def run_reset_command(command):
    terminal_text.insert(END, "Resetting USRPs clocking...\n")
    stream = os.popen(command)
    out = stream.read()
    terminal_text.insert(END, out)

def run_server_command(command):
    return os.popen(command)

antennas = ['1a']
freq = "1420.406"

def activate_antenna_clicked():
    run_reset_command("/opt/ata-flowgraphs/usrp_reset_clocking.py")
    run_server_command("python /home/vgajjar/reu-2023/Hydrogen_line/server.py")
    ant_free = str(ac.list_antenna_group('none'))
    if ant_free.find('1a') == -1:
        terminal_text.insert(END, "WARNING: Antenna 1a has already been reserved.\n")
    else:
        ac.move_ant_group(antennas, 'none', 'atagr')
        terminal_text.insert(END, "Antenna 1a has been reserved.\n")
    ac.set_freq(freq, antennas, 'd')
    ac.autotune(antennas)
    terminal_text.insert(END, "Frequency set to 1420.406 MHz and autotuned.\n")
    terminal_text.insert(END, ac.get_ascii_status()[:348] + "\n")
    time.sleep(45)
    terminal_text.insert(END, "Calibration complete!\n")

def show_ant_status_clicked():
    terminal_text.insert(END, ac.get_ascii_status()[:348] + "\n")

def shut_down_antenna_clicked():
    run_server_command("pkill -f \"python /home/vgajjar/reu-2023/Hydrogen_line/server.py\"")
    terminal_text.insert(END, "Disconnected from server.\n")
    ac.park_antennas(antennas)
    terminal_text.insert(END, "Antenna 1a has been parked.\n")
    ant_free = str(ac.list_antenna_group('none'))
    if ant_free.find('1a') != -1:
        terminal_text.insert(END, "Antenna 1a has already been released\n")
    else:
        ac.move_ant_group(antennas, 'atagr', 'none')
        terminal_text.insert(END, "Antenna 1a has been released.\n")

def ga2equ(ga):
    g1 = radians(ga[0])
    g2 = radians(ga[1])
    c1 = radians(192.85948)
    c2 = radians(27.12825)
    c3 = radians(32.93192)
    a = degrees(atan2((cos(g2) * cos(g1 - c3)), (sin(g2) * cos(c2) - cos(g2) * sin(c2) * sin(g1 - c3))) + c1)
    b = degrees(asin(cos(g2) * cos(c2) * sin(g1 - c3) + sin(g2) * sin(c2)))
    if a >= 360:
        a -= 360
    return [a, b]

def radec2alt(ga):
    ga = [ga[0], ga[1]]
    ga = ga2equ(ga)
    ata = EarthLocation(lat=40.8175 * u.deg, lon=-121.469 * u.deg, height=1007 * u.m)
    delta = 0 * u.m
    altaz = AltAz(location=ata, obstime=Time(datetime.datetime.now(), scale='utc'))
    radec = SkyCoord(ra=ga[0] * u.deg, dec=ga[1] * u.deg)
    altitude = radec.transform_to(altaz).alt.deg
    return altitude

targets = [[45,0],[50,0],[55,0],[60,0],[65,0],[70,0],[75,0],[80,0],[85,0],[90,0],[95,0],[100,0],[105,0],[110,0],[115,0],[120,0],[125,0],[130,0],[135,0],[140,0],[145,0],[150,0],[155,0],[160,0],[165,0],[170,0],[175,0],[180,0],[185,0],[190,0],[195,0],[200,0],[205,0],[210,0],[215,0]]

def list_avail_targets_clicked():
    for i in range(0, 35):
        dd_radec = ga2equ(targets[i])
        c = SkyCoord(ra=dd_radec[0] * u.deg, dec=dd_radec[1] * u.deg)
        RA = c.ra.hms
        DEC = c.dec.dms
        elevation = radec2alt(ga2equ(targets[i]))

        if elevation > 20:
            terminal_text.insert(END, "Galactic longitude " + str(targets[i][0]) + " has an elevation of " + str(elevation)[0:4] + " degrees above the horizon.\n")
            avail_long = targets[i][0]

def track_source_clicked():
    gl = int(galactic_longitude_entry.get())
    dd_radec = ga2equ([gl, 0])
    c = SkyCoord(ra=dd_radec[0] * u.deg, dec=dd_radec[1] * u.deg)
    RA = c.ra.hms
    DEC = c.dec.dms
    ac.track_source(antennas, radec=[Angle(str(int(RA[0])) + "h" + str(int(RA[1])) + "m" + str(int(RA[2])) + "s").hour,
                                     Angle(str(int(DEC[0])) + "d" + str(int(abs(DEC[1]))) + "m" + str(int(abs(DEC[2]))) + "s").deg])
    terminal_text.insert(END, "Arrived at galactic coordinate (" + str(gl) + ",0). RA " + str(int(RA[0])) + "h" + str(int(RA[1])) + "m" + str(int(RA[2])) + "s Dec " + str(int(DEC[0])) + "d" + str(int(abs(DEC[1]))) + "m" + str(int(abs(DEC[2]))) + "s\n")

activate_antenna_button = customtkinter.CTkButton(master=control_frame, text="Activate Antenna", command=activate_antenna_clicked)
activate_antenna_button.pack(padx=5, pady=5)

avail_targets_button = customtkinter.CTkButton(master=control_frame, text="Show Available Targets", command=list_avail_targets_clicked)
avail_targets_button.pack(padx=5, pady=5)

galactic_longitude_entry = customtkinter.CTkEntry(master=control_frame, placeholder_text="Galactic Longitude")
galactic_longitude_entry.pack(padx=5, pady=5)

track_source_button = customtkinter.CTkButton(master=control_frame, text="Track Source", command=track_source_clicked)
track_source_button.pack(padx=5, pady=5)

show_ant_status_button = customtkinter.CTkButton(master=control_frame, text="Show Antenna Status", command=show_ant_status_clicked)
show_ant_status_button.pack(padx=5, pady=5)

shut_down_antenna_button = customtkinter.CTkButton(master=control_frame, text="Shut Down Antenna", command=shut_down_antenna_clicked)
shut_down_antenna_button.pack(padx=5, pady=5)

root.mainloop()
