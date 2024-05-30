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
root = tk.Tk()
root.title("Easy ATA GUI")
root.geometry("1200x800")

# Create the main frame
main_frame = Frame(root)
main_frame.pack(fill=BOTH, expand=True)

# Control frame on the left
control_frame = Frame(main_frame, width=300)
control_frame.pack(side=LEFT, fill=Y, padx=10, pady=10)

# Image frame on the right
image_frame = Frame(main_frame, width=900, height=400)
image_frame.pack(side=TOP, fill=BOTH, expand=True, padx=10, pady=10)

# Terminal frame at the bottom
terminal_frame = Frame(main_frame)
terminal_frame.pack(side=BOTTOM, fill=BOTH, expand=True, padx=10, pady=10)

# Create a text widget to display the terminal output
terminal_text = Text(terminal_frame, height=10, font=("DejaVu Sans Mono", 12))
terminal_text.pack(fill=BOTH, expand=True)

# Load and resize the image
image_path = "MWimg.jpg"
img = Image.open(image_path)
img = img.resize((400, 400))
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
    #Input: [l,b] in decimal degrees
    #Returns: [ra,dec] in decimal degrees
    l = radians(ga[0])
    b = radians(ga[1])
    # North galactic pole (J2000)
    pole_ra = radians(192.859508)
    pole_dec = radians(27.128336)
    posangle = radians(122.932-90.0)   
    ra = atan2( (cos(b)*cos(l-posangle)), (sin(b)*cos(pole_dec) - cos(b)*sin(pole_dec)*sin(l-posangle)) ) + pole_ra
    dec = asin( cos(b)*cos(pole_dec)*sin(l-posangle) + sin(b)*sin(pole_dec) )
    return np.array([degrees(ra), degrees(dec)])

def radec2alt(RADEC):
    coord = SkyCoord(RADEC[0] * u.deg, RADEC[1] * u.deg)
    aa = coord.transform_to(alt_az)
    altitude = aa.alt.value
    return altitude

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

pxlib = np.array([[2800, 500],[2233, 650],[1666, 800],[1100, 950],[833, 1443],[566, 1936],[300, 2430],[350, 2903],[400, 3376],[450, 3850],
                                    [733, 4140],[1016, 4430],[1300, 4720],[1550, 4846],[1800, 4973],[2050, 5100],[2300, 5140],[2550, 5180],[2800, 5220],[3033, 5180],
                                    [3266, 5140],[3500, 5100],[3766, 4973],[4033, 4846],[4300, 4720],[4583, 4430],[4866, 4140],[5150, 3850],[5200, 3376],[5250, 2903],
                                    [5300, 2430],[5026, 1936],[4753, 1443],[4480, 950],[3920, 800],[3360, 650]])

avail_long = []
def list_avail_targets_clicked():
    global avail_long, img_tk, img
    avail_long = []  # Reset the list each time the button is clicked
    for i in range(0, 35):
        dd_radec = ga2equ(targets[i])
        c = SkyCoord(ra=dd_radec[0] * u.deg, dec=dd_radec[1] * u.deg)
        elevation = radec2alt(ga2equ(targets[i]))
        
        if elevation > 20:
            terminal_text.insert(END, "Galactic longitude " + str(targets[i][0]) + " has an elevation of " + str(elevation)[:4] + " degrees above the horizon.\n")
            avail_long.append(targets[i][0])
    
    if avail_long:
        min_long = avail_long[0]
        max_long = avail_long[-1]
        min_pix = pxlib[min_long // 10] // 14
        max_pix = pxlib[max_long // 10] // 14
        
        vis_min_x = 2800 // 14
        vis_min_y = 3850 // 14
        
        # Create a drawing context on the image
        draw = ImageDraw.Draw(img)
        
        # Draw the visible wedge
        draw.line([(vis_min_x, vis_min_y), (min_pix[0], min_pix[1])], fill="red", width=2)
        draw.line([(vis_min_x, vis_min_y), (max_pix[0], max_pix[1])], fill="red", width=2)
        
        # Draw the target points
        draw.ellipse((min_pix[0] - 5, min_pix[1] - 5, min_pix[0] + 5, min_pix[1] + 5), fill="blue")
        draw.ellipse((max_pix[0] - 5, max_pix[1] - 5, max_pix[0] + 5, max_pix[1] + 5), fill="blue")
        
        # Update the image with the new drawings
        img_tk = ImageTk.PhotoImage(img)
        image_label.config(image=img_tk)
        image_label.image = img_tk
    
    terminal_text.insert(END, "These are the available Galactic longitudes: " + str(avail_long) + "\n")
    return avail_long

# Restore original buttons and galactic longitude entry

activate_antenna_button = customtkinter.CTkButton(control_frame, text="Activate Antenna", command=activate_antenna_clicked)
activate_antenna_button.pack(pady=10)

antenna_status_button = customtkinter.CTkButton(control_frame, text="Show Antenna Status", command=show_ant_status_clicked)
antenna_status_button.pack(pady=10)

shut_down_antenna_button = customtkinter.CTkButton(control_frame, text="Shut Down Antenna", command=shut_down_antenna_clicked)
shut_down_antenna_button.pack(pady=10)

avail_targets_button = customtkinter.CTkButton(control_frame, text="List Available Targets", command=list_avail_targets_clicked)
avail_targets_button.pack(pady=10)

gal_long_label = Label(control_frame, text="Enter Galactic Longitude:")
gal_long_label.pack()

gal_long_entry = Entry(control_frame)
gal_long_entry.pack()

# Start the GUI
root.mainloop()
