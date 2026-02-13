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
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from PIL import Image
import time

# Create the main root
root = customtkinter.CTk()
root.title("Easy ATA GUI")
root.geometry("1200x800")
customtkinter.set_appearance_mode("Light")

# Control tabs
control_frame = customtkinter.CTkFrame(master=root, height=350, width=600)
control_frame.pack()

# Create a frame for the terminal output
terminal_frame = customtkinter.CTkFrame(master=root)
terminal_frame.pack(side=BOTTOM)
# Create a text widget to display the terminal output
terminal_text = customtkinter.CTkTextbox(master=terminal_frame, height=400, width=1200, font=("DejaVu Sans Mono", 12))
terminal_text.pack()

def run_test_command(command):
    terminal_text.insert(0.0, "Testing USRPs...\n")
    stream = os.popen(command)
    out = stream.read()
    terminal_text.insert(0.0, out)

def run_reset_command(command):
    terminal_text.insert(0.0, "Resetting USRPs clocking...\n")
    stream = os.popen(command)
    out = stream.read()
    terminal_text.insert(0.0, out)

def run_server_command(command):
    return os.popen(command)

antennas = ['1a']
freq = "1420.406"
def activate_antenna_clicked():
    run_reset_command("/opt/ata-flowgraphs/usrp_reset_clocking.py")
    run_server_command("python /home/vgajjar/reu-2023/Hydrogen_line/server.py")
    ant_free = str(ac.list_antenna_group('none'))
    if ant_free.find('1a') == -1:
        terminal_text.insert(0.0, "WARNING: Antenna 1a has already been reserved.\n")
    else:
        ac.move_ant_group(antennas, 'none', 'atagr')
        terminal_text.insert(0.0, "Antenna 1a has been reserved.\n")
    ac.set_freq(freq, antennas, 'd')
    ac.autotune(antennas)
    terminal_text.insert(0.0, "Frequency set to 1420.406 MHz and autotuned.\n")
    terminal_text.insert(0.0, ac.get_ascii_status()[:348]+"\n")
    time.sleep(45)
    terminal_text.insert(0.0, "Calibration complete!\n")

def show_ant_status_clicked():
    terminal_text.insert(0.0, ac.get_ascii_status()[:348]+"\n")

def shut_down_antenna_clicked():
    run_server_command("pkill -f \"python /home/vgajjar/reu-2023/Hydrogen_line/server.py\"")
    terminal_text.insert(0.0, "Disconnected from server.\n")
    ac.park_antennas(antennas)
    terminal_text.insert(0.0, "Antenna 1a has been parked.\n")
    ant_free = str(ac.list_antenna_group('none'))
    if ant_free.find('1a') != -1:
        terminal_text.insert(0.0, "Antenna 1a has already been released\n")
    else:
        ac.move_ant_group(antennas, 'atagr', 'none')
        terminal_text.insert(0.0, "Antenna 1a has been released.\n")

# target availability calculator
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
    for i in range(0,35):
        dd_radec = ga2equ(targets[i])
        c = SkyCoord(ra = dd_radec[0]*u.deg, dec = dd_radec[1] * u.deg)
        RA = c.ra.hms
        DEC = c.dec.dms
        elevation = radec2alt(ga2equ(targets[i]))

        if elevation > 20:
            terminal_text.insert(0.0, "Galactic longitude "+str(targets[i][0])+" has an elevation of "+str(elevation)[0:4]+" degrees above the horizon.\n")

def track_source_clicked():
    gl = int(galactic_longitude_entry.get())
    dd_radec = ga2equ([gl,0])
    c = SkyCoord(ra = dd_radec[0]*u.deg, dec = dd_radec[1] * u.deg)
    RA = c.ra.hms
    DEC = c.dec.dms
    ac.track_source(antennas, radec=[Angle(str(int(RA[0]))+"h"+str(int(RA[1]))+"m"+str(int(RA[2]))+"s").hour, Angle(str(int(DEC[0]))+"d"+str(int(abs(DEC[1])))+"m"+str(int(abs(DEC[2])))+"s").deg])
    terminal_text.insert(0.0, "Arrived at galactic coordinate ("+str(gl)+",0)."+" RA "+str(int(RA[0]))+"h"+str(int(RA[1]))+"m"+str(int(RA[2]))+"s"+" Dec "+str(int(DEC[0]))+"d"+str(int(abs(DEC[1])))+"m"+str(int(abs(DEC[2])))+"s\n")

def plot_galaxy(ga_min, ga_max, ga_obs):
    min_pos = int(ga_min/10)
    max_pos = int(ga_max/10)
    obs_pos = int(ga_obs/10)
    px_min = pxlib[min_pos]
    px_max = pxlib[max_pos]
    px_obs = pxlib[obs_pos]
    vis_min_x = [2800, px_min[0]]
    vis_min_y = [3850, px_min[1]]
    vis_max_x = [2800, px_max[0]]
    vis_max_y = [3850, px_max[1]]
    obs_x = [2800, px_obs[0]]
    obs_y = [3850, px_obs[1]]
    data = image.imread("C:\\Users\\jtear\\Downloads\\Artist's_impression_of_the_Milky_Way_(updated_-_annotated) (1).jpg")
    fig, ax = plt.subplots()
    ax.imshow(data)
    ax.plot(vis_min_x, vis_min_y, color="white", linewidth=2)
    ax.plot(vis_max_x, vis_max_y, color="white", linewidth=2)
    ax.plot(obs_x, obs_y, color="red", linewidth=2)
    ax.plot(2800, 3850, marker='o', color="white")
    return fig

def display_galaxy_plot():
    ga_min = int(min_longitude_entry.get())
    ga_max = int(max_longitude_entry.get())
    ga_obs = int(observed_longitude_entry.get())
    fig = plot_galaxy(ga_min, ga_max, ga_obs)
    canvas = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas.draw()
    canvas.get_tk_widget().pack()

# GUI for testing
testing = customtkinter.CTkLabel(master=control_frame, text="USRP testing controls", font=("Roboto", 24))
testing.grid(row=0, column=0, padx=10, pady=10)

reset_button = customtkinter.CTkButton(master=control_frame, text="Reset USRPs", command=lambda:run_reset_command("/opt/ata-flowgraphs/usrp_reset_clocking.py"))
reset_button.grid(row=1, column=0, padx=10, pady=10)

test_button = customtkinter.CTkButton(master=control_frame, text="Test USRPs", command=lambda:run_test_command("uhd_find_devices"))
test_button.grid(row=2, column=0, padx=10, pady=10)

# GUI for ATA
ata = customtkinter.CTkLabel(master=control_frame, text="ATA controls", font=("Roboto", 24))
ata.grid(row=0, column=1, padx=10, pady=10)

activate_button = customtkinter.CTkButton(master=control_frame, text="Activate Antenna 1a", command=activate_antenna_clicked)
activate_button.grid(row=1, column=1, padx=10, pady=10)

show_status_button = customtkinter.CTkButton(master=control_frame, text="Show Antenna Status", command=show_ant_status_clicked)
show_status_button.grid(row=2, column=1, padx=10, pady=10)

shutdown_button = customtkinter.CTkButton(master=control_frame, text="Shut Down Antenna", command=shut_down_antenna_clicked)
shutdown_button.grid(row=3, column=1, padx=10, pady=10)

list_avail_targets_button = customtkinter.CTkButton(master=control_frame, text="List Available Targets", command=list_avail_targets_clicked)
list_avail_targets_button.grid(row=4, column=1, padx=10, pady=10)

track_source_button = customtkinter.CTkButton(master=control_frame, text="Track Source", command=track_source_clicked)
track_source_button.grid(row=5, column=1, padx=10, pady=10)

galactic_longitude_label = customtkinter.CTkLabel(master=control_frame, text="Galactic Longitude:")
galactic_longitude_label.grid(row=6, column=0, padx=10, pady=10)
galactic_longitude_entry = customtkinter.CTkEntry(master=control_frame)
galactic_longitude_entry.grid(row=6, column=1, padx=10, pady=10)

min_longitude_label = customtkinter.CTkLabel(master=control_frame, text="Min Galactic Longitude:")
min_longitude_label.grid(row=7, column=0, padx=10, pady=10)
min_longitude_entry = customtkinter.CTkEntry(master=control_frame)
min_longitude_entry.grid(row=7, column=1, padx=10, pady=10)

max_longitude_label = customtkinter.CTkLabel(master=control_frame, text="Max Galactic Longitude:")
max_longitude_label.grid(row=8, column=0, padx=10, pady=10)
max_longitude_entry = customtkinter.CTkEntry(master=control_frame)
max_longitude_entry.grid(row=8, column=1, padx=10, pady=10)

observed_longitude_label = customtkinter.CTkLabel(master=control_frame, text="Observed Galactic Longitude:")
observed_longitude_label.grid(row=9, column=0, padx=10, pady=10)
observed_longitude_entry = customtkinter.CTkEntry(master=control_frame)
observed_longitude_entry.grid(row=9, column=1, padx=10, pady=10)

plot_button = customtkinter.CTkButton(master=control_frame, text="Plot Galaxy", command=display_galaxy_plot)
plot_button.grid(row=10, column=0, columnspan=2, padx=10, pady=10)

# Create a frame for the plot
plot_frame = customtkinter.CTkFrame(master=root)
plot_frame.pack()

root.mainloop()
