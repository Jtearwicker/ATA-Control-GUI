import tkinter as tk
from tkinter import *
from tkinter import ttk
import customtkinter
#import webview
import os
import subprocess
import numpy as np
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from math import *
import numpy as np
import pyautogui
from tkinterweb import HtmlFrame
customtkinter.CTkButton

server = None

# Create the main root
root = customtkinter.CTk()
root.title("ATA Control GUI")
root.geometry("1200x800")
root.configure(background='black')
customtkinter.set_appearance_mode("Dark")

#Configure frames

#Control tabs

tabview = customtkinter.CTkTabview(master=root, height=350, width=600)
tabview.pack()
tabview.add("Calibrate")  
tabview.add("Antenna Setup")
tabview.add("Observe")  
tabview.set("Calibrate") 

# Create a frame for the terminal output
terminal_frame = customtkinter.CTkFrame(master=root)
terminal_frame.pack()
# Create a text widget to display the terminal output
terminal_text =  customtkinter.CTkTextbox(master=terminal_frame, height=400, width=1200, font=("DejaVu Sans Mono", 12))
#terminal_text.configure(state="disabled")  # configure textbox to be read-only
terminal_text.pack()

#Create a frame for the ATA camera view
camera_frame = customtkinter.CTkFrame(master=root, height=200, width=300)
#webview.create_window('ATA Live View', 'http://10.3.0.30/view/view.shtml?id=342&imagepath=%2Fmjpg%2Fvideo.mjpg&size=1')
#webview.start()


def run_command(command):
  stream = os.popen(command)
  out = stream.read()
  terminal_text.insert(0.0, out)

def run_server_command(command):
    return os.popen(command)


#calibration tabs
def test_usrp_clicked():
    run_command("/opt/ata-flowgraphs/usrp_test.py")


test_usrp_button = customtkinter.CTkButton(master=tabview.tab("Calibrate"), text="Test USRP", command=test_usrp_clicked)
test_usrp_button.pack(padx=5, pady=5)

def reset_clocking_clicked():
    run_command("/opt/ata-flowgraphs/usrp_reset_clocking.py")

reset_clocking_button = customtkinter.CTkButton(master=tabview.tab("Calibrate"), text="Reset Clocking", command=reset_clocking_clicked)
reset_clocking_button.pack(padx=5, pady=5)

def server_clicked():
    server = run_server_command("python /home/vgajjar/reu-2023/Hydrogen_line/server.py &")

server_button = customtkinter.CTkButton(master=tabview.tab("Calibrate"), text="Connect to Server", command=server_clicked)
server_button.pack(padx=5, pady=5)

def disconnect_server_clicked():
    if server:
        print("Closing Server")
        server.kill()

disconnect_server_button = customtkinter.CTkButton(master=tabview.tab("Calibrate"), text="Disconnect from Server", command=disconnect_server_clicked)
disconnect_server_button.pack(padx=5, pady=5)


#antenna setup tabs
from ATATools import ata_control as ac

def show_ant_status_clicked():
  terminal_text.insert(0.0, ac.get_ascii_status())

show_ant_status_button = customtkinter.CTkButton(master=tabview.tab("Antenna Setup"), text="Show Antenna Status", command=show_ant_status_clicked)
show_ant_status_button.pack(padx=5, pady=5)

def reserve_ant_clicked():
  antennas=['1a']
  ac.move_ant_group(antennas, 'none', 'atagr')

reserve_ant_button = customtkinter.CTkButton(master=tabview.tab("Antenna Setup"), text="Reserve Antenna", command=reserve_ant_clicked)
reserve_ant_button.pack(padx=5, pady=5)

freq_text = customtkinter.CTkTextbox(master=tabview.tab("Antenna Setup"), height=100, width=210,fg_color="transparent")
freq_text.pack(padx=5, pady=5)
freq_text.insert(tk.END, "Enter frequency below (MHz)\nUse 1420.405 for Hydrogen")

freq_entry = customtkinter.CTkEntry(master=tabview.tab("Antenna Setup"), width=80)
freq_entry.pack(padx=5, pady=5)

def set_freq_and_autotune_clicked():
  ac.set_freq({freq_entry}, antennas, 'd')
  ac.autotune(antennas)

set_freq_and_autotune_button = customtkinter.CTkButton(master=tabview.tab("Antenna Setup"), text="Set frequency and \nautotune antennas", command=set_freq_and_autotune_clicked)
set_freq_and_autotune_button.pack(padx=5, pady=5)



#observation tabs

#target availability calculator
targets = [[0,0],[10,0],[20,0],[30,0],[40,0],[50,0],[60,0],[70,0],[80,0],
[90,0],[100,0],[110,0],[120,0],[130,0],[140,0],[150,0],[160,0],[170,0],
[180,0],[190,0],[200,0],[210,0],[220,0],[230,0],[240,0],[250,0],[260,0],
[270,0],[280,0],[290,0],[300,0],[310,0],[320,0],[340,0],[350,0]]

ATA_location = EarthLocation(lat=40.817 * u.deg, lon=-121.47 * u.deg, height=3235 * u.m)
obs_time = Time('2024-05-07 15:00:00')
alt_az = AltAz(location=ATA_location, obstime=obs_time)

def radec2alt(RADEC):
    coord = SkyCoord(RADEC[0] * u.deg, RADEC[1] * u.deg)
    aa = coord.transform_to(alt_az)
    #alt_str = str(aa.alt.value)[0:4]
    #print("The elevation is "+alt_str+" degrees")
    altitude = aa.alt.value
    return altitude

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

def list_avail_targets_clicked():
    for i in range(0,35):
        elevation = radec2alt(ga2equ(targets[i]))
        if elevation>20:
            print("Galactic corrdinate "+str(targets[i])+" has an elevation of "+str(elevation)[0:4])


date_time_text = customtkinter.CTkTextbox(master=tabview.tab("Observe"), height=100, width=190,fg_color="transparent")
date_time_text.pack(padx=5, pady=5)
date_time_text.insert(tk.END, "Enter the current date and time\nin format (YYYY-MM-DD HH:MM:SS)")


date_time_entry = customtkinter.CTkEntry(master=tabview.tab("Observe"), placeholder_text="Date Time")
date_time_entry.pack(padx=5, pady=5)

avail_targets_button = customtkinter.CTkButton(master=tabview.tab("Observe"), text="Show Available Targets", command=list_avail_targets_clicked)
avail_targets_button.pack(padx=5, pady=5)

def track_source_clicked():
  ac.track_source(antennas, radec=[Angle('').hour, Angle('').deg])

ra_dec_entry = customtkinter.CTkEntry(master=tabview.tab("Observe"), placeholder_text="RA Dec")
ra_dec_entry.pack(padx=5, pady=5)



root.mainloop()
