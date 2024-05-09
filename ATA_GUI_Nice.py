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
customtkinter.CTkButton

# Create the main root
root = customtkinter.CTk()
root.title("ATA Control GUI")
root.geometry("1200x700")
root.configure(background='black')

# Create a frame for the terminal output
terminal_frame = Frame(root, height=400, width=800)
terminal_frame.grid(row=1, column=0)

# Create a text widget to display the terminal output
terminal_output = Text(terminal_frame, wrap=WORD, height=20, width=100)
terminal_output.pack(side=LEFT, fill=BOTH, expand=YES)

# Create a scrollbar for the text widget
scrollbar = Scrollbar(terminal_frame, command=terminal_output.yview)
scrollbar.pack(side=RIGHT, fill=Y)
terminal_output.config(yscrollcommand=scrollbar.set)


def run_command(command):
  stream = os.popen(command)
  out = stream.read()
  pyautogui.alert(out)

tabview = customtkinter.CTkTabview(master=root)
tabview.grid(row=0,column=0)

tabview.add("Calibrate")  
tabview.add("Antenna Setup")
tabview.add("Observe")  
tabview.set("Calibrate")  


#calibration tabs
def test_usrp_clicked():
    run_command("/opt/ata-flowgraphs/usrp_test.py")

test_usrp_button = customtkinter.CTkButton(master=tabview.tab("Calibrate"), text="Test USRP", command=test_usrp_clicked)
test_usrp_button.grid(row=1, column=0, padx=5, pady=5)

def reset_clocking_clicked():
    run_command("/opt/ata-flowgraphs/usrp_reset_clocking.py")

reset_clocking_button = customtkinter.CTkButton(master=tabview.tab("Calibrate"), text="Reset Clocking", command=reset_clocking_clicked)
reset_clocking_button.grid(row=2, column=0, padx=5, pady=5)

def server_clicked():
  run_command("python /home/vgajjar/reu-2023/Hydrogen_line/server.py")

server_button = customtkinter.CTkButton(master=tabview.tab("Calibrate"), text="Connect to Server", command=server_clicked)
server_button.grid(row=3, column=0, padx=5, pady=5)



#antenna setup tabs
from ATATools import ata_control as ac

def show_ant_status_clicked():
  print(ac.get_ascii_status())

show_ant_status_button = customtkinter.CTkButton(master=tabview.tab("Antenna Setup"), text="Show Antenna Status", command=show_ant_status_clicked)
show_ant_status_button.grid(row=4, column=0, padx=5, pady=5)

def reserve_ant_clicked():
  antennas=['1a']
  ac.move_ant_group(antennas, 'none', 'atagr')

reserve_ant_button = customtkinter.CTkButton(master=tabview.tab("Antenna Setup"), text="Reserve Antenna", command=reserve_ant_clicked)
reserve_ant_button.grid(row=5, column=0, padx=5, pady=5)

freq_text = customtkinter.CTkTextbox(master=tabview.tab("Antenna Setup"), height=48, width=190,fg_color="transparent")
freq_text.grid(row=6, column=0, padx=5, pady=5)
freq_text.insert(tk.END, "Enter frequency below (MHz)\nUse 1420.405 for Hydrogen")

freq_entry = customtkinter.CTkEntry(master=tabview.tab("Antenna Setup"), width=80)
freq_entry.grid(row=7, column=0, padx=5, pady=5)

def set_freq_and_autotune_clicked():
  ac.set_freq({freq_entry}, antennas, 'd')
  ac.autotune(antennas)

set_freq_and_autotune_button = customtkinter.CTkButton(master=tabview.tab("Antenna Setup"), text="Set frequency and \nautotune antennas", command=set_freq_and_autotune_clicked)
set_freq_and_autotune_button.grid(row=8, column=0, padx=5, pady=5)



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
date_time_text.grid(row=10, column=0, padx=5, pady=5)
date_time_text.insert(tk.END, "Enter the current date and time\nin format (YYYY-MM-DD HH:MM:SS)")


date_time_entry = customtkinter.CTkEntry(master=tabview.tab("Observe"), placeholder_text="Date Time")
date_time_entry.grid(row=11, column=0)

avail_targets_button = customtkinter.CTkButton(master=tabview.tab("Observe"), text="Show Available Targets", command=list_avail_targets_clicked)
avail_targets_button.grid(row=12, column=0, padx=5, pady=5)

def track_source_clicked():
  ac.track_source(antennas, radec=[Angle('').hour, Angle('').deg])

ra_dec_entry = customtkinter.CTkEntry(master=tabview.tab("Observe"), placeholder_text="RA Dec")
ra_dec_entry.grid(row=13, column=0)

#p://10.3.0.30/view/view.shtml?id=342&imagepath=%2Fmjpg%2Fvideo.mjpg&size=1')
#webview.start()

root.mainloop()
