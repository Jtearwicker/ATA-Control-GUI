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
from PIL import Image

# Create the main root
root = customtkinter.CTk()
root.title("Easy ATA GUI")
root.geometry("1200x800")
#root.configure(background='black')
customtkinter.set_appearance_mode("Light")

#Control tabs
control_frame = customtkinter.CTkFrame(master=root, height=350, width=600)
control_frame.pack()

#Galaxy pointing image frame
#image_frame = customtkinter.CTkFrame(master=root, height=350, width=600)
#image_frame.pack()

# Create a frame for the terminal output
terminal_frame = customtkinter.CTkFrame(master=root)
terminal_frame.pack(side=BOTTOM)
# Create a text widget to display the terminal output
terminal_text =  customtkinter.CTkTextbox(master=terminal_frame, height=400, width=1200, font=("DejaVu Sans Mono", 12))
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
	terminal_text.insert(0.0, "Connecting to server...\n")
	return os.popen(command)
	time.sleep(45)
	terminal_text.insert(0.0, "Successfully connected to server.\n")

antennas = ['1a']
freq = "1420.406"
def activate_antenna_clicked():
	run_reset_command("/opt/ata-flowgraphs/usrp_reset_clocking.py")
	run_test_command("/opt/ata-flowgraphs/usrp_test.py")
	run_server_command("python /home/vgajjar/reu-2023/Hydrogen_line/server.py")
	ant_free = str(ac.list_antenna_group('none'))
	if ant_free.find('1a') == -1:
		terminal_text.insert(0.0, "WARNING: Antenna 1a has already been reserved.\n")
	else:
		ac.move_ant_group(antennas, 'none', 'atagr')
		terminal_text.insert(0.0, "You have reserved Antenna 1a!\n")
	ac.set_freq(freq, antennas, 'd')
	ac.autotune(antennas)
	terminal_text.insert(0.0, "Frequency set to 1420.406 MHz and autotuned.\n")
	terminal_text.insert(0.0, ac.get_ascii_status()[:348]+"\n")

def show_ant_status_clicked():
  terminal_text.insert(0.0, ac.get_ascii_status()[:348]+"\n")

def shut_down_antenna_clicked():
	run_server_command("pkill -f \"python /home/vgajjar/reu-2023/Hydrogen_line/server.py\"")
	terminal_text.insert(0.0, "Disconnected from server.\n")
	terminal_text.insert(0.0, "Parking Antenna 1a...\n")
	ac.park_antennas(antennas)
	terminal_text.insert(0.0, "Antenna 1a has been parked.\n")
	ant_free = str(ac.list_antenna_group('none'))
	if ant_free.find('1a') != -1:
		terminal_text.insert(0.0, "Antenna 1a has already been released\n")
	else:
		ac.move_ant_group(antennas, 'atagr', 'none')
		terminal_text.insert(0.0, "Antenna 1a has been released.\n")

#target availability calculator
targets = [[0,0],[10,0],[20,0],[30,0],[40,0],[50,0],[60,0],[70,0],[80,0],
[90,0],[100,0],[110,0],[120,0],[130,0],[140,0],[150,0],[160,0],[170,0],
[180,0],[190,0],[200,0],[210,0],[220,0],[230,0],[240,0],[250,0],[260,0],
[270,0],[280,0],[290,0],[300,0],[310,0],[320,0],[340,0],[350,0]]

ATA_location = EarthLocation(lat=40.817 * u.deg, lon=121.47 * u.deg, height=3235 * u.m)
utc = pytz.timezone('UTC')
now = utc.localize(datetime.datetime.utcnow())
la = pytz.timezone('America/Los_Angeles')
obs_time = now.astimezone(la)
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

		if elevation>20:
			terminal_text.insert(0.0, "Galactic corrdinate "+str(targets[i])+" has an elevation of "+str(elevation)[0:4] +" degrees. RA = "+str(int(RA[0]))+"h"+str(int(RA[1]))+"m"+str(int(RA[2]))+"s"+" Dec = "+str(int(DEC[0]))+"d"+str(int(abs(DEC[1])))+"m"+str(int(abs(DEC[2])))+"s"+".\n")

'''
	ga_min = min(avail_targets_long)
	ga_max = max(avail_targets_long)
	min_pos = int(ga_min/10)
	max_pos = int(ga_max/10)
	px_min = pxlib[min_pos]
	px_max = pxlib[max_pos]
	vis_min_x = [2800, px_min[0]]
	vis_min_y = [3850, px_min[1]]
	vis_max_x = [2800, px_max[0]]
	vis_max_y = [3850, px_max[1]]
	milky_way_img = image.imread(MWimg.jpg)
	plt.plot(vis_min_x, vis_min_y, color="white", linewidth=2)
	plt.plot(vis_max_x, vis_max_y, color="white", linewidth=2)
	plt.plot(2800, 3850, marker='o', color="white")
	plt.imshow(data) 
	#plt.show()
	vis_image = plt.savefig('vis_image.png')
	MW_image = customtkinter.CTkImage(light_image = Image.open("vis_image.png"))
'''
	          
def decdeg2dms(dd):
   is_positive = dd >= 0
   dd = abs(dd)
   minutes,seconds = divmod(dd*3600,60)
   degrees,minutes = divmod(minutes,60)
   degrees = degrees if is_positive else -degrees
   return str(degrees)+"d"+str(minutes)+"m"+str(seconds)+"s"

def track_source_clicked():
	gl = int(galactic_longitude_entry.get())
	dd_ra_dec = ga2equ([gl,0])
	dms_ra = decdeg2dms(dd_ra_dec[0])
	dms_dec = decdeg2dms(dd_ra_dec[1])
	#print(ra)
	#print(type(ra))
	ac.track_source(antennas, radec=[Angle(dms_ra).deg, Angle(dms_dec).deg])
	terminal_text.insert(0.0, "Arrived at RA "+dms_ra+" Dec "+dms_dec+"\n")

activate_antenna_button = customtkinter.CTkButton(master=control_frame, text="Activate Antenna", command=activate_antenna_clicked)
activate_antenna_button.pack(padx=5, pady=5)

show_ant_status_button = customtkinter.CTkButton(master=control_frame, text="Show Antenna Status", command=show_ant_status_clicked)
show_ant_status_button.pack(padx=5, pady=5)

avail_targets_button = customtkinter.CTkButton(master=control_frame, text="Show Available Targets", command=list_avail_targets_clicked)
avail_targets_button.pack(padx=5, pady=5)

'''
ra_entry = customtkinter.CTkEntry(master=control_frame, placeholder_text="RA")
dec_entry = customtkinter.CTkEntry(master=control_frame, placeholder_text="Dec")
ra_entry.pack(padx=5, pady=5)
dec_entry.pack(padx=5, pady=5)
'''
galactic_longitude_entry = customtkinter.CTkEntry(master=control_frame, placeholder_text="Galactic Longitude")
galactic_longitude_entry.pack(padx=5, pady=5)

track_source_button = customtkinter.CTkButton(master=control_frame, text="Track Source", command=track_source_clicked)
track_source_button.pack(padx=5, pady=5)

shut_down_antenna_button = customtkinter.CTkButton(master=control_frame, text="Shut Down Antenna", command=shut_down_antenna_clicked)
shut_down_antenna_button.pack(padx=5, pady=5)

#check_var = customtkinter.StringVar(value="off", command=trigger_check_box)
#checkbox = customtkinter.CTkCheckBox(master=control_frame, text="Show Galaxy Plot", variable=check_var, onvalue="on", offvalue="off")

'''
def plot_galaxy(ga_min_entry, ga_max_entry, ga_obs_entry):
		if(len(ga_min_entry.get())==0):
			ga_min = 0
		else:
			ga_min = ga_min_entry.get()
		if(len(ga_max_entry.get())==0):
			ga_max = 0
		else:	
			ga_max = ga_max_entry.get()
		if(len(ga_obs_entry.get())==0):
			ga_obs = 0
		else:
			ga_obs = ga_obs_entry.get()


		if(check_var.get()==1):
			min_pos = int(float(ga_min))
			max_pos = int(float(ga_max))
			obs_pos = int(float(ga_obs))
			px_min = pxlib[min_pos]
			px_max = pxlib[max_pos]
			px_obs = pxlib[obs_pos]
			vis_min_x = [2800, px_min[0]]
			vis_min_y = [3850, px_min[1]]
			vis_max_x = [2800, px_max[0]]
			vis_max_y = [3850, px_max[1]]
			obs_x = [2800, px_obs[0]]
			obs_y = [3850, px_obs[1]]
			milky_way_img = image.imread("MWimg.jpg")
			plt.plot(vis_min_x, vis_min_y, color="white", linewidth=2)
			plt.plot(vis_max_x, vis_max_y, color="white", linewidth=2)
			plt.plot(obs_x, obs_y, color="red", linewidth=2)
			plt.plot(2800, 3850, marker='o', color="white")
			plt.imshow(milky_way_img) 
			plt.show()

gl_min_entry = customtkinter.CTkEntry(master=control_frame, placeholder_text="GL Min")
gl_min_entry.pack(padx=5, pady=5)

gl_max_entry = customtkinter.CTkEntry(master=control_frame, placeholder_text="GL Min")
gl_max_entry.pack(padx=5, pady=5)

gl_obs_entry = customtkinter.CTkEntry(master=control_frame, placeholder_text="GL Min")
gl_obs_entry.pack(padx=5, pady=5)

plot_galaxy_button = customtkinter.CTkButton(master=control_frame, text="Plot Available Sources", command=plot_galaxy(gl_min_entry,gl_max_entry,gl_obs_entry))
plot_galaxy_button.pack(padx=5, pady=5)
'''
root.mainloop()