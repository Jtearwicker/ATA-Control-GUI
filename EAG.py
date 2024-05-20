import customtkinter
import numpy as np
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, Angle
from astropy.time import Time
import datetime
from pytz import timezone
import pytz

# Create the main root
root = customtkinter.CTk()
root.title("Easy ATA GUI")
root.geometry("1200x800")
#root.configure(background='black')
customtkinter.set_appearance_mode("Light")

#Control tabs
control_frame = customtkinter.CTkFrame(master=root, height=350, width=600)
control_frame.pack()

# Create a frame for the terminal output
terminal_frame = customtkinter.CTkFrame(master=root)
terminal_frame.pack()
# Create a text widget to display the terminal output
terminal_text =  customtkinter.CTkTextbox(master=terminal_frame, height=400, width=1200, font=("DejaVu Sans Mono", 12))
terminal_text.pack()

def run_command(command):
  stream = os.popen(command)
  out = stream.read()
  terminal_text.insert(0.0, out)

def run_server_command(command):
    return os.popen(command)

antennas = ['1a']
freq = "1420.406"
def activate_antenna_clicked():
	run_command("/opt/ata-flowgraphs/usrp_reset_clocking.py")
	terminal_text.insert(0.0, "USRP clocking reset.\n")
	run_command("/opt/ata-flowgraphs/usrp_test.py")
	run_server_command("python /home/vgajjar/reu-2023/Hydrogen_line/server.py")
	terminal_text.insert(0.0, "Connected to server.\n")
	ant_free = str(ac.list_antenna_group('none'))
    if ant_free.find('1a') == -1:
    	terminal_text.insert(0.0, "WARNING: Antenna 1a has already been reserved.\n")
  	else:
	    ac.move_ant_group(antennas, 'none', 'atagr')
	    terminal_text.insert(0.0, "You have reserved Antenna 1a!\n")
    ac.set_freq(freq, antennas, 'd')
    ac.autotune(antennas)
    terminal_text.insert(0.0, "Frequency set to 1420.406 MHz and autotuned.\n")
    terminal_text.insert(0.0, ac.get_ascii_status())

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

activate_antenna_button = customtkinter.CTkButton(master=control_frame, text="Activate Antenna", command=activate_antenna_clicked)
activate_antenna_button.pack(padx=5, pady=5)

shut_down_antenna_button = customtkinter.CTkButton(master=control_frame, text="Shut Down Antenna", command=shut_down_antenna_clicked)
shut_down_antenna_button.pack(padx=5, pady=5)