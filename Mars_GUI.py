import tkinter as tk
from tkinter import *
import customtkinter
import os
import time
import datetime
import pytz
from ATATools import ata_control as ac

# Create the main root
root = tk.Tk()
root.title("Mars Orbiter Observation GUI")
root.geometry("1000x600")  # Smaller window size

# Create the main frame
main_frame = customtkinter.CTkFrame(master=root)
main_frame.pack(fill=BOTH, expand=True)

# Control frame on the top for buttons
control_frame = customtkinter.CTkFrame(master=main_frame, height=100)
control_frame.pack(side=TOP, fill=X, padx=10, pady=10)

# Terminal frame at the bottom for feedback
terminal_frame = customtkinter.CTkFrame(master=main_frame)
terminal_frame.pack(side=BOTTOM, fill=BOTH, expand=True, padx=10, pady=10)

# Create a text widget to display terminal output
terminal_text = customtkinter.CTkTextbox(master=terminal_frame, height=400, width=800, font=("DejaVu Sans Mono", 16))
terminal_text.pack(fill=BOTH, expand=True)

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
freq = "8431.0"  # Frequency for Tianwen-1

def activate_antenna_clicked():
    run_reset_command("/opt/ata-flowgraphs/usrp_reset_clocking.py")
    run_server_command("python /home/gnuradio/Hydrogen_line/server.py")
    ant_free = str(ac.list_antenna_group('none'))
    if ant_free.find('1a') == -1:
        terminal_text.insert(0.0, "WARNING: Antenna 1a has already been reserved.\n")
    else:
        ac.move_ant_group(antennas, 'none', 'atagr')
        terminal_text.insert(0.0, "Antenna 1a has been reserved.\n")
        ac.set_freq(freq, antennas, 'd')
        ac.autotune(antennas)
        terminal_text.insert(0.0, "Frequency set to 8431 MHz and autotuned.\n")
        terminal_text.insert(0.0, ac.get_ascii_status()[:348] + "\n")
        time.sleep(45)
        terminal_text.insert(0.0, "Calibration complete!\n")

def show_ant_status_clicked():
    terminal_text.insert(0.0, ac.get_ascii_status()[:348] + "\n")

def shut_down_antenna_clicked():
    run_server_command("pkill -f \"python /home/gnuradio/Hydrogen_line/server.py\"")
    terminal_text.insert(0.0, "Disconnected from server.\n")
    ac.park_antennas(antennas)
    terminal_text.insert(0.0, "Antenna 1a has been parked.\n")
    ant_free = str(ac.list_antenna_group('none'))
    if ant_free.find('1a') != -1:
        terminal_text.insert(0.0, "Antenna 1a has already been released\n")
    else:
        ac.move_ant_group(antennas, 'atagr', 'none')
        terminal_text.insert(0.0, "Antenna 1a has been released.\n")

def go_to_Mars_clicked():
    ac.track_source(antennas, source='mars')
    terminal_text.insert(0.0, "Slewing to Mars.\n")

def go_away_from_Mars_clicked():
    ac.track_source(antennas, azel=[180.0, 40.0])
    terminal_text.insert(0.0, "Slewing away from Mars (Azimuth = 180, Elevation = 40).\n")

# Create buttons in the control frame arranged horizontally
activate_antenna_button = customtkinter.CTkButton(
    control_frame, width=150, height=40, text="Activate Antenna", font=("Arial", 16), command=activate_antenna_clicked
)
activate_antenna_button.pack(side=LEFT, padx=5, pady=5)

Mars_button = customtkinter.CTkButton(
    control_frame, width=150, height=40, text="Go To Mars", font=("Arial", 16), command=go_to_Mars_clicked
)
Mars_button.pack(side=LEFT, padx=5, pady=5)

away_from_Mars_button = customtkinter.CTkButton(
    control_frame, width=150, height=40, text="Go Away From Mars", font=("Arial", 16), command=go_away_from_Mars_clicked
)
away_from_Mars_button.pack(side=LEFT, padx=5, pady=5)

antenna_status_button = customtkinter.CTkButton(
    control_frame, width=150, height=40, text="Show Antenna Status", font=("Arial", 16), command=show_ant_status_clicked
)
antenna_status_button.pack(side=LEFT, padx=5, pady=5)

shut_down_antenna_button = customtkinter.CTkButton(
    control_frame, width=150, height=40, text="Shut Down Antenna", font=("Arial", 16), command=shut_down_antenna_clicked
)
shut_down_antenna_button.pack(side=LEFT, padx=5, pady=5)

# Start the GUI
root.mainloop()
