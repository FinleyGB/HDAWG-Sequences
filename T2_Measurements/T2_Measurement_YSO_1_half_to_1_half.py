import zhinst.utils as zu
import numpy as np
import os
from measurements.libs.QPLser.AWGmanager import HDAWG_PLser

##  Initiliase HDAWG system
device = 'dev8416'  # device ID for G14.
awgMod = HDAWG_PLser(device)
command_table=1

## Common Parameter Definition
#Clock Parameters
sampling_rate = 2.4E9 # Hz; Sampling rate, should be the same as what you set on the right panel!

# Burning Pulse Parameters
Number_of_burning_pulses = 300 # Number of burning pulse repetitions
centre_freq_burning = 250E6 #Hz; Central frequency set to drive the AOM for burning
chirpAmplitude_burning = 0.10 # V; amplitude of burning pulse
freq_sweeping_burning=3.5E6 # Hz; set the scanning frequency range of burning pulse (the actual scanning range should be 4*freq_detuning)
burning_duration=0.6e-3 # s; burning time

# Burning back Pulse Parameters
Number_of_burning_back_pulses = 200 # Number of burn-back pulse repetitions 
centre_freq_burning_back=263.75E6 # Hz; Central frequency set to drive the AOM for burn-back
chirpAmplitude_burning_back = 0.047 # V; amplitude of burn-back
freq_sweeping_burning_back=0.8e6 # Hz; set the scanning frequency range of burn-back (the actual scanning range should be 4*freq_detuning)
burning_back_duration=0.1e-3 # s; burning time burn-back

# Cleaning Pulse Parameters
Number_of_cleaning_pulses = 260 # Number of clean pulse repetitions 
centre_freq_cleaning=252.7E6 # Hz; Central frequency set to drive the AOM for cleaning
chirpAmplitude_cleaning = 0.05 #  V; amplitude of cleaning
freq_sweeping_cleaning=1E6 # Hz; set the scanning frequency range of cleaning (the actual scanning range should be 4*freq_detuning)
cleaning_duration=0.5e-3 # s; burning time cleaning

# Shuffle Pulse Parameters
centre_freq_shuffle=250E6 # Hz; Central frequency of the shuffle pulse
chirpAmplitude_shuffle = 0.21 # V; amplitude of shuffling pulse
freq_sweeping_shuffle=20E6 # Hz; set the scanning frequency range of shuffling pulse
shuffle_duration=10e-3 # s; shuffling time

# T2 Parameters
Number_of_T2_measurements = 20  # Number of T2 measuremenets
tau_start = 9e-6    # s; initial time delay between pi/2 and pi pulses
tau_stop = 49e-6    # s; final time delay between pi/2 and pi pulses

centre_freq_pi = 252.7E6 # Hz; Central frequency set to drive the AOM for pi/2 and pi pulses
pi_pulse_duration = 4.74E-6 # s; Pulse duration of pi pulse

# Pi/2 Pulse Parameters
amplitude_half_pi = 0.097 # V; amplitude of pi/2 pulse

# Pi Pulse Parameters
amplitude_pi = 0.097 # V; amplitude of pi pulse

# Load sequence file
HDAWG_filename = ('C:\Codes\HDAWG\Sequences\T2_Measurement_YSO_1_half_to_1_half.txt')

with open(HDAWG_filename, "r") as file:
    awg_string = file.read()
    awg_program = awg_string.format(

        # Clock Parameters
        sampling_rate = sampling_rate,

        # Burning Parameters
        Number_of_burning_pulses = Number_of_burning_pulses,
        centre_freq_burning = centre_freq_burning,
        chirpAmplitude_burning = chirpAmplitude_burning,
        freq_sweeping_burning = freq_sweeping_burning,
        burning_duration = burning_duration,

        # Burning Back Parameters
        Number_of_burning_back_pulses = Number_of_burning_back_pulses,
        centre_freq_burning_back = centre_freq_burning_back,
        chirpAmplitude_burning_back = chirpAmplitude_burning_back,
        freq_sweeping_burning_back = freq_sweeping_burning_back,
        burning_back_duration = burning_back_duration,

        # Cleaning Parameters
        Number_of_cleaning_pulses = Number_of_cleaning_pulses,
        centre_freq_cleaning = centre_freq_cleaning,
        chirpAmplitude_cleaning = chirpAmplitude_cleaning,
        freq_sweeping_cleaning = freq_sweeping_cleaning,
        cleaning_duration = cleaning_duration,

        # T2 Parameters
        Number_of_T2_measurements = Number_of_T2_measurements,
        tau_start = tau_start,
        tau_stop = tau_stop,

        centre_freq_pi = centre_freq_pi,
        pi_pulse_duration = pi_pulse_duration,

        # Pi/2 Parameters
        amplitude_half_pi = amplitude_half_pi,

        # Pi Parameters
        amplitude_pi = amplitude_pi,

        # Shuffle Parameters
        centre_freq_shuffle = centre_freq_shuffle,
        chirpAmplitude_shuffle = chirpAmplitude_shuffle,
        freq_sweeping_shuffle = freq_sweeping_shuffle,
        shuffle_duration = shuffle_duration
    )

awgMod.compile(device, awg_program)

# Set HDAWG parameters/settings

awgMod.set_value(f"/{device}/sines/0/enables/0", 0)

awgMod.set_value(f"/{device}/triggers/out/0/source", 4) # set up trigger, Output 1 Marker 1

# Setup output channels
awgMod.set_value(f"/{device}/sigouts/0/on", 1) # Channel 1 is ON
awgMod.set_value(f"/{device}/sigouts/1/on", 0) # Channel 2 is OFF
awgMod.set_value(f"/{device}/sigouts/2/on", 1) # Channel 3 is ON
awgMod.set_value(f"/{device}/sigouts/3/on", 0) # Channel 4 is OFF

awgMod.set_value(f"/{device}/awgs/0/single",0) #Rerun sequence