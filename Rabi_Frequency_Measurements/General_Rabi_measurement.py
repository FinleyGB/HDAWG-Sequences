'''
This codes are used to do general Rabi measurements
'''

from measurements.libs.QPLser.AWGmanager import HDAWG_PLser

##  Initiliase HDAWG system
device = 'dev8416'  # device ID for G14.
awgMod = HDAWG_PLser(device)

sampling_rate=2.4E9 # Hz; Sampling rate, should be the same as what you set on the right panel!

centre_freq=250E6 # Hz; Central frequency set to drive the AOM for the long pulse for Rabi measurement
amplitude = 0.21 # V; Amplitude of long pulse
pulse_duration = 2E-6 # s; Long pulse duration 

centre_freq_shuffle=80E6 #Hz; Central frequency of the shuffle pulse
chirpAmplitude_shuffle = 0.6 # V; amplitude of shuffling pulse
freq_sweeping_shuffle=20E6 # Hz; set the scanning frequency range of shuffling pulse
shuffle_duration=5e-3 #s; shufling time

# Load sequence file
HDAWG_filename = ('C:\HDAWG_control\HDAWG-Sequences\Rabi_Frequency_Measurements\General_Rabi_measurement.txt')

with open(HDAWG_filename, "r") as file:
    awg_string = file.read()
    awg_program = awg_string.format(

        sampling_rate=sampling_rate, # Hz; Sampling rate, should be the same as what you set on the right panel!
        
        centre_freq=centre_freq, # Hz; Central frequency set to drive the AOM for the long pulse for Rabi measurement
        amplitude = amplitude, # V; Amplitude of long pulse
        pulse_duration = pulse_duration, # s; Long pulse duration

        centre_freq_shuffle=centre_freq_shuffle, #Hz; Central frequency of the shuffle pulse
        chirpAmplitude_shuffle = chirpAmplitude_shuffle, # V; amplitude of shuffling pulse
        freq_sweeping_shuffle=freq_sweeping_shuffle, # Hz; set the scanning frequency range of shuffling pulse
        shuffle_duration=shuffle_duration #s; shufling time
    )

awgMod.compile(device, awg_program)

# Set HDAWG parameters/settings

awgMod.set_value(f"/{device}/sines/0/enables/0", 0)

awgMod.set_value(f"/{device}/system/awg/oscillatorcontrol", 0) # Disable AWG Oscillator Control
awgMod.set_value(f"/{device}/triggers/out/0/source", 4) # Set up trigger, Output 1 Marker 1

# Setup output channels
awgMod.set_value(f"/{device}/sigouts/0/on", 1) # Channel 1 is ON
awgMod.set_value(f"/{device}/sigouts/1/on", 0) # Channel 2 is OFF
awgMod.set_value(f"/{device}/sigouts/2/on", 0) # Channel 3 is OFF
awgMod.set_value(f"/{device}/sigouts/3/on", 0) # Channel 4 is OFF

awgMod.set_value(f"/{device}/awgs/0/single",0) #Rerun sequence