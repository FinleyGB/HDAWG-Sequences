import zhinst.utils as zu
import numpy as np
import os
from measurements.libs.QPLser.AWGmanager import HDAWG_PLser
from scipy.signal import chirp

##  Initiliase HDAWG system
device = 'dev8416'  # device ID for G14.
awgMod = HDAWG_PLser(device)
command_table=1

# awgMod.set_value(f"/{device}/awgs/0/enable", 1)
# awgMod.set_value(f"/{device}/awgs/1/enable", 1)

awgMod.set_value(f"/{device}/system/awg/channelgrouping", 1) # sets channel grouping 0 for 2x2, 2 for 1x4

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

# Reading Pulse Parameters
centre_freq_reading=250E6 # Hz; Central frequency set to read out the burned spectral hole
chirpAmplitude_reading = 0.03 # V; amplitude of reading-out pulse
freq_sweeping_reading=5E6 # Hz; set the scanning frequency range of reading-out pulse (the actual scanning range should be 4*freq_detuning) 1.32877326E6
startFreq_reading = (centre_freq_reading-freq_sweeping_reading) # Hz; starting frequency of the burn-back
stopFreq_reading = (centre_freq_reading+freq_sweeping_reading) # Hz; ending frequency of the burn-back
reading_duration = 4e-3 # s; reading-out time
chirpSamples_reading = round(reading_duration*sampling_rate) # Calculating sampling points 
reading_arr = np.linspace(0,reading_duration,chirpSamples_reading)

chirpSine_reading = chirpAmplitude_reading*chirp(reading_arr, startFreq_reading, reading_duration, stopFreq_reading, method='linear') # chieped sine wave for reading pulse

# Shuffle Pulse Parameters
centre_freq_shuffle=250E6 # Hz; Central frequency of the shuffle pulse
chirpAmplitude_shuffle = 0.21 # V; amplitude of shuffling pulse
freq_sweeping_shuffle=20E6 # Hz; set the scanning frequency range of shuffling pulse
shuffle_duration=10e-3 # s; shuffling time

# Stark Shift Pulse Parameters
Number_of_stark_shifts = 20 # Number of stark shift E fields applied to the sample
rect_amplitude_min = 0  # V; Minimum amplitude sent to amplifier
rect_amplitude_max = 1  # V; Maxmimum amplitude sent to amplifier
rectSamples_stark = round(reading_duration*sampling_rate)

zero_pad_duration = 1e-3 # s; Time for electric field to come on before reading pulse
zero_pad_Samples = round(zero_pad_duration*sampling_rate)

Stark_Samples_total = 2*zero_pad_Samples + chirpSamples_reading

# Stark_Shift_Pulse = np.ones((1,Stark_Samples_total), dtype= int)

chirpSine_reading_zero_pad = np.concatenate((np.zeros((1,zero_pad_Samples), dtype= int),
                            chirpSine_reading,
                            np.zeros((1,zero_pad_Samples), dtype= int)), axis = None)

## Create Command Table ## 

if command_table==1:
    # Save Location 
    save_directory = 'C:/HDAWG_control/measurements/libs/QPLser/ZI_HDAWG_Scripts/Command_Tables/'

    # file name
    file_name= 'CT_Stark_shift_YSO_1_half_to_1_half'
    # Open the file
    f = open(os.path.join(save_directory+ file_name), 'w')
    # Write the basic intro of the file
    f.write('{' + '\n' + '  '
        '"$schema": "http://docs.zhinst.com/hdawg/commandtable/v2/schema",' + '\n' + '  '
        '"header": {' + '\n' + '\t' + '"version": "0.2",' + '\n' + '\t' 
        '"UserString": "' + file_name + '",' + '\n' + '\t'  + '"partial": true,' + '\n' + '\t'  + '"description": "Command table for Stark shift measurement"' +  '\n' + '  },' + '\n'
        '  "table": [' + '\n')
    # Write the index assigning

    a = np.linspace(rect_amplitude_min,rect_amplitude_max,Number_of_stark_shifts)


    for i in range(0, len(a)): 
        f.write('\t'+'{' + '\n'
                + '\t' + '"index":'+str( i)+',' + '\n'
                + '\t' + '  "waveform": {' + '\n'
                + '\t' + '\t'  +  '  "index":'+str(10) + '\n'
                + '\t' + '  },' + '\n'
                + '\t' + '  "amplitude0": {' + '\n'
                + '\t' + '\t'  +   '"value": ' + str(1.0) + '\n'
                + '\t' + '  },' + '\n'
                + '\t' + '  "amplitude1": {' + '\n'
                + '\t' + '\t'  +   '"value": ' + str(a[i]) + '\n'
                + '\t' + '\t'  +   '}' + '\n'
                )
        if i==len(a)-1:
            f.write('\t' + '}' + '\n')
        else:
            f.write('\t' + '},' + '\n')
# end parenthesis
    f.write('  ]' + '\n' + '}')
# close the file
    f.close() 

# Load sequence file
HDAWG_filename = ('C:\HDAWG_control\HDAWG-Sequences\Stark_Shift_Measurements\Stark_Shift_YSO_1_half_to_1_half.txt')

with open(HDAWG_filename, "r") as file:
    awg_string = file.read()
    awg_program = awg_string.format(
        
        # Clock parameters
        sampling_rate = sampling_rate,

        # Burning parameters
        Number_of_burning_pulses = Number_of_burning_pulses,
        centre_freq_burning = centre_freq_burning,
        chirpAmplitude_burning = chirpAmplitude_burning,
        freq_sweeping_burning = freq_sweeping_burning,
        burning_duration = burning_duration,

        # Burning back parameters
        Number_of_burning_back_pulses = Number_of_burning_back_pulses,
        centre_freq_burning_back = centre_freq_burning_back,
        chirpAmplitude_burning_back = chirpAmplitude_burning_back,
        freq_sweeping_burning_back = freq_sweeping_burning_back,
        burning_back_duration = burning_back_duration,

        # Cleaning parameters
        Number_of_cleaning_pulses = Number_of_cleaning_pulses,
        centre_freq_cleaning = centre_freq_cleaning,
        chirpAmplitude_cleaning = chirpAmplitude_cleaning,
        freq_sweeping_cleaning = freq_sweeping_cleaning,
        cleaning_duration = cleaning_duration,

        # Reading parameters
        centre_freq_reading = centre_freq_reading,
        chirpAmplitude_reading = chirpAmplitude_reading,
        freq_sweeping_reading = freq_sweeping_reading,
        reading_duration = reading_duration,

        # Shuffle Parameters
        centre_freq_shuffle = centre_freq_shuffle,
        chirpAmplitude_shuffle = chirpAmplitude_shuffle,
        freq_sweeping_shuffle = freq_sweeping_shuffle,
        shuffle_duration = shuffle_duration,

        # Stark Shift Pulse Parameters
        Number_of_stark_shifts = Number_of_stark_shifts,
        zero_pad_duration = zero_pad_duration
    )

awgMod.compile(device, awg_program)
awgMod.upload_command_table(device, command_table='CT_Stark_shift_YSO_1_half_to_1_half')

# Convert and send them to the instrument

marker_0 = np.zeros((1,Stark_Samples_total), dtype= int)
# marker_1 = np.zeros((1,rectSamples_stark), dtype= int)

marker_0 = np.concatenate((np.zeros((1,Stark_Samples_total * 5 // 10), dtype= int),
                          np.ones((1,Stark_Samples_total * 1 // 10), dtype= int),
                          np.zeros((1,Stark_Samples_total * 4 // 10), dtype= int)), axis = None)   # Marker on core 0

# marker_1 = np.concatenate((np.zeros((1,rectSamples_stark * 5 // 10), dtype= int),
#                           np.ones((1,rectSamples_stark * 1 // 10), dtype= int),
#                           np.zeros((1,rectSamples_stark * 4 // 10), dtype= int)), axis = None)   # Marker on core 1

wave_reading_pulse = zu.convert_awg_waveform(chirpSine_reading_zero_pad, markers= marker_0)
# wave_rect_pulse = zu.convert_awg_waveform(Stark_Shift_Pulse, markers= marker_1)

# awgMod.daq.set([(f"/{device}/awgs/0/waveform/waves/10", wave_reading_pulse),
#                 (f"/{device}/awgs/1/waveform/waves/10", wave_rect_pulse)])

awgMod.daq.set(f'/{device:s}/awgs/0/waveform/waves/10', wave_reading_pulse)
# awgMod.daq.set(f'/{device}/awgs/1/waveform/waves/10', wave_rect_pulse)

# Set HDAWG parameters/settings

awgMod.set_value(f"/{device}/sines/0/enables/0", 0)

awgMod.set_value(f"/{device}/triggers/out/0/source", 4) # set up trigger, Output 1&2 Marker 1
# awgMod.set_value(f"/{device}/triggers/out/1/source", 4) # set up trigger, Output 3&4 Marker 1

# Setup output channels
awgMod.set_value(f"/{device}/sigouts/0/on", 1) # Channel 1 is ON
awgMod.set_value(f"/{device}/sigouts/1/on", 1) # Channel 2 is OFF
awgMod.set_value(f"/{device}/sigouts/2/on", 1) # Channel 3 is ON
awgMod.set_value(f"/{device}/sigouts/3/on", 1) # Channel 4 is OFF

awgMod.set_value(f"/{device}/awgs/0/single",0) #Rerun sequence