import obspy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, convolve
from obspy import read, Trace, UTCDateTime
from scipy.io.wavfile import read
from scipy.fftpack import fft
import matplotlib.gridspec as gridspec
import os
import math

def calculate_dispersion(f, H, b1, c, b2, u1, u2):
    # Calculate the left-hand side (LHS)
    y = H * math.sqrt(1 / b1**2 - 1 / c**2)
    lhs = math.tan(2 * math.pi * f * y)

    # Check if the denominator is zero
    if b1 == 1 or c == b1:
        return "No dispersion"

    # Calculate the right-hand side (RHS)
    rhs = (u2 * math.sqrt(1 - c**2 / b2**2)) / (u1 * math.sqrt(c**2 / b1**2 - 1))

    # Compare LHS and RHS to determine dispersion
    if lhs == rhs:
        return "No dispersion"
    elif lhs < rhs:
        return "Normal dispersion"
    else:
        return "Anomalous dispersion"

def calculate_depth_kernel(f, H, b1, c, b2, u1, u2):
    # Calculate the depth kernel
    y = H * math.sqrt(1 / b1**2 - 1 / c**2)
    
    # Check if the denominator is zero
    if c == b1:
        return float('inf')
    
    depth_kernel = (2 * math.pi * f * y) / (c * math.sqrt(c**2 / b1**2 - 1))

    return depth_kernel

# Define the parameters
f = 0.5  # Frequency
H = 10  # Height
b1 = 3  # Lower bound
b2 = 5  # Upper bound
d1 = 2.8 #density 1
d2 = 3.2 #density 2
u1 = 25.2 #d1* b1**2 # Parameter u1
u2 = 80 #d2 * b2**2 # Parameter u2

# Vary c in the range of b1 to b2
c_values = []
dispersion_values = []
depth_kernels = []

for c in range(b1, b2+1):
    dispersion = calculate_dispersion(f, H, b1, c, b2, u1, u2)
    depth_kernel = calculate_depth_kernel(f, H, b1, c, b2, u1, u2)
    c_values.append(c)
    dispersion_values.append(dispersion)
    depth_kernels.append(depth_kernel)

# Print the results
for c, dispersion, depth_kernel in zip(c_values, dispersion_values, depth_kernels):
    print(f"c = {c}, Dispersion: {dispersion}, Depth Kernel: {depth_kernel}")

# Calculate the values for c and LHS
c_values = np.linspace(b1, b2, 50000)  # Changed the third argument to 'num' instead of a float value
lhs_values = []
for c in c_values:
    y = H * math.sqrt(1 / b1**2 - 1 / c**2)
    lhs = math.tan(2 * np.pi * f * y)
    lhs_values.append(lhs)

#Calculate the values for c and RHS
rhs_values = []
for c in c_values:
    if c == b1:
        rhs_values.append(float('inf'))
    else:
        rhs = (u2 * math.sqrt(1 - c**2 / b2**2)) / (u1 * math.sqrt(c**2 / b1**2 - 1))
        rhs_values.append(rhs)

diff =[]
diff = np.array(rhs_values) - np.array(lhs_values)
diff_filtered = diff[np.logical_and(diff < 0.1, diff > -0.001)]
diff_filtered_index = np.where(diff == diff_filtered[0])
print(c_values[diff_filtered_index], diff_filtered_index)
diff_filtered_index_1 = np.where(diff_filtered)
print(c_values[diff_filtered_index_1][:2])
print(diff_filtered[:5])
print(rhs_values[795])
#print(c_values[diff_filtered_index])

# Plot the LHS curve
plt.plot(c_values, rhs_values, color='red',label='RHS')
plt.plot(c_values, lhs_values, color='blue',label='LHS')
plt.scatter(c_values[diff_filtered_index],rhs_values[795], color='green',label='Difference')
plt.xlabel("c")
plt.title("dispersion value")
plt.xlim(b1, b2)  # Set the x-axis limits from b1 to b2
plt.ylim(-10,30)
plt.xticks(np.arange(b1-0.1, b2+0.1, 0.2))  # Set x-axis ticks at every π/2 interval
plt.legend()
plt.grid(True)
plt.savefig("LHS_RHS_curve.png")


    
