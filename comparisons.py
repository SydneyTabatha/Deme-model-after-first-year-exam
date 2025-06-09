'''
Compare single deme and multideme models and different growth laws 
'''
  
import time
from math import sin
from numpy import array, arange, zeros
from pylab import plot, xlabel, show, legend, yscale, figure, ylabel, title
import os

#Parameters # %-------------------------------------------|
growth_rate = 'saturated'
num_wt = 1
r1 = 0.5 
r2 = 0.49 
mu = 5e-5
d =0.05
epsilon= 1.5
alpha = 0.5
h= 0.001
eta = 0.0001
Patches = 100
threshold = 5#50000

# then need to change these parameters in stochastic models and also update their saturated growth
# also - For the single patch model, multiply the patch growth rate (g in my notation) by n^(1./3). That should give you identical growth curves. 

def deme(growth_law, threshold):
    k = 10000
    K=k

    xpoints =  [] 
    ypoints = [] 
    Kpoints = [k]
    xsum = []
    ysum = []
    timepoints = [0]
    iteration = 0
    t = 0
    current_size = 0
    sizes = []
    
    xnext = num_wt
    ynext = 0
    
    #while timepoints[-1] < 70: #NEED TO USE THIS FOR LOGISITC GROWTH
    while current_size < threshold:
        #print(current_size)
        #print(K)
        
        xpoints.append(xnext)
        ypoints.append(ynext)
        
        x = xpoints[-1]
        y = ypoints[-1]

        
        if growth_law == 'logistic':
            fx_element = (1.0 - mu) * r1 * x * (1.0 - (x + y) / k) - d * x 
            fy_element = mu * r1 * x * (1.0 - (x + y) / k) + r2 * y * (1.0 - (x + y) / k) - d * y 
            
        elif growth_law == 'exponential':
            fx_element = (1.0 - mu) * r1 * x  - d * x 
            fy_element = mu * r1 * x  + r2 * y  - d * y
            
        elif growth_law == 'saturated':
            K = Kpoints[-1]
            fx_element = (1.0 - mu) * r1 * x * (1.0 - (x + y) / K) - d * x 
            fy_element = mu * r1 * x * (1.0 - (x + y) / K) + r2 * y * (1.0 - (x + y) / K) - d * y 
            g = alpha*((epsilon*((x+y)/K)) - 1)
            g= g*(Patches**(1/3))
            if g<0:
                g = 0
            K_element = g*(K**(2/3)) 
            Knext = K + h*K_element
            Kpoints.append(Knext)
            
        else:
            
            print('pick a growth law!')
        
        xnext = x + h*fx_element 
        ynext = y + h*fy_element
        
    
        #current_size = xnext + ynext
        
        iteration += 1
        t += h
        timepoints.append(t)
        
        previous_size = current_size
        current_size = xnext + ynext
        sizes.append(current_size)
        
        
        if current_size >= threshold:
            #print(current_size)
            # Linearly interpolate the mutant number at threshold
            frac = (threshold - previous_size) / (current_size - previous_size)
            estimated_mutants = y + frac * (ynext - y)
            print(f"Estimated number of mutants at threshold {threshold:.0f} is {estimated_mutants:.2f}")
            break
        

        
    # Plot
    
    figure(1)
    plot(timepoints[:-1], xpoints, label='total number of wildtype')
    plot(timepoints[:-1], ypoints, label='total number of mutants')
    legend()
    xlabel("t")
    ylabel("population abundances")
    show()
    print(f"The number of mutants at size {current_size} is ", ypoints[-1])     
    
    figure(2)
    plot(timepoints[:-1], xpoints, label='total number of wildtype')
    plot(timepoints[:-1], ypoints, label='total number of mutants')
    yscale('log')
    legend()
    xlabel("t")
    ylabel("logscale of population abundances")
    show()
    print(f"The number of mutants at size {current_size} is ", ypoints[-1])     
        
    return [timepoints[:-1], ypoints, sizes, xpoints]

# CHANGE TO JUST PRINT EVERY 100 
def multideme(growth_law, threshold, eta): # should epsilon be divided by the number of patches?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    k = 100#10**4
    K=k
    
    xpoints =  [[0] for _ in range(Patches)] 
    xpoints[0] = [num_wt] # add a wildtype to the first patch
    ypoints = [[0] for _ in range(Patches)]  
    xsum = []
    ysum = []
    timepoints = [0]
    iteration = 0
    t = 0
    current_size = 0
    xtemp = zeros(Patches)
    ytemp = zeros(Patches)
    sizes = []
    
    
    # each patch has it's own k
    Kpoints = [[k] for _ in range(Patches)]  
    
    
    while current_size < threshold: # plateaus at 79999.9
    #while timepoints[-1]<70:

        xs = 0 # each iteration start with a brand new sum for summing across populations 
        ys = 0
        
        for n in range(Patches):
            
            # new iteration
            x = xpoints[n][-1] # current x value
            y = ypoints[n][-1] # current y value

    
            if iteration==0: 
                sum_x = sum_y = 0.0
                if Patches ==0:
                    sum_x=1
            else:
                sum_x = sum_x1-x
                sum_y=sum_y1 - y 
            
            if growth_law == 'logistic':
                fx_element = (1.0 - mu) * r1 * x * (1.0 - (x + y) / k) - d * x - eta * (x - sum_x / (Patches - 1))
                fy_element = mu * r1 * x * (1.0 - (x + y) / k) + r2 * y * (1.0 - (x + y) / k) - d * y - eta * (y - sum_y / (Patches - 1))
            
            elif growth_law == 'exponential':
                fx_element = (1.0 - mu) * r1 * x - d * x - eta * (x - sum_x / (Patches - 1))
                fy_element = mu * r1 * x  + r2 * y  - d * y - eta * (y - sum_y / (Patches - 1))
                
            elif growth_law =='saturated':
                K = Kpoints[n][-1]
                fx_element = (1.0 - mu) * r1 * x * (1.0 - (x + y) / K) - d * x - eta * (x - sum_x / (Patches - 1))
                fy_element = mu * r1 * x * (1.0 - (x + y) / K) + r2 * y * (1.0 - (x + y) / K) - d * y - eta * (y - sum_y / (Patches - 1))
                g = alpha*((epsilon*((x+y)/K)) - 1) # K ISNT UPDATING , TRY PRINTING G
                #print('g = ', g)
                if g<0:
                    g = 0                
                K_element = g*(K**(2/3))
                Knext = K+ h*K_element
                Kpoints[n].append(Knext)
                
            else:
                print('pick a growth law!')
            
            # euler method
            xtemp[n] = x + h*fx_element # store xnext as temporary variable so i can update all patches simultaneously at the end
            ytemp[n] = y + h*fy_element
            
    
            xs += xpoints[n][iteration] # im adding a patch per iteration to the total sum
            ys += ypoints[n][iteration]
            
        sum_x1 = sum_y1 = 0.0
        for n in range(Patches): # adding xnext and ynext to the list of timepoints
            
            xpoints[n].append(xtemp[n])
            ypoints[n].append(ytemp[n])   
    
            sum_x1 += xpoints[n][-1] 
            sum_y1 += ypoints[n][-1]  
        
        xsum.append(xs)
        ysum.append(ys) 
        #Kpoints.append(Knext) ?
    
        iteration += 1
        t += h
        timepoints.append(t)
        
        previous_size = current_size
        current_size = xs + ys
        sizes.append(current_size)
        
        # Check if we just crossed the threshold
        
        if current_size >= threshold:
            #print(current_size)
            # Linearly interpolate the mutant number at threshold
            frac = (threshold - previous_size) / (current_size - previous_size)
            estimated_mutants = ysum[-1] + frac * (ys - ysum[-1])
            print(f"Estimated number of mutants at threshold {threshold:.0f} is {estimated_mutants:.2f}")
            break
        

        
    # Plot
    figure(3)
    plot(timepoints[:-1], xsum, label='total number of wildtype')
    plot(timepoints[:-1], ysum, label='total number of mutants')
    legend()
    xlabel("t")
    ylabel("population abundances")
    show()
    
    # Plot
    figure(4)
    plot(timepoints[:-1], xsum, label='total number of wildtype')
    plot(timepoints[:-1], ysum, label='total number of mutants')
    yscale('log')
    legend()
    xlabel("t")
    ylabel("logscale of population abundances")
    show()
    print(f"The number of mutants at size {current_size} is ", ysum[-1])
    
    return [timepoints[:-1], ysum, sizes, xsum]



# call functions
result1 = deme(growth_rate, threshold) # 700000
ypoints = result1[1]
sizes_1 = result1[2]
wt1 = result1[3]
time_single = result1[0]

print('NEXT FUNCTION------------------------------------------')

resultp = multideme(growth_rate, threshold, eta)
ysum = resultp[1]
sizes_patch = resultp[2]
wtp = resultp[3]
time_multi = resultp[0]


## NEW #################################
import glob
import numpy as np
import matplotlib.pyplot as plt

# Gather all replicate files
file_list = sorted(glob.glob("replicate_*.txt"))

# Prepare lists to collect each replicate’s time series
wt_list = []
mut_list = []
size_list = []
time_list = []

for filename in file_list:
    if os.path.getsize(filename) == 0:
        print(f"Skipping empty file: {filename}")
        continue  # skip this file because it's empty
    data = np.loadtxt(filename, skiprows=1)  # Skip header line

    time = data[:, 0]
    wt = data[:, 1]
    mut = data[:, 2]
    size = data[:, 3]

    time_list.append(time)
    wt_list.append(wt)
    mut_list.append(mut)
    size_list.append(size)

# Step 1: Find the shortest length across replicates
min_length = min(len(t) for t in time_list)

# Step 2: Truncate all lists to the shortest length
wt_array = np.array([w[:min_length] for w in wt_list])
mut_array = np.array([m[:min_length] for m in mut_list])
size_array = np.array([s[:min_length] for s in size_list])
time_array = np.array([t[:min_length] for t in time_list])

# Step 3: Compute averages safely
avg_wt = np.mean(wt_array, axis=0)
avg_mut = np.mean(mut_array, axis=0)
avg_size = np.mean(size_array, axis=0)
avg_time = np.mean(time_array, axis=0)


file_list = sorted(glob.glob("MULTIPATCHreplicate_*.txt"))
file_list.extend(sorted(glob.glob("MULTIPATCHreplicate1_*.txt")))
file_list.extend(sorted(glob.glob("MULTIPATCHreplicate2_*.txt")))
file_list.extend(sorted(glob.glob("MULTIPATCHreplicate3_*.txt")))
file_list.extend(sorted(glob.glob("MULTIPATCHreplicate4_*.txt")))
file_list.extend(sorted(glob.glob("MULTIPATCHreplicate5_*.txt")))
file_list.extend(sorted(glob.glob("MULTIPATCHreplicate6_*.txt")))
file_list.extend(sorted(glob.glob("MULTIPATCHreplicate7_*.txt")))
file_list.extend(sorted(glob.glob("MULTIPATCHreplicate8_*.txt")))
file_list.extend(sorted(glob.glob("MULTIPATCHreplicate9_*.txt")))

# Prepare lists to collect each replicate’s time series
mwt_list = []
mmut_list = []
msize_list = []
mtime_list = []

for filename in file_list:
    if os.path.getsize(filename) == 0:
        print(f"Skipping empty file: {filename}")
        continue  # skip this file because it's empty    
    data = np.loadtxt(filename, skiprows=1)  # Skip header line
    mtime = data[:, 0]
    mwt = data[:, 1]
    mmut = data[:, 2]
    msize = data[:, 3]

    mtime_list.append(mtime)
    mwt_list.append(mwt)
    mmut_list.append(mmut)
    msize_list.append(msize)

# Step 1: Find the shortest length across replicates
min_length = min(len(t) for t in mtime_list)

# Step 2: Truncate all lists to the shortest length
mwt_array = np.array([w[:min_length] for w in mwt_list])
mmut_array = np.array([m[:min_length] for m in mmut_list])
msize_array = np.array([s[:min_length] for s in msize_list])
mtime_array = np.array([t[:min_length] for t in mtime_list])

# Step 3: Compute averages safely
mavg_wt = np.mean(mwt_array, axis=0)
mavg_mut = np.mean(mmut_array, axis=0)
mavg_size = np.mean(msize_array, axis=0)
mavg_time = np.mean(mtime_array, axis=0)

figure(5)
#plot(sizes_1, ypoints, label='total number of mutants in single-deme model', linestyle='-', color = 'blue')
#plot(sizes_patch, ysum, label='total number of mutants in multi-deme model', linestyle='-', color = 'green')
plot(avg_size, avg_mut, label="num of mutants in stochastic single-deme model", linestyle='--', color = 'blue')
plot(mavg_size, mavg_mut, label='number of mutants in stochastic multi-deme model', linestyle='--', color = 'green')
legend(fontsize=14)
xlabel("size")
show()

figure(6)
#plot(time_single, wt1, label='total number of wildtypes in single-deme model', linestyle='-', color = 'blue')
#plot(time_multi, wtp, label='total number of wildtypes in multi-deme model', linestyle='-', color = 'green')
plot(avg_time, avg_wt, label="num of wildtypes in stochastic single-deme model", linestyle='--', color = 'blue')
plot(mavg_time, mavg_wt, label='number of wildtypes in stochastic multi-deme model', linestyle='--', color = 'green') 
legend(fontsize=14)
xlabel("time")
show()
