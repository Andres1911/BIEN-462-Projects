import numpy as np 
import os,csv
import matplotlib.pyplot as plt
from scipy import signal
from scipy.fft import fft, fftfreq

os.chdir("c:/Users/trivi/Documents/Bien462/csvdata") # CHANGE THIS TO YOUR CSV DIRECTORY
directory = os.getcwd()
alldata = []
start = 0
getall = True
threshold = 7
filenames = []

if getall:
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        if os.path.isfile(f):
            patient = open(filename)
            read= list(csv.reader(patient))
            data = np.array(read).astype(float)
            alldata.append(data)
            filenames.append(f[-13:-4])
else:
    patient = open("GaCo01_01.csv")
    read = list(csv.reader(patient))
    data = np.array(read).astype(float)
    alldata.append(data)
    patient = open("GaPt03_01.csv")
    read = list(csv.reader(patient))
    data = np.array(read).astype(float)
    alldata.append(data)

#Divide data into strides
    #Idea: Determine time indices for each stride
    #Afterwards, extract and place into array
tgrounded = []
for num,data in enumerate(alldata):
    tgrounded.append(np.zeros(len(data)))
    for time, force in enumerate(data):
        if force[18] > 80:
            tgrounded[num][int(force[0])] = 1
#print(tgrounded[3])
allsteps=[]
for mask in tgrounded:
    start = 0
    steps = []
    for i in range(len(mask)):
        if mask[i] == 0 and mask[i-1] == 1:
            end = i
            steps.append((start,end))
        if mask[i] == 1 and mask[i-1] == 0:
            start = i
    allsteps.append(steps)

#Validating decomposition of data into strides 
""" fig, axs = plt.subplots(2)
axs[0].plot(tgrounded)
axs[0].plot(alldata[0][:,19]/max(alldata[0][:,19]))
axs[0].set_title('Control')
axs[1].plot(alldata[1][:,19]/max(alldata[1][:,19]))
axs[1].set_title('PD')
fig.suptitle('PD vs Control')
plt.show() """

#Peak force for each stride and variability
for num,intervals in enumerate(allsteps):
    peaks = []
    peak_totals = []
    intervals.pop(0)
    for start, end in intervals:
        peaks.append(np.max(alldata[num][start:end,2:9]))        
        peak_totals.append(np.max(alldata[num][start:end,19]))
    peak_totals.append(np.max(peak_totals))
    peaks.append(np.max(peak_totals))
    np.savetxt("csvforce\MaxTotRF" + filenames[num] +".csv",peak_totals)
    #Peak total force for each stride, find variability
    np.savetxt("csvforce\MaxRF"+ filenames[num] +".csv",peaks)
    #Cross correlating between feet: Between columns 1-5,2-6,3-7,4-8, etc.

""" corr = signal.correlate(alldata[0][:,2],alldata[0][:,10])
lags = signal.correlation_lags(len(alldata[0][:,10]),len(alldata[0][:,2]))
plt.plot(lags,corr)
plt.show() """
#Get center of pressure measurements and graph displacement of cop during stride. Find variability
#FFT 
""" N=8192
yf1 = fft(alldata[0][500:8692,19]-np.mean(alldata[0][500:8692,19])) # Remove mean to get rid of strong peak at 0 hz
yf2 = fft(alldata[1][500:8692,19]-np.mean(alldata[1][500:8692,19]))
xf = fftfreq(8192,0.01)[:8192//2] """
""" x= np.linspace(0,81.92,0.01)
yf1 = np.sin() """

""" 
fig, axs = plt.subplots(2)
axs[0].plot(xf, 2.0/N * 20*np.log10(np.abs(yf1[0:N//2])))
axs[0].set_title('Control')
axs[1].plot(xf, 2.0/N * 20*np.log10(np.abs(yf2[0:N//2])))
axs[1].set_title('PD')
fig.suptitle('PD vs Control')
plt.grid()
plt.show() """

#Welch PSD

#get ratio of next harmonics