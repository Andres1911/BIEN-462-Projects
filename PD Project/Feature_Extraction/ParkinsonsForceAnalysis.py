import numpy as np 
import os,csv
import matplotlib.pyplot as plt
from scipy import signal
from scipy.fft import fft, fftfreq

os.chdir("c:/Users/trivi/Documents/Bien462/csvdata")
directory = os.getcwd()
alldata = []
start = 0
getall = False
threshold = 7

if getall:
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        if os.path.isfile(f):
            patient = open(filename)
            read= list(csv.reader(patient))
            data = np.array(read)
            alldata.append(data)
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
tgrounded = np.zeros(len(alldata[0]))
for time, force in enumerate(alldata[0]):
    if force[19] > 80:
        tgrounded[int(force[0])] = 1

start = 0
steps = []
stride = 0
for i in range(len(tgrounded)):
    if tgrounded[i] == 0 and tgrounded[i-1] == 1:
        end = i
        stride = stride + 1
        steps.append(np.array(alldata[0][start:end,0]))
    if tgrounded[i] == 1 and tgrounded[i-1] == 0:
        start = i

#Validating decomposition of data into strides 
#plt.plot(tgrounded)
fig, axs = plt.subplots(2)
axs[0].plot(alldata[0][:,19]/max(alldata[0][:,19]))
axs[0].set_title('Control')
axs[1].plot(alldata[1][:,19]/max(alldata[1][:,19]))
axs[1].set_title('PD')
fig.suptitle('PD vs Control')
plt.show()

#Peak force for each stride and variability

#Peak total force for each stride, find variability

#Cross correlating between feet: Between columns 1-5,2-6,3-7,4-8, etc.

""" corr = signal.correlate(alldata[0][:,2],alldata[0][:,10])
lags = signal.correlation_lags(len(alldata[0][:,10]),len(alldata[0][:,2]))
plt.plot(lags,corr)
plt.show() """
#Get center of pressure measurements and graph displacement of cop during stride. Find variability
#FFT 
N=8192
yf1 = fft(alldata[0][500:8692,19])
yf2 = fft(alldata[1][500:8692,19])
xf = fftfreq(8192,0.01)[:8192//2]
""" x= np.linspace(0,81.92,0.01)
yf1 = np.sin() """


fig, axs = plt.subplots(2)
axs[0].plot(xf, 2.0/N * 20*np.log10(np.abs(yf1[0:N//2])))
axs[0].set_title('Control')
axs[1].plot(xf, 2.0/N * 20*np.log10(np.abs(yf2[0:N//2])))
axs[1].set_title('PD')
fig.suptitle('PD vs Control')
plt.grid()
plt.show()

#get ratio of next harmonics