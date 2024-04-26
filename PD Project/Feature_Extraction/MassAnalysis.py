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
""" tgrounded = []
for num,data in enumerate(alldata):
    tgrounded.append(np.zeros(len(data)))
    for time, force in enumerate(data):
        if force[18] > 80:
            tgrounded[num][int(force[0])] = 1 """
#print(tgrounded[3])
""" allsteps=[]
for mask in tgrounded:
    start = 0
    steps = []
    for i in range(len(mask)):
        if mask[i] == 0 and mask[i-1] == 1:
            end = i
            steps.append((start,end))
        if mask[i] == 1 and mask[i-1] == 0:
            start = i
    allsteps.append(steps) """

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
""" for num,intervals in enumerate(allsteps):
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
    np.savetxt("csvforce\MaxRF"+ filenames[num] +".csv",peaks) """
    #Cross correlating between feet: Between columns 1-5,2-6,3-7,4-8, etc.

""" corr = signal.correlate(alldata[0][:,2],alldata[0][:,10])
lags = signal.correlation_lags(len(alldata[0][:,10]),len(alldata[0][:,2]))
plt.plot(lags,corr)
plt.show() """
#Get center of pressure measurements and graph displacement of cop during stride. Find variability
#FFT 
#Welch PSD

def find(val,arr):
    return np.argmin(np.abs(arr-val))

N=8192
allharmonics = []
for feature in alldata:
    xf, yf = signal.welch(feature[500:8692,19],100,"hamming",nperseg=N/4,nfft=N)
    #Find peaks
    peaks, properties = signal.find_peaks(20*np.log(yf[:500]), prominence = 30)
    # Get ratio of harmonics vs fundamental amplitude
    fundamental = 20*np.log(yf[peaks[0]])
    fundfreq = xf[peaks[0]]
    #print(peaks)
    harmrat = [20*np.log(yf[find(fundfreq*2,xf[peaks[0]:500])+peaks[0]])/fundamental,20*np.log(yf[find(fundfreq*3,xf[peaks[0]:500])+peaks[0]])/fundamental,20*np.log(yf[find(fundfreq*4,xf[peaks[0]:500])+peaks[0]])/fundamental,20*np.log(yf[find(fundfreq*5,xf[peaks[0]:500])+peaks[0]])/fundamental]
    allharmonics.append(harmrat)
    print((fundfreq,xf[find(fundfreq*2,xf[peaks[0]:500])+peaks[0]],xf[find(fundfreq*3,xf[peaks[0]:500])+peaks[0]],xf[find(fundfreq*4,xf[peaks[0]:500])+peaks[0]]))


with open(r"harmonics/results2.txt", 'w') as fp:
    for n,item in enumerate(allharmonics):
        # write each item on a new line
        fp.write(filenames[n]+ "%s\n" % item)
    print('Done')



