#!/home/inorton/bin/pywin
import comtypes, comtypes.client
from ctypes import *
from comtypes.automation import *
import sys
import os
import tkinter 
from tkinter import filedialog
from tkinter import *
import numpy
def parameters():
    global masses
    global ppmWin
    global minInt
    masses = uiMasses.get()
    ppmWin = uiPpmWin.get()
    minInt = uiMinInt.get()
    mGui.destroy()

def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

# GetMassListFromScanNum(long FAR* pnScanNumber, LPCTSTR szFilter,  
# long nIntensityCutoffType, long nIntensityCutoffValue,  
# long nMaxNumberOfPeaks, BOOL bCentroidResult, 
# VARIANT FAR* pvarMassList, 
# VARIANT FAR* pvarPeakFlags,  
# long FAR* pnArraySize)

#Default values for MSFileReader.GetMassListFromScanNumber call:

# create tkinter gui
mGui = Tk()
mGui.geometry('450x100+500+300')
##mGui.mainloop()
uiMasses = StringVar(mGui, value="355.06994,371.10124,371.31559,391.28429,413.26623,429.08873,445.12003,503.10752") # 356.07776504,372.10906504,372.32341504,392.29211504,414.27405504,430.09655504,446.12785504,504.11534504) + 1.00782504
uiPpmWin = StringVar(mGui, value='20')
uiMinInt = StringVar(mGui, value='100')
mGui.mLabel = Label(text='choose lock masses (comma separated no whitespace):').grid(row=0,column=0,sticky=W)
mGui.mLabel2 = Label(text='select ppm window:').grid(row=1,column=0,sticky=E)
mGui.mLabel3 = Label(text='minimum absolute intensity:').grid(row=2,column=0,sticky=E)
mGui.title('Select lock masses and ppm window')
mGui.mbutton = Button(text='OK',command=parameters).grid(row=3,column=1)
mGui.mEntry = Entry(textvariable=uiMasses).grid(row=0,column=1,columnspan=3)
mGui.mEntry2 = Entry(textvariable=uiPpmWin).grid(row=1,column=1)
mGui.mEntry3 = Entry(textvariable=uiMinInt).grid(row=2,column=1)
mGui.mainloop()
# lock mass values
# split gui text to list
masses = masses.split(',')
# convert to floating point
masses = list(map(float, masses))
# convert ppm window to float or integer is no decimal place
ppmWin = num(ppmWin)
# convert min intensity to float or integer if no decimal place
minInt = num(minInt)

# set mass tolerance
maxmz = max(masses) + 10
minmz = min(masses) - 10
scanFilter = u''
scanIntensityCutoffType = 0 # 0 = none, 1=Abs, 2=Rel. to basepk
scanIntensityCutoffValue = 0
scanMaxNumberOfPeaks = 0
scanCentroidResult = 0
pl = VARIANT() #Unused variable
# ml set up later
arsize = c_long()
# example file
#dirTmp = "J:\\Adductomics\\Bz exposed\\"
dirTmp = filedialog.askdirectory()
dirTmp = dirTmp.replace('/', '\\')
dirTmp = dirTmp + '\\'

print('Calculating mass drift for ' +  str(len(masses)) + ' lock masses identified within a ' + str(ppmWin) + ' ppm window.\n')

for fInName in os.listdir(dirTmp):
    if fInName.endswith((".raw", ".RAW")):
        print("reading file: " + fInName)
        fInName = dirTmp + fInName
        # Set up the COM object interface
        xr = comtypes.client.CreateObject('MSFileReader.XRawfile')
        xr.open(fInName)
        res = xr.SetCurrentController(0,1)
        # print("res: " + str(res))
        ns = c_long()
        xr.GetNumSpectra(ns)
        print("Num scans: " + str(ns.value))
        data = numpy.array('f')
        outStrings = []
        # testing
        ##ml = VARIANT()
        ##xr.GetMassListFromScanNum(c_long(i),scanFilter,c_long(scanIntensityCutoffType),c_long(scanIntensityCutoffValue),
        ##c_long(scanMaxNumberOfPeaks),
        ##c_long(scanCentroidResult),
        ##c_double(0),
        ##ml,pl,arsize
        ##)
        # testing
        # empty results array
        massDriftAr = numpy.zeros((3,ns.value))
        # scan numbers
        massDriftAr[0] = [i for i in range(1,ns.value+1)]
        for i in range(1,ns.value+1):
            ml = VARIANT()
            xr.GetMassListFromScanNum(
                c_long(i),scanFilter,
                c_long(scanIntensityCutoffType),
                c_long(scanIntensityCutoffValue),
                c_long(scanMaxNumberOfPeaks),
                c_long(scanCentroidResult),
                c_double(0),
                ml,pl,arsize
                )
            data = numpy.array(ml.value)
            # between mz values
            rng = [(data[0]>=minmz)&(data[0]<=maxmz)]
            data = numpy.array((data[0][rng],data[1][rng]))
            # empty array to store results
            ppmDrift = numpy.zeros((4,len(masses)))
            # lock masses in first column
            ppmDrift[0] = masses
            for j in range(len(masses)):
                # calculate ppm difference with lock masses and boolean less than ppmWin
                ppmDiff = ((data[0] - masses[j]) / masses[j]) * 1E06
                # within mass window and above zero intensity
                inMassWin = (abs(ppmDiff) <= ppmWin) & (data[1] > minInt)
                if any(inMassWin):
                    # boolean highest intensity within window
                    #minDiff = abs(ppmDiff[inMassWin]) == min(abs(ppmDiff[inMassWin]))
                    highestInt = data[1][inMassWin] == max(data[1][inMassWin])
                    # observed mass
                    #ppmDrift[1,j] = data[0][inMassWin][minDiff]
                    ppmDrift[1,j] = data[0][inMassWin][highestInt]
                    # intensity value
                    #ppmDrift[2,j] = data[1][inMassWin][minDiff]
                    ppmDrift[2,j] = data[1][inMassWin][highestInt]
                    # ppm drift
                    #ppmDrift[3,j] = ppmDiff[inMassWin][minDiff]
                    ppmDrift[3,j] = ppmDiff[inMassWin][highestInt]
            # weighted meanMassDrift
##            zeroIndx = ppmDrift[2] > 0
##            if any(zeroIndx):
##                massDriftAr[1,i - 1] = numpy.average(ppmDrift[3][zeroIndx],weights=ppmDrift[2][zeroIndx])
##            else:
##                massDriftAr[1,i - 1] = float('nan')
            massDriftAr[1,i - 1] = numpy.mean(ppmDrift[3])
            # n lock masses detected
            massDriftAr[2,i - 1] = sum(ppmDrift[1] > 0)
            del(ml)

        # transpose and write csv
        massDriftAr = massDriftAr.transpose()
        fInName = fInName.replace('.RAW', '')
        fInName = fInName.replace('.raw', '')+'.massDrift.csv'
        of = open(fInName,'w',massDriftAr.shape[0])
        of.write('scanNo,ppmDrift,nDetected of ' + str(len(masses)) + ' lock masses\n')
        for item in massDriftAr:
            of.write("%.3f,%.3f,%.0f\n" % (item[0],item[1],item[2]))

        of.flush()
        of.close()
        print('...done.\n')

#ml = c_double()
# create variable type LP_c_double dynamically with ctypes
#LP_c_double = POINTER(c_double)
#ml = c_long(0)
#xr.GetNumberOfMassCalibratorsFromScanNum(c_long(1),ml)

#ml = c_double()
#xr.GetMassCalibrationValueFromScanNum(c_long(1),c_long(7),ml) 
#minmz=600
#maxmz=1000


#stLogRt = 0.0
#arsize = 0
#varLab = VARIANT()
#arsize = c_long()
#xr.GetStatusLogForScanNum(c_long(2),c_double(stLogRt),varLab,varVal,c_long(arsize))
#xr.GetStatusLogLabelsForScanNum(c_long(1),stLogRt,varLab,arsize)
# did not work

#stLogRt = 0.02
#varLab = VARIANT()
#arsize = c_long()
#xr.GetStatusLogLabelsForRT(c_double(stLogRt),varLab,arsize)
# did not work

#varLab = VARIANT()
#varVal = VARIANT()
#arsize = c_long()
#xr.GetTrailerExtraForScanNum(c_long(1),varLab,varVal,arsize)
# did not work

#if (len(sys.argv) < 2):
#    sys.exit(-1)
#else:
#    fInName = sys.argv[1]
