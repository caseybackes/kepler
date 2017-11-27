import kplr 
import matplotlib.pyplot as plt
import scipy.stats as st
from PyAstronomy.pyasl import foldAt
import math
import numpy as np
import os
global os

from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText


##############################################
def findIndex(listToSearch, maxVal):
    '''For "listToSearch", give list or an array (which will be temporarily recast as a list). 
    For "maxVal", give the value or max value to find in listToSearch. For example, 
    "myList":
        499.90  index = 0 
        500.00  index = 1
        500.10  index = 2
    and findIndex(myList, 500.1)  will return an index of 2
    
    '''
    myList = list(listToSearch)          
    for index in range(len(myList)):
        if myList[index] > maxVal:
            return index-1
            break
##############################################
def EstablishKOI():
    '''This function will let you choose what KOI in your chosen system to investigate.
    '''
    global koi,usercall,period
    totalPlanetsInSystem=0
    KOIlist = []
    
    usercall = int(input("Enter the KOI system to lookup in MAST: "))
    client = kplr.API()
    '''
    while client.koi(usercall+.01) == ValueError:
        print "Sorry, no KOI associated with that number. Please choose another."
        usercall
    else:
        koi = client.koi(usercall+.01)
    '''
    koi = client.koi(usercall+.01)
    period= koi.koi_period
    #Determine how many KOI are in the system associated with the user's choice
    print"Item   KOI#                 KOI Name    Period (days)"
    print"_"*50
    for PlanetInSystem in range(1,15,1):

        try: 
            KOI = client.koi(float(usercall+.01*PlanetInSystem))
            totalPlanetsInSystem+=1
            KOIlist = np.append(KOIlist, KOI)
            print "(",PlanetInSystem,") ",KOI," ",KOI.kepoi_name," ",round(KOI.koi_period,1)
            
        except ValueError:
            break #exit the for loop, all planets have been found. 
               
    print '\n',totalPlanetsInSystem, " total KOI in system ",usercall
    KOIChoice = int(input("Choose a specific KOI from above (enter Item #): "))
    while KOIChoice > totalPlanetsInSystem :
        print "That is not a valid option, try again."
        KOIChoice = int(input("Choose a specific KOI from above (enter Item #): "))
    else:
        
        koi = client.koi(usercall+KOIChoice*.01)
        print "That is ", koi.kepoi_name
############################################
############################################

def Normalize_LightCurves():
           
    # based on user input, declare the operation to acquire the light curve data. 
    print"Working..."
    lcs = koi.get_light_curves(short_cadence = True, clobber = True)

    # Initiate empty arrays into which we append data from MAST
    # time array - used extensively for data plotting and model plotting
    # flux array - used to build the fluxmean, fluxmean2 arrays
    # fluxmean(2) arrays - used for data plotting and normalization. 
    # data array - a copy of the fluxmean2 array used for plotting lc data. 
    global time, flux, fluxmean, fluxmean2, data
    
    time,data = [],[]#flux,fluxmean,fluxmean2,data =[],[],[],[],[]
    for lc in lcs:
        with lc.open() as f:
            hdu_data=f[1].data
            time=np.append(time,hdu_data["time"])
            #flux=np.append(flux,(hdu_data["pdcsap_flux"]))
            #fluxmean=np.append(fluxmean, (hdu_data['pdcsap_flux'])/(np.nanmean(hdu_data["pdcsap_flux"])))
            # Fluxmean2 is for finding the base intensity of the host star. 
            #fluxmean2=np.append(fluxmean2, (hdu_data['pdcsap_flux'])/(np.nanmean(hdu_data["pdcsap_flux"])))
            data    = np.append(data,      (hdu_data['pdcsap_flux'])/(np.nanmean(hdu_data["pdcsap_flux"])))
    print "Complete."
    print "Generating plot of normalized, preconditioned data."
    ###############################################     
    '''
    Fix all the values in the time array to be numbers on cadence with the 
    time interval of the data. There are several sections in the time array that
    are huge groups of NaN's. These need to be fixed. For each index in the time
    array, if it is a nan, re-declare it to be time[i-1]+ (time interval). This 
    should fix all NaN's in the time array.   
    '''   
    for i in range(len(time)): 
        if np.isnan(time[i]) == True:
            # point = previous time  +   time interval
            time[i] =  time[i-1]     +  (time[5]-time[4])
    ###################################################   
    '''
    Establish vertical lines in the plot to show the transits by first finding the 
    time values for the centroids of each transit. This is found by by adding the
    "first observed centriod" (as given by MAST) plus an iteration times the period, 
    for each iteration that goes from 1 to the expected number of transits. The 
    expected number of transits is found by dividing the number of days that Kepler 
    operated for the primary mission, divided by the period of the KOI in question. 
    The number of expected transits will be larger for smaller-period KOI.
    At each of these "centroid" times, a vertical line is placed on the plot from 
    y = 1 to y = (lowest observed flux value, or "min(y)"). 
    '''
    centroids = []
    for i in np.arange(0,max(time)/period):
        nextCentroid = koi.koi_time0bk+i*period
        if nextCentroid < max(time):
            centroids = np.append(centroids, nextCentroid)
    print "centriod list:", centroids
    
    try:
        # state the uncertainty of the duration as what is found by MAST. 
        durationErr1 = round(koi.koi_duration_err1,1)
        if type(durationErr1) =="NoneType":
            print "none type for durration error 1"
    except TypeError:
        # if the duration error is not given by MAST, take it as the square root of the duration ( in hours).
        durationErr1 = str("$\pm$ "+str(round(math.sqrt(koi.koi_duration),1)))
    ###############################################
    '''
    The plot of the normalized light curves:
    '''
    
    plt.figure(figsize=(24,14))
    plt.plot(time,data,label ='Normalized Flux')
    
    # a horizontal line that shows the "1" value:
    plt.hlines(1.0000,min(time), max(time),'r')
    plt.title("Light Curves for KOI {}".format(koi))
    plt.xlabel("BJD minus 2455200 (days)")
    plt.ylabel("Normalized Flux Intensity")
    ax = plt.gca()
    ax.get_xaxis().get_major_formatter().set_useOffset(False)

    ###############################################
    '''
    The "Characteristics" legend: 
    '''
    KOI_Properties = AnchoredText("            Characteristics\nPeriod: {} days\nFirst Centroid: {} days\nDisposition: {}\nDuration: {} {} hours\nKepID: {}".format(round(period,1), round(koi.koi_time0bk,1),koi.koi_disposition,round(koi.koi_duration,2),durationErr1,koi.kepid),loc= 2, prop=dict(size = 10),frameon=True)
    KOI_Properties.patch.set_boxstyle("round, pad = 0.5,rounding_size= 0.2")
    plt.gca().add_artist(KOI_Properties)
    
    plt.legend()
    plt.show()
    ###############################################
def main():
    EstablishKOI()
    Normalize_LightCurves()
main()
##############################
'''
# Phase folding the light curves
centroid_list = []
timeunit = time[5] - time[4]
first_centroid = koi.koi_time0bk
period = koi.koi_period
expected_transits = int(round(max(time)/period,0))
for each in range(0,expected_transits,1):
    next_centroid = first_centroid+each*period
    centroid_list = np.append(centroid_list,next_centroid)
print centroid_list


periodIndexRange = int(round(period/timeunit ,0))
halfperiod = int(round(periodIndexRange/2,0))
firstCentIndex = findIndex(time, first_centroid)

temptime = list(time)
tempdata = list(data)

tempdata = list(tempdata)
plt.plot(tempdata[2*firstCentIndex -halfperiod: 2*firstCentIndex + halfperiod])
del tempdata[0:halfperiod+firstCentIndex+1]
plt.plot(tempdata[0:halfperiod+firstCentIndex])
'''













