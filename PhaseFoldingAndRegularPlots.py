import kplr 
import matplotlib.pyplot as plt
import scipy.stats as st
from PyAstronomy.pyasl import foldAt
from PyAstronomy.pyasl import binningx0dt
from astropy.time import Time
import math
import numpy as np
import os
global os
import sys
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText
from time import gmtime, strftime


##############################################

##############################################
def EstablishKOI():
    '''This function will let you choose what KOI in your chosen system to investigate.
    '''
    global koi,period, usercall,KOIChoice
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
    global centroids,durationErr1,lcs,time, flux, fluxmean, fluxmean2, data
    # based on user input, declare the operation to acquire the light curve data. 
    print"Working..."
    lcs = koi.get_light_curves(short_cadence = True, clobber = True)
      
  
    # Initiate empty arrays into which we append data from MAST 
    time,data = [],[]
    
    #Open each LC file and extract the data from the "pdcsap_flux" header
    for lc in lcs:
        with lc.open() as f:
            hdu_data=f[1].data
            #print f[1].header
            time=np.append(time,hdu_data["time"])
            
            #   THIS LINE NORMALIZES THE FLUX TO AVERAGE VALUE OF 1 !!
            data    = np.append(data,(hdu_data['pdcsap_flux'])/(np.nanmean(hdu_data["pdcsap_flux"])))
    print "Complete. Building plots..."
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
    #print "centroid list:", centroids
    
    try:
        # state the uncertainty of the duration as what is found by MAST. 
        durationErr1 = round(koi.koi_duration_err1,1)
        if type(durationErr1) =="NoneType":
            print "none type for durration error 1"
            durationErr1 = str("?")
    except TypeError:
        # if the duration error is not given by MAST, take it as the square root of the duration ( in hours).
        durationErr1 = str("$\pm$ "+str(round(math.sqrt(koi.koi_duration),1)))

    
def PhaseFolderOfAwesomeness():
    global phases,data2
    # Obtain the phases with respect to some
    # reference point (in this case T0=first centroid time plus a half-period)
    phases = foldAt(time, period, T0=(koi.koi_time0bk+period/2))
    # Sort with respect to phase
    #     First, get the order of indices ...
    sortIndi = np.argsort(phases)
    # ... and, second, rearrange the arrays.
    phases = phases[sortIndi]
    data2 = data[sortIndi]
    
    # Plot the result
    fig = plt.figure(figsize=(21,16))    
    ax =  fig.add_axes([0.05, 0.1, 0.8, 0.85])
    
    plt.hlines(1.000000,0,1,'r',label="Normalized Flux")
    plt.title("Phase-Folded Light Curves for {}".format(koi))
    plt.xlabel("Phase Position ( 0 to 1)")
    plt.ylabel("Normalized Flux")
    
    plt.plot(phases, data2, 'bp')
    plt.ylim(ymax = 1.007, ymin = np.nanmin(data2))
    try:
        hours =period%int(period)*24
        minutes = hours%int(hours)*60
        seconds = minutes%int(minutes)*60
        txt3 =       str(str(int(period))+" d : "+str(int(hours))+' m: '+str(int(seconds))+' s')
    except ZeroDivisionError:
        txt3 = period

    myTimeZone = str(strftime("%z", gmtime()))
    
    charText = " CHARACTERISTICS \n\nPeriod: {} days \n        {}\n\nFirst Centroid: {}d\n\nDispostion: {}\n\nNext Predicted Transit (UTC):\n {}".format(round(period,3),txt3,round(koi.koi_time0bk,3),koi.koi_disposition,FindNextEclipse().iso )
    plt.figtext(0.865, .76, charText, color='black', backgroundcolor='white',
           weight='roman', size=14)


    plt.legend()
    plt.show()

    
def RegularLC_Plot():
    fig = plt.figure(figsize=(21,16))    
    ax =  fig.add_axes([0.05, 0.1, 0.8, 0.85])    
    plt.ylim(ymax = 1.007, ymin = np.nanmin(data2))

    hours =period%int(period)*24
    minutes = hours%int(hours)*60
    seconds = minutes%int(minutes)*60
    txt3 =       str(str(int(period))+" d : "+str(int(hours))+' m: '+str(int(seconds))+' s')
    myTimeZone = str(strftime("%z", gmtime()))
    
    charText = " CHARACTERISTICS \n\nPeriod: {} days \n            {}\n\nFirst Centroid: {}d\n\nDispostion: {}\n\nNext Predicted Transit (UTC):\n {}".format(round(period,3),txt3,round(koi.koi_time0bk,3),koi.koi_disposition,FindNextEclipse().iso )
    plt.figtext(0.865, .76, charText, color='black', backgroundcolor='white',
           weight='roman', size=14)
    plt.plot(time,data,label ='Normalized Flux')
    # a horizontal line that shows the "1" value:
    plt.hlines(1.0000,min(time), max(time),'r', label="Base Flux")
    
    # Make it fancy!
    plt.title("Light Curves for {}".format(koi))
    plt.xlabel("BJD minus 2455200 (days)")
    plt.ylabel("Normalized Flux Intensity")
   

    plt.legend()

    plt.show()

    
# bin the data and conduct data analysis
def BinTheData():
    # Generate some data
    plt.figure(figsize=(24,14))    # REALLY big figure size  :)        
    x = phases
    y = data2
    
    # Bin using fixed number of bins and start at x0 = -10.
    # Use beginning of bin as starting value.
    r1, dt1 = binningx0dt(x, y, nbins=100, x0=0, useBinCenter=False)
    # Use fixed bin width. Specify another (wrong) error estimate and
    # use bin center.
    r2, dt2 = binningx0dt(x, y,  dt=dt1, \
                  x0=0, useBinCenter=True, removeNoError=True)
    #print "dt1, dt2: ", dt1, dt2
    #print "Input data points in last bin: ", r2[::,3]
    
    # Use the reducedBy flag to indicate the binning. In this case, x0
    # will be set to the lowest x value in the data, and the number of
    # bins will be calculated as: int(round(len(x)/float(reduceBy))).
    # Here, we will, thus, obtain 100 bins.
    r3, dt3 = binningx0dt(x, y, \
                        useBinCenter=True, removeNoError=True, reduceBy=10)

    print "dt3: ", dt3
    print "Number of bins in third version: ", len(r3[::,0])
    
    
    # Plot the output
    plt.plot(x,y)
    plt.errorbar(r1[::,0], r1[::,1], yerr=r1[::,2], fmt='kp--')
    plt.errorbar(r2[::,0], r2[::,1], yerr=r2[::,2], fmt='rp--')
    plt.errorbar(r3[::,0], r3[::,1], yerr=r3[::,2], fmt='gp--')
    plt.show()
            

    
    ###############################################
def FindNextEclipse():
    global star, nextCalculatedCentroid
    # kepler time is the 'then' BJD minus 2454833.0 days. This allows the kepler computers to store only
    # a few floating point numbers for the dates, eg: the first centroid of KOI5059.01 is 
    # given as 234.15043700000001. The full BJD of this occurance is the given date PLUS the offset of 
    # 2454833.0 JD = BJD:2455067.150437 JD  . So we add the offset to the given date stamp and get the 
    # full BJD using the code below. Now, if this time object is less than today's julian date, we can add
    # a period's time length of the koi and get the next date in terrestrial time. We do this process
    # iteratively until we find the JD in TT that most recently occured, and the one occurance after - the 
    # next transit in TT, in the UTC timezone. 
    ''' still need the date (iso format) of the 'most recent' transit, TT, in UTC timezone '''
    BJDoffset = 2454833.0
    todayJD = Time.now().jd
    
    # the last centroid in the kepler data:
    lastObservedCentroid = BJDoffset + centroids[-1]
    
    # calculate the next observable transit in the kepler data:
    nextCalculatedCentroid= Time(lastObservedCentroid,format = 'jd',scale = 'tdb' )
    
    while nextCalculatedCentroid.jd < Time.now().jd:
        nextCalculatedCentroid = Time(nextCalculatedCentroid.jd+period,format= 'jd', scale = 'tdb')
    
    
    star = koi.star
    #print  "\n\nNext eclipse of ", koi, "is on OR AROUND :\n",nextCalculatedCentroid.iso
    #print "\n**** LOCATION **** \n","DEC: ",star.kic_dec,"deg. \nRA:  ",star.kic_ra," deg."
    if __name__ == "__main__":
        return nextCalculatedCentroid

    ###############################################

def main():
    EstablishKOI()
    Normalize_LightCurves()
    PhaseFolderOfAwesomeness()
    RegularLC_Plot()
    #BinTheData()
    FindNextEclipse()
    print  "\n\nNext eclipse of ", koi, "is on OR AROUND :\n",nextCalculatedCentroid.iso
    print "\n**** LOCATION **** \n","DEC: ",star.kic_dec,"deg. \nRA:  ",star.kic_ra," deg."
    print "Disposition: ", koi.koi_disposition

# with the following line, functions can be imported from this file without auto-running on import
if __name__ == "__main__":
    main()
