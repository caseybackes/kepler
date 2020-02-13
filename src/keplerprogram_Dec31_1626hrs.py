#### IF THIS IS THE FIRST TIME RUNNING THIS PROGRAM: You must pip install the following libraries:
## kplr, matplotlib, scipy, math, numpy, pylab, pyfits. Some of these take a significant amount of storage space

import kplr 
import matplotlib.pyplot as plt
import scipy.stats as st
import math
import numpy as np
import os
global os


####################################################
#       Define some functions
####################################################

def PlotLC(startTime, endTime):
    '''arguments: startTime - time to begin the plot ( x-axis)
    endTime - the upper bound of the x-axis to plot.
    If you state endTime as "nlines", it will default to the end 
    of the light curves data.'''
    if startTime not in list(time):
        if startTime < time[0]:
            startTime = time[0]
            
    startTimeIndex = findIndex(time, startTime)
    endTimeIndex   = findIndex(time, endTime)
    midpoint = (endTime - startTime) / 2.0
    midTimeIndex   = findIndex(time, midpoint)
    
    print "Plotted From time = ",round(time[startTimeIndex],0), " to ",round(time[endTimeIndex],0)
    plt.figure()
    plt.hlines(Ibase, time[startTimeIndex],time[endTimeIndex], 'r', label = "Host Star's Base Intensity")
    plt.plot(time[startTimeIndex:endTimeIndex], data[startTimeIndex:endTimeIndex], "b", label = "Normalized Flux for KOI {}".format(str(usercall)))
    plt.title("Normalized KST Light Curves for KOI {}".format(str(usercall)))
    plt.text(time[midTimeIndex], max(data)*1.1, r'Status:{}'.format(koi.koi_disposition))
    plt.xlabel("BJD minus 2455200 (days)   Period: %s Terestrial Days"%(round(period,2)))
    plt.ylabel("Normalized Flux Intensity")
    plt.legend()
    plt.show()

######################
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
            return index
            break
#####################
def NanToIbase():
    '''**This function is run with no arguments. 
    
    Later, when we do the Chi Square analysis, we will measure the difference between the KOI's data points 
    and the model's data points. If there are any 'nan' values, we can eliminate their influence on the 
    Chi Square measurement by redefining the 'nan' values to be equal to the base intensity. This way 
    the Chi sqr measurement at that particular point is measured to be exactly zero ( the difference 
    between the theory curve and the actual data at any and all given points), and will therefore 
    have no influence on the Chi sqr measurement. '''
    for point in range(len(data)): 
        # if the data at any index turns out to be a 'nan', 
        if np.isnan(data[point]) == True: 
            #redefine that data point to be equal to the base intensity for data and fluxmean2 arrays
            data[point]       = Ibase
            fluxmean2[point]  = Ibase
    print "SubRoutine: NanToIbase ...complete"
######################
def EstablishCompatableKOI():
    '''** No arguments, please. The function will issue a prompt for user's choice of KOI.'''
    # declare that the koi established in this function be used as a global variable 
    global koi
    global usercall
    invalidAttempts = 0 
    totalPlanetsInSystem = 0 
    selectableList= []
    KOIlist = []
    usercall = int(input("Enter the KOI system to lookup in MAST: "))
    client = kplr.API()
    
    #Determine how many KOI are in the system associated with the user's choice
    print" KOI#     KOI Name         Period         Usable?"
    print"_"*55
    for PlanetInSystem in range(1,15,1):

        try: 
            KOI = client.koi(float(usercall+.01*PlanetInSystem))
            totalPlanetsInSystem+=1
            KOIlist = np.append(KOIlist, KOI)
            if KOI.koi_period < 15.4: 
                selectableList = np.append(selectableList, PlanetInSystem)
                print "(",PlanetInSystem,") ", KOI, round(KOI.koi_period,2), "   usable KOI"
            else: 
                invalidAttempts+=1
                print "(",PlanetInSystem,") ", KOI, round(KOI.koi_period,2), "   non-usable KOI"
        except ValueError:
            break
            
    print totalPlanetsInSystem, " total KOI in system"
    print totalPlanetsInSystem-invalidAttempts, " compatable KOI"
    
    
    #while the system has no compatable KOI for this program, require another system from the user.
    while totalPlanetsInSystem == invalidAttempts:
        usercall = int(input("Sorry, non are compatable, enter another system: "))
        # reset the total count and invalid count
        totalPlanetsInSystem = 0 
        invalidAttempts = 0 
        print" KOI#     KOI Name    Period     Usable?"
        print"_"*55
        for PlanetInSystem in range(1,15,1):
            
            try: 
                KOI = client.koi(float(usercall+.01*PlanetInSystem))
                totalPlanetsInSystem+=1
                KOIlist = np.append(KOIlist, KOI) 
                if KOI.koi_period < 15.4: 
                    selectableList = np.append(selectableList, PlanetInSystem)
                    print "(",PlanetInSystem,") ", KOI, round(KOI.koi_period,2), "   usable KOI"
                else: 
                    invalidAttempts+=1
                    print "(",PlanetInSystem,") ", KOI, round(KOI.koi_period,2), "   non-usable KOI"
            except ValueError:
                break
                
        print totalPlanetsInSystem, " total KOI in system"
        print totalPlanetsInSystem-invalidAttempts, " compatable KOI"
   
    if len(selectableList) > 1:
        KOIChoice = int(input("Choose from "+str(selectableList)+" : "))
        print "that is ", KOIlist[KOIChoice-1]
    else:
        KOIChoice = 1
    # now we have  
    #print"          PROPERTIES:"
    koi = KOIlist[KOIChoice-1]
######################    
def EstablishKOI():
    '''This function will let you choose what KOI in your chosen system to investigate.
    '''
    global koi
    global usercall
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
               
    print '\n',totalPlanetsInSystem, " total KOI in system"
    KOIChoice = int(input("Choose a specific KOI from above (enter Item #): "))
    while KOIChoice > totalPlanetsInSystem :
        print "That is not a valid option"
        KOIChoice = int(input("Choose a specific KOI from above (enter Item #): "))
    else:
        
        koi = client.koi(usercall+KOIChoice*.01)
        print "That is ", koi.kepoi_name
   
######################
def Charateristics_of_KOI():
    '''** No arguments, please. Function is fully automated.
    
    Here we are establishing (from MAST) some basic parameters of the 
    Stellar system and the KOI. '''
    # Declare global variables: name, 
    global name, HostStar, period, sma, sma_AU, AU, STeff
    
    
    AU = 149597870.7 # in km
    name = koi.kepler_name
    HostStar = koi.kepid
    period = koi.koi_period
    sma = koi.koi_sma * AU
    sma_AU = sma/AU
    STeff = koi.koi_steff    
    
    print "Planet: ", name
    print "Host Star: ", HostStar
    print "Period in days: ", round(period, 4), ' days'
    print "Period in years: ", round(period/365.4, 4),' years'
    print "Semi-major axis in AU: ", sma, ' AU'
    print "Effective Stellar Temperature (K): ", STeff
####################################################
            ##  IN PROGRESS ## 
    # This function would search the MAST database for KOI that match user-specified criteria
    
def FindKOI(period_days = 365, planet_eq_temp = 6000, semi_major_axis= 1, disposition = "candidate" ):
    # filter and display info for koi's of specific criteria:
    # period in days
    # effective planetary temperature in kelvin
    # semi major axis, in AU
    # by disposition ( confirmed, candidate, false positive)
    per = str("koi_period>"+period_days)
    pteff = planet_eq_temp
    sm_axis = semi_major_axis
    status = disposition
    temp = []
    sma = []
    dispo = []
    client = kplr.API()    
    kois = client.kois(where=str(per), sort=("koi_period", -1))
    #for i in kois:
        #temp = np.append(temp,
    return kois
    print "Not Complete yet"
####################################################
def Normalize_LightCurves():
    '''                 
    In this function, we normalize the flux of the lightcurves to where the y axis is "normalized" to the value of 1. 
    The lightcurve data is downloaded to the current directory and the path I found my downloaded light
    curves in was "C:\Users\Casey\.kplr\data\lightcurves" and I found folders for each object I have investigated
    and in each folder there were about 16 files - the lightcurves. With short cadence, this is 
    a SHIT TON of data! Watch, download a LC with short cadence and long cadence. Then look at the 
    size of each file, and you'll see the difference. Not only that but the amount of time it takes 
    to download the short cadence data is FUCKING retarded! You'll grow old. Trust me. 
    Be sure of your choice if you want short cadence data!  
    
    
    In this function, we normalize the flux of the light curves to the value of 1 and plot it on the y axis as 
    a time series. The light curve data is downloaded to the current working directory, a folder for each KOI 
    investigated.
    
    Running this function will invoke user input for the "cadence" and the "fetch" options. Selecting short
    candence light curve data will download the raw, un-normalized light curve files to the current working
    directory. Short cadence files will tend to be very data-heavy and take several moments to download and process
    in this program.   
    
    '''
    # prompt for user input parameters for the acquisition command to get LC data. 
    cadence = str(raw_input("Choose (S)hort cadence or (L)ong cadence data : "))
    if "l" in cadence or "L" in cadence:
        cadence = False # set short cadence option to False
    elif "s" in cadence or "S" in cadence:
        conf = str(raw_input("Are you sure you want short cadence data? This may take a while. (y)es / (n)o:"))
        if "y" in conf or "Y" in conf or "yes" in conf or "Yes" in conf or "YES" in conf or "just fucking do it" in conf:
            print "Downloading short cadence data for light curve data, but this WILL take a few moments to download.","\n"
            cadence = True
    # prompt for saving light curves to the current directory
    fetch   = str(raw_input("Once downloaded , save light curve data to current working directory? (Y)es / (N)o: "))
    if "Y" in fetch or "y" in fetch or "yes" in fetch or "Yes" in fetch or "YES" in fetch: 
        ftch = True; print"Ok, saving downloaded data to directory: {}".format(os.getcwd())
    else:
        ftch = False 
        
    # based on user input, declare the operation to acquire the light curve data. 
    print"Working..."
    lcs = koi.get_light_curves(short_cadence = cadence, fetch = ftch, clobber = True)

    # Initiate empty arrays into which we append data from MAST
    # time array - used extensively for data plotting and model plotting
    # flux array - used to build the fluxmean, fluxmean2 arrays
    # fluxmean(2) arrays - used for data plotting and normalization. 
    # data array - a copy of the fluxmean2 array used for plotting lc data. 
    global time, flux, fluxmean, fluxmean2, data, nlines, fmax, fmin, tmin, tmax
    
    time,flux,fluxmean,fluxmean2,data =[],[],[],[],[]
    for lc in lcs:
        with lc.open() as f:
            hdu_data=f[1].data
            time=np.append(time,hdu_data["time"])
            flux=np.append(flux,(hdu_data["pdcsap_flux"]))
            fluxmean=np.append(fluxmean, (hdu_data['pdcsap_flux'])/(np.nanmean(hdu_data["pdcsap_flux"])))
            # Fluxmean2 is for finding the base intensity of the host star. 
            fluxmean2=np.append(fluxmean2, (hdu_data['pdcsap_flux'])/(np.nanmean(hdu_data["pdcsap_flux"])))
            data    = np.append(data,      (hdu_data['pdcsap_flux'])/(np.nanmean(hdu_data["pdcsap_flux"])))
    wait =input("Press enter to continue")        
    #For reference later in other functions, passed to main function through global declaration
    nlines = len(fluxmean2)
    fmax=np.nanmax(flux)
    fmin=np.nanmin(flux)
    tmin=np.nanmin(time)
    tmax=np.nanmax(time)
    
    print "Complete."
    print "Generating plot of normalized, preconditioned data."
    plt.plot(time, data)
    plt.title("Normalized KST Light Curves for KOI {}".format(koi))
    plt.xlabel("BJD minus 2455200 (days)   Period: %s Terestrial Days"%(round(koi.koi_period,2)))
    plt.ylabel("Normalized Flux Intensity")
    plt.legend()
    plt.show()

    #print"Complete. You now have access to the arrays:"
    #print"time, flux, fluxmean, fluxmean2, data"
####################################################
def Base_Stellar_Intensity():
    '''No arguments, please. Function is fully automated.'''
    global Ibase
    # Establish the signal from the star alone, while the KOI is not in transit. This 
    # will be the star's baseline intensity - without a KOI blocking any of its light. 
    # We do this by removing the the 2-sigma and higher data points from the mean. 
    # doing this twice is sufficient to remove the great majority of data that doesnt 
    # reflect the most reasonalbly accurate baseline intensity of the host star. 
    # First, we get the initial "average" intensity , or flux, of the hos star...
    Ibase=np.nanmean(fluxmean2[0:nlines])
    sigma= np.nanstd(fluxmean2[0:nlines])  # define one sigma deviation, without the 'nan' data points. 
    # Now remove the 2 sigma data points by redefining them to be equal to the base intensity
    # and take another average. Do this twice (the for loop) to get a really accurate base intensity. 
    for _ in range(5):
        for i in range(0,nlines):
            if fluxmean2[i] <(Ibase -(2*sigma) ) or fluxmean2[i] > (Ibase + (2*sigma)):
                fluxmean2[i]=Ibase
        sigma=np.nanstd(fluxmean2[0:nlines])
        Ibase=np.nanmean(fluxmean2[0:nlines])
    print "One NORMALIZED standard deviation = ", sigma
    print "Base Stellar Intensity, 'Ibase' = ",Ibase
####################################################
def GaussianModel( int ):
    global duration, d, ic, w, centroids, tunit, phalf,gmaster
    ModelLength = int
    #some necessary primary definitions#
    duration=koi.koi_duration/24.  # used in the width parameter, width=duration/2 .   actually more like 2.3, also, this is given by MAST in units of hours, so divide by 24 to get back into (used) units of days 
    d=koi.koi_depth/10e+5    # d parameter in the gaussian function
    ic=koi.koi_time0bk   # first observed centroid... MAST already has this 'figured out'... more or less. 
    w=duration/2.3     # w parameter in the gaussian function. 
    centroids=[]       # an array of the time location of each centroid (empty right now, but next step fills the array. :)

    totalTransitWindow=tmax-ic   # final observing time minus first transit gives total time that transits occur in
    expectedTransits=totalTransitWindow/period  # span of time to see transits divided by period gives expected number of periods...approximately...more or less. 
    
    #create an array with the locations of all centroids, starting from the first observed transit.
    for i in np.arange(0,expectedTransits):
        c=ic+i*period
        centroids=np.append(centroids, c)   # this is where all the centriods should lay on the plot. 
    
    # The time between data recordings ( the "candence" from earlier) declared as "tunit"
    tunit=time[5]-time[4]  # just geting the delta t here, i think it was about three minutes between readings on the telescope, but the indecies of 4 and 5 are completely arbitrary. 
    
    phalf=period/2.0
      
    #remove any (artificial or authentic) centroids before time_0 of the lightcurve data
    for cent in centroids:
        if cent < time[0]:
            centroids = np.delete(centroids,0)
        if cent == time[0] or cent > time[0]:
            break
                        
    #************************************************************************
    #    The Mathematical Transit Model: The Inverted Gaussian "Bell Curve":
    #   
    #      GaussModel=Ibase-Ibase*d*math.exp(-(theoryFlux-Centroid)**2/w**2)
    #************************************************************************
    
    # initialize the array of theoretical flux values for the model
    gmaster =[]
    
    # The first section of the model extends a stright line from time[0] 
    # to half a period before the first centroid (transit center) 
    initialTimeRange = np.arange(time[0], centroids[0] - phalf , tunit) 
    
    # Sometimes the first transit occurs immediately after the initial recording time, 
    # and therefore the data doesnt start at the "half-period before first centroid" mark, 
    # but starts afterwards by a small amount
    if len(initialTimeRange) == 0 : 
        pass
    else:    
        # use the initial time range to build identical sections of the Gaussian Model
        for timepoint in initialTimeRange: 
            GaussModel = Ibase 
            gmaster = np.append(gmaster, GaussModel)
    
    #establish first theoretical event
    timeRange = np.arange(time[0], centroids[0] + phalf, tunit) 
    for timepoint in timeRange:
        GaussModel = Ibase - Ibase*d*math.exp(-(timepoint - centroids[0])**2/w**2)
        gmaster = np.append(gmaster, GaussModel)
    while len(gmaster) > len(data[0:findIndex(time, centroids[0]+phalf)]):
        gmaster = np.delete(gmaster, 0)
    #print"Length of first theory section: ", len(gmaster)
    #print"Length of data from t_0 to half past first event", len(data[0:findIndex(time, centroids[0]+phalf)])
    
    # This gives us a general time range to use, does not need to be centriod specific, 
    # any centroid will be fine timeRange = np.arange( centroids[0]-phalf +tunit, centroids[0]+phalf, tunit)
    timeRange = np.arange(centroids[0]-phalf + tunit, centroids[0] + phalf , tunit) 
    
    for _ in range(1,ModelLength):
        for timepoint in timeRange: 
            GaussModel=Ibase-Ibase*d*math.exp(-(timepoint-centroids[0])**2/w**2)  
            gmaster = np.append(gmaster, GaussModel)
    plt.figure()
    plt.plot(time[0:len(gmaster)], data[0:len(gmaster)], "b")
    plt.plot(time[0:len(gmaster)], gmaster, 'r')
    #plt.plot(gmaster, "r")
    #plt.plot(time[0:len(gmaster)], gmaster,"r")
    #plt.plot(data[0:len(gmaster)], "b")
    plt.show()
####################################################
def Chi_Square_Gaussian():
    '''No arugments needed. Run this sub after GaussianModel function'''
    global X2,SuperData, SuperModel, FluxAverage
    
    print '\n               Chi Square Analysis\n                 "Gaussian Model"            '
    print 'This is the chi squared measurement of how well the model fits the data.\nThis is given as:'
    print '    X2= sum (  (observed - model)^2/ model   ).\nWhere "observed" is the normalized flux data multiplied by the mean of the flux, and "expected" is the\nvalue of the model at each point.\nIf this value is less then the length of the model, the model is a good fit.\nIf X2 > len(gmaster), then the model isnt a good fit, and the model doesnt represent this planet.'
    print '\nLength of model array: ',len(gmaster),'elements.'
    
    # rescale the values in the data and gmaster arrays by average of original mean of flux to get nearly the 
    # original flux values back, applying the model to the "life size" flux data. 
    FluxAverage= np.nanmean(flux[0:len(gmaster)])
    print"Flux Average: ",FluxAverage
    SuperData = data*FluxAverage
    SuperModel = gmaster * FluxAverage
    # Initialize the Chi Square to zero, adding the ratio of (1) squared differences between model and data to (2) the model
    X2 =0
    for point in range(len(SuperModel)):
        x2_point = (SuperData[point] - SuperModel[point] )**2 / (SuperModel[point])
        x2_point = round(x2_point,5)
        X2 = X2+x2_point
        
        #print "( ", SuperData[point], " - ", SuperModel[point], ")^2 / ", SuperModel[point], " = ", x2_point
        #print x2_point
    print "Chi square = ",X2,"with ",len(SuperModel)," elements." 
    
    
    # Initialize the Chi Square to be zero
    X2_norm = 0
    for point in range(len(gmaster)):
        x2_point = (data[point] - gmaster[point])**2 / (gmaster[point])
        #x2_point = round(x2_point,5)
        X2_norm = X2_norm + x2_point
        #print x2_point
    print "Chi square of normalized data and model = ",X2_norm, " with ",len(gmaster)," elements."
        
    # Plot the overlay of the SuperData with the SuperModel and annotate X2 on the plot in the title 
    plt.plot(time[0:len(SuperModel)],SuperData[0:len(SuperModel)],'b',label='KOI Light Curve Data')
    xgmaster=np.arange(time[0],time[0]+len(SuperModel)*tunit,tunit)
    plt.plot(xgmaster,SuperModel,'r', label='Gaussian Model')
    plt.legend()
    plt.xlabel('Time (days)')
    plt.ylabel('Un-Normalized Flux Intensity')
    plt.title(r"Gaussian Model Fitting to {}, $\chi^2$ = {}".format(name,round(X2,3)))
    plt.show()
    
    
    # Plot the overaly of the normalized data and model and annotate X2 on the plot in the title
    plt.figure()
    plt.plot(time[0:len(gmaster)], data[0:len(gmaster)], 'b', label = 'KOI Light Curve Data')
    xgmaster = np.arange(time[0], time[0]+len(gmaster)*tunit, tunit)
    plt.plot(xgmaster, gmaster, 'r', label = 'Gaussian Model')
    plt.legend()
    plt.xlabel('Time (days)')
    plt.ylabel('Normalized Flux Intensity')
    plt.title(r"Gaussian Model Fitting to {}, $\chi^2$ = {}".format(name,round(X2_norm,8)))
    plt.show()
    
####################################################
def BoxModel(B_limit):
    global boxmaster, modeltimerange
    #new figure
    plt.figure()

    # the time values for the boxmodel's duration
    modeltimerange=[]     
    # the theory curve values
    boxmaster=[]           
    counter = 0 # counts how many times the model is built from for-loop below
    # Construct the model for eight transits
    for ic in centroids[0:B_limit]:
        x1 = np.arange(time[0], ic-w, tunit)
        modeltimerange=np.append(modeltimerange,x1)
        x2 = np.arange(max(x1)+tunit, ic+w, tunit)
        modeltimerange=np.append(modeltimerange,x2)
        x3 = np.arange(max(x2)+tunit, ic+phalf, tunit)
        modeltimerange=np.append(modeltimerange,x3)
        
        y1=np.ones(len(x1))*Ibase 
        boxmaster=np.append(boxmaster,y1)
        y2=np.ones(len(x2))*(Ibase-d)
        boxmaster=np.append(boxmaster,y2)
        y3=np.ones(len(x3))*Ibase
        boxmaster=np.append(boxmaster,y3)
        counter= counter+1
        plt.plot(x1,y1,'r'),plt.plot(x2,y2,'r')
        plt.plot(x3,y3,'r',)
        plt.vlines(max(x1),Ibase-d,Ibase,'r')
        if counter == B_limit:
            plt.vlines(max(x2), Ibase-d,Ibase,'r', label = "Box Model")
        else:
            plt.vlines(max(x2), Ibase-d,Ibase,'r')
        del x1,x2,x3,y1,y2,y3
        
        
    plt.plot(time[0:len(boxmaster)],data[0:len(boxmaster)],'b', label= "KOI Data (Normalized)")
    plt.legend()
    plt.title(r'Box Model Algorithm over KOI {}'.format(name))
    plt.xlabel('Time (days)')
    plt.ylabel('Normalized Flux of Host Star')
    
        
    plt.legend(),plt.show()
####################################################
    
####################################################
def Chi_Square_Box():
    print '\n               Chi Square Analysis\n                 "Box Model"            '
    print 'This is the chi squared measurement of how well the box model fits the data.'
    print 'The chi squared measurement of how well the model fits the data is given as:'
    print 'X2= sum (  (observed - model)^2/ model   ) for each measurement. If this value is close to zero, the model is a good fit.\nIf X2 > len(boxmaster), then the model isnt a good fit, and the data isnt likely of a planet transit.'
    # the cboxhai square mesurement needs to not count 'nan' values
    X2=np.nansum ( ((fluxmean[0:len(boxmaster)]-boxmaster)**2.)/(boxmaster))
    print 'The X2 measurement for the box model over eight transits is ', round(X2,6)
    FluxAverage = np.nanmean(flux[0:len(boxmaster)])
    print 'The Chi Square before rescaling by average flux: X2 = ', round(X2,3)
    Rescaled_X2=X2*FluxAverage
    print 'X2 x (mean of flux) = ', X2,' x ' ,np.nanmean(flux),' is the properly rescaled Chi Square value.'
    print 'Bins = ', len(boxmaster)
    print 'Reduced Chi Square = X2/ d, d=n-c, n=len(gmaster), c= 3 (width, depth, and centroid)'
    Reduced_X2=Rescaled_X2/(len(boxmaster)-3)
    print 'Reduced Chi Square of rescaled Chi Square =', Reduced_X2
    print 'A reduced Chi Square is ideal when on the order of or at a value of 1. '
    print 'Original flux average for Q1: ', round(np.nanmean(flux[0:len(gmaster)]),2),'electrons/sec' 
    # Begin the plot for the box model
    plt.title(r"Box Model Fitting to {}, Reduced $\chi^2$ = {}".format(name,Reduced_X2,3))
    plt.xlabel('Time (days)')
    plt.ylabel('Normalized Flux Intensity')
    plt.show()
    print 'Degress of freedom equals (# bins)-parameters. \nParameters used are depth, width, and centroid, totals 3 parameters.\nDegrees of freedom   =  ', len(boxmaster),' -  3 =  ', len(boxmaster)-3,'.'
    print 'The measurement of the reduced chi-squared measurement is chi-squared divided by degrees of freedom.'
    print 'X2/d = ',X2/(len(boxmaster)-3)
####################################################
'''
def FindSing():
    # If any data point is greater than two sigma from the 
    for datapoint in data:
        if datapoint > 1.002:
            datapoint = Ibase
    # now establish one standard deviation of the mean of all data points
    oneSTD = np.std(data)

    for datapoint in data:
        break
'''
####################################################
def PTeff():
    global Stefan_Bolotzman_Constant, SRad, L, errL, PTeff_Bond
    
    ############################################
    ##### Plot the possible values for     #####
    ##### the effective planetary temps    #####
    ##### given different albedo values    #####
    #####           (really cool!!)        #####
    ############################################
    
    # Establish the luminosity of the parent star
    Stefan_Boltzman_Constant=5.6704e-08
    # call for the stellar radius from MAST and convert from Solar Masses to kg
    SRad=koi.koi_srad*695500# km    Stellar Radius in kilometers
    # Calculate luminosity 
    L=4*math.pi*SRad**2*Stefan_Boltzman_Constant*STeff**4 
    
    # Establish mass of the parent star
    Mstar_sm=koi.koi_smass # mass of parent star in solar masses. 
    Mstar_kg=Mstar_sm*1.9891e+30 # mass of parent star converted to kg.
    Mstar_kg = round(Mstar_kg, 3) 
    
    
    plt.figure()
    # Calculate the effective planetary temperature based on varying albedo
    A=np.arange(0,1,.01) #where the .01 is the resolution of the curve.
    PTeff=STeff*((1-A)/4)**(.25)*(SRad/(sma))**.5; #print PTeff,"K :Effective Planetary Temp, Albedo of",A
    PTeff_Bond=STeff*((1-.3)/4)**(.25)*(SRad/(sma))**.5

    plt.plot(A,PTeff,label="Theoretical Temperature");#plot shows how planetary eq temp relates to albedo
    plt.axhline(y=PTeff_Bond,xmin=0,xmax=.3,color='r')
    plt.axvline(x=.3, ymin=0, ymax=.843, color='r',label="Temp at BOND Albedo")
    plt.xlabel('Varying Values of Planetary Albedo')
    plt.ylabel('Effective Planetary Temperature (K)')
    plt.title("Varying Teff for Kepler Object of Interest (KOI) %s"%(name))
    plt.legend()
    plt.show();#plt.figure()#
    
    
    print "Mass of parent star in solar masses: {} SM".format(Mstar_sm)
    print "Mass of parent star in kg:",Mstar_kg," kg"
    print "Luminosity of parent star: {} W/m^2".format(L)
    
    #####################################################
    # planetary effective temperature = 
    #     stellar temperature x ( (1- albedo) / 4 )^(1/4) x sqrt(stellar radius / semi-major axis)
    print '\n'*2
    print "     Effective Planetary Temperature as a function of albedo",'\n',""*15,"(see generated plot,'Verying Teff for Kepler Object of Interest')",'\n'
    for A in range(0,11,1):
        if A!=3:
            # calculate effective planetary temp
            PTeff=STeff*((1-A*.1)/4)**(.25)*(SRad/(sma))**.5; 
            PTeffRounded = int(round(PTeff,-1));print PTeffRounded,"K :Effective Planetary Temp, Albedo of",A
            # print "effective planetary temp at (current albedo)
        else:
            # calculate planetary temperature anyway, but identify this as the temperature at the "bond albedo" - the albedo which is similar to Earth. 
            print" "*22,"***"
            PTeff=STeff*((1-A*.1)/4)**(.25)*(SRad/(sma))**.5
            PTeffRounded = int(round(PTeff,-1));print PTeffRounded,"K :Effective Planetary Temp, BOND Albedo of",A
            print" "*22,"***"
    
    # Automated analysis of Effective Planetary temperature as a function of Albedo:
    print'\n'
    print "According to the Kepler Database, {} has an equilibruim temperature of {} Kelvin. With this program I have a calculated difference of about {}%, likely due to rounding differences between this program and the Kelpler Science Team's methods.".format(name,koi.koi_teq,math.fabs(round(100-100*koi.koi_teq/PTeff_Bond,2)))
    if 273<=PTeff_Bond<=373.15:
        print "Based on the calculated equilibrium temperature of the planet with an estimaed albedo similar to Earth's, THIS COULD BE A HABITABLE PLANET!!!"
    elif PTeff_Bond<273:
        print "Based on an equilibrium temperature of {} Kelvin, this is not a habitable planet because the temperature is to low for liquid water to exist.".format(round(PTeff_Bond),2)
    elif PTeff_Bond>373.15:
        print "Based on the equilibrium temperature of {} Kelvin, this is not a habitable planet because the temperature is to high for liquid water to exist.".format(round(PTeff_Bond),2)
####################################################

# main program
def KeplerProgram():
    # Execute function that establishes the user's choice of KOI and the link with MAST. 
    EstablishKOI()
    ####################################################
    # Execute funtion that establishes the KOI's name, period, host star designation, 
    # semi-major axis in meters and AU, Effective stellar temperature of host star
    Charateristics_of_KOI()

    
    ###########################################################
    ##   NORMALIZE THE LIGHT CURVES!!! ########################  this was seriously a real bastard to figure out how to do. 
    ###########################################################
    Normalize_LightCurves()

    # For reference in analysis
    
    
    # Establish the base stellar flux without the transits. 
    Base_Stellar_Intensity()

    # Replace in the 'data' array any flux values that are 'nan's with the value of base stellar intensity, "Ibase"
    NanToIbase()
    
    # Can only allow a few transits in the Gaussian model:
    maxTransits = int((130-time[0])/period)+1
    print "index in time array where t=130 occurs: ",findIndex(time, 130), 
    sentence = str("Number of transits, up to ")
    sentence2 = sentence + str(maxTransits) +", " + str("to cover in the Gaussian Model: ")
    try:
        G_limit = int(raw_input(sentence2))
        while G_limit > int(maxTransits):
            print" Please read the instructions above."
            G_limit = int(raw_input(sentence2))
    except ValueError as shitHitTheFan:
        print"...um, I think something went wrong there..."
        G_limit = int(raw_input(sentence2))
    
    
    
        
    GaussianModel(G_limit)
    Chi_Square_Gaussian()
    
    B_limit = int(raw_input("Number of transits to cover in Box model: "))  
    BoxModel(B_limit)
    
    #Chi_Square_Box()
    
    PTeff()
#####################################################################
# run main KeplerProgram
#KeplerProgram()
    



     
'''
############################################
## Chi - Square Analysis of Box Model ##
############################################
print '\n               Chi Square Analysis\n                 "Box Model"            '
print 'This is the chi squared measurement of how well the box model fits the data.'
print 'The chi squared measurement of how well the model fits the data is given as:'
print 'X2= sum (  (observed - model)^2/ model   ) for each measurement. If this value is close to zero, the model is a good fit.\nIf X2 > len(boxmaster), then the model isnt a good fit, and the data isnt likely of a planet transit.'
# the Chi square mesurement needs to not count 'nan' values
X2=np.nansum ( ((fluxmean[0:len(boxmaster)]-boxmaster)**2.)/(boxmaster))
print 'The X2 measurement for the box model over eight transits is ', round(X2,6)

# Begin the plot for the box model
plt.title(r'Box Model Fitting to {}, X2={}'.format(koi,round(X2,6)))
plt.xlabel('Time (days)')
plt.ylabel('Normalized Flux Intensity')
plt.show()
print 'Degress of freedom equals (# bins)-parameters. \nParameters used are depth, width, and centroid, totals 3 parameters.\nDegrees of freedom   =  ', len(boxmaster),' -  3 =  ', len(boxmaster)-3,'.'
print 'The measurement of the reduced chi-squared measurement is chi-squared divided by degrees of freedom.'
print 'X2/d = ',X2/(len(boxmaster)-3)
############################################

'''