# Developed independantly by Casey Backes at the University of Colorado Boulder for the "Scientific Computing and Astrophysical Statistical Data Analysis" course with Dr. Webster Cash (2014)
# Also, this is my time learning computer programming and first programming project. 





#### IF THIS IS THE FIRST TIME RUNNING THIS PROGRAM: You must pip install the following libraries:
## kplr, matplotlib, scipy, math, numpy, pylab, pyfits. Some of these take a significant amount of storage space

## !!! If running in Oracle VM: Several plots will be produced (like 4 or 5) but you must close each one in order to allow the next to open. Oracle Sucks. Use Enthought's Canopy iPython instead. Its a thousand times better. 
import kplr # to install 'kplr', go to cmd prompt and type "pip install kplr"  - no quotes. 
import matplotlib.pyplot as plt
import scipy.stats as st
import math
import numpy as np
import os

####################################################
#       Define some functions

def PlotLC(startTime, endTime):
    '''arguments: startTime - time to begin the plot ( x-axis)
    endTime - the upper bound of the x-axis to plot.
    If you state endTime as "nlines", it will default to the end 
    of the light curves data.'''
    IndexForStartTime = list(time).index(startTime)
    start = time[IndexForStartTime]
    IndexForEndTime = list(time).index(endTime)
    end = time[IndexForEndTime]

    if endTime == "":
        endTime == nlines
    plt.hlines(Ibase, start, end, 'r', label = "Host Star's Base Intensity")
    plt.plot(time[0:nlines], data[0:nlines], "b", label = "Normalized Flux {}".format(str(usercall)))
    plt.title("Normalized KST Light Curves for KOI {}".format(str(usercall)))
    plt.text(151, (Ibase + 8*sigma), r'Status:{}'.format(koi.koi_disposition))
    plt.xlabel("BJD minus 2455200 (days)   Period: %s Terestrial Days"%(round(period,2)))
    plt.ylabel("Normalized Flux Intensity")

    
    '''
    plt.figure()
    plt.hlines(Ibase,start,end,'r',label='Normalized Stellar Flux')
    plt.plot(time[findIndex(time, startTime):findIndex(time, endTime)],data[findIndex(time, startTime):findIndex(time, endTime)]/Ibase,'b',label='KOI Light Curve Data')
    plt.title("Kepler Space Telescope Data")
    plt.text(151,(Ibase+4*sigma),r'Status:{}'.format(koi.koi_disposition))
    plt.xlabel("BJD - 2455200 (days)   Period: %s"%(round(period,2)))
    plt.ylabel("Normalized Flux Intensity")
    plt.title('KST Light Curves for KOI  %s'%(name)),plt.show() 
    '''
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
####################################################
# good objects to look at:
# 1,2,10,
# Here, we use 'kplr' third party library to interface with the MAST database.
usercall=int(input("Enter the number of koi you wish to lookup: "))
print"Analyzing object..."
print"This may take a few moments. "
usercall=float(int(usercall) + .01)
client=kplr.API()
koi =client.koi(usercall) #picking a kepler object of interest from the database, always includes the .01 at the end
############################################

# these two lines call the name of the koi as 
# designated by the MAST database
name=koi.kepoi_name
print "Name of KOI is:",name
print "Designation of Host Star:",koi.kepid
print ''
###############################

# here we fetch the period data from MAST database 
# to get the period in days of the koi
period=koi.koi_period
print "Period in terrestrial days:",round(period,3) #  self explanitory.
print 'Period in terrestrial years:',round(period/365,4) # self explanitorypl.
print ''# blank space in python shell
#######################################################

sma=(period/365.)**(2/3.)# semi major axis in AU as calculated by NVK3L
au=149597870700.#meters = 1 AU. This converts AU to meters
SMA=sma*au #converts the semi major axis to AU
print "Semi major axis in AU:",sma
print "Semi major axis in meters:", '%s'%float('%.3g'%SMA)
print ''
###########################################################

# here, we ask the MAST database for the surface temperature 
# of the parent star of the koi, and use it later 
# for calculating the effective ( or "equilibrium" ) planetary temperature 

STeff=koi.koi_steff
print "Effective Stellar Surface Temperature:",STeff,"Kelvin."
print ""# just prints a blank space in the python window

    
###########################################################
##   NORMALIZE THE LIGHT CURVES!!! ########################  this was seriously a real bastard to figure out how to do. 
###########################################################
#using koi xx, known as xx.01,
# we dont want the short cadence data, and the option "fetch=true/false" means to download the data or not.
#elect fetch=False if you dont want to download the light curves to your machine. 
## ** WARNING: Light Curve data is VERY costly in mahcine memory. The more objects you investigate, 
## **          the more the data stacks up in the machine's memory if you download the data. 
## **           Highly recommended that you dont download. 

lcs=koi.get_light_curves(short_cadence=False,fetch=False )
time,flux,fluxmean,fluxmean2,data=[],[],[],[],[]
for lc in lcs:
    with lc.open() as f:
        hdu_data=f[1].data
        time=np.append(time,hdu_data["time"])
        flux=np.append(flux,(hdu_data["pdcsap_flux"]))
        # calculate the averge line without the 3sigma y-values
        # plot red line at that average, call Io (I <not>) which just the intensity of the star itself. 
        # fluxmean = mean of values <= 3sigma
        fluxmean=np.append(fluxmean, (hdu_data['pdcsap_flux'])/(np.nanmean(hdu_data["pdcsap_flux"])))
        fluxmean2=np.append(fluxmean2, (hdu_data['pdcsap_flux'])/(np.nanmean(hdu_data["pdcsap_flux"])))
        data    = np.append(data,      (hdu_data['pdcsap_flux'])/(np.nanmean(hdu_data["pdcsap_flux"])))

# for analysis and data plotting: 
nlines=len(flux)
fmax=np.nanmax(flux)
fmin=np.nanmin(flux)
tmin=np.nanmin(time)
tmax=np.nanmax(time)


# Plot the signal from the star alone, while the KOI is not in transit. This 
# will be the star's baseline intensity - without a KOI blocking any of its light. 
# We do this by removing the the 2-sigma and higher data points from the mean. 
# doing this twice is sufficient to remove the great majority of data that doesnt 
# reflect the true baseline intensity of the star. 

Ibase=np.nanmean(fluxmean[0:nlines])
sigma= np.nanstd(fluxmean[0:nlines])  # define one sigma deviation, without the 'nan' data points. 
for _ in range(6):
    for i in range(0,nlines):
        if fluxmean[i] <(Ibase -(2*sigma) ) or fluxmean[i] > (Ibase - (2*sigma)):
            fluxmean[i]=Ibase
    sigma=np.nanstd(fluxmean[0:nlines])
    Ibase=np.nanmean(fluxmean[0:nlines])
#for i in range(0,nlines):
 #    if fluxmean[i] <(Ibase-(2*sigma) ):
  #      fluxmean[i]=Ibase     
#Ibase= np.nanmean(fluxmean[0:nlines])   
 
         #########################
# Re-value the 'nan' values to the same as the average
for point in range(len(data)): 
    if np.isnan(data[point]) == True: 
        data[point]       = Ibase
        fluxmean2[point]  = Ibase
        ##########################
# Plot the normalized lightcurves with function defined above
PlotLC(0,nlines)

#################################################################
#################################
# Model Fitting: Gaussian Model #
#################################


##################################
#some necessary primary definitions#
duration=koi.koi_duration/24.  # used in the width parameter, width=duration/2 .   actually more like 2.3, also, this is given by MAST in units of hours, so divide by 24 to get back into (used) units of days 
d=koi.koi_depth/10e+5    # d parameter in the gaussian function
ic=koi.koi_time0bk   # first observed centroid... MAST already has this 'figured out'... more or less. 
w=duration/2.3     # w parameter in the gaussian function. 
centroids=[]       # an array of the time location of each centroid (empty right now, but next step fills the array. :)

a=tmax-ic   # final observing time minus first transit gives total time that transits occur in
b=a/period  # span of time to see transits divided by period gives expected number of periods...approximately...more or less. 
#create an array with the locations of all centroids, starting from the first observed transit.
for i in np.arange(0,b):
    c=ic+i*period
    centroids=np.append(centroids, c)   # this is where all the centriods should lay on the plot. 
tunit=time[5]-time[4]  # just geting the delta t here, i think it was about three minutes between readings on the telescope, but the indecies of 4 and 5 are completely arbitrary. 


phalf=period/tunit/2.

############################
#DOWN AND DIRTY WITH PLOTS #
############################
#remove any (artificial or authentic) centroids before time_0 of the lightcurve data
for cent in centroids:
    if cent < time[0]:
        centroids = np.delete(centroids,0)
    if cent == time[0] or cent > time[0]:
        break
                    # THE MODEL
#GaussModel=Ibase-Ibase*d*math.exp(-(theoryFlux-Centroid)**2/w**2)


'''
create one section ( after the first centroid) 
gmaster = []

time range from cent - phalf +tunit to cent + phalf
########### ####**#### ####**
          # #        # # 
          # #        # #
          # #        # #
          ###        ###
          
'''

# initialize the array of theoretical flux values for the model
gmaster =[]

# The first section of the model extends a stright line from time[0] 
# to half a period before the first centroid (transit center) 
initialTimeRange = np.arange(time[0], centroids[0] - phalf , tunit) 

# Sometimes the first transit occurs immediately after the initial recording time, 
# and therefore the data doesnt start at the "half-period before first centroid" mark, 
# but starts afterwards by a small amount
if len(initialTimeRange) == 0 : 
    # then the first transit event occured right after the initial recording time. 
    print " len of initial Time range is zero"
    #initialTimeRange # = np.arange( time[0], centroids[0] - phalf , tunit) 
    
# use the initial time range to build identical sections of the Gaussian Model
for timepoint in initialTimeRange: 
    GaussModel = Ibase - Ibase*d*math.exp(-(timepoint - centroids[0])**2/w**2)
    gmaster = np.append(gmaster, GaussModel)


# This gives us a general time range to use, does not need to be centriod specific, 
# any centroid will be fine timeRange = np.arange( centroids[0]-phalf +tunit, centroids[0]+phalf, tunit)
timeRange = np.arange(centroids[0]-phalf + tunit, centroids[0] + phalf , tunit) 

for _ in range(1,8):
    for timepoint in timeRange: 
        GaussModel=Ibase-Ibase*d*math.exp(-(timepoint-centroids[0])**2/w**2)  
        gmaster = np.append(gmaster, GaussModel)
plt.figure()
plt.plot(gmaster, "r")
plt.plot(time[0:len(gmaster)], data[0:len(gmaster)], "b+")
#plt.plot(time[0:len(gmaster)], gmaster,"r")
#plt.plot(data[0:len(gmaster)], "b")
plt.show()


########
# find the time of " cent 1 - phalf + tunit " 
# what is the time of centriod 8 ?

###################################################################
####   Chi - Square Measurement of the Gaussian Model Fit #########
###################################################################
print '\n               Chi Square Analysis\n                 "Gaussian Model"            '
print 'This is the chi squared measurement of how well the model fits the data.\nThis is given as:'
print '    X2= sum (  (observed - model)^2/ model   ).\nWhere "observed" is the normalized flux data, and "expected" is the\nvalue of the model at each point.\nIf this value is close to zero, the model is a good fit.\nIf X2 > len(gmaster), then the model isnt a good fit, and the data isnt likely of a planet transit.'
print '\nLength of model array: ',len(gmaster),'elements.'
X2=np.nansum ( ((data[0:len(gmaster)]-gmaster)**2.)/(gmaster))
#######
#Because flux was reduced by a factor of the mean of flux 
#( for a previously user declared range of time), we must scale the 
# X2 of the normalized flux data by the same factor of the mean of flux. 
#######
print 'X2 = ', round(X2,3)
b=X2*np.nanmean(flux[0:len(gmaster)])
print 'X2 x (mean of flux) = ', b,' is the properly rescaled Chi Square value.'
print 'bins = ', len(gmaster)
print 'Reduced Chi Square = X2/ d, d=n-c, n=len(gmaster), c= 3 (width, depth, and centroid)'
x=b/(len(gmaster)-3)
print 'Reduced Chi Square of rescaled Chi Square =', x
print 'A reduced Chi Square is ideal at a value of (1), on the order of one. '
print 'Original flux average for Q1: ', round(np.nanmean(flux[0:len(gmaster)]),2),'electrons/sec'


plt.plot(time[0:len(gmaster)],data[0:len(gmaster)],'b',label='KOI Data')
xgmaster=np.arange(time[0],time[0]+len(gmaster)*tunit,tunit)
plt.plot(xgmaster,gmaster,'r', label='Gaussian Model')
plt.legend()
plt.xlabel('Time (days)')
plt.ylabel('Normalized Flux Intensity')
plt.title(r'Gaussian Model Fitting to {}, X2={}'.format(koi,round(b,6)))
plt.show()
###################################################################
####   Lets try another model, a 'box' model  #####################
###################################################################
# Initiate a new figure to plot, a "model time range" for the box model, and 
# the array to put the model's theoretical flux values. 

#new figure
plt.figure()

# the time values for the boxmodel's duration
modeltimerange=[]     
# the theory curve values
boxmaster=[]           

# Construct the model for eight transits
for ic in centroids[0:8]:
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
    
    plt.plot(x1,y1,'r'),plt.plot(x2,y2,'r'),plt.plot(x3,y3,'r')
    plt.vlines(max(x1),Ibase-d,Ibase,'r')
    plt.vlines(max(x2), Ibase-d,Ibase,'r')

#plt.plot(time[0:len(boxmaster)],fluxmean2[0:len(modeltimerange)],'b',label='KOI Data')
plt.legend(),plt.show()
############################################
## Chi - Square Analysis of Box Model ##
############################################
print '\n               Chi Square Analysis\n                 "Box Model"            '
print 'This is the chi squared measurement of how well the box model fits the data.'
print 'The chi squared measurement of how well the model fits the data is given as:'
print 'X2= sum (  (observed - model)^2/ model   ) for each measurement. If this value is close to zero, the model is a good fit.\nIf X2 > len(boxmaster), then the model isnt a good fit, and the data isnt likely of a planet transit.'
# the chai square mesurement needs to not count 'nan' values
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



#############################################
##  calculate luminosity of host star  ######
#############################################
print "Luminosity of Parent Star:"
sigma=5.6704e-08
# call for the stellar radius from MAST. 
SRad=koi.koi_srad*695500# km
#equation for luminosity:
L=4*math.pi*SRad**2*sigma*STeff**4 
#uncertainty of luminosity is given by error in teff,
#which is stated as +- 200 kelvin.  
errL=L*math.sqrt((4*200./STeff)**2)
### state SIGFIGS in dL, use same SIGFIGS in L


############################################
##### Calculate the possible values for#####
##### the effective planetary temps    #####
##### given different albedo values    #####
#####           (really cool!!)        #####
############################################

# the SMA is the semi major axis, used as the distance in 
AU=149597870.700 # in km
au = AU
Mstar=koi.koi_smass  # mass of parent star in solar masses. 
Mstar_kg=Mstar*1.9891e+30 # mass of parent star converted to kg.
Mstar='%s'%float('%.3g'%Mstar) # mass of star in Solar masses
Mstar_kg='%s'%float('%.3g'%Mstar_kg)# returns value with sigfigs called for in paren's. 

print "Mass of parent star : {} solar masses, {} kg".format(Mstar, Mstar_kg)


SMA=koi.koi_sma;print SMA,"AU: Semi Major Axis"
sma=SMA*au; print sma,"km: Semi Major Axis"
print ''
print "Luminosity of star",L,"plus or minus",errL,"W/m^2"
plt.figure()
#The following algorithm calculates the values
## sought above, and plots them on a smooth curve.
A=np.arange(0,1,.01) #where the .01 is the resolution of the curve.
PTeff=STeff*(1- A)**.25*(SRad/(2*sma))**.5; #print PTeff,"K :Effective Planetary Temp, Albedo of",A
pmax= STeff*(1-.3)**.25*(SRad/(2*sma))**.5;
plt.plot(A,PTeff,label="Theoretical Temperature");#plot shows how planetary eq temp relates to albedo
plt.axhline(y=pmax,xmin=0,xmax=.3,color='r')
plt.axvline(x=.3, ymin=0, ymax=.843, color='r',label="Temp at BOND Albedo")
plt.xlabel('Varying Values of Planetary Albedo')
plt.ylabel('Effective Planetary Temperature (K)')
plt.title("Varying Teff for Kepler Object of Interest (KOI) %s"%(name))
plt.legend()
plt.show();#plt.figure()#
#####################################################
# planetary effective temperature = stellar temperature x ( 1- albedo)^(1/4) x sqrt(stellar radius / major axis)

for A in range(0,9,1):
    if A!=3:
        # calculate effective planetary temp
        PTeff=STeff*((1-A*.1)/4)**(.25)*(SRad/(sma))**.5; 
        PTeffRounded = int(round(PTeff,-1));print PTeffRounded,"K :Effective Planetary Temp, Albedo of",A
        # print "effective planetary temp at (current albedo)
    else:
        # calculate planetary temperature anyway, but identify this as the temperature at the "bond albedo" - the albedo which is similar to Earth. 
        print"***"
        PTeff=STeff*((1-A*.1)/4)**(.25)*(SRad/(sma))**.5
        PTeffRounded = int(round(PTeff,-1));print PTeffRounded,"K :Effective Planetary Temp, BOND Albedo of",A
        print"***"
        # this time print out that this is the same albedo as earth - the bond albedo
        

print'\n'
print "According to the Kepler Database, {} has an equilibruim temperature of {} Kelvin. With this program I have a calculated difference of about {}%, likely due to rounding differences between this program and the Kelpler Science Team's methods.".format(name,koi.koi_teq,math.fabs(round(100-100*koi.koi_teq/pmax,2)))
if 273<=pmax<=373.15:
    print "Based on the calculated equilibrium temperature of the planet with an estimaed albedo similar to Earth's, THIS COULD BE A HABITABLE PLANET!!!"
elif pmax<273:
    print "Based on an equilibrium temperature of {} Kelvin, this is not a habitable planet because the temperature is to low for liquid water to exist.".format(round(pmax),2)
elif pmax>373.15:
    print "Based on the equilibrium temperature of {} Kelvin, this is not a habitable planet because the temperature is to high for liquid water to exist.".format(round(pmax),2)
