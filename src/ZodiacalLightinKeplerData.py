#### IF THIS IS THE FIRST TIME RUNNING THIS PROGRAM: You must pip install the following libraries:
## kplr, matplotlib, scipy, math, numpy, pylab, pyfits. Some of these take a significant amount of storage space

import kplr # to install 'kplr', go to cmd prompt and type "pip install kplr"  - no quotes. 
import matplotlib.pyplot as plt
#import scipy
import numpy as np
import pylab
import sys
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText

'''
Current Status: 
    working on conceptual description of science process
            
        
Interesting Objects:
    5073.01
    5059.01
    422.01
    523.01
    771.01
    1089.01
    
'''

'''
OBJECTIVE
    What we want to do is make a measurement of how much detectable zodiacal 
light is in a given system. The nature of the zodiacal light will be to scatter 
light that passes through through the dust/debris cloud, leaving a small 
transit signiture in the Kepler Light Curve fits file.  

METHOD
    One way to make this measurement is to take the light curve data file, extract
the flux values into an array, and plot this data iteratively from -period/2 to 
period/2. This will 'fold' each light curve file over the course of one orbit 
of the given KOI. At this point, there are N transits overlaid on each other in 
one common interval, where N is necessarily the number of transits in the entire
light curve file (except in periods where no data is recorded). In each "time bin"
there will then be N data points of flux for that time in the orbit. The average of
these can be taken as the best estimate of the true flux value at that point in time,
from which a model of the true flux vs time can be descirbed and plotted as an overlaid
red line. 
    The base intensity of the host star is taken as the best estimate of flux values
that are not associated with the transit event. Because we will be mainly looking
at KOI with periods greater than 300 days(**needs justification - larger objects 
farther away??), we can determine the average base intensity of the star by taking
the average of flux that lies between 
    For each time bin, the best estimate (model value) of flux can be subtracted 
from the base intensity, from which the maxmium value outside'''

# Establish the search criteria for kepler objects of interest. 
# here, we search for koi that have period greater than 500 days, producing 
# a list 'kois' from which we can work with. 
kois = kplr.API().kois(where = "koi_period<800", sort=("koi_period",1))
print "Found {} KOI list with specified criteria.".format(len(kois))

# here, we further filter the results to those who are not false positives, 
# leaving only the candidates and the confirmed planets. 

good_kois =[]
for i in kois:
    if i.koi_disposition != 'CONFIRMED':
        kois.remove(i)
    else:
        good_kois = np.append(good_kois, i)
        print i.kepoi_name, i.koi_disposition, "Period:", i.koi_period
count = len(good_kois)
print count, " Confirmed Planets."

# here we ask the user how many of the above objects to plot light curves for. 
ask_continue = str(raw_input("Do you want to generate plot figures for (a)ll of these, (s)ome, or (n)one? : "))
while ask_continue != "a" and ask_continue != "n" and ask_continue != "s":
    print "Sorry, choose again."
    ask_continue = str(raw_input("Do you want to generate plot figures for (a)ll of these, (s)ome, or (n)one? : "))


def generate_figures(limit):
    #For each koi in the the search results, to the user-specified limit, 
    for i in range(0,limit,1):
        
        # if the koi is not a false positive, and its period is less than 800 days
        #     so that we can see at least two transits
        if good_kois[i].koi_disposition != "FALSE POSITIVE":
            if good_kois[i].koi_period<800.0:
                
                # get the light curve flux, normalize it by dividing by the mean, and the time for each flux data point
                lcs = good_kois[i].get_light_curves(short_cadence = False, fetch = False, clobber = False)
                time, data = [],[]
                # for each light curve file, fill time and data (flux/mean_flux) arrays
                for lc in lcs: 
                    with lc.open() as f:
                        hdu_data = f[1].data
                        time = np.append(time, hdu_data["time"])
                        data = np.append(data, (hdu_data['pdcsap_flux'])/(np.nanmean(hdu_data["pdcsap_flux"])))
                # plot a new figure for plotting the data for this koi
                plt.figure(figsize = (24,14))
                plt.plot(time,data,label ='Normalized Flux')
                plt.legend()
                plt.xlabel('Time (days)')
                plt.ylabel('Normalized Flux for {}'.format(good_kois[i].kepoi_name))
                plt.title("Light Curve {}".format(good_kois[i].kepoi_name))
            
                # State some properties of the koi in a caption box
                KOI_Properties = AnchoredText("Period: {} days\nFirst Centroid: {} days\nDisposition: {}".format(round(good_kois[i].koi_period,1), round(good_kois[i].koi_time0bk,1),good_kois[i].koi_disposition),loc= 2, prop=dict(size = 10),frameon=True)
                KOI_Properties.patch.set_boxstyle("round, pad = 0.,rounding_size= 0.2")
                
                plt.gca().add_artist(KOI_Properties)
                
                #plt.show()      # ONLY UNCOMMENT THIS LINE IF YOU ARE SURE THERE ARE LESS THAN 20 GENERATED FIGURES. 
                            
                figureName = str('Kepler LC (long cadence) for ') + str(good_kois[i].kepoi_name)+str('.png')
                #print "(",i+1,")",figureName
                pylab.savefig('{}.png'.format(figureName))
                plt.close()
                print "Saved ",figureName," to working directory."
    print "Finished. Check working directory for saved figures"#, '{}'"#.format(figureName)

if ask_continue == "a":
    generate_figures(len(good_kois))
elif ask_continue == "s":
    how_many = int(raw_input("How many (integer) figures do you want to plot? : "))
    while how_many > count:
        print "There aren't that many KOI that matched your search criteria. There are only",count,"objects. Choose a number of plots equal or less than",count,"."
        how_many = int(raw_input("How many (integer) figures do you want to plot? : "))
    generate_figures(how_many)
elif ask_continue == "n":print "\nGoodbye", sys.exit()
        