import kplr 
import scipy.stats as st
from astropy.time import Time
import math
import numpy as np

import os
global os

from PhaseFoldingAndRegularPlots import EstablishKOI,Normalize_LightCurves, FindNextEclipse

EstablishKOI()
Normalize_LightCurves()
FindNextEclipse()