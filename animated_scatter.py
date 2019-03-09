import csv
import astropy.units as u
import gzip
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from scipy.stats import kde
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import time

# Suppress warnings. Comment this out if you wish to see the warning messages
import warnings
warnings.filterwarnings('ignore')
#  I have made an array of the Messier cluster objects for ease of use later.
targets = []
sizes = []
with open("C:/Users/User/GITHUB/Clusters/clusters.csv") as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
        targets.append(row[0])
        sizes.append(row[1])
