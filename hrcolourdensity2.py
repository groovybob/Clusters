import astropy.units as u
import gzip
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from scipy.stats import kde
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from scipy.ndimage.filters import gaussian_filter
import pandas as pd
from scipy.interpolate import interpn
# Suppress warnings. Comment this out if you wish to see the warning messages
import warnings
warnings.filterwarnings('ignore')

target = input("""Here is a list of clusters for you to examine:\nm2    m3    m4    m5\nm6    m7    m9    m10
m11    m12    m13    m14\nm15    m18    m19    m21
m22    m23\n""")
print("Please wait while I collect your data.")

#print(r['source_id'])
def plotradec (r):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.scatter(r['ra'], r['dec'], s=2 * (22 - r['phot_g_mean_mag']), color='w', alpha=0.3, marker='.')
    # plt.xlim(-80,80)
    # plt.ylim(-120,120)
    ax.set_facecolor("k")
    plt.show()
    plt.title("HR diagram for "+target)
    ax.set_facecolor("k")
    plt.xlim(-2,2)
    plt.xlabel("bp_rp")
    plt.ylabel("phot_g_mean_mag")
    ax.invert_yaxis()
    #ax.invert_xaxis()
    plt.show()
def colour_plot(x,y):
    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
    nbins = 300
    k = kde.gaussian_kde([x, y])
    xi, yi = np.mgrid[x.min():x.max():nbins * 1j, y.min():y.max():nbins * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    # Make the plot
    ax.set_title(f'H-R Diagram for {target}.\n (Gaia Dataset II)')
    ax.title.set_fontsize(20)
    ax.set_xlabel('Color index B-R')
    ax.xaxis.label.set_fontsize(20)
    ax.set_ylabel('Mean magnitude in the G band')
    ax.yaxis.label.set_fontsize(20)
    ax.tick_params(axis='both', labelsize=14)
    plt.pcolormesh(xi, yi, zi.reshape(xi.shape))

    #ax.invert_yaxis()
    ax.set_xlim(-0.5, 2)
    ax.set_ylim(22, 10)
    ax.grid()
    plt.colorbar()
    plt.show()

    # Change color palette
    #plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.Greens_r)
    #plt.show()


object = Simbad.query_object(target)
ra = object["RA"]
RA = ra.data[0]
dec = object["DEC"]
DEC = dec.data[0]
string = str(RA+ ' '+DEC)
c = SkyCoord(string, unit=(u.hourangle, u.deg))

Ra = str(c.ra.degree)
Dec = str(c.dec.degree)

string2 = str(Ra+","+Dec)
bigstring = str("SELECT all source_id,ra,ra_error,dec,dec_error,parallax,parallax_error,phot_g_mean\
_mag,bp_rp,radial_velocity,radial_velocity_error,phot_variable_flag,teff_val,a_g_val FROM gaiadr2.gaia_source  \
WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',"+string2+",0.5))=1")
job = Gaia.launch_job_async(bigstring, dump_to_file=True)
r = job.get_results()
r = np.asarray(r)
r = pd.DataFrame(r)
x = r['bp_rp']
x = x.dropna()
mask = pd.isna(r['bp_rp'])
#x = x[np.logical_not(np.isnan(x))]
y = (r['phot_g_mean_mag'])[-mask]
#y = y[np.logical_not(np.isnan(y))]
action = input(f"""What would you like to look at?\n a) The RA DEC plot of {target}\n b) The HR diagram for {target}\n 
 Type a or b.""")
if action =='b':
    colour_plot(x, y)
if action =='a':
    plotradec(r)

plt.show()
#plt.hist(r['parallax'],bins = 10000)