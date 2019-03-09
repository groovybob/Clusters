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


def plotradec(r):
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


def density_scatter(x, y, ax=None, sort=True, bins=20, **kwargs):
    """
    Scatter plot colored by 2d histogram
    """
    if ax is None :
        fig , ax = plt.subplots()
    data, x_e, y_e = np.histogram2d(x, y, bins=bins)
    z = interpn((0.5*(x_e[1:] + x_e[:-1]), 0.5*(y_e[1:]+y_e[:-1])), data, np.vstack([x, y]).T, method="splinef2d",
                bounds_error=False)

    # Sort the points by density, so that the densest points are plotted last
    if sort:
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax.scatter(x, y, c=z, **kwargs)
    ax.invert_yaxis()
    return ax


object = Simbad.query_object(target)
ra = object["RA"]
RA = ra.data[0]
dec = object["DEC"]
DEC = dec.data[0]
string = str(RA + ' '+DEC)
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
y = (r['phot_g_mean_mag'])[-mask]
action = input(f"""What would you like to look at?\n a) The RA DEC plot of {target}\n b) The HR diagram for {target}\n 
 Type a or b.""")
if action == 'b':
    density_scatter(x, y, bins=[300, 300])
if action == 'a':
    plotradec(r)

plt.show()
