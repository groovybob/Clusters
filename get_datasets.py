import astropy.units as u
import gzip
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from scipy.stats import kde
import matplotlib.pyplot as plt
import numpy as np
# Suppress warnings. Comment this out if you wish to see the warning messages
import warnings
warnings.filterwarnings('ignore')

target = input("""Here is a list of clusters for you to examine:\nm2    m3    m4    m5\nm6    m7    m9    m10
m11    m12    m13    m14\nm15    m18    m19    m21
m22    m23    m25    m26\nm28    m29    m30    m34
m35    m36    m37    m38\nm39    m41    m45    m46
m47    m48    m50    m52\nm53    m54    m55    m56
m67    m68    m69    m70\nm70    m71    m72    m75
m79    m80    m92    m93\nm103   m107""")
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

def plotHRdiag(r):
    #plt.scatter(r['bp_rp'], r['phot_g_mean_mag'], color='w', alpha=0.3, marker='.')
    fig, ax = plt.subplots(figsize=(8, 10))

    ax.set_xlim(-0.5, 2)
    ax.set_ylim(22, 10)
    #ax.invert_yaxis()
    ax.grid()
    ax.set_title(f'H-R Diagram for {target}.\n (Gaia Dataset II)')
    ax.title.set_fontsize(20)
    ax.set_xlabel('Color index B-R')
    ax.xaxis.label.set_fontsize(20)
    ax.set_ylabel('Mean magnitude in the G band')
    ax.yaxis.label.set_fontsize(20)

    ax.scatter(r['bp_rp'], r['phot_g_mean_mag'],
               #           s=50, edgecolors='none', alpha=0.015, c='k')
               s=1, edgecolors='none', c='k')
    plt.legend(f"Number of stars {len(r['bp_rp'])}")
    ax.tick_params(axis='both', labelsize=14)




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

action = input(f"""What would you like to look at?\n a) The RA DEC plot of {target}\n b) The HR diagram for {target}\n 
 Type a or b.""")
if action =='b':
    plotHRdiag(r)
if action =='a':
    plotradec(r)
plt.show()