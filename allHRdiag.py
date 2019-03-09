#  this script creates HR diagrams for all stars in 0.5geg of each cluster including foreground stars

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
targets = ['m2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm9', 'm10', 'm11', 'm12', 'm13', 'm14', 'm15', 'm18', 'm19', 'm21',\
           'm22', 'm23', 'm25', 'm26', 'm28', 'm29', 'm30', 'm34', 'm35', 'm36', 'm37', 'm38', 'm39', 'm41', 'm45',\
           'm46', 'm47', 'm48', 'm50', 'm52', 'm53', 'm54', 'm55', 'm56', 'm67', 'm68', 'm69', 'm70', 'm71', 'm72',\
           'm75', 'm79', 'm80', 'm92', 'm93', 'm103', 'm107']

#  creates a function that plots BP-RP against Magnitude in G band

def plotHRdiag(r):
    fig, ax = plt.subplots(figsize=(8, 10))
    ax.set_xlim(-0.5, 4)
    ax.set_ylim(21.2, 10)
    ax.grid()
    ax.set_title(f'H-R Diagram for {target}.\n (Gaia Dataset II)')
    ax.title.set_fontsize(20)
    ax.set_xlabel('Color index B-R')
    ax.xaxis.label.set_fontsize(20)
    ax.set_ylabel('Mean magnitude in the G band')
    ax.yaxis.label.set_fontsize(20)
    textstr = f"Number of stars {len(r['bp_rp'])}"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    t = r['bp_rp']
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    ax.scatter(r['bp_rp'], r['phot_g_mean_mag'],s=1, edgecolors='none',vmin=-1., vmax=4., c=t, cmap = 'jet')
    fig.patch.set_facecolor('dimgray')
    ax.set_facecolor('dimgray')
    ax.tick_params(axis='both', labelsize=14)
    fig.savefig(f'C:/Users/User/GITHUB/Clusters/HRdiags/H-R Diagram for {target}', dpi=fig.dpi)
    #plt.show()

#  function to get target RA DEC for each target
#  uses the RA DEC (conerted to decimal degrees) to get the desired data from Gaia database

def generate_r(target):
    target = str(target)
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
    bigstring = str("SELECT all gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error,gaia_source.phot_g_mean_mag,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.bp_rp,gaia_source.bp_g,gaia_source.g_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.teff_val,gaia_source.a_g_val,gaia_source.lum_val,gaia_source.lum_percentile_lower,gaia_source.lum_percentile_upper  FROM gaiadr2.gaia_source  WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',"+string2+",0.5))=1")
    job = Gaia.launch_job_async(bigstring, dump_to_file=True)
    r = job.get_results()
    return r

#  loops over each target to download data, plot and save image. I recommend you use your own directory adresses.

for target in targets:
    r = generate_r(target)
    plotHRdiag(r)
    dir_name = "C:/Users/User/GITHUB/Clusters"
    test = os.listdir(dir_name)
    for item in test:
        if item.endswith(".vot"):
            os.remove(os.path.join(dir_name, item))



