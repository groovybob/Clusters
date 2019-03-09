#  this script creates HR diagrams for all stars in 0.5geg of each cluster and tries to
#  discard stars outside of cluster by means of proper motion matching
import csv
import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
import matplotlib.pyplot as plt
import os


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

#  creates a function that plots BP-RP against Magnitude in G band


def plotHRdiag(r):
    targetx = r['pmdec']
    meanx = targetx.mean()
    sdx = targetx.std()
    targety = r['pmra']
    meany = targety.mean()
    sdy = targety.std()
    selected_r = r[(targetx > meanx - sdx) & (targetx < meanx + sdx) & (targety > meany - sdy) & (targety < meany + sdx)]
    variables = r[(r['phot_variable_flag'] == b'VARIABLE')]
    fig, ax = plt.subplots(figsize=(8, 10))
    ax.set_xlim(-0.5, 4)
    ax.set_ylim(21.2, 10)
    ax.grid()
    ax.set_title(f'H-R Diagram for {str(targets[i])}.\n (Gaia Dataset II)')
    ax.title.set_fontsize(20)
    ax.set_xlabel('Color index B-R')
    ax.xaxis.label.set_fontsize(20)
    ax.set_ylabel('Mean magnitude in the G band')
    ax.yaxis.label.set_fontsize(20)
    textstr = f"Number of stars {len(selected_r['bp_rp'])}\nBlack + = variable stars"
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    t = selected_r['bp_rp']
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    ax.scatter(selected_r['bp_rp'], selected_r['phot_g_mean_mag'],s=1, edgecolors='none',vmin=-1., vmax=4., c=t, cmap='jet')
    ax.scatter(variables['bp_rp'], variables['phot_g_mean_mag'],s=5, marker = 'x',color = 'k')
    fig.patch.set_facecolor('dimgray')
    ax.set_facecolor('dimgray')
    ax.tick_params(axis='both', labelsize=14)
    fig.savefig(f'C:/Users/User/GITHUB/Clusters/HRdiags/H-R Diagram for {str(targets[i])}', dpi=fig.dpi)
    #plt.show()

#  function to get target RA DEC for each target
#  uses the RA DEC (conerted to decimal degrees) to get the desired data from Gaia database


def generate_r(target):
    target = str(targets[i])
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
    bigstring = str("SELECT all gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,gaia_source.pmra,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error,gaia_source.phot_g_mean_mag,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.bp_rp,gaia_source.bp_g,gaia_source.g_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.teff_val,gaia_source.a_g_val,gaia_source.lum_val,gaia_source.phot_variable_flag,gaia_source.lum_percentile_lower,gaia_source.lum_percentile_upper  FROM gaiadr2.gaia_source  WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',"+string2+","+size+"))=1")
#  previous line kept failing if I broke it up to fit in view so I left it long.
    job = Gaia.launch_job_async(bigstring, dump_to_file=True)
    r = job.get_results()
    return r

#  loops over each target to download data, plot and save image. I recommend you use your own directory adresses.


for i in range(len(targets)):
    size = str(sizes[i])
    if not os.path.isfile(f'C:/Users/User/GITHUB/Clusters/HRdiags/H-R Diagram for {str(targets[i])}.png'):
        print("processing " + str(targets[i]))
        r = generate_r(targets[i])
        plotHRdiag(r)
        dir_name = "C:/Users/User/GITHUB/Clusters"
        test = os.listdir(dir_name)
        for item in test:
            if item.endswith(".vot"):
                os.remove(os.path.join(dir_name, item))

    else:
        print("skipping " + str(targets[i]) )

