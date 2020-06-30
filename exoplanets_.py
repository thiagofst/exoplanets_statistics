# 16 Jun 2020
# Different catalogues: EXO; OEC; KOI

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os, statistics
import cmocean as oc
import seaborn as sns
import pandas as pd

from astropy.io import ascii

plt.rcdefaults()
plt.rc('text', usetex = True)
plt.rc('font', family = 'serif')

os.chdir('/home/volkoff/VVV/exoplanets/plots')

# EXO: exoplanet.eu
exo_ = ascii.read('/home/volkoff/VVV/exoplanets/exoplanet.eu_catalog.csv')
m = exo_['planet_status'] == 'Confirmed'
exo = exo_[m]
# Flags for detection methods
T = exo['detection_type'] == 'Primary Transit';
R = exo['detection_type'] == 'Radial Velocity';
G = exo['detection_type'] == 'Microlensing';
I = exo['detection_type'] == 'Imaging'

# OEC: Open Exoplanet Catalogue (!wget https://raw.githubusercontent.com/OpenExoplanetCatalogue/oec_tables/master/comma_separated/open_exoplanet_catalogue.txt)
names = ['id', 'binary', 'mass', 'radius', 'period', 'a', 'e', 'per', 'long', 'asc', 'incl',
         'teq', 'age', 'method', 'year', 'updated', 'ra', 'dec', 'dist', 's_mass', 's_radius',
         's_feh', 's_teff', 's_age']
oec = pd.read_csv('/home/volkoff/VVV/exoplanets/open_exoplanet_catalogue.txt', skiprows = 30, names = names)

# KOI: Kepler Objects of Interest (!wget 'http://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&select=*' -O kois.csv)
# https://exoplanetarchive.ipac.caltech.edu/docs/API_kepcandidate_columns.html Documentation
kois = pd.read_csv("/home/volkoff/VVV/exoplanets/kois.csv")

m = kois.koi_disposition == "CONFIRMED" # Flag for confirmed exoplanets only

(x, y) = kois[m].koi_period, kois[m].koi_prad
(x_, y_) = exo['orbital_period'], exo['radius'] # From EXO
(_x_, _y_) = oec['period'], oec['radius'] # From OEC

# Plot I
# Orbital period vs. planet radius (KOI: Confirmed exoplanets)
plt.figure()
plt.xlim(1.1, 450.); plt.ylim(0.3, 30) # For KOI and OEC
plt.grid(True, alpha = 0.3)
plt.xlabel(r'Orbital period/days', fontsize = 10)
plt.ylabel(r'Radius/(R/R$\oplus$)', fontsize = 10)
plt.loglog(x, y, '.k', ms = 4) # Main plot

x0 = np.exp(np.linspace(np.log(0.1), np.log(9e4), 500))
for f in np.exp(np.linspace(np.log(0.001), np.log(100), 8)):
    y0 = np.sqrt(f * np.sqrt(x0))
    # y0[x0 > 365 * 4. / 3.] = np.sqrt(f * 4.5 / 3.)
    plt.plot(x0, y0, "k", lw=0.5, alpha=0.5)

fmt = matplotlib.ticker.FormatStrFormatter("%.0f")
plt.gca().xaxis.set_major_formatter(fmt)
plt.gca().yaxis.set_major_formatter(fmt)

''' Solar system planets '''
radius = np.array(['2439.7','6051.8','6371.00','3389.5','69911','58232','25362','24622'], dtype = float)
period = np.array(['0.2408467','0.61519726','1.0000174','1.8808476','11.862615','29.447498','84.016846','164.79132'], dtype = float)
mass   = np.array(['3.30e23', '4.86e24', '5.97e24', '6.41e23', '1.89e27', '5.68e26', '8.68e25', '1.024e26'], dtype = float)
a      = np.array(['0.3871', '0.7233', '1', '1.5237', '5.2034', '9.5371', '19.1913', '30.0690'], dtype = float)

plt.plot(period*365, radius/radius[2], 'o', color = 'red')
plt.xlim(0.5, max(x0))

''' Unconfirmed/Candidate KOIs in the background'''
m = kois.koi_disposition == "CANDIDATE"
x, y = kois[m].koi_period, kois[m].koi_prad
plt.plot(x, y, ".", color = "#6baed6", ms = 4, alpha = 0.3, zorder = -1) 
plt.tight_layout(); plt.savefig('orbital_radius.png', dpi = 200); 
plt.close()

# Plot II
# Planet mass vs. Star mass with EXO
plt.figure()
plt.grid(True, alpha = 0.3)

plt.xlabel(r'Planet mass (M/M$_{Jup}$)', fontsize = 10)
plt.ylabel(r'Star mass (M/M$\odot$)', fontsize = 10)

plt.loglog(exo['mass'], exo['star_mass'], '.', color = 'black', ms = 4)

''' Unconfirmed/Candidate EXO in the background'''
j = exo_['planet_status'] != 'Confirmed'
(u, v) = exo_[j]['mass'], exo_[j]['star_mass']
plt.plot(u, v, ".", color = "#6baed6", ms = 4, alpha = 0.3, zorder = -1) 

#plt.loglog(exo[T]['mass'], exo[T]['star_mass'], 'o', color = 'black', ms = 4, label = r'Primary Transit ({0})'.format(len(exo[T])))
#plt.loglog(exo[R]['mass'], exo[R]['star_mass'], 'o', color = 'red', ms = 4, label = r'Radial Velocity ({0})'.format(len(exo[R])))
#plt.loglog(exo[G]['mass'], exo[G]['star_mass'], 'o', color = 'navy', ms = 4, label = r'Gravitational Microlensing ({0})'.format(len(exo[G])))
#plt.loglog(exo[I]['mass'], exo[I]['star_mass'], 'o', color = 'yellow', ms = 4, label = r'Direct Imaging ({0})'.format(len(exo[I])))

#plt.legend(fontsize = 10, markerscale = 1, shadow = 'True', loc = 0);
plt.tight_layout(); plt.savefig('exo_star_mass', dpi = 200) 
plt.close()

# Plot III
# Planet mass (Mjup) vs. Planet radius with EXO
j = exo_['planet_status'] != 'Confirmed' # Unconfirmed then
(x, y) = exo_[j]['mass'], exo_[j]['radius']
(xT, yT)	= exo[T]['mass'], exo[T]['radius']
(xR, yR) 	= exo[R]['mass'], exo[R]['radius']
(xG, yG) 	= exo[G]['mass'], exo[G]['radius']
(xI, yI) 	= exo[I]['mass'], exo[I]['radius']

plt.figure()
plt.grid(True, alpha = 0.4)

plt.xlim(1e-4, 1e2)
plt.ylim(1e-3, 1e2)

plt.xlabel(r'Planet mass (M/M$_{Jup}$)', fontsize = 10)
plt.ylabel(r'Planet radius (R/R$\oplus$)', fontsize = 10)

plt.loglog(y, x, ".", color = "#6baed6", ms = 4, alpha = 0.3, zorder = -1, label = r'__nolabel__') 
plt.loglog(xT, yT, '.', color = 'black', ms = 4, label = r'Transit')
plt.loglog(xR, yR, '.', color = 'darkmagenta', ms = 4, label = r'Radial velocity')
plt.loglog(xG, yG, '.', color = 'yellow', ms = 4, label = r'Gravitational Microlensing')
plt.loglog(xI, yI, '.', color = 'navy', ms = 4, label = r'Imaging')

plt.legend(fontsize = 10, markerscale = 1, shadow = 'True', loc = 0);

plt.tight_layout(); plt.savefig('exo_mass_radius.png', dpi = 200)
plt.close()

# Plot IIII
# Relationship between semi major axis and masses
plt.figure(); 
plt.grid(True, alpha = 0.4);

plt.xlabel(r'log$_{10}$~a/AU', fontsize = 10);
plt.ylabel(r'log$_{10}$~M/M$_{Jup}$', fontsize = 10)

plt.loglog((exo[T]['semi_major_axis']), (exo[T]['mass']), '.', color = 'black', label = r'Transit ({0})'.format(len(exo[T])))
plt.loglog((exo[R]['semi_major_axis']), (exo[R]['mass']), '.', color = 'red', label = r'R. Velocity ({0})'.format(len(exo[R])))
plt.loglog((exo[G]['semi_major_axis']), (exo[G]['mass']), '.', color = 'navy', label = r'Microlensing ({0})'.format(len(exo[G])))
plt.loglog((exo[I]['semi_major_axis']), (exo[I]['mass']), '.', color = 'yellow', label = r'Imaging ({0})'.format(len(exo[I])))

plt.loglog(a, mass/mass[4], 'og', ms = 4)

snow_line = 2.7 # Snow line 170 K at 2.7 AU (Hayashi, 1981)[

plt.axvline(snow_line, lw = 1.5, color = 'blue', label = '__nolabel__')
plt.legend(fontsize = 10, markerscale = 1, shadow = 'True', loc = 0);
plt.tight_layout(); plt.savefig('a-Mass.png', dpi = 200)
plt.close()

# Plot V
# Planet mass vs. Planet radius
plt.figure(); 
plt.grid(True, alpha = 0.4);
plt.title(r'N = {0} exoplanets'.format(len(exo['mass'])), fontsize = 10, loc = 'right')

plt.xlabel(r'log$_{10}$~M/M$_{Jup}$', fontsize = 10);
plt.ylabel(r'log$_{10}$~R/R$_{Jup}$', fontsize = 10)

plt.errorbar(np.log10(0.00314558), np.log10(6371./69911.), fmt = 'o', color = 'red', label = r'Earth')
plt.axvline(np.mean(np.log10(exo['mass'])), ls = 'dashdot', lw = 1., c = 'black', alpha = 0.5)
plt.axhline(np.mean(np.log10(exo['radius'])), ls = 'dashdot', lw = 1., c = 'black', alpha = 0.5)

#plt.yscale('log'); plt.xscale('log')

points = plt.scatter(np.log10(exo['mass']), np.log10(exo['radius']), c = (exo['temp_calculated']), s = 6.0, cmap = oc.cm.ice, label = '__nolabel__')
plt.colorbar(points, label = r'T$_{star}$/K')

plt.legend(fontsize = 10, markerscale = 1, shadow = 'True', loc = 0);
plt.tight_layout(); plt.savefig('mass_radius.png', dpi = 200)
plt.close()

# Plot VI
# Colour-magnitude and colour-colour diagrams (KOI)
m = kois.koi_disposition == "CONFIRMED" # Flag for confirmed exoplanets only
(gmag, rmag, imag, zmag, jmag, hmag, kmag) = kois[m].koi_gmag, kois[m].koi_rmag, kois[m].koi_imag, kois[m].koi_zmag, kois[m].koi_jmag, kois[m].koi_hmag, kois[m].koi_kmag
(gmag_, rmag_, imag_, zmag_, jmag_, hmag_, kmag_) = kois[~m].koi_gmag, kois[~m].koi_rmag, kois[~m].koi_imag, kois[~m].koi_zmag, kois[~m].koi_jmag, kois[~m].koi_hmag, kois[~m].koi_kmag

# J-K vs. K
plt.figure()
plt.grid(True, alpha = 0.4)
plt.gca().invert_yaxis();
plt.xlim(-0.25, 1.5); plt.ylim(16, 6)
plt.xlabel(r'(J$-$K)/mag', fontsize = 10)
plt.ylabel(r'K/mag', fontsize = 10)
plt.plot(jmag - kmag, kmag, '.k', ms = 4)
plt.plot(jmag_ - kmag_, kmag_, ".", color = "#6baed6", ms = 4, alpha = 0.3, zorder = -1) 
plt.tight_layout(); plt.savefig('jks.png', dpi = 200)

# J-H vs. H-K
plt.figure()
plt.grid(True, alpha = 0.4)
plt.xlim(-0.15, 1.); plt.ylim(0.5, -0.4)
plt.xlabel(r'(J$-$H)/mag', fontsize = 10)
plt.ylabel(r'(H$-$K)/mag', fontsize = 10)
plt.plot((jmag - hmag),(hmag - kmag), '.k', ms = 4)
plt.plot((jmag_ - hmag_), (hmag_ - kmag_), ".", color = "#6baed6", ms = 4, alpha = 0.3, zorder = -1) 
plt.tight_layout(); plt.savefig('jhhk.png', dpi = 200)

plt.close()

'''# g-i vs g
plt.figure()
plt.grid(True, alpha = 0.4)
plt.gca().invert_yaxis();
plt.xlabel(r'(g$-$i)/mag', fontsize = 10)
plt.ylabel(r'g/mag', fontsize = 10)
plt.plot(gmag - imag, imag, '.k', ms = 4)
plt.plot(gmag_ - imag_, imag_, ".", color = "#6baed6", ms = 4, alpha = 0.3, zorder = -1) 
plt.tight_layout(); plt.show()
'''

# Plot VII
# log Distance vs. K
plt.figure(); 
plt.grid(True, alpha = 0.4); 
plt.gca().invert_yaxis()

plt.xlabel(r'Distance/kpc', fontsize = 10);
plt.ylabel(r'K/mag', fontsize = 10);

points = plt.scatter(np.log10(exo[T]['star_distance']), exo[T]['mag_k'], c = np.log10(exo[T]['star_radius']), s = 10, cmap = oc.cm.tempo, label = '%s exoplanets'%(len(exo[T])))

plt.colorbar(points, label = r'Host star radius/log$_{10}$R$\odot$')
plt.legend(fontsize = 10, markerscale = 1, shadow = 'True', loc = 0);
plt.tight_layout(); plt.savefig('distance_Kmag.png', dpi = 200)

# Plot VIII
# Planet surface equilibrium temperature and star effective temperature
# OEC
(tp, ts) = oec['teq'], oec['s_teff']
plt.figure(); 

plt.xlabel(r'Planet surface temperature [K]', fontsize = 10); 
plt.ylabel(r'Star effective temperature [K]', fontsize = 10); 

#plt.loglog(tp, ts, '.k', ms = 4); 

plt.tight_layout(); #plt.savefig('teff_oec.png', dpi = 200)
plt.close()

# KOI
m = kois.koi_disposition == "CONFIRMED" # Flag for confirmed exoplanets only
(tp, ts) = kois[m].koi_teq, kois[m].koi_steff
(tp_u, ts_u) = kois[~m].koi_teq, kois[~m].koi_steff # Unconfirmed/Candidate exoplanets

plt.figure()
plt.grid(True, alpha = 0.4)
plt.xlim(0.9e2, 1e4); plt.ylim(2.6e3, 1e4)

sun_teff = 5778
earth_teff = 287.15 # K. https://www.space.com/17816-earth-temperature.html#:~:text=GISS%20data%20show%20global%20average,57%20F%20(14%20C).

plt.xlabel(r'Planet surface temperature [K]', fontsize = 10); 
plt.ylabel(r'Star effective temperature [K]', fontsize = 10); 

plt.loglog(tp, ts, '.k', ms = 4); 
plt.loglog(earth_teff, sun_teff, '.r', ms = 10); 
plt.loglog(tp_u, ts_u, ".", color = "#6baed6", ms = 4, alpha = 0.3, zorder = -1) 

plt.tight_layout(); plt.savefig('teff.png', dpi = 200); 
plt.show()