import sys
gitpaths=['/Users/remy/lustre/naasc/users/rindebet/github/quickclump/',
          '/Users/remy/lustre/naasc/users/rindebet/github/pyprops/']
for gitpath in gitpaths:
    if not gitpath in sys.path:
        sys.path.insert(0,gitpath)





        
# C0:
fitsfile="12CO.combined.20130113.nopbcor.Jybm.kms.fits"

# assignment file for comparison
assignfile2=None


# C023 rimc
dir='/Users/remy/lustre/naasc/users/rindebet/30dor/chevance/'
fluxfile=dir+'dor12coS_r1.5.pb.fapex.ch0.trim.fits'
datafile=dir+"dor12coS_r1.5.fapex.trim.pbcor.fits"

cutoff=0.04
pkmin=0.1
dTleaf=0.004
dTleaf=0.008





# C023 rimc
dir='/Users/remy/lustre/naasc/users/rindebet/30dor/chevance/'
fluxfile=dir+'Dor13co_sT1.0_smoT2.0_cutT0.01_noiseT5.0.pb.ch0.trim.fits'
datafile=dir+"Dor13co_sT1.0_smoT2.0_cutT0.01_noiseT5.0.trim.fits"
N13file="Dor_peak_n13cube.fits"
tau13file="Dor_peak_tau13.fits"

# 13 ~ 0.25*12 often.
cutoff=0.01
pkmin=0.025
dTleaf=0.002
dTleaf=0.004
pkmin=0.03

# match dendro settings:
sigma=0.00252
cutoff=3.*sigma
dTleaf=2.5*sigma
pkmin=10*sigma

# increase the floor and see if we can get closer to dendro:

cutoff=20.*sigma








parmstr="c"+str(cutoff)+"_dT"+str(dTleaf)+"_min"+str(pkmin)
clumpfile=datafile[:-4]+"clumps_"+parmstr+".fits"
import os
if not os.path.exists(clumpfile):
    import qc,pdb
    qc.main(argv=['--Tcutoff',str(cutoff),'--Tpkmin',str(pkmin),
                  '--dTleaf',str(dTleaf),'--verbose',datafile])
    
    import shutil
    shutil.move(datafile[:-4]+"clumps.fits",clumpfile)
    clumptxt=datafile[:-4]+"clumps_"+parmstr+".txt"
    shutil.move(datafile[:-4]+"clumps.txt",clumptxt)

from pyprops import pyprops
pyprops(datafile,fluxfile,clumpfile,datafile[:-4]+"clumps_"+parmstr,assignfile2=assignfile2,montecarlo=0,doplots=True)

import numpy as np
# calculate masses
dist=48000.

from astropy.io import fits
cubehdr=fits.getheader(datafile)
dx=cubehdr['cdelt2']*3600/206265*dist*3.09e18 # cm
nu0=cubehdr['restfrq']
dnu=cubehdr['cdelt3']
dv=2.99792458e5 *np.absolute(dnu)/nu0
K2jybm = 8.185e-7 * cubehdr['bmaj']*cubehdr['bmin']*3600**2 * (nu0/1e9)**2

abundance=1e6 # n(H2)/n(13CO)

from astropy.io import fits
masscube=fits.getdata(N13file)
taucube=fits.getdata(tau13file)
factor=dv*dx**2 *abundance *1.67e-24 /1.99e33 # Npix > M
    
assignfile=clumpfile

import pickle
mclfile=datafile[:-4]+"clumps_"+parmstr+".mcl.pkl"
if os.path.exists(mclfile):
    mcl=pickle.load(open(mclfile,'rb'))
else:
    assigncube=fits.getdata(assignfile)
    ncl=assigncube.max()
    mcl=np.zeros(ncl)

    for i in range(ncl):
        z=np.where(assigncube==(i+1))
        mcl[i]=np.nansum(masscube[z])*factor
    pickle.dump(mcl,open(mclfile,'wb'))

taufile=datafile[:-4]+"clumps_"+parmstr+".taupk.pkl"
if os.path.exists(taufile):
    taupk=pickle.load(open(taufile,'rb'))
else:
    assigncube=fits.getdata(assignfile)
    ncl=assigncube.max()
    taupk=np.zeros(ncl)

    for i in range(ncl):
        z=np.where(assigncube==(i+1))
        taupk[i]=np.nanmax(taucube[z])
    pickle.dump(taupk,open(taufile,'wb'))




