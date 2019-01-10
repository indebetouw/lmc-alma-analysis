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
root="dor.rimc.12co"

cutoff=0.04
pkmin=0.1
dTleaf=0.004
dTleaf=0.008



# C023 rimc
dir='/Users/remy/lustre/naasc/users/rindebet/30dor/chevance/'
fluxfile=dir+'Dor13co_sT1.0_smoT2.0_cutT0.01_noiseT5.0.pb.ch0.trim.fits'
datafile=dir+"Dor13co_sT1.0_smoT2.0_cutT0.01_noiseT5.0.trim.fits"
root="dor.rimc.13co"
# 13 ~ 0.25*12 often.
cutoff=0.01
pkmin=0.025
dTleaf=0.002
dTleaf=0.004
pkmin=0.03

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
pyprops(datafile,fluxfile,clumpfile,root+"_"+parmstr,assignfile2=assignfile2,montecarlo=0,doplots=True)




