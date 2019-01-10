import sys
gitpaths=['/Users/remy/lustre/naasc/users/rindebet/github/quickclump/',
          '/Users/remy/lustre/naasc/users/rindebet/github/pyprops/']
for gitpath in gitpaths:
    if not gitpath in sys.path:
        sys.path.insert(0,gitpath)


dir='/Users/remy/lustre/naasc/users/rindebet/30dor/chevance/'
fluxfile=dir+'dor12coS_r1.5.pb.fapex.ch0.trim.fits'
datafile=dir+"dor12coS_r1.5.fapex.trim.pbcor.fits"
root="dor.rimc.12co"
