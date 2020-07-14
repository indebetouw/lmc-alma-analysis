nsigma=3.
min_delta=2.5
min_bms=2.
redo=False
#---------

region="PerA"
line="13co"
distpc=250.
nsigma=5. # raise this for COMPLETE data
min_delta=3.
min_bms=4.

# region="OphA"
# line="13co"
# distpc=131.
# min_bms=6. # Oph is closer, eliminate tiny things




# only used for input data files:
root="trim"
completename="FCRAO_F_xyv."+root

#-------------
gitdirs=["/Users/remy/local/github/astrodendro","/Users/remy/local/github/lmc-alma-analysis"]

import sys,os
for adir in gitdirs:
    if adir not in sys.path:
        sys.path.append(adir)

from astropy.io import fits
from astropy import stats
import numpy as np


newcubefile=region+'_'+line+completename+'.masked.fits'
if not os.path.exists(newcubefile):
    newcubefile=region+'_'+line+completename+'.fits'

#-------------
# mom0 and mom8 files based on unmaskedcubefile
cubefile=region+'_'+line+completename+'.fits'

if not (os.path.exists(region+'_'+line+completename+'.mom8.fits') and
        os.path.exists(region+'_'+line+completename+'.mom0.fits')):
        
    cubehdu=fits.open(cubefile)[0]
    cubedata=cubehdu.data
    sigma=stats.mad_std(cubedata[~np.isnan(cubedata)])
    mom8=np.nanmax(cubedata*(cubedata>sigma),axis=0)
    mom0=np.nansum(cubedata*(cubedata>sigma),axis=0)
    dv=cubehdu.header['cdelt3']
    if cubehdu.header['cunit3'] == 'm/s':
        dv=dv/1000
    elif cubehdu.header['cunit3'] != 'km/s':
        print("Error can't handle cube vel unit "+cubehdu.header['cunit3'])
        print("Units of mom0,8 images will be wrong")
    cubehdu.data=mom8
    cubehdu.writeto(region+'_'+line+completename+'.mom8.fits')
    cubehdu.header['bunit'] += '*km/s'
    cubehdu.data=mom0*dv
    cubehdu.writeto(region+'_'+line+completename+'.mom0.fits')


if "masked" in newcubefile and not os.path.exists(region+'_'+line+completename+'.masked.mom8.fits'):
        
    cubehdu=fits.open(newcubefile)[0]
    cubedata=cubehdu.data
    sigma=stats.mad_std(cubedata[~np.isnan(cubedata)])
    mom8=np.nanmax(cubedata*(cubedata>sigma),axis=0)
    mom0=np.nansum(cubedata*(cubedata>sigma),axis=0)
    dv=cubehdu.header['cdelt3']
    if cubehdu.header['cunit3'] == 'm/s':
        dv=dv/1000
    elif cubehdu.header['cunit3'] != 'km/s':
        print("Error can't handle cube vel unit "+cubehdu.header['cunit3'])
        print("Units of mom0,8 images will be wrong")
    cubehdu.data=mom8
    cubehdu.writeto(region+'_'+line+completename+'.masked.mom8.fits')
    cubehdu.header['bunit'] += '*km/s'
    cubehdu.data=mom0*dv
    cubehdu.writeto(region+'_'+line+completename+'.masked.mom0.fits')


#-------------

dendrofile=region+"_"+line+"_dendrogram_rmsmap.hdf5"
try:
    if len(mydendro)>1 and dloaded==dendrofile:
        dloaded=dendrofile
except:
    dloaded=None
if os.path.exists(dendrofile) and dloaded!=dendrofile and not redo:
    from astrodendro import Dendrogram
    print("loading "+dendrofile)
    mydendro = Dendrogram.load_from(dendrofile)
    dloaded=dendrofile

if redo:
    os.system("rm "+dendrofile)

print()
if redo or not (os.path.exists(region+"_"+line+'_full_catalog.txt') and os.path.exists(region+"_"+line+'_full_catalog_clipped.txt')):
    from run_dendro import *
    print("running dendro and ppvcat > "+region+"_"+line+'_full_catalog.txt')
    if dloaded==dendrofile:
        run_dendro(label=region+"_"+line,
                   cubefile=newcubefile,
                   flatfile=region+"_"+line+completename+".mom8.fits",
                   position_dependent_noise=True,
                   nsigma=nsigma, min_delta=min_delta,min_bms=min_bms,
                   dendro_in=mydendro)
    else:
        run_dendro(label=region+"_"+line,
                   cubefile=newcubefile,
                   flatfile=region+"_"+line+completename+".mom8.fits",
                   nsigma=nsigma, min_delta=min_delta,min_bms=min_bms,
                   position_dependent_noise=True)
    


#-------------
print()
n13file=region+"_peak_n13cube.fits"

if not os.path.exists(n13file):
    if os.path.exists(n13file+".gz"):
        os.system("gunzip *gz")

for xline in ["12co","13co"]:
    if not (os.path.exists(region+"_"+xline+completename+".rms.fits") and redo==False):
        cubedata=fits.getdata(region+"_"+xline+completename+".fits")
        sigma = stats.mad_std(cubedata[~np.isnan(cubedata)])
        mask3d = cubedata<3*sigma
        rmsmap = np.nanstd(cubedata*mask3d,axis=0) # assumes spectral 1st
        mom8hdu=fits.open(region+"_"+xline+completename+".mom8.fits")[0]
        mom8hdu.data=rmsmap
        mom8hdu.writeto(region+"_"+xline+completename+".rms.fits")

if not (os.path.exists(region+"_12co"+completename+".mask.fits") and redo==False):
    cubehdu = fits.open(region+"_12co"+completename+".fits")[0]
    cubedata = cubehdu.data
    sigma = stats.mad_std(cubedata[~np.isnan(cubedata)])
    mask3d = cubedata>1.*sigma
    cubehdu.data=np.int16(mask3d)
    cubehdu.writeto(region+"_12co"+completename+".mask.fits")

if not (os.path.exists(n13file) and redo==False):
    from lte import *
    print("calculating "+n13file)
    lte(files=[region+"_12co"+completename+".fits",
               region+"_13co"+completename+".fits",
               region+"_12co"+completename+".rms.fits",
               region+"_13co"+completename+".rms.fits",
               region+"_12co"+completename+".mask.fits"],datainfo=region)
    os.system("gunzip *gz")

    
#--------------------
# calc_phys_props
print()

clipstr=""
ptab=region+"_"+line+'_physprop'+clipstr+'.txt'
boot_iter=10 # test
#boot_iter=100
boot_iter=400 # production

if not (os.path.exists(ptab) and redo==False):
    from calc_phys_props import *
    print("running calc_phys_props > "+ptab)
    calc_phys_props(label=region+"_"+line,
                    cubefile=region+"_"+line+completename+".fits",
                    dendrofile=region+"_"+line+"_dendrogram_rmsmap.hdf5",
                    boot_iter=boot_iter, efloor=0,
                    alphascale=1, distpc=distpc,
                    copbcor=None, conoise=None, ancfile=None, anclabel=None,
                    verbose=False, clipping=False, co13toh2=1e6,
                    n13cube=region+"_peak_n13cube.fits",
                    n13errcube=region+"_peak_n13cubeerr.fits",
                    dendro_in=mydendro)

#--------------------
# filaments

if not (os.path.exists(region+".filprops.pkl") and redo==False):
    from run_filfinder import *
    fils=run_filfinder(label=region,
                  mom8file=region+"_"+line+completename+".masked.mom8.fits",
                  mom0file=region+"_"+line+completename+".masked.mom0.fits",
                  redo=redo,distpc=distpc,xra=[0,1100],yra=[0,800],glob_thresh=0.2, adapt_beams=25)
    # how make pickle of phys props?

# TODO for complete, lower glob_thresh from 0.1 increase adapt_beams from 14
# TODO in run_filfinder make xra,yra size of image
