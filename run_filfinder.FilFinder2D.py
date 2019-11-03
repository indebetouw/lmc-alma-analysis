from astropy.io import fits
from astropy import units as u
import sys, os, pickle

# this version switches to filfinder2D (from filfind_class::fil_finder_2D)

gitpaths=['/Users/remy/lustre/naasc/users/rindebet/github/lmc-alma-analysis/']
for gitpath in gitpaths:
    if not gitpath in sys.path:
        sys.path.insert(0,gitpath)

from fil_finder import FilFinder2D
import pylab as pl
pl.ion()
import numpy as np
from scipy.signal import savgol_filter 
import pdb 

label="GMC1_12CO"
mom8file = "GMC1_12CO_12m.maximum.fits"
mom8file = "GMC1_12CO_12m7mTPF.maximum.fits"
cubefile = "GMC1_12CO_12m.image.fits"
#mom0file = "GMC1_12CO_12m.integrated.fits"  # use mom0pm5 instead?
distpc=4.8e4

#def run_filfinder(label='mycloud', cubefile=None, mom8file=None, mom0file=None, redo=False,
#distpc=4.8e4):
    
plotfiledir = label+".plots" # NO TRAILING / !
filpropfile = label+".filprops.pkl" # phys props

# initial fil finding in mom8 image:
fits_hdu = fits.open(mom8file)[0]
bmwidth=np.sqrt(fits_hdu.header['BMIN']*fits_hdu.header['BMAJ']) *3600 # arcsec

# setting beamwidth will have no effect if a fits header is given
# beamwidth=bmwidth*u.arcsec,
# it turns out it uses sigma not FWHM so later, nbeam=#sigma=1.6pix here
fils=FilFinder2D(fits_hdu, distance=distpc*u.pc)

# 30Dor: glob=72, flatten=100
# GMC1 1st time: 50, 85
fils.save_name=plotfiledir
if not os.path.exists(plotfiledir):
    os.mkdir(plotfiledir)
    redo=True

# give radial profiles more breathing room.
fils.skeleton_pad_size=5

fils.preprocess_image(flatten_percent=90)

# 30dor: adapt=75, size=800
# GMC1 first time: 50,500
# can lower glob_thresh, get more small things, then lower adapt_thresh to
# make them break up more, and finally let size_thresh kill them
fils.create_mask(verbose=True, border_masking=False, adapt_thresh=55*u.pix, size_thresh=300*u.pix**2, use_existing_mask=False,save_png=True,smooth_size=2*u.pix, glob_thresh=50)

#Adaptive thresholding patch is larger than 40pixels. Regridding has been disabld.
#  warnings.warn("Adaptive thresholding patch is larger than 40"

print "calling medskel"
fils.medskel(verbose=True,save_png=True)

print "getting cube"
cube=fits.getdata(cubefile)[0]
cube[np.isnan(cube)]=0

print "calling analyze_skeletons"
# 30dor:
# fils.branch_thresh=200 # typical lengths are 1-100
# raise from 10 to 30 to elminate more little things
# fils.analyze_skeletons(verbose=True, save_png=True, skel_thresh=10,cubefile=cubefile)

# mc3: double filter instead, with branch_thresh=30; branch_thresh says it overrides all previous settings
fils.analyze_skeletons(verbose=True, save_png=True, prune_criteria="length",max_iter=1,
                       #cubefile=cubefile,
                       branch_nbeam_lengths=5,nbeam_lengths=3, skel_thresh=10)
#relintens_thresh=0.1, branch_thresh=30,

pl.clf()
pl.imshow(fils.skeleton,origin="bottom")
pl.xlim(200,850)
pl.ylim(100,650)
pl.draw()
import pdb
pdb.set_trace()

fils.analyze_skeletons(verbose=True, save_png=True, cubefile=cubefile, branch_thresh=30)

# the labeled filament arrays are left separated even though the branches
# have been pruned, so we need to recalculate the labelled_filament_arrays
# pix_identify is supposed to do that;
# fils.filament_arrays comes from make_final_skeletons

if False:
    if redo:
        pl.clf()
        vmin = np.percentile(fils.flat_img[np.isfinite(fils.flat_img)], 20)
        vmax = np.percentile(fils.flat_img[np.isfinite(fils.flat_img)], 99)
        pl.imshow(fils.flat_img, interpolation=None, origin="lower",
                 vmin=vmin, vmax=vmax*1.5, cmap='jet')
        pl.contour(fils.skeleton_nopad, colors='r')
        offs=fils.array_offsets
        for i in range(len(offs)):
            pl.text(0.5*(offs[i][0][1]+offs[i][1][1]),0.5*(offs[i][0][0]+offs[i][1][0]),str(i),color='orange')
        pl.savefig(plotfiledir+"/"+plotfiledir+".skel.png")
    
        #if rerun:
        pl.clf()
        for n in range(len(fils.branch_properties['length_2d'])):
            z=np.where(np.array(fils.branch_properties['length'][n])>0.5)[0]    
            pl.plot(np.array(fils.branch_properties['length_2d'][n])[z],
                   np.array(fils.branch_properties['length'][n])[z]/
                   np.array(fils.branch_properties['length_2d'][n])[z],'o',label=str(n))
    
        pl.xlabel("2D length")
        pl.ylabel("3D/2D length")
        pl.legend(prop={"size":8})
        pl.xscale("log")
        pl.yscale("log")
        pl.savefig(plotfiledir+"/"+plotfiledir+".lengths_3d_rat2d.png")

        # find_widths calculates the 2d distance transform, and uses self.image, 
        # self.imgscale
        fils.image=fits.getdata(mom0file)
    
        fils.save_name=plotfiledir+".mom0wid"
        fils.find_widths(verbose=True,save_png=True,pad_to_distance=0.25,max_distance=0.5)
        mom0width=fils.width_fits.copy()

        pl.clf()
        vmin = np.percentile(fils.image[np.isfinite(fils.image)], 20)
        vmax = np.percentile(fils.image[np.isfinite(fils.image)], 99)
        pl.imshow(fils.image, interpolation=None, origin="lower",
                 vmin=vmin, vmax=vmax*1.5, cmap='jet')
        pl.contour(fils.skeleton_nopad, colors='r')
        offs=fils.array_offsets
        for i in range(len(offs)):
            pl.text(0.5*(offs[i][0][1]+offs[i][1][1]),0.5*(offs[i][0][0]+offs[i][1][0]),str(i),color='orange')
        pl.savefig(plotfiledir+"/"+plotfiledir+".int.skel.png")

        # widths and orientations on the mom8 directly:
        fils.image=fits.getdata(mom8file)
        fils.save_name=plotfiledir+".mom8wid"
        fils.find_widths(verbose=True,save_png=True,pad_to_distance=0.25,max_distance=0.5)
    
        mom8width=fils.width_fits.copy()
    
        pl.clf()
        x=mom0width['Parameters'][:,1]
        dx=mom0width['Errors'][:,1]
        y=mom8width['Parameters'][:,1]
        dy=mom8width['Errors'][:,1]
        pl.plot(x,y,'yo')
    
        #ok=  [0,4,6,12,13,15,18,20,21,22,25,27,28,30,31,34]
        #good=[1,5,9,10,14,16,26,32,33,35]
        #pl.plot(x[good],y[good],'ro')
        #pl.plot(x[ok],y[ok],'mo')    
        #for i in good:
        #    if not (np.isnan(x[i]) or np.isnan(y[i])):
        #        pl.text(x[i],y[i],str(i))
    
        pl.savefig(plotfiledir+"/"+plotfiledir+".widths.mom8.mom0.png")
    
        pl.clf()
        x=mom0width['Parameters'][:,1]
        dx=mom0width['Errors'][:,1]
        y=mom8width['Parameters'][:,1]/x
        dy=y*np.sqrt( (mom8width['Errors'][:,1]/mom8width['Parameters'][:,1])**2+
                      (mom0width['Errors'][:,1]/mom0width['Parameters'][:,1])**2 )
        
        pl.errorbar(x,y,xerr=dx,yerr=dy,fmt='o')
        pl.xlabel("width of mom0")
        pl.ylabel("mom8 width / mom0 width")
        pl.plot(pl.xlim(),[1,1],'k',linestyle='dotted')
    
        pl.savefig(plotfiledir+"/"+plotfiledir+".widtherrs.mom8.mom0.png")

#    return fils


if __name__ == "__main__":
    run_filfinder(label=sys.argv[1], cubefile=sys.argv[2], mom0file=sys.argv[3], mom8file=sys.pargv[4])
