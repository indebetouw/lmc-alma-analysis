from astropy.io import fits
from astropy import units as u
import sys, os, pickle

# uses the new filfinder2D (from filfind_class::fil_finder_2D)
# also uses my version of filFinder, which has some bugs fixed:

gitpaths=['/Users/remy/local/github/lmc-alma-analysis/',
          '/Users/remy/local/github/FilFinder/',
          '/Users/remy/local/github/FilFinder/fil_finder']
for gitpath in gitpaths:
    if not gitpath in sys.path:
        sys.path.insert(0,gitpath)

from fil_finder import FilFinder2D
import pylab as pl
pl.ion()
import numpy as np
from scipy.signal import savgol_filter 
import pdb 

def run_filfinder(label='mycloud', mom8file=None, mom0file=None, redo=False,
                  distpc=4.8e4,xra=[0,2000],yra=[0,2000],glob_thresh=0.1, adapt_beams=14):
    
    plotfiledir = "fil_plots" # no trailing /
    # filpropfile = label+".filprops.pkl" # phys props
    
    fits_hdu = fits.open(mom8file)[0]    
    # setting beamwidth will have no effect if a fits header is given
    fils=FilFinder2D(fits_hdu, distance=distpc*u.pc)
    
    fils.save_name=plotfiledir+"/"+label
    if not os.path.exists(plotfiledir):
        os.mkdir(plotfiledir)
        redo=True
    
    fils.preprocess_image(flatten_percent=90)
    
    # GMC1 adapt=50pix, size=300pix2, smooth=8pix    ppbm=1.76"/0.5=3.53
    # so in bms, that's adapt=14, size=24, smooth=2.3
    # for gmc1, using adapt=10, size=15 gives 23 fils instead of 8 and too much chaff,
    # size=24, adapt=14, smooth=2.0: 9 fils instead of 8, more branches, and more vel perpendicular to fils.
    # size=24, adapt=14, smooth=3: 11fils, a bit fewer small fils.
    
    # 30dor: glob=72, flatten=100, adapt=75pix, size=800pix2  ppbm=0.29/0.032=9.1
    # so in bms, adapt=8.2, size=9.7 but bm elongated so narrow dir adapt=10
    
    # one can lower glob_thresh, get more small things, then lower adapt_thresh to
    # make them break up more, and finally let size_thresh kill them
    # in new FilFinder, glob_thresh is no longer a percentage its a value;
    # setting to 50% of the image max is too high; 20% is closer, 10% quite good

    bmwidth = pl.sqrt(fits_hdu.header['bmaj']*fits_hdu.header['bmin'])*3600*u.arcsec
    # 30dor tuned parameters
    #fils.create_mask(verbose=True, border_masking=False, 
    #                 use_existing_mask=False,save_png=True,
    #                 adapt_thresh=14*bmwidth, size_thresh=24*bmwidth**2,
    #                 smooth_size=3*bmwidth, glob_thresh=pl.nanmax(fits_hdu.data)*glob_thresh)
    fils.create_mask(verbose=True, border_masking=False, 
                     use_existing_mask=False,save_png=True,
                     adapt_thresh=adapt_beams*bmwidth, size_thresh=24*bmwidth**2,
                     smooth_size=3*bmwidth, glob_thresh=pl.nanmax(fits_hdu.data)*glob_thresh)
    # GMC104: decreasing adapt_thresh makes more filaments but also more arcs
    #fils.create_mask(verbose=True, border_masking=False, 
    #                 use_existing_mask=False,save_png=True,
    #                 adapt_thresh=5*bmwidth, size_thresh=24*bmwidth**2,
    #                 smooth_size=3*bmwidth, glob_thresh=pl.nanmax(fits_hdu.data)*glob_thresh)
    
    fils.medskel(verbose=True,save_png=True)
    
    print("calling analyze_skeletons")
    # 30dor:
    # fils.branch_thresh=200 # typical lengths are 1-100
    # raise from 10 to 30 to elminate more little things
    # fils.analyze_skeletons(verbose=True, save_png=True, skel_thresh=10,cubefile=cubefile)
    
    # mc3: old version was doing double filter instead, with branch_thresh=30;
    # branch_thresh says it overrides all previous settings
    fils.analyze_skeletons(verbose=True, save_png=True, prune_criteria="length",
                           branch_nbeam_lengths=5,nbeam_lengths=5, skel_thresh=10*u.pix)
    #relintens_thresh=0.1, branch_thresh=30,
    
    # the labeled filament arrays are left separated even though the branches
    # have been pruned, so we need to recalculate the labelled_filament_arrays
    # pix_identify is supposed to do that;
    # fils.filament_arrays comes from make_final_skeletons
    
    if redo:
        pl.clf()
        vmin = np.percentile(fils.flat_img[np.isfinite(fils.flat_img)], 20)
        vmax = np.percentile(fils.flat_img[np.isfinite(fils.flat_img)], 99)
        pl.imshow(fils.flat_img.value, interpolation=None, origin="lower",vmin=vmin, vmax=vmax*1.5, cmap='jet')
        pl.contour(fils.skeleton, colors='r')
        offs=fils.array_offsets
        for i in range(len(offs)):
            pl.text(0.5*(offs[i][0][1]+offs[i][1][1]),0.5*(offs[i][0][0]+offs[i][1][0]),str(i),color='orange')
        pl.xlim(xra)
        pl.ylim(yra)
        pl.savefig(plotfiledir+"/"+plotfiledir+".skel.png")
        

        # find_widths calculates the 2d distance transform, and uses self.image, self.imgscale
        # first widths on mom0
        fils.set_image(fits.getdata(mom0file))
        
        fils.save_name=plotfiledir+".mom0wid"
        fils.find_widths(verbose=True,save_png=True,max_dist=2.*u.pc)
        mom0width=fils.width_fits().copy()
        
        pl.clf()
        vmin = np.percentile(fils.image[np.isfinite(fils.image)], 20)
        vmax = np.percentile(fils.image[np.isfinite(fils.image)], 99)
        pl.imshow(fils.image.value, interpolation=None, origin="lower",vmin=vmin, vmax=vmax*1.5, cmap='jet')
        pl.contour(fils.skeleton, colors='r')
        offs=fils.array_offsets
        for i in range(len(offs)):
            pl.text(0.5*(offs[i][0][1]+offs[i][1][1]),0.5*(offs[i][0][0]+offs[i][1][0]),str(i),color='orange')
        pl.xlim(xra)
        pl.ylim(yra)
        pl.savefig(plotfiledir+"/"+plotfiledir+".int.skel.png")

        
        # now widths and orientations on the mom8 directly:
        fils.set_image(fits.getdata(mom8file))

        fils.save_name=plotfiledir+".mom8wid"
        fils.find_widths(verbose=True,save_png=True,max_dist=2.*u.pc) # CHANGE FOR 30DOR
        mom8width=fils.width_fits().copy()
        
        pl.clf()
        x =fils.converter.from_pixel(mom8width['fwhm'],u.pc).value
        dx=fils.converter.from_pixel(mom8width['fwhm_err'],u.pc).value
        y =fils.converter.from_pixel(mom0width['fwhm'],u.pc).value
        dy=fils.converter.from_pixel(mom0width['fwhm_err'],u.pc).value
        pl.errorbar(x,y,xerr=dx,yerr=dy,fmt='.')
        xmax=pl.nanmax([x,y])
        pl.xlim(0,xmax*1.1)
        pl.ylim(0,xmax*1.1)
        pl.xlabel("FWHM fit to mom8 [pc]")
        pl.ylabel("FWHM fit to mom0 [pc]")
        pl.plot([0,xmax],[0,xmax],':k')
        pl.savefig(plotfiledir+"/"+plotfiledir+".widths.mom8.mom0.png")
     
    
    return fils


if __name__ == "__main__":
    distpc=5.0e4
    #label="GMC1_12CO_12m7mT"
    #mom8file = "GMC1_12CO_12m7mTPF.maximum.fits"
    #mom0file = "GMC1_12CO_12m7mTPF.mom0.fits"
    #xra=[180,850]
    #yra=[120,650]
    #glob_thresh=0.1
    #run_filfinder(label=label, mom0file=mom0file, mom8file=mom8file, xra=xra, yra=yra, glob_thresh=glob_thresh, distpc=distpc, redo=True)

    label="GMC104_12CO_12m7mT"
    mom8file = "GMC104_12CO_12m7mTPF.maximum.fits"
    mom0file = "GMC104_12CO_12m7mTPF.integrated.fits"
    xra=[200,620]
    yra=[200,620]
    glob_thresh=0.15
    run_filfinder(label=label, mom0file=mom0file, mom8file=mom8file, xra=xra, yra=yra, glob_thresh=glob_thresh, distpc=distpc, redo=True)
