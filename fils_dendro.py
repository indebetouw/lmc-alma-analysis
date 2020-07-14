distpc=4.8e4
label="GMC1_12CO_12m7mT"
mom8file = label+"PF.maximum.fits"
mom0file = label+"PF.mom0.fits"
mom1file = label+"PF_dil.mom1.fits"
xra=[180,850]
yra=[120,650]
glob_thresh=0.1


label="Dor_13CO"
cubefile="both.13co.dv02.p05.multiscale.cubic.cbm.trim.0.375x0.25.hanning.gridto12.rotate2.fits"
mom0file="both.13co.dv02.p05.multiscale.cubic.cbm.trim.0.375x0.25.hanning.gridto12.rotate2.integrated.fits"
mom8file="both.13co.dv02.p05.multiscale.cubic.cbm.trim.0.375x0.25.hanning.gridto12.rotate2.maximum.fits"
mom1file="both.13co.dv02.p05.multiscale.cubic.cbm.trim.0.375x0.25.hanning.gridto12.rotate2.mom1.gt0.02.fits"
xra=[0,2100]
yra=[0,800]

label="allDor_13CO"
cubefile="Dor13co_sT1.0_smoT2.0_cutT0.01_noiseT5.0.trim.fits"
mom0file="dor12coS_r1.5.fapex.gtr3mJy.integrated.trim.fits"
mom8file="dor12coS_r1.5.fapex.gtr3mJy.maximum.trim.fits"
mom1file="both.13co.dv02.p05.multiscale.cubic.cbm.trim.0.375x0.25.hanning.gridto12.rotate2.mom1.gt0.02.fits"
xra=[0,2100]
yra=[0,800]


#imsubimage("dor12coS_r1.5.fapex.fits",region='box[[250pix,350pix],[1700pix,1750pix]]',outfile="dor12coS_r1.5.fapex.trim")
#imsubimage("Dor13co_sT1.0_smoT2.0_cutT0.01_noiseT5.0.fits",region='box[[250pix,350pix],[1700pix,1750pix]]',outfile="Dor13co_sT1.0_smoT2.0_cutT0.01_noiseT5.0.trim")


#===============================
label="allDor_12CO.new"
cubefile="dor12coS_r1.5.fapex.trim.pbgt0.3.fits"
mom0file="dor12coS_r1.5.fapex.trim.pbgt0.3.gtr3mJy.integrated.fits"
mom8file="dor12coS_r1.5.fapex.trim.pbgt0.3.gtr3mJy.maximum.fits"
mom1file="dor12coS_r1.5.fapex.trim.pbgt0.3.gtr10mJy.mom1.fits"
xra=[150,1400]
yra=[50,1150]
glob_thresh=0.17
#0.15 still gets a little chaff on the top/bottom; 0.2 misses fainter fils
zelong_tmax=15
zelong_area=[.02,.5]  # area_pc
boot_iter=10  # 10 for quick and dirty, 400 for tony/proper
 

##===============================
#label="allDor_13CO.new"
#cubefile="Dor13co_sT1.0_smoT2.0_cutT0.01_noiseT5.0.trim.fits"
#mom0file="Dor13co_sT1.0_smoT2.0_cutT0.01_noiseT5.0.trim.gtr3mJy.integrated.fits"
#mom8file="Dor13co_sT1.0_smoT2.0_cutT0.01_noiseT5.0.trim.gtr3mJy.maximum.fits"
#mom1file="Dor13co_sT1.0_smoT2.0_cutT0.01_noiseT5.0.trim.gtr10mJy.mom1.fits"
#xra=[50,1400]
#yra=[50,1300]
#glob_thresh=0.07
#zelong_tmax=1.
#zelong_area=[.02,1]
#boot_iter=10




# can use fils.filaments[i].median_brightness(image) to get brtness in new image
# ff.ridge_profile(fils.image)

import pylab as pl
pl.ion()
pl.clf()

from astropy import units as u
from astropy import constants as const
from astropy.table import Table, Column
from copy import copy

import sys,os
gitpaths=['/Users/remy/lustre/naasc/users/rindebet/github-local/lmc-alma-analysis/','/Users/remy/lustre/naasc/users/rindebet/github-local/FilFinder/','/Users/remy/lustre/naasc/users/rindebet/github-local/FilFinder/fil_finder/']
for gitpath in gitpaths:
    if not gitpath in sys.path:
        sys.path.insert(0,gitpath)

# this uses my version of filfinder
from run_filfinder import run_filfinder
from fil_finder.pixel_ident import pix_identify

# unfortunately I haven't successfully save/restored fils, so have to regenerate
try:
    dir(fils)
except:
    fils=run_filfinder(label=label, mom0file=mom0file, mom8file=mom8file, xra=xra, yra=yra, glob_thresh=glob_thresh, distpc=distpc, redo=True)




# read in dendro props
from astrodendro import Dendrogram, ppv_catalog, analysis

try:
    x=len(d)
except:
    x=0
if x<=0:
    import os
    if not os.path.exists(label+'_dendrogram.hdf5'):
        from run_dendro import run_dendro
        run_dendro(criteria=['volume'], label=label, cubefile=cubefile, mom0file=mom0file, nsigma=5.)
    d = Dendrogram.load_from(label+'_dendrogram.hdf5')


cat = Table.read(label+'_full_catalog.txt', format='ascii.ecsv')

if not os.path.exists(label+'_physprop.txt'):
    from calc_phys_props import *
#    calc_phys_props(label=label, cubefile=cubefile, boot_iter=400, efloor=0,
# quick path: use 10 iterations for uncerts
    calc_phys_props(label=label, cubefile=cubefile, boot_iter=boot_iter, efloor=0,
                    alphascale=1, distpc=4.8e4, copbcor=None, conoise=None, ancfile=None, anclabel=None, verbose=True)

pcat = Table.read(label+'_physprop.txt', format='ascii.ecsv')




if False:
    def plotleaves(s,cat,xfield,yfield,color='k'):
        x0=cat[xfield][s.idx]
        y0=cat[yfield][s.idx]
        for c in s.children:
            pl.plot([x0,cat[xfield][c.idx]],[y0,cat[yfield][c.idx]],color=color)
            plotleaves(c,cat,xfield,yfield,color=color)
    
    for t in d.trunk:
        xfield='area_pc2'
        yfield='axratio'
        x0=pcat[xfield][t.idx]
        y0=pcat[yfield][t.idx]
        myplot,=pl.plot(x0,y0,'o')
        plotleaves(t,pcat,xfield,yfield,color=myplot.get_color())
    pl.xlabel(xfield)
    pl.ylabel(yfield)
    pl.savefig("plots/"+label+"treeplot."+xfield+"_"+yfield+".png")
    
    #pl.xscale("log")


    
from astropy.io import fits 
mom0=fits.getdata(mom0file)

# select elongated dendro branches - these are indices in the dendro structure.
#zelong=pl.where((pcat['axratio']<0.3)*
#               (pcat['area_pc2']>3)*(pcat['area_pc2']<40))[0]
zelong=pl.where((pcat['axratio']<0.4)*
               (pcat['area_pc2']>3)*(pcat['area_pc2']<40))[0]


# for 30dor this needs to be a lot smaller dendros; 
# for 12CO also raise the peak to 7 sigma = .035 Jy/bm
# Jy/bm = .00027 *K
zelong=pl.where((pcat['axratio']<0.4)*(cat['tmax']>zelong_tmax)*
               (pcat['area_pc2']>zelong_area[0])*(pcat['area_pc2']<zelong_area[1]))[0]


# keep only lowest level of each selected branch
nz=len(zelong)
keep=pl.ones(nz,dtype=bool)
for i in range(nz):
    if keep[i]:
        t=d[zelong[i]]
        for k in t.descendants:
            if k.idx in zelong:
                keep[i]=False

zelong=zelong[keep]
n_elong_dendro=len(zelong)

# this will be for the list of fil,branch pairs that overlap this dendro, and by how many pixels:
overlapping_fil_branches=pl.repeat({},n_elong_dendro)
for i in range(n_elong_dendro):
    overlapping_fil_branches[i]={"dendro_id":0,"filbranches":[],"npixoverlap":[]}.copy()
                                   




#================================================================
mom8=fits.getdata(mom8file)
pl.clf()
pl.imshow((pl.nanmax(mom8)-mom8),origin="bottom",cmap="gray")

pl.subplots_adjust(left=0.05,right=0.98,bottom=0.05,top=0.98)
pl.xlim(xra)
pl.ylim(yra)
pl.xticks([])
pl.yticks([])
pl.contour(fils.skeleton, colors='c',linewidths=1,vmin=0,vmax=1,levels=[0.5])

for ff in fils.filaments:
    for b in ff.end_pts:
        pl.plot(b[1],b[0],'ob',markersize=2)
    if len(ff.intersec_pts)>0:
        for b in pl.concatenate(ff.intersec_pts):
            pl.plot(b[1],b[0],'ob',markersize=2)


pl.savefig(label+".plots/"+label+".fils_branches.png")




            


    
#================================================================

mom1=fits.getdata(mom1file)
ly, lx = mom1.shape
x, y = range(0, lx), range(0,ly)
xi, yi = pl.meshgrid(x, y)


wid=9
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
kernel = Gaussian2DKernel(stddev=wid/2.354)
vsmooth = convolve(mom1, kernel)

vgrad = pl.gradient(vsmooth)
vsmooth[pl.where(pl.isnan(mom1))]=pl.nan

for i in [0,1]:
    vgrad[i][pl.where(pl.isnan(mom1))]=pl.nan
    vgrad[i][pl.where(pl.absolute(vgrad[i])>0.1)]=pl.nan

vv=(vgrad[0]**2+vgrad[1]**2)**0.6

#pl.streamplot(xi, yi, vgrad[0], vgrad[1])
s=7 # for MC
s=10 # for dor - finer pixellation
skip = (slice(None, None, s), slice(None, None, s))
pl.clf()
pl.quiver(xi[skip], yi[skip], vgrad[0][skip], vgrad[1][skip], vv[skip],scale=3, cmap="jet",pivot="middle")
pl.subplots_adjust(left=0.05,right=0.98,bottom=0.05,top=0.98)
pl.xlim(xra)
pl.ylim(yra)
pl.contour(fils.skeleton, colors='k',linewidths=1,vmin=0,vmax=1,levels=[0.5])
pl.xticks([])
pl.yticks([])
#pl.contour(mom8,levels=[0.1],colors='gray',linewidths=1)



pl.savefig("fil_plots/"+label+".fils_velfield.png")





#================================================================
# show the image with the elongated dendros and the filaments, and then the overlapping filaments as we find them below


pl.figure(1)
pl.clf()
pl.imshow((pl.nanmax(mom0)-mom0)**3,origin="bottom",cmap="gray")
pl.subplots_adjust(left=0.05,right=0.98,bottom=0.05,top=0.98)
pl.xlim(xra)
pl.ylim(yra)

pl.figure(2)
pl.clf()
pl.imshow((pl.nanmax(mom0)-mom0)**3,origin="bottom",cmap="gray")
pl.subplots_adjust(left=0.05,right=0.98,bottom=0.05,top=0.98)
pl.xlim(xra)
pl.ylim(yra)

import pdb

plotted=pl.zeros(len(pcat),dtype=bool)
for i in range(n_elong_dendro):
    overlapping_fil_branches[i]["dendro_id"]=zelong[i]
    if not plotted[zelong[i]]:
        t=d[zelong[i]]
        mask2d=pl.byte(t.get_mask().max(axis=0))
        #mycont=pl.contour(mask2d,1,levels=[0.1+0.2*t.level],vmin=0,vmax=1,cmap="plasma")
        #mycont=pl.contour(mask2d,1,levels=[0.1],vmin=0,vmax=1,cmap="plasma")
        pl.figure(1)
        mycont=pl.contour(mask2d,1,levels=[0.1],vmin=0,vmax=1,colors='g')
        pl.figure(2)
        mycont=pl.contour(mask2d,1,levels=[0.1],vmin=0,vmax=1,colors='g')
        plotted[zelong[i]]=True
        zyx=t.indices()
        # only on plot 2:
        pl.text(zyx[2].mean(),zyx[1].mean(),t.idx)
        #print t.idx,t.level

pl.figure(1)
pl.contour(fils.skeleton, colors='c',linewidths=1,vmin=0,vmax=1,levels=[0.5])
pl.xticks([])
pl.yticks([])
pl.figure(2)
pl.contour(fils.skeleton, colors='c',linewidths=1,vmin=0,vmax=1,levels=[0.5])
pl.xticks([])
pl.yticks([])

# need explicit X,Y coords to contour subarrays on top of full image,
# and to match to dendro structures
s=mom0.shape
XX=pl.arange(s[1])
YY=pl.arange(s[0])



# match each fil to (could be several) elongated dendro structs
associated_dendros=[]

labeled_filament_arrays=[ff._labeled_mask for ff in fils.filaments]
# this will be the skeleton truncated within dendros
overlapping_arrays=[]
# and this the skeleton with all branches that go into a dendro
overlapping_xarrays=[]

for ifil in range(len(labeled_filament_arrays)):
    off=fils.array_offsets[ifil]
    lflarr=labeled_filament_arrays[ifil]
    farr=fils.filaments[ifil].skeleton() # not padded apparently
    sfarr=farr.shape

    # this will be the skeleton truncated within dendros
    farr_overlap=pl.zeros(sfarr,dtype=int)
    # and this the skeleton with all branches that go into a dendro
    xfarr_overlap=pl.zeros(sfarr,dtype=int)

    # first determine how much each elongated dendro struct overlaps with this filament
    pixoverlap=pl.zeros(n_elong_dendro,dtype=int)

    # and how much each branch of the fil overlaps each elongated dendro
    nbranches=fils.filaments[ifil].branch_properties['number']
    pixoverlapbranches=pl.zeros([nbranches,n_elong_dendro],dtype=int)

    for iz in range(n_elong_dendro):
        # dendro indices relative to fil. subarray
        drelx=d[zelong[iz]].indices()[2]-off[0][1]
        drely=d[zelong[iz]].indices()[1]-off[0][0]
        zz=pl.where( (drelx>=0)*(drely>=0)*(drelx<sfarr[1])*(drely<sfarr[0]) )[0]
        if len(zz)>0:
            pixoverlap[iz]=farr[drely[zz],drelx[zz]].sum()
            farr_overlap[drely[zz],drelx[zz]]=farr[drely[zz],drelx[zz]]

    zorder=pl.argsort(pixoverlap)[::-1]
    # now if any dendros overlap this entire fil, 
    if pixoverlap[zorder[0]]>0:
        zoverlap=pl.where(pixoverlap[zorder]>0)[0]
        # remember which dendros overlap:
        associated_dendros.append({"indices":zelong[zorder[zoverlap]],"pixoverlap":pixoverlap[zorder[zoverlap]]})
        pl.figure(2)
        pl.text(off[0][1],off[0][0],ifil,color="m")

        # go through overlapping dendros and analyze overlap w/subbranches from _labeled_mask
        thisbranchlist=[]

        # TODO add longest path from subbranch at least to a hub also?


        # go through each overlapping dendro:
        overlapthreshold=20 # >1 to only count the significantly overlapping branches
        for izo in zorder[zoverlap]:

            t=d[zelong[izo]]
            mask2d=pl.byte(t.get_mask().max(axis=0))
            pl.figure(1)
            mycont=pl.contour(mask2d,1,levels=[0.1],vmin=0,vmax=1,cmap="plasma")
            pl.figure(2)
            mycont=pl.contour(mask2d,1,levels=[0.1],vmin=0,vmax=1,cmap="plasma")
            # 0.6: salmon-colored;  0.1=purple

            
            # again the dendro indices relative to the fil subarray:
            drelx=d[zelong[izo]].indices()[2]-off[0][1]
            drely=d[zelong[izo]].indices()[1]-off[0][0]
            zz=pl.where( (drelx>=0)*(drely>=0)*(drelx<sfarr[1])*(drely<sfarr[0]) )[0]
            for ibr in pl.arange(nbranches)+1:
                mask = lflarr==ibr
                pixoverlapbranches[ibr-1,izo] = mask[drely[zz],drelx[zz]].sum()
                if pixoverlapbranches[ibr-1,izo]>overlapthreshold:  
                    thisbranchlist.append(ibr)
                    xfarr_overlap[pl.where(mask)]=1  # could keep labels by setting these pixels of xfarr_overlap to ibr
            # for each dendro, keep a list of all branches which overlap >thresh, and what the actual overlap is; later, can filter based on the ratio of the overlap to the fil length
            zz=pl.where(pixoverlapbranches[:,izo]>overlapthreshold)[0]
            if len(zz)>0:
                # keep the actual branch index, labeled starting from 1 not 0:
                #overlapping_fil_branches[izo]["filbranches"].append([[ifil,iz+1] for iz in zz])
#                overlapping_fil_branches[izo]["npixoverlap"].append(pixoverlapbranches[zz,izo])
                if len(overlapping_fil_branches[izo]["filbranches"])>0:
                    overlapping_fil_branches[izo]["filbranches"]=pl.append(overlapping_fil_branches[izo]["filbranches"],[[ifil,iz+1] for iz in zz],axis=0)
                else:
                    overlapping_fil_branches[izo]["filbranches"]=[[ifil,iz+1] for iz in zz]
                overlapping_fil_branches[izo]["npixoverlap"]=pl.append(overlapping_fil_branches[izo]["npixoverlap"],pixoverlapbranches[zz,izo])
        
        # TODO: enumerate which branches overlap which dendro, so that below, I can associate one dendro with each branch;  maybe decide based on whcih branch has the greatest # overlap pixels?  does the associateion have beo unique or can several branches have the same associated dendro?  i think it'd be better to be unique.

        associated_dendros[-1]["overlapbranches"]=thisbranchlist

        # interpts, hubs, ends, filbranches, labeled_fil_arrays =  pix_identify([farr_overlap], 1)

        # now we could rerun the fil dissection and graph analysis
        # on these overlapping arrays - including separation into
        # separated skels.

    else:
        associated_dendros.append({"indices":[],"npixoverlap":[],
                                   "overlapbranches":[]})

    pl.figure(1)
    pl.contour(XX[off[0][1]:off[1][1]+1],YY[off[0][0]:off[1][0]+1],xfarr_overlap,1,levels=[0.5],colors='r',vmin=0,vmax=1,linewidths=1)
    pl.figure(2)
    pl.contour(XX[off[0][1]:off[1][1]+1],YY[off[0][0]:off[1][0]+1],xfarr_overlap,1,levels=[0.5],colors='r',vmin=0,vmax=1,linewidths=1)
    overlapping_arrays.append(farr_overlap)
    overlapping_xarrays.append(xfarr_overlap)



pl.figure(2)
pl.savefig(label+".plots/"+label+".fils_dendro_labels.png")
pl.figure(1)
pl.savefig(label+".plots/"+label+".fils_dendro.png")





#-------------
# orientation of velocity gradient, relative to the fil. direction, for each branch

mom1=fits.getdata(mom1file)
fils.exec_rht(verbose=True,branches=True,gradimage=mom1)

goodorients=[]
goodmedmom1=[]
goodrmsmom1=[]
pl.clf()
pl.imshow(mom1-pl.nanmean(mom1),origin="bottom",cmap="jet")
pl.xticks([])
pl.yticks([])
pl.subplots_adjust(left=0.05,right=0.98,bottom=0.05,top=0.98)
pl.xlim(xra)
pl.ylim(yra)
cbar=pl.colorbar(ticks=5*pl.frange(5)+240-pl.nanmean(mom1))
cbar.ax.set_yticklabels(["%i km/s"%i for i in 240+5*pl.frange(5)])

pl.contour(fils.skeleton, colors='c',linewidths=1,vmin=0,vmax=1,levels=[0.5])

pad_size = 1
from fil_finder.utilities import pad_image

plotlabels=False
for i in range(len(overlapping_fil_branches)):
    for fbranch in overlapping_fil_branches[i]['filbranches']:
        off=fils.array_offsets[fbranch[0]]
        pl.contour(XX[off[0][1]:off[1][1]+1],YY[off[0][0]:off[1][0]+1],overlapping_xarrays[fbranch[0]],1,levels=[0.5],colors='r',vmin=0,vmax=1,linewidths=1)
        xy=fils.filaments[fbranch[0]].branch_properties['pixels'][fbranch[1]-1]
        if plotlabels: pl.text(XX[off[0][1]]+xy[:,1].mean(),YY[off[0][0]]+xy[:,0].mean(),"%i"%(fils.orientation_branches[fbranch[0]][fbranch[1]-1].value*180/pl.pi))
        #pl.text(XX[off[0][1]]+xy[0,1],YY[off[0][0]]+xy[0,0],"%i"%(fils.orientation_branches[fbranch[0]][fbranch[1]-1].value*180/pl.pi))
        goodorients.append(fils.orientation_branches[fbranch[0]][fbranch[1]-1].value*180/pl.pi)
        goodmedmom1.append(fils.filaments[fbranch[0]].median_brightness(mom1,branchid=fbranch[1]))

        mom1pad = pad_image(mom1, fils.filaments[fbranch[0]].pixel_extents, pad_size)
        skels = fils.filaments[fbranch[0]].skeleton(pad_size=pad_size, out_type="branch", branchid=fbranch[1])
        if mom1pad.shape != skels.shape:
            mom1pad = fils.filaments[fbranch[0]].image_slicer(mom1pad, skels.shape,pad_size=pad_size)
        assert mom1pad.shape == skels.shape
        goodrmsmom1.append(pl.nanstd(mom1pad[skels]))


pl.savefig(label+".plots/"+label+".fils_dendro_angle_mom1.png")


#>>>>> TODO: do some kind of sorting by the mean value of the gradient image,
# or better yet the RMS of the mom1, to highlight those with *large* vel gradients

rmsmom1=[]
for ff in fils.filaments:
    mom1pad = pad_image(mom1, ff.pixel_extents, pad_size)
    for i in range(len(ff.branch_properties['length'])):
        skels = ff.skeleton(pad_size=pad_size, out_type="branch", branchid=i)
        if mom1pad.shape != skels.shape:
            mom1pad = ff.image_slicer(mom1pad, skels.shape,pad_size=pad_size)
        assert mom1pad.shape == skels.shape
        rmsmom1.append(pl.nanstd(mom1pad[skels]))





orients=pl.concatenate([k.value*180/pl.pi for k in fils.orientation_branches])
z=pl.where(pl.isnan(orients)==False)[0]
orients=orients[z]
rmsmom1=pl.array(rmsmom1)[z]
goodrmsmom1=pl.array(goodrmsmom1)
goodorients=pl.array(goodorients)

pl.subplots_adjust(left=0.1,bottom=0.1)
pl.clf()
z=pl.where(rmsmom1>1)[0]
n,bb=pl.histogram(pl.absolute(orients),bins=10*pl.frange(9))
pl.plot(0.5*(bb[:-1]+bb[1:]),n,label="all branches")
#n,bb=pl.histogram(pl.absolute(goodorients),bins=10*pl.frange(9))
#pl.plot(0.5*(bb[:-1]+bb[1:]),n,label="matched to dendros")
n,bb=pl.histogram(pl.absolute(orients[z]),bins=10*pl.frange(9))
pl.plot(0.5*(bb[:-1]+bb[1:]),n,label="vel. rms > 1 km/s")
z=pl.where(goodrmsmom1>1)[0]
#n,bb=pl.histogram(pl.absolute(goodorients[z]),bins=10*pl.frange(9))
#pl.plot(0.5*(bb[:-1]+bb[1:]),n,label="with large vel rms and match dendro")
pl.xlabel(" <---- aligned ---- perpendicular ---->")
pl.legend(loc="best",prop={"size":10})

pl.savefig(label+".plots/"+label+".fils_dendro_angle_mom1_histogram.png")

#associated_dendros - for each filament, 
#"indices":zelong[zorder[zoverlap]],
#"pixoverlap":pixoverlap[zorder[zoverlap]]})        
#"overlapbranches":which branches overlap with each dendro in the indices list





        
#-------------
# now calculate properties of filaments, but highlight the ones associated with dendros,
# since they're more "real"

# for each branch of each fil,
Imean = [] # mean intensity
fwhm = [] # fitted FWHM
dfwhm = [] # fitted FWHM
length = [] # length
assdendro = []
#branchprops

for ifil in range(len(fils.filaments)):
    ff = fils.filaments[ifil]
    for ibr in range(ff.branch_properties['number']):
        length.append(ff.branch_properties['length'][ibr])
        Imean.append(ff.branch_properties['intensity'][ibr])
        ff.width_analysis(fils.image,single_branch=True,branchid=ibr+1,deconvolve_width=False)
        fwhm.append( ff.radprof_fwhm(u.pix)[0].value)
        dfwhm.append(ff.radprof_fwhm(u.pix)[1].value)
#        if ibr+1 in associated_dendros[ifil]['overlapbranches']
            
    
        #-------
#        self.skeleton_longpath = \
#            recombine_skeletons(self.filament_arrays["long path"],
#
#        length_output = main_length(max_path, edge_list, labeled_fil_arrays,
#                                    interpts,
#                                    self.branch_properties["length"],
#                                    self.imgscale,
#                                    verbose=verbose, save_png=save_png,
#                                    save_name=self.save_name,
#                                    vskel=self.vskeleton, array_offsets=self.array_offsets)
#        self.lengths, self.filament_arrays["long path"] = length_output

# pre_graph: calculate graphx
# longest_path - take output of pre_graph
# main_length

# need something that puts the intersection points back into the labeled_fil_array that's been decimated to only contain the overlapping fils
# may want to even consider truncating the overlapping fils so they _only_
# have the partial branches that overlap the dendro area - that would help
# avoid the long branches going off outside of the dendro area

# once the skeleton is reattached, feed it to pix_identify and then
# pre_graph and longest_path

                                
# not sure how to do alignment that way. I guess longest path starting with the
# >> run longest path algorithm with just the subbranches identified!!

# subfilament - go to each hub and choose the longer path from there?
# look at how he does longest_path

    
import pickle
pickle.dump({"associated_dendros":associated_dendros,
             "overlapping_fil_branches":overlapping_fil_branches,
             "overlapping_arrays":overlapping_arrays, # just the part of the skel that overlaps
             "overlapping_xarrays":overlapping_xarrays}, # extended to end of branch
            open(label+".fils_dendro.pkl","wb"))








