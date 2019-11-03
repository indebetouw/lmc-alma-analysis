# run after fils_dendro

pl.clf()
pl.imshow((pl.nanmax(mom0)-mom0)**2,origin="bottom",cmap="gray")
pl.subplots_adjust(left=0.05,right=0.98,bottom=0.05,top=0.98)
pl.xlim(xra)
pl.ylim(yra)

plotted=pl.zeros(len(pcat),dtype=bool)
for i in range(n_elong_dendro):
    overlapping_fil_branches[i]["dendro_id"]=zelong[i]
    if not plotted[zelong[i]]:
        t=d[zelong[i]]
        mask2d=pl.byte(t.get_mask().max(axis=0))
        #mycont=pl.contour(mask2d,1,levels=[0.1+0.2*t.level],vmin=0,vmax=1,cmap="plasma")
        #mycont=pl.contour(mask2d,1,levels=[0.1],vmin=0,vmax=1,cmap="plasma")
        #mycont=pl.contour(mask2d,1,levels=[0.1],vmin=0,vmax=1,colors='g')
        mycont=pl.contourf(mask2d,levels=[0.1,10],colors='yellow',alpha=0.7)
        plotted[zelong[i]]=True

pl.contour(fils.skeleton, colors='g',linewidths=1,vmin=0,vmax=1,levels=[0.5],linestyles=':')
pl.xticks([])
pl.yticks([])

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

        # go through each overlapping dendro:
        overlapthreshold=20 # >1 to only count the significantly overlapping branches
        for izo in zorder[zoverlap]:

            t=d[zelong[izo]]
            mask2d=pl.byte(t.get_mask().max(axis=0))
            #mycont=pl.contour(mask2d,1,levels=[0.1],vmin=0,vmax=1,cmap="plasma")
            # 0.6: salmon-colored;  0.1=purple
            mycont=pl.contourf(mask2d,levels=[0.1,10],colors='red')

            # again the dendro indices relative to the fil subarray:
            drelx=d[zelong[izo]].indices()[2]-off[0][1]
            drely=d[zelong[izo]].indices()[1]-off[0][0]
            zz=pl.where( (drelx>=0)*(drely>=0)*(drelx<sfarr[1])*(drely<sfarr[0]) )[0]
            for ibr in pl.arange(nbranches)+1:
                mask = lflarr==ibr
                pixoverlapbranches[ibr-1,izo] = mask[drely[zz],drelx[zz]].sum()
                if pixoverlapbranches[ibr-1,izo]>overlapthreshold:  
                    thisbranchlist.append(ibr)
                    xfarr_overlap[pl.where(mask)]=1  # could keep labels by setting these pi

                    
    pl.contour(XX[off[0][1]:off[1][1]+1],YY[off[0][0]:off[1][0]+1],xfarr_overlap,1,levels=[0.5],colors='blue',vmin=0,vmax=1,linewidths=2)


pl.savefig(label+".plots/"+label+".fils_dendro_2.png")


