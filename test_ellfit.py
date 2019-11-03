import sys,os
gitpaths=['/Users/remy/lustre/naasc/users/rindebet/github/lmc-alma-analysis/']
for gitpath in gitpaths:
    if not gitpath in sys.path:
        sys.path.insert(0,gitpath)

from ellfit import ellfit
import pylab as pl
pl.ion()

rmaj=4.
rmin=3.

x=[]
y=[]

d=rmaj/30 # delta
ny=int(rmin*2./d)
yvals=(pl.frange(ny)-0.5*ny)*2/ny*rmin
nx=int(rmaj*2./d)
xx=(pl.frange(nx)-0.5*nx)*2/nx*rmaj

for yval in yvals:
    z=pl.where((xx/rmaj)**2 <= 1.-(yval/rmin)**2)[0]
    if len(x)>0:
        x=pl.concatenate([x,xx[z]])
    else:
        x=xx[z]
    if len(y)>0:
        y=pl.concatenate([y,yval+pl.zeros(len(z))])
    else:
        y=yval+pl.zeros(len(z))


pl.clf()
pl.plot(x,y,'.')

efit=ellfit(x,y)
print "major:",efit[0]/2,rmaj
print "minor:",efit[1]/2,rmin
