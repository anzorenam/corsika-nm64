#!/usr/bin/env python3.7
# -*- coding: utf8 -*-

import matplotlib as mat
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.colors as colors
import matplotlib.colorbar as cbar

sns.set(rc={"figure.figsize":(8,4)})
sns.set_context('paper',font_scale=2.0,rc={'lines.linewidth':1.5})
sns.set_style('ticks')
mat.rc('text',usetex=True)
mat.rcParams['text.latex.preamble']=[r'\usepackage[utf8]{inputenc}',r'\usepackage[T1]{fontenc}',r'\usepackage[spanish]{babel}',r'\usepackage[scaled]{helvet}',r'\renewcommand\familydefault{\sfdefault}',r'\usepackage{amsmath,amsfonts,amssymb}',r'\usepackage{siunitx}']

#muon0=np.loadtxt('muon_n.dat')[:,0:5588]
#muon1=np.loadtxt('muon_p.dat')[:,0:2543]
#muons=np.hstack([muon0,muon1])
neutrons=np.loadtxt('neutron_pos.dat')[:,0:468601]
xl,xh=-10000,10000
yl,yh=-10000,10000
fig,ax=plt.subplots(nrows=1,ncols=1,sharex=False,sharey=False)
h=ax.hexbin(neutrons[0,:],neutrons[1,:],gridsize=(100,100),bins='log',
            mincnt=1,cmap='viridis',vmin=1e0,vmax=1e2,extent=(xl,xh,yl,yh))
ax.set(xlim=(xl,xh),ylim=(yl,yh))
ax.set_xlabel(r'Coordenada X $\left[\si{\m}\right]$')
ax.set_ylabel(r'Coordenada Y $\left[\si{\m}\right]$')
fig.subplots_adjust(bottom=0.15,right=0.9,wspace=0.1)
cax=plt.axes([0.93,0.15,0.025,0.75])
fig.colorbar(h,cax=cax)
plt.show()
