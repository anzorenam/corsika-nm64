#!/usr/bin/env python3.7
# -*- coding: utf8 -*-

import matplotlib.pyplot as plt
import corsikaio
import numpy as np
import seaborn as sns
import scipy.stats as stat
import os

sns.set(rc={"figure.figsize":(8,4)})
sns.set_context('paper',font_scale=1.5,rc={'lines.linewidth':1.5})
sns.set_style('ticks')
plt.rc('text',usetex=True)
plt.rc('text.latex',preamble=r'\usepackage[utf8]{inputenc} \usepackage[T1]{fontenc} \usepackage[spanish]{babel} \usepackage{amsmath,amsfonts,amssymb} \usepackage{siunitx}')

dir='/run/media/shirokuma/neutron-data/neutrons-cdmx'
ebins=np.ravel(np.outer(10**np.arange(1,4),np.arange(1,10,0.5)))
ps_ground=np.zeros([400,54])
m=np.ones(400)
m[0:18]=np.array([1.0,0,5.10998e-4,5.10998e-4,1.0,0.105658,0.105658,
                  1.0,1.0,1.0,1.0,1.0,1.0,0.939565,0.938272,1.0,1.0,1.0])
jevent=0
jmp,jmn,jnt=0,0,0
muonp_pos=np.zeros([2,10000])
muonn_pos=np.zeros([2,10000])
neutron_pos=np.zeros([2,600000])

pos=True
spectra=True
for k in range(0,4):
  name='{0}/DAT00000{1}'.format(dir,k+1)
  f=corsikaio.CorsikaParticleFile(name)
  for shower in f:
    jevent+=1
    parts=shower.particles
    if np.size(parts)!=0.0:
      a=np.array(0.001*parts['particle_description'],dtype=np.int)
      p=parts['px']**2.0+parts['py']**2.0+parts['pz']**2.0
      ek=1e3*(np.sqrt(p+m[a])-m[a])
      eindex=np.digitize(ek,bins=ebins)
      ps_ground[a,eindex]+=1.0
      muon_p,muon_n,neut=(a==5.0),(a==6.0),(a==13.0)
      if np.any(muon_p) or np.any(muon_n) or np.any(neut):
        x,y=(1.0/100.0)*parts['x'],(1.0/100.0)*parts['y']
        muonp_pos[:,jmp]=np.sum(x*muon_p),np.sum(y*muon_p)
        muonn_pos[:,jmn]=np.sum(x*muon_n),np.sum(y*muon_n)
        neutron_pos[:,jnt]=np.sum(x*neut),np.sum(y*neut)
        jmp+=np.any(muon_p)
        jmn+=np.any(muon_n)
        jnt+=np.any(neut)
  f.close()

if pos==True:
  np.savetxt('muon_p.dat',muonp_pos,fmt='%1.4f')
  np.savetxt('muon_n.dat',muonn_pos,fmt='%1.4f')
  np.savetxt('neutron_pos.dat',neutron_pos,fmt='%1.4f')
  np.savetxt('spectrum.dat',ps_ground[0:18,:],fmt='%1.4f')
if spectra==True:
  norm=np.amax(ps_ground,axis=1)
  print(np.sum(ps_ground[0:18,:],axis=1))
  norm=np.divide(1.0,norm,where=norm!=0)
  ps_dense=np.transpose(norm[np.newaxis])*ps_ground
  fig,ax=plt.subplots(nrows=1,ncols=1,sharex=False,sharey=True)
  sns.heatmap(ps_dense[0:18,:],vmin=0,vmax=1.0,xticklabels=7,rasterized=True,cmap='viridis')
  ax.set_xticklabels(ebins[::7])
  plt.xlabel(r'Energy $\left(\si{\mega\electronvolt}\right)$',x=0.9,ha='right')
  plt.ylabel(r'Particles')
  plt.tight_layout(pad=1.0)
  plt.savefig('neutron-prod.pdf')
