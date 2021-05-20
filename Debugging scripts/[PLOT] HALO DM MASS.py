# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 10:04:19 2021

@author: Andrew
"""

import numpy as np
import matplotlib.pyplot as plt

import h5py as h5


plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-large')
plt.rc('ytick', labelsize='x-large')
plt.style.use('classic')
# plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True) 
plt.rc('axes', labelsize=18)  
plt.rc('legend', fontsize=13)
plt.rcParams['legend.title_fontsize'] = 'x-large'
hdf_filepath = '../../Input Data/mega_graph_mini_09-11.hdf5'

the_file = h5.File(hdf_filepath,'r')


header_data = the_file['Header']
part_mass = header_data.attrs['part_mass']
list_of_nparts = []

haloes = False

subhaloes = False

both = True
if haloes:

    nprogs = np.empty(0)
    redshifts = np.empty(0)
    for i in range(0,976):
        
        
        graph = the_file[str(i)]
        
    
        nprogs = np.concatenate((nprogs, graph['nparts'][:]*part_mass))
        redshifts = np.concatenate((redshifts, graph['redshifts'][:]))
        
        
    
    fig, ax1 = plt.subplots(figsize=(7,4))
    ax1.hist(nprogs[redshifts <= 0.001],label='z = 0')
    ax1.hist(nprogs[((redshifts > 0.001) & (redshifts <= 1))],zorder=-1,label='$0 < z \leq 1$')
    ax1.hist(nprogs[((redshifts > 1) & (redshifts <= 2))],label='$1 < z \leq 2$')  
    ax1.hist(nprogs[(redshifts > 2)],label='$ z > 2$')  
    ax1.set_yscale('log')
    ax1.set_ylabel('N')
    ax1.set_xlabel('Mass/$M_\odot$')
    ax1.legend()
    plt.tight_layout()
    plt.savefig('output graphs/Halo Mass Distribution.png',dpi=400)
    plt.savefig('output graphs/Halo_dist.eps')
    plt.show()


if subhaloes:
    
    nprogs = np.empty(0)
    redshifts = np.empty(0)
    for i in range(0,976):
        
        
        graph = the_file[str(i)]
        
        try:
            nprogs = np.concatenate((nprogs, graph['sub_nparts'][:]*part_mass))
            redshifts = np.concatenate((redshifts, graph['sub_redshifts'][:]))
        except KeyError:
            print('No Suhaloes')
        
    
    fig, ax1 = plt.subplots(figsize=(7,4))
    ax1.hist(nprogs[redshifts <= 0.001],label='z = 0')
    ax1.hist(nprogs[((redshifts > 0.001) & (redshifts <= 1))],zorder=-1,label='$0 < z \leq 1$')
    ax1.hist(nprogs[((redshifts > 1) & (redshifts <= 2))],label='$1 < z \leq 2$')  
    ax1.hist(nprogs[(redshifts > 2)],label='$ z > 2$')  
    ax1.set_yscale('log')
    ax1.set_ylabel('N')
    ax1.set_xlabel('Mass/$M_\odot$')
    ax1.legend()
    plt.tight_layout()
    plt.savefig('output graphs/Subhalo Mass Distribution.png',dpi=400)
    plt.savefig('output graphs/Subhalo_dist.eps')
    plt.show()


if both:
    
    nprogs = np.empty(0)
    redshifts = np.empty(0)
    for i in range(0,976):
        
        
        graph = the_file[str(i)]
        
    
        nprogs = np.concatenate((nprogs, graph['nparts'][:]*part_mass))
        redshifts = np.concatenate((redshifts, graph['redshifts'][:]))
        
        
    
    sub_nprogs = np.empty(0)
    sub_redshifts = np.empty(0)
    for i in range(0,976):
        
        
        graph = the_file[str(i)]
        
        try:
            sub_nprogs = np.concatenate((sub_nprogs, graph['sub_nparts'][:]*part_mass))
            sub_redshifts = np.concatenate((sub_redshifts, graph['sub_redshifts'][:]))
        except KeyError:
            print('No Suhaloes')
            
            
    print(nprogs.max()/part_mass,nprogs.min()/part_mass)
            



    fig, (ax2,ax1) = plt.subplots(2,1,figsize=(7,6),sharex=True,gridspec_kw={'hspace':0.05})
    
    
    sub_unique_redshifts = np.unique(sub_redshifts)
    
    
    
    sub_redshift_one = sub_unique_redshifts[0]
    sub_redshift_two = sub_unique_redshifts[36]
    sub_redshift_three = sub_unique_redshifts[45]
    
    # bins=np.histogram(np.hstack((np.log18sub_nprogs[sub_redshifts == sub_redshift_one],
    #                              sub_nprogs[sub_redshifts == sub_redshift_two],
    #                              sub_nprogs[sub_redshifts == sub_redshift_three])), bins=20)[1]
    
    
    bins = np.logspace(10,14,25)
    
    ax1.hist(sub_nprogs[sub_redshifts == sub_redshift_one],zorder=1,label='Snapshot ${}; z = 0.0$'.format(61),bins=bins)
    ax1.hist(sub_nprogs[sub_redshifts == sub_redshift_two],zorder=2,
             label='Snapshot ${}; z = {}$'.format(25,np.round(sub_redshift_two,1)),bins=bins)
    
    
    
    ax1.hist(sub_nprogs[sub_redshifts == sub_redshift_three],zorder=3,
             label='Snapshot ${}; z = {}$'.format(16,np.round(sub_redshift_three,1)),bins=bins)  
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_ylabel('N')
    ax1.set_xlabel('M$_\mathrm{DM}$/M$_\odot$')
    ax1.grid(True,linestyle='--')
    ax1.set_axisbelow(True)
    ax1.legend(title='Subhaloes')
    ax1.set_xlim(7*10**9,2*10**14)
    ax2.set_xlim(7*10**9,2*10**14)
    ax1.set_ylim(1,4*10**3)
    ax2.set_ylim(1,4*10**3)
    unique_redshifts = np.unique(redshifts)
    
    redshift_one = unique_redshifts[0]
    redshift_two = unique_redshifts[36]
    redshift_three = unique_redshifts[45]

    
    ax2.hist(nprogs[redshifts == redshift_one],zorder=1,label='Snapshot ${};z = 0.0$'.format(61),bins=bins)
    ax2.hist(nprogs[redshifts == redshift_two],zorder=2,label='Snapshot ${};z = {}$'.format(25,np.round(redshift_two,1)),bins=bins)
    ax2.hist(nprogs[redshifts == redshift_three],zorder=3,label='Snapshot ${};z = {}$'.format(16,np.round(redshift_three,1)),bins=bins)  
    ax2.set_yscale('log')
    # ax2.set_xscale('log')
    ax2.set_ylabel('N')
    # ax2.set_ylim(10**0, 2*10**5)
    # ax1.set_ylim(10**0, 2*10**5)
    ax2.grid(True,linestyle='--')
    ax2.set_axisbelow(True)
    # ax2.set_xlabel('Mass/$M_\odot$')
    ax2.legend(title='Haloes')
    plt.subplots_adjust(hspace=.0)
    plt.tight_layout()
    plt.savefig('output graphs/Halo_subhalo_dist.png',dpi=400,bbox_inches='tight')
    plt.savefig('output graphs/Halo_subhalo_dist.eps')#,bbox_inches='tight')
    plt.show()
    