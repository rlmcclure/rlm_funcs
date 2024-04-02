#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 01 2024

@author: rlmcclure
"""
#%%
import numpy as np
import h5py
import numpy as np
from glob import glob
import multiprocessing as mp

def namestr(obj):
    try:
        return([name for name in globals() if globals()[name] is obj][0])
    except IndexError:
        return([name for name in locals() if locals()[name] is obj][0])

#%%
# first we'll set where we want it to look for the snaps
snapdirpath = '/usr/data/scylla/rmcclure/cca/sims/bulgeless/runs/lvl4-MB000/output/'#'/usr/data/scylla/rmcclure/cca/sims/lvl3/output/'
nested = 0#1#
if nested:
    snapdirpath += 'snapdir_'
else:
    snapdirpath += 'snapshot_'
newdtype = [('t','d'),('id','Q'),('x','d'),('y','d'),('z','d'),('vx','d'),('vy','d'),('vz','d')]
nsnap = len(glob(snapdirpath+'*'))
verbose = 3

# here i'm going to make a dictionary so that I can give a key and from there select the part of the sim to load
# if you don't know these, you can explore the file with an with f= h5py.File(flist[randfilen], 'r'); print(help(f)); print(f.keys()); to learn about your file
aspectdict = {'gas':"/PartType0",'halo':"/PartType1",'disk':"/PartType2",'bulge':"/PartType3",'newstars':"/PartType4"}

#%% helper functions

def strsnapn(snapn):
    if snapn <10:
        snapstr = '00'+str(snapn)
    elif snapn <100:
        snapstr = '0'+str(snapn)
    else:
        snapstr = str(snapn)
    return(snapstr)
if nested:
    # use glob to get the list of all the snap files for that timestep
    getflist = lambda snapn: glob(snapdirpath+strsnapn(snapn)+'/*.hdf5')
else:
    getflist = lambda snapn: [snapdirpath+strsnapn(snapn)+'.hdf5']
    
# a quick funcitons to get the data out of the hdf5 files
getfeature = lambda closedgroup,featstr: np.array(closedgroup.__getitem__(featstr))

#to extract attributes of each 
def char(aspect):
    if verbose>3:
        print('\nFor aspect labeled '+namestr(aspect)+' keys are:')
        for kk in aspect.keys():
            print(kk)
    idds = getfeature(aspect,"ParticleIDs")
    if verbose>3:
        print('N particles:'+str(len(idds)))
    locs = getfeature(aspect,"Coordinates")
    vels = getfeature(aspect,"Velocities")
    return(idds,locs,vels)

def makesubstruct(tind,idds,locs,vels):
    #now we'll make a structured array so we can use the data efficiently 
    snp = np.empty((len(idds)),dtype=newdtype)

    #once we find the time in the header we'll fix this
    snp['t'] = np.ones(len(idds))*int(tind)
    
    #now we can fill in the info 
    snp['id'] = idds
    snp['x'] = locs[:,0]
    snp['y'] = locs[:,1]
    snp['z'] = locs[:,2]
    snp['vx'] = vels[:,0]
    snp['vy'] = vels[:,1]
    snp['vz'] = vels[:,2]
    return(snp)


#redefine for use in mp
def rtstructs(filen,flist):
    '''
    For use on a GADGET/GIZMO sim output for hdf5 files. 
    Creates a stuctured array with the locations and velocities, and particle Ids
    '''
    with h5py.File(flist[filen], 'r') as f:
        tind = f.__str__().split('_')[1].split('.')[0]
        if verbose >3:
            print(f.keys())
        hdr = f.__getitem__("/Header")

        disk = f.__getitem__(aspectdict['disk'])
        #disk particles
        diskidds,disklocs,diskvels = char(disk)
        diskstruct = makesubstruct(tind,diskidds,disklocs,diskvels)
        
        if str(f.keys()).count(aspectdict['gas']):
            gas = f.__getitem__(aspectdict['gas'])
            #gas particles
            gasidds,gaslocs,gasvels = char(gas)
            gasstruct = makesubstruct(tind,gasidds,gaslocs,gasvels)
        else:
            gasstruct = np.zeros_like(diskstruct[0:1])


        if str(f.keys()).count(aspectdict['halo']):
            halo = f.__getitem__(aspectdict['halo'])
            #halo particles
            haloidds,halolocs,halovels = char(halo)
            halostruct = makesubstruct(tind,haloidds,halolocs,halovels)
        else:
            halostruct = np.zeros_like(diskstruct[0:1])


        if str(f.keys()).count(aspectdict['bulge']):
            bulge = f.__getitem__(aspectdict['bulge'])
            #bulge particles
            bulgeidds,bulgelocs,bulgevels = char(bulge)
            bulgestruct = makesubstruct(tind,bulgeidds,bulgelocs,bulgevels)
        else:
            bulgestruct = np.zeros_like(diskstruct[0:1])

        if str(f.keys()).count(aspectdict['newstars']):
            newstars = f.__getitem__(aspectdict['newstars'])
            #newstars particles
            newstaridds,newstarlocs,newstarvels = char(newstars)
            newstarstruct = makesubstruct(tind,newstaridds,newstarlocs,newstarvels)
        else:
            newstarstruct = np.zeros_like(diskstruct[0:1])

    outputdict = {'gas':gasstruct,'halo':halostruct,'disk':diskstruct,'bulge':bulgestruct,'newstars':newstarstruct}
    return(outputdict)

#%% create multiprocessing loader for each snap, initalize things by loading the first and last to get array sizes
#do first snap
snapn = 0
n0typedict = {'gas':0,'halo':0,'disk':0,'bulge':0,'newstars':0}

if verbose>3:
    print('snapn is '+str(snapn))
flist = getflist(snapn)
    
#redef for mp 
def rtsubstruct(filen):
    return(rtstructs(filen,flist))

if __name__ == "__main__":
    with mp.Pool(processes=mp.cpu_count()) as pool:#
        for fln, output in enumerate(pool.imap(rtsubstruct, np.arange(len(flist)))):
            for kk in output.keys():
                n0typedict[kk]+=output[kk].shape[0]
print('snapn:',snapn,n0typedict)

#and final snap
snapn = nsnap-1
nftypedict = {'gas':0,'halo':0,'disk':0,'bulge':0,'newstars':0} #get final n counts

if verbose>3:
    print('snapn is '+str(snapn))
flist = getflist(snapn)

#redef for mp 
def rtsubstruct(filen):
    return(rtstructs(filen,flist))
    
if __name__ == "__main__":
    with mp.Pool(processes=mp.cpu_count()) as pool:#
        for fln, output in enumerate(pool.imap(rtsubstruct, np.arange(len(flist)))):
            for kk in output.keys():
                nftypedict[kk]+=output[kk].shape[0]
print('snapn:',snapn,nftypedict)

#%%
snapn = 3
tempntypedict = {'gas':0,'halo':0,'disk':0,'bulge':0,'newstars':0}
gassnap = np.empty((np.amax([n0typedict['gas'],nftypedict['gas']])),dtype=newdtype)
halosnap = np.empty((np.amax([n0typedict['halo'],nftypedict['halo']])),dtype=newdtype)
disksnap = np.empty((np.amax([n0typedict['disk'],nftypedict['disk']])),dtype=newdtype)
bulgesnap = np.empty((np.amax([n0typedict['bulge'],nftypedict['bulge']])),dtype=newdtype)
newstarssnap = np.empty((np.amax([n0typedict['newstars'],nftypedict['newstars']])),dtype=newdtype)
outputsnapdict = {'gas':gassnap,'halo':halosnap,'disk':disksnap,'bulge':bulgesnap,'newstars':newstarssnap}

flist = getflist(snapn)

#redef for mp 
def rtsubstruct(filen):
    return(rtstructs(filen,flist))
    
if __name__ == "__main__":
    with mp.Pool(processes=mp.cpu_count()) as pool:#
        for fln, output in enumerate(pool.imap(rtsubstruct, np.arange(len(flist)))):
            predict = tempntypedict.copy()
            for kk in output.keys():
                tempntypedict[kk]+=output[kk].shape[0]
                outputsnapdict[kk][predict[kk]:tempntypedict[kk]] = output[kk]
if verbose>2:
    print('snapn:',snapn,tempntypedict)