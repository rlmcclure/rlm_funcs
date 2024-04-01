#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 16:01:10 2020

@author: rlm
"""
#%%
import os
import numpy as np
from astropy.io import fits
import pandas as pd
from rlm_funcs.helpers import *
from astropy.wcs import WCS
from spectral_cube import SpectralCube as sc
import h5py
from glob import glob

#%%
def checkforlineerr(hdul,verbose=None):
    '''
    a function to specifically handle the GILDAS fits header error
        if 'ORIGIN' key is present and contains GILDAS 
            try 'LINE' and except to delete if that fails
        if the 'LINE' keyword is present ...
    '''
    try:
        if hdul[0].header['ORIGIN'].find('GILDAS')>-1:
            if verbose is not None:
                print('File appears to be from GILDAS, trying to remove LINE key error.')
            try:
                hdul[0].header['LINE']
            except:
                del hdul[0].header['LINE']

            # for kw in "CTYPE", "CRVAL", "CRPIX", "CDELT", "CUNIT":
            #     hdul[0].header.remove(f"{kw}{3}")
            #     print('dropped '+kw,flush=True)
            # fits.writeto(outname, hdu.data[0,0], hdu.header, clobber=True)
                return(hdul)
    except:
        return(hdul)

def getcube(cubestr,verbose=None):
    with fits.open(cubestr, ignore_missing_end=True) as hdu:
        hdu=checkforlineerr(hdu,verbose=verbose)
        cube = sc.read(hdu)
    return(cube)

def infoprint(cubestr,verbose=None,varindx=None):
    '''
    a function to print out the header info for the stuff in these fits files
    --------
    params:
        cubestr = filestring for cube
        verbose = bool, defaults to true for printing the info for all headers
        
    returns:
        fileinfo: list of string names of variables in file 
    --------
    '''
    #open file
    with fits.open(cubestr) as hdul:
        #check for GILDAS err
        hdul = checkforlineerr(hdul,verbose=verbose)

        #if there's a specified variable, print only that
        if varindx != None:
            fileinfo = hdul[varindx].header
            
            #if output on then print header info for each variable
            if verbose != 0:
                print(' ')
                print(fileinfo)
                print(' ')
                        
        #otherwise do all
        else:
            #get the string names of variables in HDUList
            fileinfo = [n[1] for n in hdul.info(output=False)]
            
            #step thorugh and print headers
            for ii in np.arange(np.size(hdul)):
                    hdr = hdul[ii].header
                    #if output on then print header info for each variable
                    if verbose != 0:
                        print(' ')
                        print(fileinfo[ii])
                        print(hdr)
                        print(' ')
                    
#    #close file
#    hdul.close()
    return(fileinfo)
    
#%%
def getvar(cubestr,varindx,verbose=None,neststr=None):
    '''
    a function to snag one variable's cube from a fits file
    --------
    params:
        cubestr = filestring for cube
        varindx = variable's number or string index
        verbose = bool, defaults to true for printing the info for all headers
        neststr = nested variable's string index (optional)
        
    returns:
        varstr: list of string names of variables in file 
        vardata: the variable file data, hdul[var_indx].data
    --------
    '''
    #open file
    with fits.open(cubestr) as hdul:
    
        #check for GILDAS err
        hdul = checkforlineerr(hdul,verbose=verbose)
        
        #get the string names of variables in HDUList
        fileinfo = [n[1] for n in hdul.info(output=False)]
        
        #get variable data
        vardata = hdul[varindx].data
        hdr = hdul[varindx].header
        
        #if nested index is present
        if neststr != None:
    #        parent_vardata = np.copy(vardata)
            vardata = hdr[neststr]
            
            
        #if output on then print variable string
        if verbose != None:
            print(' ')
            if type(varindx) != str:
                print(fileinfo[varindx])
            print(hdr)
            print(' ')
                
    return(hdr,vardata)

#%%=============================================================
#=====open different type of files, from M. Soares =============
#===============================================================

def readcolumn(var,col,infile,datformat='float',div=None,fill=False):
	fin=open(infile,mode='r')
	data=fin.readline().split(div)
	if (datformat=='float'):
		while((not len(data)==0) and (not data == [''])):
			if(not (data[0].startswith('#'))):
				try:
					var.append(float(data[col-1]))
				except ValueError:
					if(fill):
						var.append(np.nan)
					else:
						raise
				except IndexError:
					if(fill):
						var.append(np.nan)
					else:
						raise
			data=fin.readline().split(div)

	if (datformat=='str'):
		while((not len(data)==0) and (not data == [''])):
			if(not (data[0].startswith('#'))):
				try:
					var.append(data[col-1])
				except IndexError:
					if(fill):
						var.append(None)
					else:
						raise IndexError

			data=fin.readline().split(div)
	fin.close()
	return

#%%
def fistar2df(filesstr):

    '''
    Takes a string file location and outputs a pandas data frame with X&Y
    Index the data frame to get the X(Y) loc by df.loc[TARGETID]['X(Y)']
    _________________
    Params:(string) file path
    
    Returns:(pandas.DataFrame) with X, Y, & FWHM col with index from star ID in fistar
    _________________
    '''
    with open(filesstr) as f:    
        # first we allocate memory
        siz = 0
        for l in f:
            if l[0] != '#':
                siz+=1
        ID = np.empty(siz)
        X = np.empty(siz)
        Y = np.empty(siz)
        FWHM = np.empty(siz)
        
    with open(filesstr) as f:  
        ii = 0
        for l in f:
            if l[0] != '#':
                ID[ii] = l[:7]
                
                X[ii] = float(l[7:16])
                
                Y[ii] = float(l[16:25])
                
                FWHM[ii] = float(l[25:31])
                
                ii+=1
                    
    df = pd.DataFrame(data={'X':X,'Y':Y,'FWHM':FWHM},index=ID)
        
#            etc cols:  
#            NPix 
#            sigma  
#            delta  
#            Flux      
#            SN    
#            Bg      
#            Amp      
#            S      
#            D      
#            K    
#            CMax   

    return(df)
    
#%%
def lcs2df(filesstr,aps,apN):
    '''
    Takes a string file location and outputs a pandas data frame with X&Y
    Index the data frame to get the X(Y) loc by df.loc[TARGETID]['X(Y)']
    _________________
    Params:(string) file path
    
    Returns:(pandas.DataFrame) with X, Y, & FWHM col with index from star ID in fistar
    _________________
    '''
    with open(filesstr) as f:    
        # first we allocate memory
        naps = len(aps.split(','))
        siz = 0
        for l in f:
            if l[0] != '#':
                siz+=1
        ID = []
        X = np.empty(siz)
        Y = np.empty(siz)
            
        vals = np.empty([siz,naps])
        errs = np.empty([siz,naps])
    with open(filesstr) as f:  
        ii = 0
        for l in f:
            if l[0] != '#':
                ID.append(str(l[:19]))
                X[ii] = float(l[24:31])#safe to assume these may need to be bumped to the right one value
                Y[ii] = float(l[36:43])
                
                st = 64 #start parsing the apertures
                spc = 5
                dig = 8
                for aa in np.arange(naps):
                    if l[st]=='G':
                        st = st+spc
                        en = st+dig
                        vals[ii,aa] = float(l[st:en])
                        errs[ii,aa] = float(l[st+dig+spc:en+dig+spc])
                        st = en+dig+spc
                    else:
                        st = st+spc
                        en = st+dig
                        vals[ii,aa] = 0
                        errs[ii,aa] = 0
                        st = en+dig+spc
                    
                
                ii+=1
                
        ap = aps.split(',')[apN]
    #build this out to get all the apertures
    df = pd.DataFrame(data={'X':X,'Y':Y,ap:vals[:,apN],'errs':errs[:,apN]},index=ID)
    return(df)

#%%
def loadf(path,verbose = 3, commentchar = '#', spacerchar = ' ', defaultdtype='D',typedict = {'region':'i','i_r':'i','i_z':'i'}):
    with open(path) as f:
        siz = 0 
        datastart = 0
        columnsind = 0
        for ii,l in enumerate(f):
            if verbose >2:
                if l[0] == commentchar:
                    print('info line at ', ii,':', l)
            #get columns as headers
            if l.count('columns')>0:
                columnsind = ii+1
            if ii == columnsind:
                # they're not regularly spaced colums 
                colhd = l.split(' ')
                cols = [] #for the acutal headers
                for jj,c in enumerate(colhd):
                    if c!= '':
                        #so this is the complicated part. you need to make a dictionary that sets the dtype for each column. you can have a default type and then a set where you use the dictionary to set it
                        try:
                            cols.append((c,typedict[c]))
                        except KeyError:
                            cols.append((c,defaultdtype))
                            
            #get data
            if l.count('data')>0:
                datastart = ii+1
            if (datastart != 0) and (l[0] != '#'):
                    siz+=1
            

        #now with final line set the data col widths
        colinds = [0] #where to cut the columns
        emptychar = l[0]
        for jj,c in enumerate(l[:-1]):
            if l[jj]!=emptychar and l[jj+1]==emptychar:
                colinds.append(jj+1)
        colinds = np.array(colinds,dtype='int')
        
        #create your structured array now
        arr = np.empty(siz,dtype=cols)

    with open('/Users/kookoo2052/Downloads/out.dat') as f:
        #fill it
        for ii,l in enumerate(f):
            if ii > datastart:
                colind = 0
                for key in cols[:-1]:
                    # try:
                    # print(l[colinds[colind]:colinds[colind+1]])
                    arr[key[0]][ii-datastart] = l[colinds[colind]:colinds[colind+1]]
                    colind += 1
    return(arr)