#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 15:42:49 2020

@author: rlm
"""

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# try:
#     from scipy import interpolate
# except ImportError:
#     raise('Install scipy or don\'t use interpolation in doublescplot.')
import numpy as np
from rlm_funcs.helpers import *

params = {'axes.labelsize': 20,
          'axes.titlesize': 20,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16}
plt.rcParams.update(params)


import json
from matplotlib.colors import LinearSegmentedColormap
myspringcs = json.loads('["#000040","#000041","#070041","#0e0042","#150043","#1c0044","#230044","#290045","#300146","#370247","#3d0347","#440448","#4a0449","#50054a","#56064a","#5c074b","#62084c","#68094c","#6e0a4d","#740b4e","#790c4f","#7f0d4f","#840e50","#8a0f51","#8f1152","#941252","#991353","#9e1454","#a31655","#a81755","#ad1856","#b11957","#b61b58","#bb1c58","#bf1e59","#c31f5a","#c7215b","#cc225b","#d0245c","#d4255d","#d7275e","#db285e","#df2a5f","#e22c60","#e62d61","#e92f61","#ed3162","#f03363","#f33464","#f63664","#f93865","#fc3a66","#ff3c67","#ff3e67","#ff4068","#ff4269","#ff446a","#ff466a","#ff486b","#ff4a6c","#ff4c6d","#ff4e6d","#ff506e","#ff536f","#ff5570","#ff5770","#ff5971","#ff5c72","#ff5e73","#ff6073","#ff6374","#ff6575","#ff6776","#ff6a76","#ff6c77","#ff6f78","#ff7179","#ff7479","#ff777a","#ff797b","#ff7c7c","#ff7f7c","#ff817d","#ff847e","#ff877f","#ff897f","#ff8c80","#ff8f81","#ff9282","#ff9582","#ff9883","#ff9b84","#ff9e85","#ffa185","#ffa486","#ffa787","#ffaa88","#ffad88","#ffb089","#ffb38a","#ffb68b","#ffb98b","#ffbd8c","#ffc08d","#ffc38e","#ffc68f","#ffca8f","#ffcd90","#ffd091","#ffd492","#ffd792","#ffdb93","#ffde94","#ffe295","#ffe595","#fce996","#f9ec97","#f6f098","#f3f498","#f0f799","#edfb9a","#e9ff9b","#e6ff9b","#e2ff9c","#dfff9d","#dbff9e","#d7ff9f","#d4ff9f"]')
myfunnyvalentinecs = json.loads('["#410400","#450600","#490700","#4d0800","#510900","#550a00","#590b00","#5d0c00","#610d00","#650e00","#680f00","#6c1000","#701200","#741300","#781400","#7c1500","#7f1600","#831700","#871900","#8a1a00","#8e1b00","#921c00","#961d00","#991f00","#9d2000","#a02100","#a42200","#a82300","#ab2500","#af2600","#b22700","#b62800","#b92a00","#bd2b00","#c02c00","#c42e00","#c72f00","#ca3000","#ce3100","#d13300","#d43400","#d83500","#db3700","#de3800","#e23900","#e53b00","#e83c00","#eb3e04","#ef3f09","#f2400d","#f54212","#f84317","#fb451b","#fe4620","#ff4725","#ff4929","#ff4a2e","#ff4c33","#ff4d37","#ff4f3c","#ff5040","#ff5145","#ff5349","#ff544e","#ff5652","#ff5757","#ff595b","#ff5a60","#ff5c64","#ff5d69","#ff5f6d","#ff6071","#ff6276","#ff647a","#ff657f","#ff6783","#ff6887","#ff6a8b","#ff6b90","#ff6d94","#ff6f98","#ff709d","#ff72a1","#ff73a5","#ff75a9","#ff77ad","#ff78b2","#ff7ab6","#ff7cba","#ff7dbe","#ff7fc2","#ff81c6","#ff82ca","#ff84ce","#ff86d2","#ff87d6","#ff89da","#ff8bde","#ff8ce2","#ff8ee6","#ff90ea","#ff92ee","#ff93f2","#ff95f6","#ff97fa","#ff99fe","#ff9aff","#ff9cff","#ff9eff","#ffa0ff","#ffa1ff","#ffa3ff","#ffa5ff","#ffa7ff","#ffa9ff","#ffabff","#ffacff","#ffaeff","#ffb0ff","#ffb2ff","#ffb4ff","#ffb6ff","#ffb8ff","#ffb9ff","#ffbbff","#ffbdff","#ffbfff","#ffc1ff"]')
mymosscs = json.loads('["#001f0f","#00200f","#00210f","#00230f","#002410","#002510","#002710","#002810","#002910","#002b10","#002c10","#002d10","#002f11","#003011","#003211","#003311","#003511","#003611","#003811","#003912","#003b12","#003c12","#003e12","#003f12","#034112","#074212","#0b4412","#0f4612","#134713","#174913","#1b4a13","#1e4c13","#224e13","#264f13","#2a5113","#2d5313","#315513","#345614","#385814","#3b5a14","#3f5c14","#425d14","#455f14","#496114","#4c6314","#4f6514","#526714","#556915","#586a15","#5b6c15","#5e6e15","#617015","#647215","#677415","#6a7615","#6c7815","#6f7a15","#727c15","#747e15","#778016","#798216","#7c8416","#7e8616","#818816","#838b16","#858d16","#888f16","#8a9116","#8c9316","#8e9516","#909716","#929a16","#949c16","#969e17","#98a017","#9aa317","#9ba517","#9da717","#9fa917","#a0ac17","#a2ae17","#a4b017","#a5b317","#a7b517","#a8b817","#a9ba17","#abbc17","#acbf17","#adc117","#aec417","#b0c617","#b1c918","#b2cb18","#b3ce18","#b4d018","#b5d318","#b5d518","#b6d818","#b7da18","#b8dd18","#b9df18","#b9e218","#bae518","#bae718","#bbea18","#bbed18","#bcef18","#bcf218","#bcf518","#bdf718","#bdfa18","#bdfd18","#bdff18","#bdff18","#beff18","#beff18","#beff18","#bdff18","#bdff18","#bdff18","#bdff18","#bdff18","#bcff18","#bcff18","#bcff18","#bbff18","#bbff18","#baff18"]')
myterrecs = json.loads('["#002300","#002400","#002500","#002600","#002700","#002800","#002900","#002a00","#002b00","#002c00","#002d00","#002e00","#002f00","#003100","#003200","#003300","#003400","#003500","#003600","#043700","#0a3800","#0f3a00","#143b00","#193c00","#1f3d00","#243f00","#294000","#2e4100","#324200","#374400","#3c4500","#404600","#454800","#4a4900","#4e4a00","#524c00","#574d00","#5b4e00","#5f5000","#635100","#675300","#6b5400","#6f5500","#725700","#765800","#7a5a02","#7d5b05","#815d07","#845e0a","#88600d","#8b6210","#8e6312","#916515","#946619","#97681c","#9a691f","#9d6b22","#a06d26","#a36e29","#a5702d","#a87230","#aa7334","#ad7538","#af773c","#b17940","#b37a44","#b67c48","#b87e4d","#ba8051","#bc8155","#bd835a","#bf855f","#c18763","#c38968","#c48a6d","#c68c72","#c78e77","#c8907c","#ca9281","#cb9487","#cc968c","#cd9891","#ce9a97","#cf9c9d","#d09ea2","#d1a0a8","#d1a2ae","#d2a4b4","#d3a6ba","#d3a8c0","#d3aac7","#d4accd","#d4aed3","#d4b0da","#d4b2e0","#d5b4e7","#d5b6ee","#d4b8f5","#d4bbfc","#d4bdff","#d4bfff","#d3c1ff","#d3c3ff","#d3c6ff","#d2c8ff","#d1caff","#d1ccff","#d0cfff","#cfd1ff","#ced3ff","#cdd5ff","#ccd8ff","#cbdaff","#cadcff","#c8dfff","#c7e1ff","#c6e3ff","#c4e6ff","#c2e8ff","#c1ebff","#bfedff","#bdf0ff","#bbf2ff","#baf4ff","#b8f7ff","#b6f9ff","#b3fcff","#b1feff"]')
myprettycs = json.loads('["#001b00","#001a00","#001900","#001800","#001700","#00170b","#001633","#001558","#00147a","#001499","#0013b6","#0013d1","#0012e9","#0012ff","#0012ff","#0011ff","#0011ff","#0011ff","#0011ff","#0411ff","#0a11ff","#1012ff","#1712ff","#1d13ff","#2413ff","#2a14ff","#3115ff","#3816ff","#3f17ff","#4518ff","#4c19ff","#531aff","#591cff","#5f1dff","#651fff","#6b21ff","#7123ff","#7725ff","#7c27ff","#8129ff","#862cff","#8b2eff","#8f31ff","#9334fd","#9736f2","#9a39e6","#9d3cdb","#9f40cf","#a143c4","#a346b9","#a54aad","#a64ea2","#a65197","#a7558c","#a75982","#a65d77","#a5616d","#a46563","#a26a5a","#a06e51","#9d7248","#9a7740","#977b38","#938031","#8f852a","#8b8a23","#868e1e","#819318","#7c9813","#769d0f","#70a20b","#6aa708","#63ac06","#5cb104","#55b603","#4ebb02","#46c002","#3fc503","#37ca04","#2fcf05","#27d408","#1fd90b","#16de0f","#0ee313","#06e818","#00ec1d","#00f123","#00f62a","#00fa31","#00fe39","#00ff41","#00ff4a","#00ff53","#00ff5d","#00ff67","#00ff72","#00ff7d","#00ff88","#00ff94","#00ffa0","#00ffad","#00ffb9","#00ffc6","#00ffd4","#00ffe1","#00ffef","#00fffc","#00ffff","#00ffff","#00ffff","#00ffff","#00ffff","#00ffff","#00ffff","#00ffff","#00ffff","#00ffff","#00ffff","#00ffff","#06ffff","#19ffff","#2dffff","#43ffff","#5cffff","#76ffff","#92ffff","#b0fdff","#d0f5ff"]')
mypretty2cs = json.loads('["#001b00","#002100","#002600","#002b00","#003000","#003400","#003800","#003b00","#003e00","#004100","#004300","#004500","#004700","#004800","#004a00","#004b00","#004c00","#004c00","#004d00","#004d00","#004d00","#024d00","#094d00","#114d00","#184c00","#204c00","#274b00","#2f4b00","#374a00","#3f4900","#474900","#4f4800","#574800","#5f4700","#674600","#6e4600","#764500","#7d4500","#844400","#8b4402","#92440e","#984319","#9e4324","#a4432f","#a9433a","#ae4444","#b3444f","#b7445a","#bb4564","#bf466e","#c24778","#c54882","#c7498c","#c94a95","#cb4b9e","#cc4da7","#cd4fb0","#cd50b8","#cd53c0","#cc55c8","#cb57d0","#ca59d7","#c85cde","#c65fe5","#c362eb","#c065f2","#bd68f7","#b96bfd","#b56fff","#b173ff","#ac76ff","#a67aff","#a17eff","#9b83ff","#9587ff","#8f8bff","#8890ff","#8194ff","#7a99ff","#739eff","#6ba3ff","#63a8ff","#5cadff","#54b2ff","#4cb7ff","#44bcff","#3cc1ff","#34c6ff","#2cccff","#25d1ff","#1dd6ff","#16dbff","#0ee0ff","#07e5ff","#01eaff","#00efff","#00f4ff","#00f9ff","#00fdff","#00ffff","#00ffff","#00ffff","#00ffff","#00ffff","#00ffff","#00ffff","#00ffff","#00ffff","#00ffff","#00ffff","#00fffe","#00fffc","#00fffa","#00fff8","#00fff6","#00fff4","#0cfff3","#18fff2","#26fff1","#36fff1","#47fff1","#5afff1","#6ffff2","#86fff3","#9efff4","#b9fff6","#d5fff9","#f4fffc"]')
mypretty3cs = json.loads('["#001b00","#001a00","#091a00","#121900","#1b1900","#231900","#2b1900","#321900","#391900","#401a00","#461a00","#4c1b00","#521b00","#571c00","#5c1d00","#601e00","#641f00","#682000","#6c2100","#6f2200","#722400","#752500","#772700","#792800","#7b2a00","#7c2c00","#7d2d00","#7e2f00","#7f3100","#7f3300","#803500","#803700","#803900","#7f3c00","#7f3e00","#7e4000","#7d4300","#7c4500","#7a4800","#794a00","#774d00","#754f07","#73520e","#715515","#6f581c","#6d5a23","#6b5d2a","#686031","#656338","#63663f","#606946","#5d6c4d","#5a6f54","#57725a","#547561","#517868","#4e7b6f","#4b7e75","#48817c","#458483","#428789","#3e8a90","#3b8d96","#38909c","#3593a2","#3297a9","#2f9aaf","#2c9db5","#29a0bb","#26a3c0","#24a6c6","#21a9cc","#1eacd1","#1cafd7","#1ab2dc","#17b5e1","#15b8e6","#13bbeb","#11bef0","#0fc1f5","#0ec3f9","#0cc6fe","#0bc9ff","#0accff","#09ceff","#08d1ff","#08d4ff","#07d6ff","#07d9ff","#07dbff","#07ddff","#08e0ff","#09e2ff","#0ae4ff","#0be7ff","#0ce9ff","#0eebff","#10edff","#13efff","#15f1ff","#18f2ff","#1bf4ff","#1ff6ff","#23f7ff","#27f9ff","#2cfaff","#30fcff","#36fdff","#3bfeff","#41ffff","#48ffff","#4effff","#55ffff","#5dffff","#65ffff","#6dffff","#76ffff","#7fffff","#89ffff","#93ffff","#9dffff","#a8ffff","#b4ffff","#c0ffff","#ccffff","#d9ffff","#e6ffff","#f4fffc"]')
lavhaze = json.loads('["#37002d","#380031","#3a0035","#3b0039","#3d003d","#3f0041","#400044","#420047","#44004b","#46004d","#470050","#490153","#4b0355","#4d0557","#4f0759","#50095b","#520b5d","#540d5e","#560f60","#581161","#5a1462","#5c1663","#5d1864","#5f1a65","#611c65","#631e66","#652166","#672367","#692567","#6b2767","#6d2967","#6f2c67","#712e67","#733067","#753267","#763567","#783766","#7a3966","#7c3b66","#7e3e65","#804065","#824264","#844564","#864763","#884963","#894c62","#8b4e62","#8d5061","#8f5361","#915560","#925760","#94595f","#965c5f","#985e5e","#99605e","#9b635e","#9d655d","#9f685d","#a06a5d","#a26c5d","#a36f5d","#a5715d","#a7735d","#a8765d","#aa785d","#ab7a5d","#ac7d5e","#ae7f5e","#af815f","#b1845f","#b28660","#b38861","#b58b62","#b68d64","#b78f65","#b89166","#b99468","#ba966a","#bb986c","#bc9b6e","#bd9d70","#be9f72","#bfa275","#c0a478","#c1a67b","#c2a87e","#c2ab81","#c3ad85","#c4af88","#c4b18c","#c5b490","#c5b695","#c6b899","#c6ba9e","#c7bca3","#c7bfa8","#c7c1ae","#c7c3b3","#c8c5b9","#c8c7c0","#c8c9c6","#c8cccd","#c8ced4","#c7d0db","#c7d2e3","#c7d4eb","#c7d6f3","#c6d8fb","#c6daff","#c5dcff","#c5deff","#c4e0ff","#c4e2ff","#c3e4ff","#c2e6ff","#c1e8ff","#c0eaff","#bfecff","#beeeff","#bdf0ff","#bcf2ff","#baf4ff","#b9f5ff","#b8f7ff","#b6f9ff","#b5fbff","#b3fdff","#b1feff"]')
def cmapsegmenter(strorclist,nseg=256):
    try:
        if type(strorclist) == str:
            cmapp = mpl.cm.get_cmap(strorclist,lut=nseg)
        if type(strorclist) == list:
            cmapp = LinearSegmentedColormap.from_list('custom', strorclist,nseg)
    except:
        print('Issue somewhere, trouble shoot for now just giving back input.')
        cmapp = strorclist
    return(cmapp)

def wbkg(fs=None,rtax=0,alph=1,fc='white'):
    '''give a figure a white background'''
    if rtax:
        f,ax = plt.subplots(figsize=fs,tight_layout=True)
    else:

        f = plt.figure(figsize=fs,tight_layout=True)

    f.patch.set_facecolor(fc)
    f.patch.set_alpha(alph)
    if rtax:
        return(f,ax)
    else:
        return f

#%%
def quick_imshow(twod_field,field_title=None,cmapp='viridis',cmin=None,cmax=None,rtnfig=0):
    '''
    plots a two-dimensional field using imshow with option to return the figure object
    ---------
    params:
        twod_field = a two dimensional array to plot
        field_title = plot title
        cmapp = color map, defaults to viridis
        cmin = color range minimum, defaults to the field minimum
        cmax = color range maximum, defaults to the field maximum
        rtnfig = boolean to return the figure as an object, defaults to False
        
    returns:
        fig (optional) = object of the figure, for itteration over
        figure 
    '''
    if cmin==None:
        cmin=np.amin(twod_field)
        cmax=np.amax(twod_field)
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.imshow(twod_field,origin='lower',cmap=cmapp,clim=(cmin,cmax))
    if field_title is None:
        try:
            field_title = namestr(twod_field)
        except:
            field_title = 'Quick Plot'
    fig.suptitle(field_title)
    sm = plt.cm.ScalarMappable(cmap=cmapp,norm=plt.Normalize(cmin,cmax))
    sm._A=[]
    cb = plt.colorbar(sm,fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=8)
    if rtnfig==1:
        return(fig)
    else:
        plt.show()
#%%
def quick_imshow_savepngmv(twod_field,field_title=None,save_str='quickplot',n=0,cmapp='viridis',show=0,save=1,cmin=None,cmax=None,rtnfig=0):
    '''
    plots a two-dimensional field using imshow with option to return the figure object
    ---------
    params:
        twod_field = a two dimensional array to plot
        field_title = plot title
        save_str = string for the file save name
        n = for itteration over many files, defaults to 0
        cmapp = color map, defaults to viridis
        show = boolean to show the figure, defaults to False
        save = boolean to save the figure, defaults to True
        cmin = color range minimum, defaults to the field minimum
        cmax = color range maximum, defaults to the field maximum
        rtnfig = boolean to return the figure as an object, defaults to False
        
    returns:
        fig (optional) = object of the figure, for itteration over
        figure 
    '''
    if cmin==None:
        cmin=np.amin(twod_field)
        cmax=np.amax(twod_field)
    fig = plt.figure()
    ax1 = plt.subplot(111)
    ax1.imshow(twod_field,origin='lower',cmap=cmapp,interpolation='none',clim=(cmin,cmax))
    if field_title is None:
        try:
            field_title = namestr(twod_field)
        except:
            field_title = 'Quick Plot'
    fig.suptitle(field_title)
    sm = plt.cm.ScalarMappable(cmap=cmapp,norm=plt.Normalize(cmin,cmax))
    sm._A=[]
    cb = plt.colorbar(sm,fraction=0.046, pad=0.04)
    cb.ax.tick_params(labelsize=8)
    if save==1:
        plt.savefig('%s%s.png'%(save_str,n),dpi=200,bbox_inches='tight')
    if rtnfig==1:
        return(fig)
    elif show ==1:
        plt.show()
    else:
        plt.close()
        
#%% 
def quick_hist(data, field_title = None,nbins=100,rng=None,nstd=1,rtnfig=0):
    if rng != None:
        rng = (np.amin(data.flatten()),np.amax(data.flatten()))
    plt.hist(data.flatten(),bins=nbins,range=rng)
    med = np.median(data.flatten())
    plt.axvline(x=med,c='k',label='Median = %.2e'%(med))
    mNstd = med-nstd*np.std(data.flatten())
    pNstd = med+nstd*np.std(data.flatten())
    plt.axvline(x=mNstd,c='k',ls='--',label=r'$\pm%i\sigma$ = $\pm$ %.2e'%(nstd,np.std(data.flatten())))
    plt.axvline(x=pNstd,c='k',ls='--')
    plt.legend()
    if field_title is None:
        try:
            field_title = namestr(data)
        except:
            field_title = 'Quick Hist'
    plt.title(field_title)
    if rtnfig==1:
        return(fig)
    else:
        plt.show()

        #%% pixel by pixel plot from eleanor.Visualization tess package
def pixel_by_pixel(box, colrange=None, rowrange=None, cmap='BrBG', mask=None,
                    xlim=None,ylim=None, color_by_pixel=True, 
                    save=1,savestr=None):
    
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    """
    Creates a pixel-by-pixel light curve using the corrected flux.
    Contribution from Oliver Hall.
    Parameters
    ----------
    colrange : np.array, optional
         A list of start column and end column you're interested in
         zooming in on.
    rowrange : np.array, optional
         A list of start row and end row you're interested in zooming
         in on.
    cmap : str, optional
         Name of a matplotlib colormap. Default is 'viridis'.
    mask : np.array, optional
         Specifies the cadences used in the light curve. If not, default
         set to good quality cadences.
    xlim : np.array, optional
         Specifies the xlim on the subplots. If not, default is set to 
         the entire light curve.
    ylim : np.array, optional
         Specifies the ylim on the subplots, If not, default is set to 
         the entire light curve flux range.
    color_by_pixel : bool, optional
         Colors the light curve given the color of the pixel. If not,
         default is set to True.
    """
    if colrange is None:
        colrange = [0, np.shape(box)[-1]]
    
    if rowrange is None:
        rowrange = [0, np.shape(box)[-2]]
        
    nrows = int(np.round(colrange[1]-colrange[0]))
    ncols = int(np.round(rowrange[1]-rowrange[0]))
    
      
    fig = plt.figure(figsize=(20,8))
    outer = gridspec.GridSpec(1,2, width_ratios=[1,4])
    
    inner = gridspec.GridSpecFromSubplotSpec(ncols, nrows, hspace=0.1, wspace=0.1,
                                             subplot_spec=outer[1])
    
    i, j = rowrange[0], colrange[0]
    
    ## PLOTS TARGET PIXEL FILE ##
    ax = plt.subplot(outer[0])
    
    c = ax.imshow(np.nanmedian(box[:,rowrange[0]:rowrange[1],colrange[0]:colrange[1]],0),origin='lower',cmap=cmap)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.15)
    plt.colorbar(c, cax=cax, orientation='vertical')
    
    ## PLOTS PIXEL LIGHT CURVES ##
    for ind in range( int(nrows * ncols) ):
        ax = plt.Subplot(fig, inner[ind])
    
        flux = box[:,i,j]
        y = flux
        x = np.arange(len(y))
        
        if color_by_pixel is False:
            color = 'k'
        else:
            rgb = c.cmap(c.norm(np.nanmedian(box,0)[i,j]))
            color = mpl.colors.rgb2hex(rgb)
    
        ax.plot(x, y, c=color)
    
        j += 1
        if j == colrange[1]:
            i += 1
            j  = colrange[0]
    
        if ylim is None:
            ax.set_ylim(np.percentile(y, 1), np.percentile(y, 99))
        else:
            ax.set_ylim(ylim[0], ylim[1])
    
        if xlim is None:
            ax.set_xlim(np.min(x)-0.1, np.max(x)+0.1)
        else:
            ax.set_xlim(xlim[0], xlim[1])
    
        ax.set_xticks([])
        ax.set_yticks([])
    
        fig.add_subplot(ax)
    
    if save == 1:
        plt.savefig('%s_pixbypix_X%iby%i_Y%iby%i.png'%(savestr,rowrange[0],rowrange[1],colrange[0],colrange[1]),dpi=400,bbox_inches='tight')
    return fig
    

def linplot(fig,ax,m,b,x,ylims,label,ls='--',alpha=.8,c='grey',lw=3,lbl=1):
    '''
    creat poly func (y=mx+b) and then give it two x vals (array) to draw line between on existing plot and axis
    '''
    if m == 'inf': 
        #set b as the y limits for the vertical line
        ax.axvline(x,b[0],b[1],ls=ls,alpha=alpha,c=c,lw=lw) 
        if lbl:
            ax.text(x[-1],y[-1],label,alpha=alpha,c=c)
    else:
        p = np.poly1d([m,b])
        if len(x)>0:
            y = p(x)
            ax.plot(x,y,ls=ls,alpha=alpha,c=c,lw=lw)
            if lbl:
                if y[-1]>ylims[-1]:
                    ax.text(x[-1],ylims[-1],'('+label+')',alpha=alpha,c=c)
                
                elif y[-1]<ylims[0]:
                    ax.text(x[-1],ylims[0],'('+label+')',alpha=alpha,c=c)
                else:
                    ax.text(x[-1],y[-1],label,alpha=alpha,c=c)
        else:
            print('x has no length.')
    return(fig,ax)



def simovert(frame,simstr='',nestdir='plots/',xylim=20,xycen=0,zlim=5,zcen=0,sv=1,alph=.7,lw=.5,c='k',htch=None,ls='--',cmapp='magma_r',physmult=1,pustr='',binsn=2000,tind=None):
    '''
    frame has assumption of a structured array with 't' 'x' 'y' and 'z' values.
    '''
    # frame = {}
    try:
        t = '%.3f'%np.unique(frame['t'])
        if tind is None:
            tind = np.unique(frame['tind'])
    except ValueError:
        try:
            t = '%i'%np.unique(frame['t'])
            if tind is None:
                tind = np.unique(frame['t'])
        except KeyError:
            t = ''
            if tind is None:
                tind = ''
    f = wbkg((6,8))
    fig_mos = """
    AA
    BB
    BB
    """
    ax = f.subplot_mosaic(fig_mos)#,subplot_kw={'sharex':1,'sharey':1})

    x='x';y='z'
    _,bbx,bby,_=ax['A'].hist2d(frame[x],frame[y],bins=[np.linspace((-1*xylim-xycen)*physmult,(xylim-xycen)*physmult,binsn),np.linspace((-1*zlim-zcen)*physmult,(zlim-zcen)*physmult,binsn)],cmin=1,alpha=1,norm=mpl.colors.LogNorm(),cmap=cmapp)
    
    ax['A'].set_xlim((-1*xylim-xycen)*physmult,(xylim-xycen)*physmult)
    ax['A'].set_ylim((-1*zlim-zcen)*physmult,(zlim-zcen)*physmult)
    ax['A'].set_ylabel(y+pustr)
    # ax['A'].set_xlabel(x+' ['+SIM.spatialunitlabel+']')
    ax['A'].grid(lw=.5,c='grey',alpha=.5)
    # ax['A'].legend(prop={'size':10},loc=3)


    x='x';y='y'
    _,bbx,bby,_=ax['B'].hist2d(frame[x],frame[y],bins=np.linspace((-1*xylim-xycen)*physmult,(xylim-xycen)*physmult,binsn),cmin=1,alpha=1,norm=mpl.colors.LogNorm(),cmap=cmapp)
    # px = bbx[(bbx>SIM.xinnerlim*physmult) & (bbx<barlims[tind]['barmax']*physmult)]
    pl = 0 
    # pu = SIM.youtlim*physmult
    # ax['B'].fill_between(px,-1.0*pu,pu,fc='None',hatch=htch,alpha=alph,edgecolor=c,lw=lw,ls=ls)#,label='B/PX Selection Region')
    # ax['B'].fill_between(-1.0*px,-1.0*pu,pu,fc='None',hatch=htch,alpha=alph,edgecolor=c,lw=lw,ls=ls)
    # ax['B'].fill_between([xylim*physmult+1,xylim*physmult+2],[xylim*physmult+1,xylim*physmult+2],[xylim*physmult+3,xylim*physmult+4],fc='w',hatch='//',alpha=.5,edgecolor=c,lw=lw,ls=ls,label='B/PX Selection Regions')

    # if barlims[tind]['barmax']>barlims[tind]['barmid']:
    #     circle1 = Circle((0,0),barlims[tind]['barmid']*physmult,facecolor='None',edgecolor='cyan',ls='-',label='Bar Peak')
    #     ax['B'].add_patch(circle1)
    #     circle2 = Circle((0,0),barlims[tind]['barmax']*physmult,facecolor='None',edgecolor='limegreen',ls='-.',label='Bar Outer Limit')
    #     ax['B'].add_patch(circle2)
    #     ax['B'].legend(prop={'size':10},loc=2)
    # # cts,xb,yb,img=ax['B'].hist2d(xaposatframe[x],xaposatframe[y],bins=(bbx,bby),alpha=0.1,norm=mpl.colors.LogNorm(),cmap='binary')

    ax['B'].set_xlim((-1*xylim-xycen)*physmult,(xylim-xycen)*physmult)
    ax['B'].set_ylim((-1*xylim-xycen)*physmult,(xylim-xycen)*physmult)
    ax['B'].set_ylabel(y+pustr)
    ax['B'].set_xlabel(x+pustr)
    ax['B'].grid(lw=.5,c='grey',alpha=.5)
    ax['B'].set_title('%s, t %s'%(simstr,t),size=24)
    if sv:
        plt.savefig(nestdir+simstr+'/%i_simovertime.png'%tind,dpi=300,bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def doubleplot(x,y1,y2,structarr=None,c1='darkred',s2=':',c2='darkgreen',title='quick plot',fs=(9,6)):
    f = wbkg(fs)
    ax = f.add_subplot()
    if structarr is not None:
        ax.plot(structarr[x],structarr[y1],c=c1)
        ax.set_ylabel(y1)
        ax.set_xlabel(x)
        ax2 = ax.twinx()
        ax2.set_ylabel(y2)
        ax2.plot(structarr[x],structarr[y2],c=c2,ls=s2)
        if title == 'quick plot':
            try:
                plt.suptitle(namestr(structarr),fontsize=14,y=.85)
            except IndexError:
                plt.suptitle(title,fontsize=14,y=.85)
    else:
        ax.plot(x,y1,c=c1)
        try:
            ax.set_ylabel(namestr(y1))
            ax.set_xlabel(namestr(x))
        except IndexError:
            ax.set_ylabel('')
            ax.set_xlabel('')
        ax2 = ax.twinx()
        try:
            ax2.set_ylabel(namestr(y2))
        except IndexError:
            ax2.set_ylabel('')
        ax2.plot(x,y2,c=c2,ls=s2)
        plt.suptitle(title,fontsize=14,y=.85)
    ax.yaxis.label.set_color(c1)
    ax.tick_params(axis='y', colors=c1)
    ax2.yaxis.label.set_color(c2)
    ax2.tick_params(axis='y', colors=c2)
    ax.grid(c=c1,alpha=.4)
    ax2.grid(ls=':',c=c2,alpha=.3)
    plt.show()

def doublescplot(x,y1,c1,y2,c2,x2=None,title='quick plot',k1=None,k2=None,kx=None,kc1=None,kc2=None,cm1='inferno',ct1='darkred',cm2='cividis',ct2='darkblue',fs=(11,5),vlines=None,tratio=1):
    f = wbkg(fs)
    ax = f.add_subplot()
    if x2 is None:
        x2 = np.copy(x)
        
    #set title label strings
    try:
        if kx is None:
            kx = namestr(x)
        if k1 is None:
            k1 = namestr(y1)
        if kc1 is None:
            k1 = namestr(c1)
        if k2 is None:
            k2 = namestr(y2)
        if kc2 is None:
            kc2 = namestr(c2)
    except IndexError:
        kx = '';k1 = '';k2 = ''

    #plots
    o1 = ax.scatter(x,y1,c=c1,cmap=cm1,marker='o',alpha=.7)
    if vlines is not None:
        ax.axvline(x=vlines[0],ls='--',c='k',alpha=.5)
        ax.axvline(x=vlines[1],ls='--',c='k',alpha=.5)
    ax.set_ylabel(k1,fontsize=12)
    ax.set_xlabel(kx,fontsize=12)
    ax2 = ax.twinx()
    ax2.set_ylabel(k2,fontsize=12)
    o2 = ax2.scatter(x2,y2,c=c2,cmap=cm2,marker='D')

    # # add ratio through interpolating, too specific 
    # if tratio == True:

    #     #colorbars, has to padding adjusted for extra colorbar ofc
    #     cb = f.colorbar(o1,pad=-.05)#%total')#.1)
    #     cb.set_label(label=kc1,size=12,c=ct1)
    #     cb.ax.tick_params(labelsize=10,color=ct1)
    #     cb = f.colorbar(o2,pad=-.05)#%total')#.1)
    #     cb.set_label(label=kc2,size=12,c=ct2)
    #     cb.ax.tick_params(labelsize=10,color=ct2)

    #     tvals = np.linspace(0,SIM.tlensim*2,len(SIM.times)) #hard coded cause i'm tired    
    #     fn1 = interpolate.interp1d(x,c1,kind='nearest',fill_value='extrapolate')
    #     fn2 = interpolate.interp1d(x2,c2,kind='nearest',fill_value='extrapolate')
    #     o3 = ax2.scatter(tvals,np.zeros(len(tvals)),c=fn2(tvals)/fn1(tvals),s=30,marker='|',cmap='Paired_r')#,vmin=0.25,vmax=4.25)
    #     cb = f.colorbar(o3,pad=.12)
    #     cb.set_label(label='ratio of times, %s/%s'%(k2,k1),size=12,c=ct2)
    #     cb.ax.tick_params(labelsize=10,color=ct2)

    #colorbars
    cb = f.colorbar(o1,pad=-.02)#%total')#.1)
    cb.set_label(label=kc1,size=12,c=ct1)
    cb.ax.tick_params(labelsize=10,color=ct1)
    cb = f.colorbar(o2,pad=.12)#%total')#.1)
    cb.set_label(label=kc2,size=12,c=ct2)
    cb.ax.tick_params(labelsize=10,color=ct2)

    #grids
    ax.tick_params(axis='x',labelsize=10)
    ax.yaxis.label.set_color(ct1)
    ax.tick_params(axis='y', colors=ct1,labelsize=10)
    ax2.yaxis.label.set_color(ct2)
    ax2.tick_params(axis='y', colors=ct2,labelsize=10)
    ax.grid(c=ct1,alpha=.4)
    ax2.grid(ls=':',c=ct2)
        
    plt.suptitle(title,fontsize=14,y=.9)
    plt.show()


def histhexplt(plttype,xvals,yvals,cvals=None,hexfn=np.sum,f=None,ax=None,fs=None,ptitle='Quick Plot',xlim=(None,None),ylim=(None,None),vlim=(None,None),xbins=None,ybins=None,gsxdiv=100,gsydiv=100,logtog=0,densitytog=0,cmapp='viridis',xlabel=None,ylabel=None,rtax=0):
    if f is None:
        f,ax = wbkg(fs,rtax=1)
    elif ax is None:
        ax = f.add_subplot()

    #set gridsize from the values being plotted, use gridsize dividers (gsydiv, gsxdiv) to set size dynamcially
    if ybins is None:
        if gsydiv is not None:
            ybins = int(np.floor(len(np.unique(yvals))/gsydiv))
        else: 
            ybins = int(len(np.unique(yvals)))
    if xbins is None:
        if gsxdiv is not None:
            xbins = int(np.floor(len(np.unique(xvals))/gsxdiv))
        else: 
            xbins = int(len(np.unique(xvals)))
    gs_bs = (xbins,ybins) #number of bins in x and y

    #set extent and range by xlim and ylim if given, otherwise it is set by inar params later
    if xlim[0] is None:
        xv = xvals[np.isfinite(xvals)]
        xmin = np.amin(np.unique(xv))
    else:
        xmin = xlim[0]
    if xlim[1] is None:
        xv = xvals[np.isfinite(xvals)]
        xmax = np.amax(np.unique(xv))
    else:
        xmax = xlim[1]
    
    if ylim[0] is None:
        yv = yvals[np.isfinite(yvals)]
        ymin = np.amin(np.unique(yv))
    else:
        ymin = ylim[0]
    if ylim[1] is None:
        yv = yvals[np.isfinite(yvals)]
        ymax = np.amax(np.unique(yv))
    else:
        ymax = ylim[1]
    ex_rg = ([xmin,xmax],[ymin,ymax])

    if plttype == 'hist':
        if logtog == 1:
            normtog = mpl.colors.LogNorm()
        else:
            normtog=None
        if densitytog == 1:## are these indexings are wrong
            todump=ax.hist2d(xvals.flatten(),yvals.flatten(),bins=gs_bs,range=ex_rg,norm=normtog,vmin=vlim[0],vmax=vlim[1],cmap=cmapp)
            
            ax.pcolormesh(todump[1], todump[2], (todump[0]/np.sum(todump[0])).T,vmin=vlim[0],vmax=vlim[1],cmap=cmapp)        
        else:
            histout=ax.hist2d(xvals.flatten(),yvals.flatten(),bins=gs_bs,range=ex_rg,norm=normtog,vmin=vlim[0],vmax=vlim[1],cmap=cmapp)

        if densitytog == 1:
            f.colorbar(histout,ax=ax)#,label='relative density')
        else:
            f.colorbar(histout[3],ax=ax)#,label='N')

    elif plttype == 'hex':
        if logtog == 1:
            binstog ='log'
        else:
            binstog=None
        hexout=ax.hexbin(xvals,yvals,C=cvals,reduce_C_function=hexfn,gridsize=gs_bs,bins=binstog,extent=list(np.array(ex_rg).flat),vmin=vlim[0],vmax=vlim[1],cmap=cmapp)
        if densitytog == 1:
            outar = hexout.get_array()
            hexout.set_array(outar/np.sum(outar))
            if vlim[0] is None:
                hexout.set_clim(np.amin(hexout.get_array()),np.amax(hexout.get_array()))
            else:
                hexout.set_clim(vlim[0],vlim[1])
        if densitytog == 1:
            f.colorbar(hexout,ax=ax)#,label='relative density')
        else:
            f.colorbar(hexout,ax=ax)#,label='N')

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    # ax.set_xticklabels([])
    # ax.set_yticklabels([])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.grid(linestyle=':',linewidth=1)
    ax.set_title(ptitle)

    if rtax:
        return(f,ax)
    else:
        plt.show()
