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
import numpy as np
from rlm_funcs.helpers import *

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
    creat poly func (y=mx+b) and then give it two x vals (array) to draw line between on existing plot adn axis
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


def simovert(frame,simstr='',xylim = 20,xycen=0,zlim=3,zcen=0,sv=1,alph=.7,lw=.5,c='k',htch=None,ls='--',cmapp='magma_r',physmult=1,pustr='',binsn=2000):
    # frame = {}
    try:
        tind = np.unique(frame['t'])
    except KeyError:
        tind = 0
    f = wbkg((6,8))
    fig_mos = """
    AA
    BB
    BB
    """
    ax = f.subplot_mosaic(fig_mos)#,subplot_kw={'sharex':1,'sharey':1})

    x='x';y='z'
    _,bbx,bby,_=ax['A'].hist2d(frame[x],frame[y],bins=np.linspace((-1*zlim-zcen)*physmult,(zlim-zcen)*physmult,binsn),cmin=1,alpha=1,norm=mpl.colors.LogNorm(),cmap=cmapp)
    
    ax['A'].set_xlim((-1*xylim-xycen)*physmult,(xylim-xycen)*physmult)
    ax['A'].set_ylim((-1*zlim-zcen)*physmult,(zlim-zcen)*physmult)
    ax['A'].set_ylabel(y+pustr)
    # ax['A'].set_xlabel(x+' ['+SIM.spatialunitlabel+']')
    ax['A'].grid(lw=.5,c='grey',alpha=.5)
    # ax['A'].legend(prop={'size':10},loc=3)


    x='x';y='y'
    _,bbx,bby,_=ax['B'].hist2d(frame[x],frame[y],bins=np.linspace((-1*zlim-zcen)*physmult,(zlim-zcen)*physmult,binsn),cmin=1,alpha=1,norm=mpl.colors.LogNorm(),cmap=cmapp)
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
    ax['B'].set_title('%s, tind %i'%(simstr,tind),size=24)
    if sv:
        plt.savefig('plots/'+simstr+'/%i_simovertime.png'%tind,dpi=300,bbox_inches='tight')
        plt.close()
    else:
        plt.show()