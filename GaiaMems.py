#python 3.8
import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astroquery.gaia import Gaia
import hdbscan
import numpy as np
from functools import reduce
from scipy import stats
from scipy import optimize
from astropy.modeling import models, fitting
from astropy.table import Table
from astropy.table import join as Tjoin
from astropy.table import vstack as Tvstack
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.cm as cm
# Suppress warnings. Comment this out if you wish to see the warning messages
import warnings
warnings.filterwarnings('ignore')

class GaiaClusterMembers(object):
    '''
    This was dev by Aaron Geller, built on by R L McClure with a goal of fitting further clusters.
    
    This Class will grab data from the Gaia archive, and attempt to determine members using the
    proper motions, radial velocities and parallaxes.
    
    The user must provide the RA and Dec values, and the Class will return the full catalog and
    the indices of the members.
    
    Note: in this code, each membership check builds on the previous one, using only the stars that
    have passed the previous membership selection.
    
    A note about Binaries: The GeDR3 parallaxes are inferred assuming the source is single, and our
    distance estimates assume this too. If the parallax is corrupted by any binarity, our distance
    estimates will be too. Whether this is the case depends on whether we talking about an unresolved
    or a resolved binary. If resolved, and the period is much longer than the GeDR3 baseline (3 years),
    then the parallaxes might be okay. They might also be okay if the sources are unresolved and the
    period is much shorter than the GeDR3 baseline.
    '''
    
    def __init__(self, *args,**kwargs):
        
        self.titlestr = None
        
        #required inputs
        self.RA = None
        self.Dec = None
        
        #must include in either plx or dist, built for dist rn
        self.knowndist = None #parsec
        self.knowndisterr = None #parsec
        self.knownplx = None #mas
        self.knownplxerr = None #mas

        #parameters that the user could change
        self.verbosity = 1
        self.showPlots = False
        self.dmin = 0 #parsec
        self.dmax = 10000 #parsec
        self.dbins = 100
        
        #gaia pull
        self.radius = 1 #in degrees
        self.minPMmembership = 0.5
        self.minPMmembershipHDB = 0.2
        
        #general bounds within gaia pull
        self.minMag = None #set this to do magnitude cuts
        self.maxMag = None #set this to do magnitude cuts
        self.distsig = 3.0
        self.distcutP = 0.00 #cause 3 sigma gives 99.7
        self.cdistCut = None #[inner,outer] as a portion of 1 out to radius
        
        #proper motion fit params
        self.hdbPMM = False
        self.bothPMM = False
        self.PMxmin = -20 #mas/yr
        self.PMxmax = 20 #mas/yr
        self.PMxbins = 200
        self.PMymin = -20 #mas/yr
        self.PMymax = 20 #mas/yr
        self.PMybins = 200
        self.PMxwdth = 1
        self.PMywdth = 1
        self.PMfieldx = 0
        self.PMfieldy = 0
        self.PMfieldxwdth = 5
        self.PMfieldywdth = 5
        self.hdbminN = 30
        self.GaussParams = None
        self.histlogscale = False
        self.fitmembonly = False
        self.PMmean = [None, None]
        
        self.limitedC = False
        
        #options to assign HDBSCAN IDs
        self.hdbscanCID = None
        self.hdbscanFID = None
        self.hdbscanUID = None
        self.hdbminN = 30
        
        #cmd plot params
        self.CMDxmin = 0.3
        self.CMDxmax = 2.5
        self.CMDymin = 18
        self.CMDymax = 8
        
        #easy axcess params
        self.geodist = None
        self.geodlo = None
        self.geodhi = None
        self.plx = None
        self.plxerr = None
        self.Qgeodist = False #toggle for using the photogeometric Bailer-Jones dist
        '''QG prior may be inappropriate; (2) The source may be in a crowded region,\
        in which case the BP/RP spectrum may be blended and the colours incorrect.
        '''
        self.plxdist = False #toggle to use only the inverse parallax distances
        self.useNBins = False
        self.parseClusterBins = False #toggle to run scotts rule binning around the cluster specific region separately
        self.sigamClusterBins = 5
        
        #outputs
        self.catalog = None
        self.members = None
        self.distcut = None
        self.dist = None
        self.hdbscanP = None
        self.membershipP = None
        self.MagCutmembers = None #keeps all things surviving the mag and radial cuts
        #note, depending on which is run first, PLX or PM will be the same as memebrs
        self.distmembers = None
        self.PMmembers = None

    def getGaiaData(self):
        if (self.verbosity > 0):
            print("Retrieving catalog ... ")
        cmd = f"SELECT * FROM gaiaedr3.gaia_source \
        WHERE CONTAINS(POINT('ICRS',gaiaedr3.gaia_source.ra, gaiaedr3.gaia_source.dec),\
        CIRCLE('ICRS', {self.RA}, {self.Dec}, {self.radius}))=1\
         AND pmra IS NOT NULL AND pmdec IS NOT NULL;"
        if (self.verbosity > 1):
            print(cmd)
        job = Gaia.launch_job_async(cmd, dump_to_file=False) #could save this to a file
        self.catalog = job.get_results()
        
        
        if (self.verbosity >0):
            print('\nUnique length: '+str(len(np.unique(self.catalog['source_id']))))
            print('\nNon-Unique length: '+str(len(self.catalog['source_id'])))
        
        cmd = f"SELECT source_id, r_med_geo, r_lo_geo, r_hi_geo, r_med_photogeo, r_lo_photogeo, r_hi_photogeo, \
        phot_bp_mean_mag-phot_rp_mean_mag AS bp_rp_bj, phot_g_mean_mag - 5 * LOG10(r_med_geo) + 5 AS qg_geo, \
        phot_g_mean_mag - 5 * LOG10(r_med_photogeo) + 5 AS gq_photogeo \
        FROM (SELECT * FROM gaiaedr3.gaia_source WHERE 1 = CONTAINS(POINT('ICRS',ra, dec), \
        CIRCLE('ICRS', {self.RA}, {self.Dec}, {self.radius}))) AS edr3 \
        JOIN external.gaiaedr3_distance using(source_id) WHERE ruwe< 200 "
        if (self.verbosity > 1):
            print('\n'+cmd)
            
        job = Gaia.launch_job_async(cmd, dump_to_file=False)
        bj = job.get_results()
        
        self.catalog = Tjoin(self.catalog, bj,keys='source_id',join_type="left")
        
        if (self.verbosity >0):
            print('\nUnique length: '+str(len(np.unique(self.catalog['source_id']))))
            print('\nNon-Unique length: '+str(len(self.catalog['source_id'])))
        
        # idlist = '('+', '.join(map(str, self.catalog['source_id'].data.data))+')'
        
        # cmd = f"SELECT * FROM gaiaedr3.dr2_neighbourhood \
        # WHERE dr3_source_id IN {idlist}"
        
        # if (self.verbosity > 1):
        #     print('\nAnd now adding in the dr2 source match from gaia.\n')
            
        # job = Gaia.launch_job_async(cmd, dump_to_file=False)
        # matchdr2 = job.get_results()
        # matchdr2 = matchdr2[np.where(matchdr2['angular_distance']<100)]
        
        # matchdr2['dr3_source_id'].name = 'source_id'
        # print(matchdr2)
        
        # self.catalog = Tjoin(self.catalog, matchdr2,keys='source_id',join_type="left")
        
        # if (self.verbosity >0):
        #     print('\nUnique length: '+str(len(np.unique(self.catalog['source_id']))))
        #     print('\nNon-Unique length: '+str(len(self.catalog['source_id'])))
        
        
        
        if self.Qgeodist == True: #photogeometric Bailer-Jones dist
            self.geodist = self.catalog['r_med_photogeo']
            #1 sigma (ish) see Bailer-Jones 2020(1) edr3 dist sect 2.6
            self.geodlo = self.catalog['r_lo_photogeo']
            self.geodhi = self.catalog['r_hi_photogeo']
        else:
            self.geodist = self.catalog['r_med_geo']
            #1 sigma (ish) see Bailer-Jones 2020(1) edr3 dist sect 2.6
            self.geodlo = self.catalog['r_lo_geo']
            self.geodhi = self.catalog['r_hi_geo']
        
    def getMembers(self):
        self.members = np.arange(len(self.catalog)) #if fitmembonly then probably need to redo the gaia pull here
        if (self.verbosity > 0):
            print("Retrieving members ... ")
            print("Gaia pulls "+str(len(self.members))+" possible sources.")
        
        #check that we need to add units
        if type(self.knowndist) != u.quantity.Quantity:
            self.knowndist = self.knowndist*u.parsec
            self.knowndisterr = self.knowndisterr*u.parsec
        
        self.knownplx = self.knowndist.to(u.mas, equivalencies=u.parallax())
        self.knownplxerr = (self.knowndisterr/self.knowndist)*self.knownplx
      
        self.plx = self.catalog['parallax'].to(u.mas)
        self.plxerr = self.catalog['parallax_error'].to(u.mas)
        self.dist = (self.catalog['parallax']).to(u.parsec, equivalencies=u.parallax())
        cdist = np.sqrt(np.abs(self.catalog['dec']-self.Dec)**2 + np.abs(self.catalog['ra']-self.RA)**2)
        self.catalog['cdist'] = cdist
        self.catalog['cdist'].__setattr__('description', 'Cluster center distance')
        
        
    def cutMembers(self):
        # is this where doubles are coming from? is it the gaia pull itself??
        if self.minMag is not None:
            if (self.verbosity > 0):
                print("Making a hard cut brighter than %.2f mag ... " % self.minMag)
                
            #check that we need to add units
            if type(self.minMag) != u.quantity.Quantity:
                self.minMag = self.minMag*u.mag
                
            
            #pull the Gaia G Mag
            Gmag = self.catalog['phot_g_mean_mag']
            "rig"
            #members brighter than
            membersbtMag = np.where(Gmag <= self.minMag)
            
            if self.fitmembonly == True: #this actually cuts the catalog down
                self.catalog = self.catalog[membersbtMag]
            else: #this just cuts what are members
                self.members = np.intersect1d(self.members, membersbtMag)
            
            if (self.verbosity > 0):
                print("This leaves us with "+str(len(self.members))+" possible sources.")
         
        if self.maxMag is not None:
            if (self.verbosity > 0):
                print("Making a hard cut fainter than %.2f mag ... " % self.maxMag)
  
            #check that we need to add units
            if type(self.maxMag) != u.quantity.Quantity:
                self.maxMag = self.maxMag*u.mag
            
            #pull the Gaia G Mag
            Gmag = self.catalog['phot_g_mean_mag']
            
            #members fainter than
            membersftMag = np.where(Gmag > self.maxMag)
            
            if self.fitmembonly == True: #this actually cuts the catalog down
                self.catalog = self.catalog[membersftMag].copy()
            else: #this just cuts what are members
                self.members = np.intersect1d(self.members, membersftMag)
                
                
            if (self.verbosity > 0):
                print("This leaves us with "+str(len(self.members))+" possible sources.")
        
        self.MagCutmembers = self.members
        
        if self.cdistCut is not None:
            self.innerCut = self.cdistCut[0]*self.catalog['cdist'].max()
            self.outerCut = self.cdistCut[1]*self.catalog['cdist'].max() #the enclosed radius is larger than set because of the uncertainty in ra & dec brings in some larger than the end.

            if (self.verbosity > 0):
                print("Making a radial cut at  %.2f to %.2f deg... " % (float(self.innerCut),float(self.outerCut)))
            
            if self.fitmembonly == True: #this actually cuts the catalog down
                self.catalog = self.catalog[np.intersect1d(np.where(self.catalog['cdist']<=self.outerCut),np.where(self.catalog['cdist']>self.innerCut))]
            else: #this just cuts what are members
                self.ringMems = np.intersect1d(np.where(self.catalog['cdist']<=self.outerCut),np.where(self.catalog['cdist']>self.innerCut))
                self.members = np.intersect1d(self.members, self.ringMems)
                self.MagCutmembers = np.intersect1d(self.MagCutmembers, self.ringMems)
            
        
    def getBJDistMembers(self):
        if (self.verbosity > 0):
            print("Removing Bailer-Jones distance non-members ...")

        gdist = self.geodist*u.parsec
        
        loerr = np.abs(self.geodist-self.geodlo)
        hierr = np.abs(self.geodhi-self.geodist)
        
        gdistmerr = (self.geodist - self.distsig*loerr)
        gdistperr = (self.geodist + self.distsig*hierr)
        
        dmin = self.knowndist-self.distsig*self.knowndisterr
        dmax = self.knowndist+self.distsig*self.knowndisterr
        
        ind = np.where(np.logical_and(dmin <= gdistperr, dmax > gdistmerr))
        indcut = np.where(~np.logical_and(dmin <= gdistperr, dmax > gdistmerr))
        mem = np.intersect1d(self.members, ind)
        distm = gdist[mem]
        
        #if verbosity = 3 include these extra funky plots?
        
        #1D histogram (use the members here)
        hpa, bpa = np.histogram(distm, bins = self.dbins, range=(self.dmin, self.dmax))

        if (self.showPlots):
            f = plt.figure()
            hpa, bpa, im = plt.hist(gdist, bins = self.dbins, histtype='step', fill=False, range=(self.dmin, self.dmax), linewidth=2,label='all',color='steelblue')
            
            #pop in an extra plot when cutting by mag
            if self.minMag is not None:
                hpa, bpa, im = plt.hist(gdist[self.members], bins = self.dbins, histtype='step', fill=False, range=(self.dmin, self.dmax), linewidth=2,label='after mag cut',color='rosybrown')

            hpa, bpa, im = plt.hist(distm, bins = self.dbins, histtype='step', fill=False, range=(self.dmin, self.dmax), linewidth=2,label='after dist cut',color='darkgoldenrod')
           
            plt.xlabel('distance (pc)', fontsize = 16)
            if self.titlestr is not None:
                plt.title(str(self.titlestr))
            plt.legend()
            f.patch.set_facecolor('white')
            f.patch.set_alpha(0.6)
            if self.titlestr is not None:
                plt.title(str(self.titlestr))
            plt.show()
        
        self.distmembers = mem
        self.distcut = indcut
        self.members = np.intersect1d(self.members, mem)
        
        if (self.verbosity > 0):
            print("This leaves us with "+str(len(self.members))+" possible sources.")
        
        
    def getParallaxMembers(self):
        if (self.verbosity > 0):
            print("Removing parallax non-members ... ")
            
        dist = self.dist.to(u.parsec).value
        
        plxmerr = (self.plx - self.distsig*self.plxerr)
        plxperr = (self.plx + self.distsig*self.plxerr)
        
        pmin = self.knownplx-self.distsig*self.knownplxerr
        pmax = self.knownplx+self.distsig*self.knownplxerr
        
        ind = np.where(np.logical_and(pmin <= plxperr, pmax > plxmerr))
        indcut = np.where(~np.logical_and(pmin <= plxperr, pmax > plxmerr))
        mem = np.intersect1d(self.members, ind)
        distm = dist[mem]
        
        #if verbosity = 3 include these extra funky plots
        
        #1D histogram (use the members here)
        hpa, bpa = np.histogram(distm, bins = self.dbins, range=(self.dmin, self.dmax))

        if (self.showPlots):
            f = plt.figure()
            hpa, bpa, im = plt.hist(dist, bins = self.dbins, histtype='step', fill=False, range=(self.dmin, self.dmax), linewidth=2,label='all',color='steelblue')
            
            #pop in an extra plot when cutting by mag
            if self.minMag is not None:
                hpa, bpa, im = plt.hist(dist[self.members], bins = self.dbins, histtype='step', fill=False, range=(self.dmin, self.dmax), linewidth=2,label='after mag cut',color='rosybrown')

            hpa, bpa, im = plt.hist(distm, bins = self.dbins, histtype='step', fill=False, range=(self.dmin, self.dmax), linewidth=2,label='after plx cut',color='darkgoldenrod')
           
            plt.xlabel('distance (pc)', fontsize = 16)
            if self.titlestr is not None:
                plt.title(str(self.titlestr))
            plt.legend()
            f.patch.set_facecolor('white')
            f.patch.set_alpha(0.6)
            if self.titlestr is not None:
                plt.title(str(self.titlestr))
            plt.show()
        
        self.distmembers = mem
        self.distcut = indcut
        self.members = np.intersect1d(self.members, mem)
        
        if (self.verbosity > 0):
            print("This leaves us with "+str(len(self.members))+" possible sources.")
        
    def getPMMembers(self,cmarkers=None):
        if (self.verbosity > 0):
            print("Finding proper-motion members with 2D-Gaussian method ...")
        
        self.catalog['membership_P']=0
        
        x = self.catalog['pmra']
        xW = 1./self.catalog['pmra_error']
        y = self.catalog['pmdec']
        yW = 1./self.catalog['pmdec_error']
                                                                           
        xm = x[self.members]
        xmW = xW[self.members]
        ym = y[self.members]
        ymW = yW[self.members]
        mW = 1./(np.sqrt((1./ymW)**2+(1./xmW)**2))
        #mW /= np.sum(mW)
        
        #1D histograms (use the members here)
        
        #initial bins
        if self.useNBins == True:
            pmRAbins = np.linspace(self.PMxmin, self.PMxmax, self.PMxbins)
            pmDecbins = np.linspace(self.PMymin, self.PMymax, self.PMybins)
        else: #by Scotts rule 1/(N^1/3)*sigma*3.5 with sigma as the median proper motion error in that bin
            #ra
            rng = self.PMxmax-self.PMxmin
            self.PMxbins = int(np.round(rng/(1/(len(self.members)**(1/3))*np.median(self.catalog['pmra_error'][self.members])*3.5)))
            pmRAbins = np.linspace(self.PMxmin,self.PMxmax,self.PMxbins)
            
            #dec
            rng = self.PMymax-self.PMymin
            self.PMybins = int(np.round(rng/(1/(len(self.members)**(1/3))*np.median(self.catalog['pmdec_error'][self.members])*3.5)))
            pmDecbins = np.linspace(self.PMymin,self.PMymax,self.PMybins)
            
                
    
        hx1D, x1D = np.histogram(xm, bins=pmRAbins)
        hy1D, y1D = np.histogram(ym, bins=pmDecbins)

        #2D histogram on the members so far
        h2D, x2D, y2D = np.histogram2d(xm, ym, bins=[pmRAbins, pmDecbins], \
                                       range=[[self.PMxmin, self.PMxmax], [self.PMymin, self.PMymax]])
                                       
        
        #fit on the members so far
        PMxguess = x1D[np.argmax(hx1D)]
        PMyguess = y1D[np.argmax(hy1D)]
        if (self.PMmean[0] != None):
            PMxguess = self.PMmean[0]
        if (self.PMmean[1] != None):
            PMyguess = self.PMmean[1]
       
        p_init = models.Gaussian2D(np.max(h2D.flatten()), PMxguess, PMyguess, self.PMxwdth, self.PMywdth)\
                + models.Gaussian2D(np.max(h2D.flatten()), self.PMfieldx, self.PMfieldy, self.PMfieldxwdth, self.PMfieldywdth)
        fit_p = fitting.LevMarLSQFitter()
        xf, yf = np.meshgrid(x2D, y2D, indexing='ij')
        xf = xf[:-1,:-1]
        yf = yf[:-1,:-1]
        
        pmG2D = fit_p(p_init, xf, yf, h2D,)
        if (self.verbosity > 1):
            print(pmG2D)
            print(pmG2D.parameters)
            
        self.GaussParams = pmG2D.parameters
        
        Fc = models.Gaussian2D()
        Fc.parameters = pmG2D.parameters[0:6]

        #Get initial fit params to pull out the cluster
        xFcMean, yFcMean = Fc.parameters[1], Fc.parameters[2]
        xFcStd, yFcStd = Fc.parameters[3], Fc.parameters[4]
        
        #here we compare the median of the measurment errors in the current selection with the fit of the cluster standard deviation
        med_decerr = np.median(self.catalog['pmdec_error'][self.members])
        med_raerr = np.median(self.catalog['pmra_error'][self.members])
        
        #if the fit std is > the observed error then use the fit std, not the observed error cause in this case we are fitting the cluster dispersion
        if ((yFcStd > 1.5*med_raerr) or (xFcStd > 1.5*med_decerr)):#(xFcStd+yFcStd)/2.0 > (med_decerr+med_raerr)/2.0:
            usefitstd = True
            print('Using the initial cluster fit standard deviation since',str(xFcStd),str(yFcStd),'is greater than',str(med_decerr),str(med_raerr))
            print('Using the initial cluster fit standard deviation since',str((xFcStd+yFcStd)/2.0),'is greater than',str((med_decerr+med_raerr)/2.0))
        else:
            usefitstd = False
            

if self.parseClusterBins == True:
#have to get the bin edges that are closest to the cluster range so we dont double up on members
try:
    #RA
    clusterRAregion = (pmRAbins>(xFcMean-self.sigamClusterBins*xFcStd))&(pmRAbins<=(xFcMean+self.sigamClusterBins*xFcStd))
    clusterRAregionIndx = np.where(clusterRAregion)[0]
    lowboundRA = pmRAbins[clusterRAregionIndx[0]]
    lowRABins = pmRAbins[:clusterRAregionIndx[0]]
    upRABins = pmRAbins[clusterRAregionIndx[-1]+1:]
    upboundRA = pmRAbins[clusterRAregionIndx[-1]+1]
    
    #Dec
    clusterDecregion =(pmDecbins>(yFcMean-self.sigamClusterBins*yFcStd))&(pmDecbins<=(yFcMean+self.sigamClusterBins*yFcStd))
    clusterDecregionIndx = np.where(clusterDecregion)[0]
    lowboundDec = pmDecbins[clusterDecregionIndx[0]]
    lowDecBins = pmDecbins[:clusterDecregionIndx[0]]
    upDecBins = pmDecbins[clusterDecregionIndx[-1]+1:]
    upboundDec = pmDecbins[clusterDecregionIndx[-1]+1]
    
    #subset of members in the bounds
    m = self.catalog[self.members]
    mC = (xm>lowboundRA)&(xm<=upboundRA)&(ym<=upboundDec)&(ym>lowboundDec)
    xmC = xm[mC]
    ymC = ym[mC]
    
    #ra by Scotts rule 1/(N^1/3)*sigma*3.5 with sigma as the median proper motion error in that bin
    rng = upboundRA-lowboundRA
    if usefitstd == True:
        self.PMxbinsC = int(np.round(rng/(1/(len(mC)**(1/3))*xFcStd*3.5)))
    else:
        self.PMxbinsC = int(np.round(rng/(1/(len(mC)**(1/3))*np.median(self.catalog['pmra_error'][self.members][mC])*3.5)))
    pmRAbinsC = np.linspace(lowboundRA,upboundRA,self.PMxbinsC)
    
    #dec by Scotts rule 1/(N^1/3)*sigma*3.5 with sigma as the median proper motion error in that bin
    rng = upboundDec-lowboundDec
    
    if usefitstd == True:
        self.PMybinsC = int(np.round(rng/(1/(len(mC)**(1/3))*yFcStd*3.5)))
    else:
        self.PMybinsC = int(np.round(rng/(1/(len(mC)**(1/3))*np.median(self.catalog['pmdec_error'][self.members][mC])*3.5)))
    pmDecbinsC = np.linspace(lowboundDec,upboundDec,self.PMybinsC)
    
    
    #make new bins
    pmRAbins = np.concatenate((lowRABins,pmRAbinsC,upRABins[1:]),axis=0)
    self.PMxbins = pmRAbins
    pmDecbins = np.concatenate((lowDecBins,pmDecbinsC,upDecBins[1:]),axis=0)
    self.PMybins = pmDecbins
    
    hx1D, x1D = np.histogram(xm, bins=pmRAbins)
    hy1D, y1D = np.histogram(ym, bins=pmDecbins)

    #2D histogram on the members so far
    h2D, x2D, y2D = np.histogram2d(xm, ym, bins=[self.PMxbins, self.PMybins], \
                                   range=[[self.PMxmin, self.PMxmax], [self.PMymin, self.PMymax]])

    #fit on the members so far
    PMxguess = x1D[np.argmax(hx1D)]
    PMyguess = y1D[np.argmax(hy1D)]
    if (self.PMmean[0] != None):
        PMxguess = self.PMmean[0]
    if (self.PMmean[1] != None):
        PMyguess = self.PMmean[1]
   
    p_init = models.Gaussian2D(np.max(h2D.flatten()), PMxguess, PMyguess, self.PMxwdth, self.PMywdth)\
            + models.Gaussian2D(np.max(h2D.flatten()), self.PMfieldx, self.PMfieldy, self.PMfieldxwdth, self.PMfieldywdth)
    fit_p = fitting.LevMarLSQFitter()
    xf, yf = np.meshgrid(x2D[:-1], y2D[:-1], indexing='ij')
    pmG2D = fit_p(p_init, xf, yf, h2D,)
    if (self.verbosity > 1):
        print(pmG2D)
        print(pmG2D.parameters)

    self.GaussParams = pmG2D.parameters
    
    Fc = models.Gaussian2D()
    Fc.parameters = pmG2D.parameters[0:6]
except:
    print('ATTN: dynamic binning failed, using orignal')

        #membership calculation, apply the fit here to the whole catalog
        PPM = Fc(x,y)/pmG2D(x,y)
        self.pmg2d = pmG2D
        self.fc = Fc
            
        PPM[self.distcut] = PPM[self.distcut]*self.distcutP
        self.catalog['membership_P']=PPM
        
        

        if (self.showPlots):
            f = plt.figure(figsize=(8, 8))
            f.patch.set_facecolor('white')
            f.patch.set_alpha(0.6)
            gs = gridspec.GridSpec(2, 2, height_ratios = [1, 3], width_ratios = [3, 1])
            ax1 = plt.subplot(gs[0])
            ax2 = plt.subplot(gs[2])
            ax3 = plt.subplot(gs[3])

            #histograms
        
            hx1D, x1D = np.histogram(xm, bins=pmRAbins)
            hy1D, y1D = np.histogram(ym, bins=pmDecbins)
            ax1.step(x1D[:-1], hx1D)
            ax1.plot(x2D[:-1], np.sum(pmG2D(xf, yf), axis=1), color='red')
            ax3.step(hy1D, y1D[:-1])
            ax3.plot(np.sum(pmG2D(xf, yf), axis=0), y2D[:-1], color='red')

            #heatmap
            if self.limitedC == True:
                h2D, x2D, y2D, im = ax2.hist2d(xm, ym, bins=[self.PMxbins, self.PMybins],\
                                               range=[[self.PMxmin, self.PMxmax], [self.PMymin, self.PMymax]], \
                                               norm = mpl.colors.LogNorm(), cmap = 'viridis', vmin=2, vmax=6)
                cb = f.colorbar(im)
                cb.set_label("N in Bin", rotation=270)
            else:
                 h2D, x2D, y2D, im = ax2.hist2d(xm, ym, bins=[self.PMxbins, self.PMybins],\
                                           range=[[self.PMxmin, self.PMxmax], [self.PMymin, self.PMymax]], \
                                           norm = mpl.colors.LogNorm(), cmap = cm.Blues)
            if cmarkers is not None:
                ax2.scatter(self.catalog['pmra'][mems], self.catalog['pmdec'][mems], c=self.catalog[cmarkers][mems], cmap='magma',s=1,alpha=.2,label='Members')
            try:
                ax2.contourf(x2D[:-1], y2D[:-1], pmG2D(xf, yf).T, cmap=cm.Reds, bins = 20, \
                             norm=mpl.colors.LogNorm(), alpha = 0.3)
            except:
                print('contour error')
#
            ax1.set_xlim(self.PMxmin, self.PMxmax)
            ax2.set_xlim(self.PMxmin, self.PMxmax)
            ax2.set_ylim(self.PMymin, self.PMymax)
            ax3.set_ylim(self.PMymin, self.PMymax)
            ax1.set_ylim(1, 2*max(hx1D))
            if self.histlogscale == True:
                ax1.set_yscale("log")
                ax3.set_xscale("log")
            ax3.set_xlim(1, 2*max(hy1D))
            ax2.set_xlabel(r'$\mu_\alpha$ (mas yr$^{-1}$)', fontsize=16)
            ax2.set_ylabel(r'$\mu_\delta$ (mas yr$^{-1}$)', fontsize=16)
            plt.setp(ax1.get_yticklabels()[0], visible=False)
            plt.setp(ax1.get_xticklabels(), visible=False)
            plt.setp(ax3.get_yticklabels(), visible=False)
            plt.setp(ax3.get_xticklabels()[0], visible=False)
            f.subplots_adjust(hspace=0., wspace=0.)
            f.patch.set_facecolor('white')
            f.patch.set_alpha(0.6)
            if self.titlestr is not None:
                plt.title(str(self.titlestr))
            plt.show()
#
            
            #plot
            f = plt.figure()
            b, h, im = plt.hist(PPM, bins = 100, histtype='step', fill=False, range=(0,1), linewidth=2)
            plt.yscale('log')
            plt.xlabel('probability')
            f.patch.set_facecolor('white')
            f.patch.set_alpha(0.6)
            if self.titlestr is not None:
                plt.title(str(self.titlestr))
            plt.show()


        self.PMmembers = np.where(np.logical_and(PPM > self.minPMmembership, self.catalog['pmra'].mask == False))
        membersPMAll = np.where(PPM > self.minPMmembership)
        self.members = np.intersect1d(self.members, membersPMAll)
        
    def getPMMembersHDB(self,cmarkers=None):
        if (self.verbosity > 0):
            print("Finding proper-motion members with HDBSCAN...")
        
        self.catalog['hdbscanP'] = None
        self.catalog['hdbscanLabel'] = None
        self.catalog['hdbscanID'] = None
       # self.catalog['memb_P_HDBSCAN']=0
        #start with the magnitude range
        mems=self.MagCutmembers
        
        pdcat = self.catalog[mems].to_pandas(index='source_id')
        
        #here is the toggles to include different aspects of the data for the hdbscan clustering
#        if self.plxdist == True:
        data = pdcat[['pmra','pmdec','ra','dec']].to_numpy()
        #data = pdcat[['pmra','pmdec','parallax']].to_numpy()
#        else:
#            if self.Qgeodist == True: #photogeometric Bailer-Jones dist
#                data = pdcat[['pmra','pmdec','r_med_photogeo']].to_numpy()
#            else:
#                data = pdcat[['pmra','pmdec','r_med_geo']].to_numpy()

        clusterer=hdbscan.HDBSCAN(min_cluster_size=self.hdbminN, metric='euclidean')
        clusterer.fit_predict(data)
        
        if (self.showPlots):
            f = plt.figure(figsize=(10,8))
            plt.scatter(pdcat.pmra,pdcat.pmdec,c=clusterer.labels_,cmap='coolwarm',s=1)#parallax in and RA & dec
            if self.GaussParams is not None:
                plt.scatter(self.GaussParams[1],self.GaussParams[2],marker='+',color='yellow',s=100)
            plt.xlim(self.PMxmin, self.PMxmax)
            plt.ylim(self.PMymin, self.PMymax)
            plt.xlabel('pmra')
            plt.ylabel('pmdec')
            f.patch.set_facecolor('white')
            f.patch.set_alpha(0.6)
            if self.titlestr is not None:
                plt.title(str(self.titlestr))
            plt.show()

        npops = [len(pdcat.ra[clusterer.labels_==-1]),len(pdcat.ra[clusterer.labels_==0]),len(pdcat.ra[clusterer.labels_==1])]
    
        pdcat['hdbscanP']=clusterer.probabilities_
        pdcat['hdbscanLabel']=clusterer.labels_
        pdcat['hdbscanID'] = ' '


        
        
        if len(clusterer.labels_>2):
            labelids = [' ',' ',' ']
            
            if self.hdbscanCID is None:
                no = np.unique(clusterer.labels_)[int(np.argwhere(npops==np.amax(npops)))]
                uk = np.unique(clusterer.labels_)[int(np.argwhere(npops==np.amin(npops)))]
                cl = np.unique(clusterer.labels_)[int(np.argwhere((npops>np.amin(npops))*(npops != np.amax(npops))))]
            else:
                cl = self.hdbscanCID
                no = self.hdbscanFID
                uk = self.hdbscanUID
                
            labelids[uk+1] = str('Unknown, %i'%uk)
            pdcat['hdbscanID'][pdcat['hdbscanLabel']==uk] = 'Unknown'
            labelids[no+1] = str('Field, %i'%no)
            pdcat['hdbscanID'][pdcat['hdbscanLabel']==no] = 'Field'
            labelids[cl+1] = str('Cluster, %i'%cl)
            pdcat['hdbscanID'][pdcat['hdbscanLabel']==cl] = 'Cluster'

        else:
            labelids = [' ',' ']
            if self.hdbscanCID is None:
                no = np.unique(clusterer.labels_)[int(np.argwhere(npops==np.amax(npops)))]
                cl = np.unique(clusterer.labels_)[int(np.argwhere(npops<np.amax(npops))[-1])]
            else:
                cl = self.hdbscanCID
                no = self.hdbscanFID
                uk = self.hdbscanUID
            labelids[no+1] = str('Field, %i'%no)
            pdcat['hdbscanID'][pdcat['hdbscanLabel']==no] = 'Field'
            labelids[cl+1] = str('Cluster, %i'%cl)
            pdcat['hdbscanID'][pdcat['hdbscanLabel']==cl] = 'Cluster'

        self.HDBParams = (labelids,npops)

        if (self.verbosity > 1):
            print('HDBSCAN Populations: ',self.HDBParams)
        
       # self.catalog['memb_P_HDBSCAN'][np.where(pdcat['hdbscanID']=='Cluster')] = pdcat[pdcat['hdbscanID']=='Cluster'].hdbscanP.values
        #print('wtf')
        #print(clusterer.probabilities_)
       # print(pdcat[pdcat['hdbscanID']=='Cluster'].hdbscanP.values)
        #print('wtf')
        #PPM  = self.catalog['memb_P_HDBSCAN']
        
        #PPM[self.distcut] = PPM[self.distcut]*self.distcutP
       # self.catalog['memb_P_HDBSCAN']=PPM
        
        if (self.showPlots):
            #plot
            f = plt.figure()
            b, h, im = plt.hist(clusterer.probabilities_, bins = 100, histtype='step', fill=False, range=(0,1), linewidth=2)
            plt.yscale('log')
            plt.xlabel('probability')
            f.patch.set_facecolor('white')
            f.patch.set_alpha(0.6)
            if self.titlestr is not None:
                plt.title(str(self.titlestr))
            plt.show()

        self.PMmembersHDB = np.where(np.logical_and(clusterer.probabilities_ > self.minPMmembershipHDB, self.catalog[mems]['pmra'].mask == False))
        #membersPMAll = np.where(PPM > self.minPMmembershipHDB)
        #self.members = np.intersect1d(self.members, membersPMAll)
        self.catalog['hdbscanP'][mems]=pdcat.hdbscanP.values
        self.catalog['hdbscanLabel'][mems]=pdcat.hdbscanLabel.values
        self.catalog['hdbscanID'][mems] = pdcat.hdbscanID.values
        
        
    def plotPMOneMemb(self,sourceid,cmarkers=None):
        #plot the proper motion plot with one member marked with X
        x = self.catalog['pmra']
        y = self.catalog['pmdec']
        print("Gaia pull gave "+str(len(self.catalog))+" possible sources.")
        
        source = self.catalog[self.catalog['source_id']==sourceid]
        xsource = source['pmra']
        ysource = source['pmdec']
               
        print("We fit with "+str(len(self.distmembers))+" possible sources after parallax and mag cut.")
        
        ind = np.where(self.catalog['parallax'].to(u.parsec, equivalencies=u.parallax()).value < self.dmax)
        mem = np.intersect1d(self.distmembers, ind)
        xm = x[mem]
        ym = y[mem]
        
        #1D histograms (use the members here)
        pmRAbins = np.linspace(self.PMxmin, self.PMxmax, self.PMxbins)
        pmDecbins = np.linspace(self.PMymin, self.PMymax, self.PMybins)
        hx1D, x1D = np.histogram(xm, bins=pmRAbins)
        hy1D, y1D = np.histogram(ym, bins=pmDecbins)

        #2D histogram (use the members here)
        h2D, x2D, y2D = np.histogram2d(xm, ym, bins=[self.PMxbins, self.PMybins], \
                                       range=[[self.PMxmin, self.PMxmax], [self.PMymin, self.PMymax]])
            
        #fit on the members
        PMxguess = x1D[np.argmax(hx1D)]
        PMyguess = y1D[np.argmax(hy1D)]
        if (self.PMmean[0] != None):
            PMxguess = self.PMmean[0]
        if (self.PMmean[1] != None):
            PMyguess = self.PMmean[1]
       
        p_init = models.Gaussian2D(np.max(h2D.flatten()), PMxguess, PMyguess, self.PMxwdth, self.PMywdth)\
                + models.Gaussian2D(np.max(h2D.flatten()), self.PMfieldx, self.PMfieldy, self.PMfieldxwdth, self.PMfieldywdth)
        fit_p = fitting.LevMarLSQFitter()
        xf, yf = np.meshgrid(x2D[:-1], y2D[:-1], indexing='ij')
        pmG2D = fit_p(p_init, xf, yf, h2D)
        if (self.verbosity > 1):
            print(pmG2D)
            print(pmG2D.parameters)
            
        if (self.showPlots):
            f = plt.figure(figsize=(8,8))
            gs = gridspec.GridSpec(2, 2, height_ratios = [1, 3], width_ratios = [3, 1])
            ax1 = plt.subplot(gs[0])
            ax2 = plt.subplot(gs[2])
            ax3 = plt.subplot(gs[3])

            #histograms of full pull
            hx1D, x1D = np.histogram(xm, bins=pmRAbins)
            hy1D, y1D = np.histogram(ym, bins=pmDecbins)
            ax1.step(x1D[:-1], hx1D)
            ax1.plot(x2D[:-1], np.sum(pmG2D(xf, yf), axis=1), color='red')
            ax3.step(hy1D, y1D[:-1])
            ax3.plot(np.sum(pmG2D(xf, yf), axis=0), y2D[:-1], color='red')

            #heatmap
            if self.limitedC == True:
                h2D, x2D, y2D, im = ax2.hist2d(xm, ym, bins=[self.PMxbins, self.PMybins],\
                                               range=[[self.PMxmin, self.PMxmax], [self.PMymin, self.PMymax]], \
                                               norm = mpl.colors.LogNorm(), cmap = 'viridis', vmin=2, vmax=6)
               
                cb = f.colorbar(im)
                cb.set_label("N in Bin", rotation=270)
            else:
                h2D, x2D, y2D, im = ax2.hist2d(xm, ym, bins=[self.PMxbins, self.PMybins],\
                                               range=[[self.PMxmin, self.PMxmax], [self.PMymin, self.PMymax]], \
                                               norm = mpl.colors.LogNorm(), cmap = cm.Blues)
            if cmarkers is not None:
                ax2.scatter(self.catalog['pmra'][mems], self.catalog['pmdec'][mems], c=self.catalog[cmarkers][mems], cmap='magma',s=1,alpha=.2,label='Members')
            else:
                ax2.scatter(xsource,ysource,marker='o',edgecolor='k',facecolor='None',s=20)
            ax2.contourf(x2D[:-1], y2D[:-1], pmG2D(xf, yf).T, cmap=cm.Reds, bins = 20, \
                         norm=mpl.colors.LogNorm(), alpha = 0.3)

            ax1.set_xlim(self.PMxmin, self.PMxmax)
            ax2.set_xlim(self.PMxmin, self.PMxmax)
            ax2.set_ylim(self.PMymin, self.PMymax)
            ax3.set_ylim(self.PMymin, self.PMymax)
            ax1.set_ylim(1, 2*max(hx1D))
            if self.histlogscale == True:
                ax1.set_yscale("log")
                ax3.set_xscale("log")
            ax3.set_xlim(1, 2*max(hy1D))
            ax2.set_xlabel(r'$\mu_\alpha$ (mas yr$^{-1}$)', fontsize=16)
            ax2.set_ylabel(r'$\mu_\delta$ (mas yr$^{-1}$)', fontsize=16)
            plt.setp(ax1.get_yticklabels()[0], visible=False)
            plt.setp(ax1.get_xticklabels(), visible=False)
            plt.setp(ax3.get_yticklabels(), visible=False)
            plt.setp(ax3.get_xticklabels()[0], visible=False)
            f.subplots_adjust(hspace=0., wspace=0.)
            f.patch.set_facecolor('white')
            f.patch.set_alpha(0.6)
            if self.titlestr is not None:
                plt.title(str(self.titlestr))
            plt.show()
        
    def plotCMD(self,sid=None,cmarkers=None):
        f = plt.figure(figsize=(5,6))
        plt.scatter(self.catalog['g_rp'], self.catalog['phot_g_mean_mag'], s = 1,  color='steelblue', alpha = 0.3, label='All')
        if cmarkers is not None:
            plt.scatter(self.catalog['g_rp'][self.members], self.catalog['phot_g_mean_mag'][self.members], s = 2, c=self.catalog[cmarkers][self.members], cmap='magma',label='Members')
            plt.colorbar()
        else:
            plt.scatter(self.catalog['g_rp'][self.members], self.catalog['phot_g_mean_mag'][self.members], s = 2, color='red',  label='Members')
        if sid is not None:
            source = self.catalog[np.isin(self.catalog['source_id'].tolist(),sid)]
            plt.scatter(source['g_rp'],source['phot_g_mean_mag'] ,marker='o',edgecolor='k',facecolor='None',s=20,label='Source')
        plt.legend()
        plt.xlim(self.CMDxmin, self.CMDxmax)
        plt.ylim(self.CMDymin, self.CMDymax)
        plt.xlabel('G - RP (mag)', fontsize=16)
        plt.ylabel('G (mag)', fontsize=16)
        f.patch.set_facecolor('white')
        f.patch.set_alpha(0.6)
        if self.titlestr is not None:
            plt.title(str(self.titlestr))
        plt.show()
          
    def printSourceP(self,sourceid):
        source = self.catalog[self.catalog['source_id']==sourceid]
        print("Source %s has P = %e.3"%(str(sourceid),source['membership_P']))
        
    def printSource(self,sourceid,gaiaparams=None):
        source = self.catalog[self.catalog['source_id']==sourceid]
        print("Source %s has P = %.4f"%(str(sourceid),source['membership_P']))
        if gaiaparams is not None:
            for kw in gaiaparams:
                print("%s"%str(source[kw]))
        if (self.showPlots):
            self.plotPMOneMemb(sourceid)
            self.plotCMD(sid=sourceid)
        
    def runAll(self):
        self.getGaiaData()
        self.getMembers()
        self.cutMembers()
        
        if self.plxdist == True:
            self.getParallaxMembers()
        else:
            self.getBJDistMembers()
        
        if self.hdbPMM == True:
            self.getPMMembersHDB()
        elif self.bothPMM == True:
            self.getPMMembers()
            self.getPMMembersHDB()
        else:
            self.getPMMembers()
            
        if (self.showPlots):
            self.plotCMD()
        if (self.verbosity > 1):
            print("We end up with "+str(len(self.members))+" sources.")
            print("We're all done here pal, how's she look?")
        elif (self.verbosity > 0):
            print("We end up with "+str(len(self.members))+" sources.")
        
    def runMembs(self):
        self.getMembers()
        self.cutMembers()
        
        if self.plxdist == True:
            self.getParallaxMembers()
        else:
            self.getBJDistMembers()
        
        if self.hdbPMM == True:
            self.getPMMembersHDB()
        elif self.bothPMM == True:
            self.getPMMembers()
            self.getPMMembersHDB()
        else:
            self.getPMMembers()
            
        if (self.showPlots):
            self.plotCMD()
        if (self.verbosity > 1):
            print("We end up with "+str(len(self.members))+" sources.")
            print("Is she looking better or worse? Be honest.")
        elif (self.verbosity > 0):
            print("We end up with "+str(len(self.members))+" sources.")
