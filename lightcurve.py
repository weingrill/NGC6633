#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Aug 11, 2017

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import config
import numpy as np
import logging
#import matplotlib.pyplot as plt
from astropy.stats import LombScargle
import pickle

logging.basicConfig(filename=config.projectpath+'lightcurve.log', 
                    format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('NGC6633 lightcurve')

class LightCurve(object):
    """
    Lightcurve object for NGC6633
    fetch lightcurve from file.
    ability to store lightcurve to file
    """
    def __init__(self, starid):
        self.starid = starid
        self.fromfile()
        
        self.pdm_period = np.NAN
        self.pdm_error = np.NAN
        self.pdm_theta = np.NAN

    def fromfile(self, filename=None):
        """
        loads the lightcurve from a file.
        if the filename is not given it is assembled from the starid
        """
        
        if filename is None:
            filename = config.lightcurvespath+self.starid+'.dat'
        logger.info('load file %s' % filename)
        self.hjd, self.mag, self.err = np.loadtxt(filename, unpack = True)
        
        logger.info('%d datapoints' % len(self.hjd))
        
        return (self.hjd, self.mag, self.err)
    
    def normalize(self):
        try:
            self.hjd -= np.min(self.hjd)
            self.mag -= np.mean(self.mag)
        except ValueError:
            return
    
    def detrend(self):
        try:
            par = np.polyfit(self.hjd, self.mag, 1)
            self.mag -= np.polyval(par, self.hjd)
        except TypeError:
            return

    def sigma_clip(self, sigmas=2.0):
        from functions import sigma_clip
        self.hjd, self.mag, self.err = sigma_clip(self.hjd, self.mag, self.err, sigmas=sigmas)
    
    def clip(self, limit = 0.05):
        """
        clip lightcurve at given limit
        """
        m = np.mean(self.mag)
        valid = abs(self.mag-m)<limit
        self.hjd = np.compress(valid, self.hjd)
        self.mag = np.compress(valid, self.mag)
        self.err = np.compress(valid, self.err)
        
    def __len__(self):
        return len(self.hjd)
    
    def psd(self, minperiod, maxperiod):
        from psd import ppsd
        # perform a power spectrum analysis
        tpsa, mpsa = self.hjd, self.mag
        n = len(tpsa)
        # zero padded lightcurves
        t_padded = np.zeros(4*n)
        t_padded[:n] = tpsa
        t_padded[n:] = np.linspace(max(tpsa),4*max(tpsa),3*n)
        m_padded = np.zeros(4*n)
        m_padded[:n] = mpsa
        
        px, f = ppsd(t_padded, 
                     m_padded, 
                     lower=1./maxperiod, 
                     upper=1./minperiod,
                     num= 2000)
        px = np.sqrt(px)
        period = 1./f[np.argmax(px)]
        return period
    
    def pdm(self, minperiod, maxperiod):
        """
        calculate the phase dispersion minimization with uncertainties.
        returns period, uncertainty of period and theta value of period
        taken from pwkit
        see also: http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyTimingDoc/pyPDMDoc/pdm.html
        """
        from pwkit import pdm
        # look at 20 days or at most at the length of dataset
        import os
        filename = os.path.join(config.lightcurvespath,'pdm',self.starid+'.pdm')
        try:
            result = pickle.load(open(filename,'r'))
        except IOError:
            # estimate the number of periods to be tested
            n = np.round(len(self.hjd)+50,-2)
            #TODO: replace by numpy.geomspace
            periods = np.linspace(1.0,21.0, n)
            result = pdm.pdm(self.hjd, self.mag, self.err, periods, 3)
            
            pickle.dump(result, open(filename,'w')) 
        self.pdm_period = result.pmin
        self.pdm_error = result.mc_puncert
        self.pdm_theta = result.thetas[result.imin]           
        return periods, result.thetas
    
    def lomb_scargle(self, minfrequency=0.05, maxfrequency=0.8, numfrequencies=10000):
        """
        calculate the Lomb-Scargle periodogram
        returns period, period uncertainty, power
        """
        import os
        import functions as fx
        filename = os.path.join(config.lightcurvespath,'ls',self.starid+'.ls')
        try:
            frequencies, power, window = np.genfromtxt(filename, unpack=True)
        except IOError:
            # complute the window
            win = np.ones(np.shape(self.mag))
            
            frequencies = np.linspace(minfrequency, maxfrequency, numfrequencies)
            power = LombScargle(self.hjd, self.mag, self.err).power(frequencies)
            window = LombScargle(self.hjd, win, center_data=False).power(frequencies)
            
            a = np.column_stack((frequencies, power, window))
            np.savetxt(filename, a)
        best_frequency = frequencies[np.argmax(power)]
        period = 1./ best_frequency
        
        periods = 1.0 / frequencies    
        i = np.argmax(power)
        i1 = fx.largmin(power, i)
        i2 = fx.largmin(power, i,'right')
        
        _, _, ls_sigma = fx.gauss_fit(periods[i1:i2], power[i1:i2],  power[i], period, 1.0)    
        self.ls_period, self.ls_error, self.ls_power = period, ls_sigma, np.max(power)      
        return frequencies, power
        
    def clean(self,minperiod=1.0, maxperiod=20.0, gain=0.1):
        """
        calculates the CLEANed spectrum
        returns period at maximum amplitude, period uncertainty, maximum amplitude 
        """
        
        from clean import clean  # @UnresolvedImport
        import os
        import functions as fx
        filename = os.path.join(config.lightcurvespath,'clean',self.starid+'.clean')
        t = self.hjd
        x = self.mag
        try:
            p, a = np.genfromtxt(filename, unpack=True)
        except IOError:
            f, a, _ = clean(t, x, gain=gain, threshold=1e-3)
            p = 1.0/f[f>0]
            a = a[f>0]
            
            save_array = np.column_stack((p, a))
            np.savetxt(filename, save_array)
        # sort the periods and amplitudes, because, they are in reversed order
        periods = p[::-1]
        assert(periods[0]<periods[-1])
        amplitudes = a[::-1]
        j = np.where((periods>=minperiod) & (periods<maxperiod))
        periods = periods[j]
        amplitudes = amplitudes[j]
        
        i = np.argmax(amplitudes)
        period = periods[i]
        i1 = fx.largmin(amplitudes, i)
        i2 = fx.largmin(amplitudes, i,'right')
        sigma0 = np.abs(periods[i1]-periods[i2])/2.0
        try:
            _, _, sigma = fx.gauss_fit(periods[i1:i2], amplitudes[i1:i2],  amplitudes[i], period, sigma0)
        except RuntimeError:
            logger.warn('cannot determine sigma by Gauss fit in CLEAN')
            sigma = sigma0
        self.clean_period, self.clean_error, self.clean_amplitude = period, sigma, np.max(amplitudes)        
        return periods, amplitudes
            
    def phased(self, period):
        from functions import phase
        tp, yp = phase(self.hjd, self.mag, period)

            
        s1 = np.sin(2*np.pi*tp/period)
        c1 = np.cos(2*np.pi*tp/period)
        s2 = np.sin(4*np.pi*tp/period)
        c2 = np.cos(4*np.pi*tp/period)
        
        A = np.column_stack((np.ones(tp.size), s1, c1, s2, c2))
        c, resid,_,_ = np.linalg.lstsq(A,yp)
        amp_err = resid[0]

        amp = max(c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2)-\
              min(c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2)

        tp1 = np.linspace(0.0, period, 100)
        s1 = np.sin(2*np.pi*tp1/period)
        c1 = np.cos(2*np.pi*tp1/period)
        s2 = np.sin(4*np.pi*tp1/period)
        c2 = np.cos(4*np.pi*tp1/period)
        return amp, amp_err

    def rebin(self, interval = 0.5, method='mean'):
        """
        rebin to new interval using the mean of each bin
        interval determines the length of each bin
        medianbins calculates the median in each bin otherwise the mean is taken
        https://people.ucsc.edu/~ianc/python/_modules/tools.html#errxy
        """
        
        
        if method not in ['mean','median','weighted']:
            raise ValueError
        
        data = self.hjd
        # ...+interval so that the last bin includes the last epoch
        bins = np.arange(self.hjd[0]-0.5*interval, self.hjd[-1]+0.5*interval, interval)
        nbins = len(bins)-1
        t = np.zeros(nbins)
        f = np.zeros(nbins)
        f_mean = np.zeros(nbins)
        f_median = np.zeros(nbins)
        f_weight = np.zeros(nbins)
        e = np.zeros(nbins)
        # adopted from Ian's Astro-Python Code v0.3
        # https://people.ucsc.edu/~ianc/python/_modules/tools.html#errxy
        idx = [[data.searchsorted(bins[i]), \
                data.searchsorted(bins[i+1])] for i in range(nbins)]
        np.seterr(invalid='ignore')
        for i in range(nbins):
            f_mean[i] = np.mean(self.mag[idx[i][0]:idx[i][1]])
            #f_median[i] = np.median(self.mag[idx[i][0]:idx[i][1]])
            #f_weight[i] = np.average(self.mag[idx[i][0]:idx[i][1]], weights=1./self.err[idx[i][0]:idx[i][1]])
            t[i] = np.mean(self.hjd[idx[i][0]:idx[i][1]])
            e[i] = np.std(self.hjd[idx[i][0]:idx[i][1]])

        
        if method=='mean':
            f = f_mean
        elif method=='median':
            f = f_median
        elif method=='weighted':
            f = f_weight
        
                    
        np.seterr(invalid='warn')
        valid = ~np.isnan(t)
        self.mag = np.compress(valid,f)
        self.hjd = np.compress(valid,t)
        self.err = np.compress(valid,e)

