#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 15, 2017

@author: Joerg Weingrill <jweingrill@aip.de>
'''

import config
import numpy as np
import logging
import matplotlib.pyplot as plt
from astropy.stats import LombScargle

logging.basicConfig(filename=config.projectpath+'analysis.log', 
                    format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('NGC6633 analysis')

class LightCurve(object):
    """
    Lightcurve object for NGC6633
    fetch lightcurve from file.
    ability to store lightcurve to file
    """
    def __init__(self, starid):
        self.starid = starid
        self.fromfile()

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
        from pdm import pdm
        # look at 20 days or at most at the length of dataset
        import os
        filename = os.path.join(config.lightcurvespath,'pdm',self.starid+'.pdm')
        try:
            
            periods, thetas = np.genfromtxt(filename, unpack=True)
        except IOError:
            periods, thetas = pdm(self.hjd, self.mag, minperiod, maxperiod, 0.5/24)
            a = np.column_stack((periods, thetas))
            np.savetxt(filename, a)
        return periods, thetas
    
    def lomb_scargle(self, minfrequency=0.05, maxfrequency=0.8):
        """
        calculate the Lomb-Scargle periodogram
        """
        import os
        filename = os.path.join(config.lightcurvespath,'ls',self.starid+'.ls')
        try:
            frequencies, power, window = np.genfromtxt(filename, unpack=True)
        except IOError:
            t = self.hjd
            t -= t[0]
            y1 = np.ones(np.shape(self.mag))
            
            frequencies = np.linspace(minfrequency, maxfrequency, 10000)
            power = LombScargle(self.hjd, self.mag, self.err).power(frequencies)
            window = LombScargle(self.hjd, y1, center_data=False).power(frequencies)
            
            a = np.column_stack((frequencies, power, window))
            np.savetxt(filename, a)
        return frequencies, power, window
        
    def clean(self, minperiod, maxperiod):
        from clean import clean  # @UnresolvedImport
        import os
        filename = os.path.join(config.lightcurvespath,'clean',self.starid+'.clean')
        try:
            p, cf = np.genfromtxt(filename, unpack=True)
        except IOError:
            t = self.hjd
            x = self.mag
            f, cleaned, _ = clean(t, x, threshold=1e-3)
            n2 = len(f) /2
            cf = cleaned[n2+1:]/(2.0*np.var(x))
            p = 1./f[n2+1:]
            cf = cf[(p>=minperiod) & (p<maxperiod)]
            p = p[(p>=minperiod) & (p<maxperiod)]
            #i = np.argmax(cf)
            #period = p[i]
            a = np.column_stack((p, cf))
            np.savetxt(filename, a)
        except ValueError:
            t = self.hjd
            x = self.mag
            f, cleaned, _ = clean(t, x, threshold=1e-3)
            n2 = len(f) /2
            cf = cleaned[n2+1:]/(2.0*np.var(x))
            p = 1./f[n2+1:]
            cf = cf[(p>=minperiod) & (p<maxperiod)]
            p = p[(p>=minperiod) & (p<maxperiod)]
            #i = np.argmax(cf)
            #period = p[i]
            a = np.column_stack((p, cf))
            np.savetxt(filename, a)
        return p, cf
            
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
        # http://www.mpia-hd.mpg.de/homes/ianc/python/_modules/tools.html
        # def errxy()
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

class Analysis(object):
    '''
    Analysis class for NGC 6633
    '''


    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource
        self.wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost)
        self._gyroperiods()
        
    def _gyroperiods(self):
        from gyroage import gyroperiod
        
        self.bv_array = np.arange(0.4, 1.6, 0.01)
        
        self.P11_array = gyroperiod(self.bv_array, 600)
        self.P34_array = gyroperiod(self.bv_array, 600, P0=3.4)
        self.P01_array = gyroperiod(self.bv_array, 600, P0=0.1)
        
    def _p01(self,bv):
        i = np.abs(self.bv_array-bv).argmin()
        return self.P01_array[i]
        
    def _p11(self,bv):
        i = np.abs(self.bv_array-bv).argmin()
        return self.P11_array[i]
        
    def _p34(self,bv):
        i = np.abs(self.bv_array-bv).argmin()
        return self.P34_array[i]
        
    def load_candidates(self, candidatesfile):
        self.candidates = []
        
        self.candidates = np.genfromtxt(candidatesfile,
                                        dtype=('S25', np.float16, np.float16, np.float32, np.float32),
                                         names=True, comments=';')
    
    def save_candidates(self, candidatesfile):
        np.savetxt(candidatesfile, 
                   self.good_candidates, 
                   fmt = '%s %.3f %.3f %.5f %.5f',
                   header='starid bv vmag ra dec', 
                   comments=';')
        
    def make_lightcurves(self):
        from estimator import Estimator
        import os
        
        e = Estimator(0,len(self.candidates))
        good_candidates = []
        for k,candidate in enumerate(self.candidates):
            starid = candidate['starid']
            # if file exists, we dont need to make he effort
            if os.path.isfile(config.lightcurvespath+'%s.dat' % starid):
                good_candidates.append(tuple(candidate))
                continue
            
            # apply mid-exposure correction
            query = "SELECT hjd+(expt/43200.), mag_isocor-corr, magerr_isocor" + \
                    " FROM frames, phot" + \
                    " WHERE frames.object LIKE 'NGC 6633 rot%%'" + \
                    " AND phot.objid=frames.objid " + \
                    " AND circle(point(%(ra)f,%(dec)f),0.6/3600.0)@>circle(phot.coord,0)" %candidate + \
                    " AND filter='V' AND flags < 4" + \
                    " AND magerr_isocor < 0.01" + \
                    " AND mag_isocor+corr>0.0" + \
                    " ORDER BY frames.objid"
            result = self.wifsip.query(query)
            
            if len(result)>=50: 
                np.savetxt(config.lightcurvespath+'%s.dat' % starid, 
                           result, fmt='%.5f %.4f %.4f', header = 'hjd mag magerr') 
                #candidate['num'] = len(result)
                good_candidates.append(tuple(candidate))  
                
            e.estimate(k, comment='%-25s %d' % (starid, len(result)))
        columns = ['starid', 'bv', 'vmag', 'ra', 'dec']
        data_types = ['S25', np.float16, np.float16, np.float32, np.float32]
        self.good_candidates = np.array(good_candidates, dtype = zip(columns, data_types))
        

    def getstars(self, candidatesfile = config.datapath+'candidates_good.txt'):
        """
        build up a list of stars, where we do not have periods yet
        """
        self.load_candidates(candidatesfile)
        self.stars = []
        self.vmag = []
        self.bv = []
        for candidate in self.candidates:
            self.stars.append(candidate['starid'])
            self.vmag.append(candidate['vmag'])
            self.bv.append(candidate['bv'])
        print '... %d stars found' % len(self.candidates)
    
    def setperiod(self, starid, period):
        self.starid = starid
        params = {'starid': starid, 'period': period}
        if np.isfinite(period):
            query = """UPDATE ngc6633 SET period = %(period)f WHERE starid='%(starid)s';""" % params
        else: 
            query = """UPDATE ngc6633 SET period = NULL WHERE starid='%(starid)s';""" % params
        self.wifsip.execute(query) 

    def setamp(self, starid, amplitude):
        self.starid = starid
        params = {'starid': starid, 'amp': amplitude}
        if np.isfinite(amplitude):
            query = """UPDATE ngc6633 SET amp = %(amp)f WHERE starid='%(starid)s';""" % params
        else: 
            query = """UPDATE ngc6633 SET amp = NULL WHERE starid='%(starid)s';""" % params
        self.wifsip.execute(query) 


    def setbad(self, starid):
        self.starid = starid
        query = """UPDATE ngc6633 SET good = False WHERE starid='%s';""" % starid
        self.wifsip.execute(query) 
        
    #def __setattr__(self, name, value):
    #    params = {'starid': self.starid, 'name': name, 'value': str(value)}
    #    
    #    query = """UPDATE ngc6633 SET %(name)s = %(value)s WHERE starid='%(starid)s';""" % params
    #    self.wifsip.execute(query)
    
    #def __getattribute__(self, name):
    #    params = {'starid': self.starid, 'name': name}
    #    query = "SELECT %(name)s from ngc6633 WHERE starid='%(starid)s';" % params
    #    result = self.wifsip.query(query)
    #    return result[0]

    def plot_lightcurve(self):
        """
        plot the lightcurve for a given star
        """
        mean = np.mean(self.mag)
        plt.hlines(mean,min(self.hjd),max(self.hjd),linestyle='--')
        plt.xlim(min(self.hjd),max(self.hjd))
        plt.grid()
        plt.errorbar(self.hjd, self.mag, yerr=self.err*0.5, fmt='o')
        ylim=plt.ylim()
        plt.ylim(ylim[1],ylim[0])

    def plot_clean(self, periods, amplitudes, bv):
        plt.plot(periods, amplitudes, 'k')
        
        plt.axvspan(self._p01(bv), self._p34(bv), alpha=0.2)
        plt.axvline(self._p11(bv), ls='-.')
        i = np.argmax(amplitudes)
        period = periods[i]
        plt.axvline(x = period, color='red', alpha=0.5)
        plt.axhline(np.mean(amplitudes), color='b', ls='--')
        plt.axhline(5.*np.mean(amplitudes), color='g', ls='--')
        plt.xlim(1.0, 15.0)
        plt.minorticks_on()
    
    def plot_pdm(self, periods, thetas, bv):
        from scipy import signal
        
        kernel = signal.gaussian(101, 2)
        n = len(thetas)/2
        padded_thetas = np.lib.pad(thetas, n, mode='constant', constant_values=(1.0,1.0))
        smoothed = signal.fftconvolve(padded_thetas, kernel, mode='same')[n:-n]/5
        plt.plot(periods, smoothed, 'k')
        i = np.argmin(smoothed)
        period = periods[i]
        
        plt.axvspan(self._p01(bv), self._p34(bv), alpha=0.2)
        plt.axvline(self._p11(bv), ls='-.')
        
        
        plt.axvline(x = period, color='red', alpha=0.5)
        plt.axhline(0.8, color='b', ls='--')
        
        plt.xlim(1.0, 15.0)
        plt.ylim(0.0,1.0)
        plt.minorticks_on()
        
    def plot_phase(self, period):    
        from functions import phase
        tp, yp = phase(self.hjd, self.mag, period)

            
        s1 = np.sin(2*np.pi*tp/period)
        c1 = np.cos(2*np.pi*tp/period)
        s2 = np.sin(4*np.pi*tp/period)
        c2 = np.cos(4*np.pi*tp/period)
        
        A = np.column_stack((np.ones(tp.size), s1, c1, s2, c2))
        c, _,_,_ = np.linalg.lstsq(A,yp)

        tp1 = np.linspace(0.0, 2*period, 100)
        s1 = np.sin(2*np.pi*tp1/period)
        c1 = np.cos(2*np.pi*tp1/period)
        s2 = np.sin(4*np.pi*tp1/period)
        c2 = np.cos(4*np.pi*tp1/period)

        plt.scatter(tp, yp, edgecolor='none', alpha=0.75)
        plt.scatter(tp+period, yp, edgecolor='none', alpha=0.75)
        plt.plot(tp1,c[1]*s1+c[2]*c1+c[3]*s2+c[4]*c2, 'k', 
                  linestyle='--', linewidth=2)
        plt.axvline(period,ls='-.')
        plt.xlim(0.0, 2*period)
        plt.ylim(max(yp-np.mean(yp)),min(yp-np.mean(yp)))
        plt.xlabel('P = %.4f' % period)

    def _analyzestar(self, starid, vmag, bv, minperiod = 1.0, maxperiod = 20):
        print '%-24s '% starid,
        try:
            lc = LightCurve(starid)
            
        except IOError:
            logger.error("Can't load lightcurve %s" % starid)
            print 'no lightcurve'
            self.setbad(starid)
            return
        
        if len(lc)<50:
            logger.warn("%s: not enough datapoints" % starid)
            print 'not enough datapoints'
            return                    
        lc.rebin(0.1)
        lc.normalize()
        #lc.clip(0.1)
        lc.detrend()
        lc.sigma_clip()
        lc.detrend()
        lc.normalize()
        
        self.hjd = lc.hjd
        self.mag = lc.mag
        self.err = lc.err
        
        clean_periods, clean_amplitudes = lc.clean(minperiod, maxperiod)
        pdm_periods, pdm_thetas = lc.pdm(minperiod, maxperiod)
       
        #period = clean_periods[np.argmax(clean_amplitudes)]

        from scipy import interpolate
        
        i = np.argsort(clean_periods)
        c_periods = clean_periods[i]
        c_amps = clean_amplitudes[i]
        c_periods = np.insert(c_periods, 0 , 0.0)
        c_amps = np.insert(c_amps, 0 , 0.0)
        c_periods = np.insert(c_periods, -1 , maxperiod)
        c_amps = np.insert(c_amps, -1 , 0.0)
        c_int = interpolate.interp1d(c_periods, c_amps)
        
        # use interpolation function returned by `interp1d`
        sum_amp = c_int(pdm_periods)*(1.-pdm_thetas)  
        sum_amp /= max(sum_amp) 
        i = np.argmax(sum_amp)
        period = pdm_periods[i] 
        theta = pdm_thetas[i]
        
        import functions as fx
        
        
        ci = np.argmax(c_amps)
       
        
            
        try:
            j = np.where((pdm_periods>period-2) & (pdm_periods<period+2))
            _, _, pgf_sigma = fx.gauss_fit(pdm_periods[j], 1.-pdm_thetas[j],  1.0-theta, period, 1.0)      
        except RuntimeError:
            logger.warning('%s: unable to fit pgf gaussian'  % starid)
            pgf_sigma = 'NULL'
        cgf_sigma = 0.0    
        try:
            j = np.where((c_periods>c_periods[ci]-2) & (c_periods<c_periods[ci]+2))
            _, _, cgf_sigma = fx.gauss_fit(c_periods[j], c_amps[j],  c_amps[ci], c_periods[ci], 1.0)
        except RuntimeError:
            logger.warning('%s: RuntimeError, unable to fit cgf gaussian' % starid)
        except TypeError:
            logger.warning('%s: TypeError, unable to fit cgf gaussian' % starid)
            cgf_sigma= 'NULL'
        
        amp, rms = 0,0
        try:
            amp, rms = lc.phased(period)
            snr = np.sqrt(len(self.hjd))*amp/rms
        except IndexError:
            logger.warning('%s: unable to calculate phased function' % starid)
            amp = 'NULL'
            rms = 'NULL'
        
        sigma_limit = 3.0
        theta_limit = 0.875
        try:
            ci = np.argmax(c_amps)
            
            
            params = {'period': period,
                      'period_err': pgf_sigma,
                      'clean_period': c_periods[ci],
                      'clean_amp': c_amps[ci],
                      'clean_sigma': cgf_sigma,
                      'theta': theta,
                      'freq': period,
                      'amp': amp,
                      'rms': rms,
                      'ni': len(self.hjd)
                      }
            
            keys = ', '.join(params.keys())
            values = ', '.join([str(v) for v in params.values()])
            
            query = """UPDATE ngc6633 SET (%s) = (%s) WHERE starid='%s';""" % (keys, values, starid)
            #logger.info(query)
            self.wifsip.execute(query)
        except ValueError:
            print 'Cannot store params for starid %s' % starid
            logger.error("Cannot store params for starid %s: %s" % (starid, query))
        
        try:    
            mamp = np.mean(clean_amplitudes)
        except ValueError:
            print 'No Clean amplitudes for starid %s' % starid
            logger.error("No Clean amplitudes for starid %s" % starid)
            return
        try:
            if max(clean_amplitudes)> sigma_limit*mamp and pdm_thetas[i] < theta_limit:
                print "%.2f %.1f %.2f" % (period,   max(clean_amplitudes)/mamp, cgf_sigma)
                self.setperiod(starid, period)
            else:
                print '< %.1f sigmas or  theta %.2f > %.2f' % (sigma_limit, pdm_thetas[i], theta_limit)
                self.setperiod(starid, np.nan)
                return
        except ValueError:
            print 'No maximum clean amplitude for starid %s' % starid
            logger.error("No maximum clean amplitude for starid %s" % starid)
            return           
            
        #self.setamp(starid, amp)
            
        star = {'tab':0, 'bv':bv}
        
        plt.subplot(411) ##################################################
        plt.title('%s (%d) V = %.2f B-V=%.2f S/N = %.1f' % (starid, star['tab'], vmag, bv,snr))
        self.plot_lightcurve()
        
        plt.subplot(412) ##################################################
        self.plot_clean(clean_periods, clean_amplitudes, bv)
        

        plt.subplot(413) ##################################################
        plt.plot(pdm_periods, sum_amp,'g')
        plt.axvline(period, color='b')
        self.plot_pdm(pdm_periods, pdm_thetas, bv)
        
        plt.subplot(414) ##################################################
        self.plot_phase(period)
        
        plt.savefig(config.plotpath+'%s(%d).pdf' % (starid,star['tab']))
        plt.close()

    def analysis(self):
        """
        perform various pariod and frequency analysis on each star with a 
        lightcurve
        """
        from matplotlib import rcParams
        
        print 'Analysis'

        fig_width = 18.3/2.54  # width in inches, was 7.48in
        fig_height = 23.3/2.54  # height in inches, was 25.5
        fig_size =  [fig_width,fig_height]
        #set plot attributes
        params = {'backend': 'Agg',
          'axes.labelsize': 12,
          'axes.titlesize': 12,
          'font.size': 12,
          'xtick.labelsize': 12,
          'ytick.labelsize': 12,
          'figure.figsize': fig_size,
          'savefig.dpi' : 300,
          'font.family': 'sans-serif',
          'axes.linewidth' : 0.5,
          'xtick.major.size' : 2,
          'ytick.major.size' : 2,
          }
        rcParams.update(params)
        
        for starid,vmag,bv in zip(self.stars, self.vmag, self.bv):
            self._analyzestar(starid, vmag, bv, minperiod = 1.3, maxperiod = 15)

if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description='NGC 6633 analysis')
    parser.add_argument('-a', '--analysis', action='store_true', help='perform analysis')
    parser.add_argument('--make', action='store_true', help='make lightcurves')
#    parser.add_argument('filter', default='V', type=str, help='filter color to process')

    args = parser.parse_args()

    if args.make:
        analysis = Analysis()
        analysis.load_candidates(config.datapath+'candidates.txt')
        analysis.make_lightcurves()
        analysis.save_candidates(config.datapath+'candidates_good.txt')
    
    if args.analysis: 
        analysis = Analysis()
        analysis.getstars()
        analysis.analysis()

        


    