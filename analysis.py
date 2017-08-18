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
from lightcurve import LightCurve

logging.basicConfig(filename=config.projectpath+'analysis.log', 
                    format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('NGC6633 analysis')
logging.captureWarnings(True)

class Analysis(object):
    '''
    Analysis class for NGC 6633
    '''
    def __init__(self, minperiod = 1.25, maxperiod = 20):
        '''
        Constructor
        '''
        from datasource import DataSource
        
        self.minperiod = minperiod
        self.maxperiod = maxperiod
        self.sigma_limit = 3.5
        self.theta_limit = 0.97
        
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
            
        # restart from 20140721A-0002-0014#1504 ?
            
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
        
    def isbad(self, starid):
        query = """SELECT good FROM ngc6633 WHERE starid='%s';""" % starid
        result = self.wifsip.query(query)
        if result[0][0] is None or result[0][0] is True:
            return False
        else:
            return True
        
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
        std = np.std(self.mag)
        t = self.hjd - min(self.hjd)
        plt.axhline(mean,linestyle='--')
        plt.axhline(mean+2.0*std,linestyle=':')
        plt.axhline(mean-2.0*std,linestyle=':')
        
        plt.xlim(min(t),max(t))
        
        plt.ylabel('mag')
        
        plt.minorticks_on()
        plt.errorbar(t, self.mag, yerr=self.err, fmt='.k',
                     capsize=None, capthick=0.0, elinewidth=0.0)
        ylim=plt.ylim()
        plt.ylim(ylim[1],ylim[0])

    def plot_clean(self, periods, amplitudes, residuals, bv, clean_period, clean_error):
        mean_amps = np.mean(amplitudes)
        # we plot normalized by sigmas
        plt.plot(periods, amplitudes/mean_amps, 'k')
        plt.plot(periods, residuals/mean_amps, 'b')
        
        plt.axvspan(self._p01(bv), self._p34(bv), alpha=0.2)
        plt.axvline(self._p11(bv), ls='-.')
        
        plt.axvline(x = clean_period, color='red', alpha=0.5)
        plt.axvline(x = clean_period+clean_error, color='red', alpha=0.5, ls=':')
        plt.axvline(x = clean_period-clean_error, color='red', alpha=0.5, ls=':')
        
        plt.axhline(self.sigma_limit, color='k', ls='--')
        plt.xlim(self.minperiod, self.maxperiod)
        plt.ylabel('clean $\sigma$')
        plt.minorticks_on()
        
    def plot_lombscargle(self, periods, powers, window, bv, ls_period, ls_error):
        mean_pwrs = np.mean(powers)
        plt.plot(periods, powers/mean_pwrs, 'k')
        plt.plot(periods, window/np.mean(window), 'b')

        plt.axvspan(self._p01(bv), self._p34(bv), alpha=0.2)
        plt.axvline(self._p11(bv), ls='-.')

        plt.axvline(x = ls_period, color='blue', alpha=0.5)
        plt.axvline(x = ls_period+ls_error, color='blue', alpha=0.5, ls=':')
        plt.axvline(x = ls_period-ls_error, color='blue', alpha=0.5, ls=':')

        plt.axhline(self.sigma_limit, color='k', ls='--')
        plt.xlim(self.minperiod, self.maxperiod)
        plt.ylabel('lomb $\sigma$')
        plt.minorticks_on()

    
    def plot_pdm(self, periods, thetas, bv, pdm_period, pdm_error):
        plt.plot(periods, thetas, 'k')
        
        plt.axvspan(self._p01(bv), self._p34(bv), alpha=0.2)
        plt.axvline(self._p11(bv), ls='-.')
        
        plt.axhline(self.theta_limit, color='k', ls='--')
        
        plt.axvline(x = pdm_period, color='green', alpha=0.5)
        
        plt.xlim(self.minperiod, self.maxperiod)
        plt.ylabel('$\Theta$')
        
        plt.ylim(min(thetas),1.0)
        plt.minorticks_on()
        
    def plot_phase(self, period):    
        from functions import phase
        tp, yp = phase(self.hjd-np.min(self.hjd), self.mag-np.mean(self.mag), period)

            
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
        plt.ylabel('mag')

    def _analyzestar(self, starid, vmag, bv):
        logger.info('%-24s: Vmag=%.2f B-V=%.2f' % (starid, vmag, bv))
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
        
        if self.isbad(starid):
            logger.info("%s: already marked bad" % starid)
            print 'already marked bad'
            return  
                          
        #lc.rebin(0.1)
        #lc.normalize()
        #lc.clip(0.1)
        lc.detrend()
        lc.sigma_clip(2.0)
        # since sigma clipping might change the trend of the lightcurve, we are detrending again 
        lc.detrend()
        #lc.normalize()
        
        self.hjd = lc.hjd
        self.mag = lc.mag
        self.err = lc.err
        
        clean_periods, clean_amplitudes, clean_residuals = lc.clean(self.minperiod, self.maxperiod, gain=0.1)
        ls_periods, ls_powers, ls_window = lc.lomb_scargle(minfrequency=1./self.maxperiod, maxfrequency=1./self.minperiod)
        pdm_periods, pdm_thetas = lc.pdm(self.minperiod, self.maxperiod)
       
        def issimilar(value1, error1, value2, error2):
            if error1<0.0 or error2<0.0:
                raise(ValueError,"error values must be positive")
            if np.abs(value1 - value2) < error1 + error2:
                return True
            else:
                return False
            
            
        period = lc.pdm_period
        period_err = lc.pdm_error
            
        
        amp, rms = 0,0
        # try to determine the rms of the residual lightcurve
        try:
            amp, rms = lc.phased(period)
            snr = np.sqrt(len(self.hjd))*amp/rms
        except IndexError:
            logger.warning('%s: unable to calculate phased function' % starid)
            amp = 'NULL'
            rms = 'NULL'
        
        
        
        # try to store the results
        try:
            params = {'pdm_period':     lc.pdm_period,
                      'pdm_error':      lc.pdm_error,
                      'pdm_theta':      lc.pdm_theta,
                      'clean_period':   lc.clean_period,
                      'clean_amp':      lc.clean_amplitude,
                      'clean_sigma':    lc.clean_error,
                      'ls_period':      lc.ls_period, 
                      'ls_error':       lc.ls_error, 
                      'ls_power':       lc.ls_power, 
                      'period':         period,
                      'period_err':     period_err,
                      'amp':            amp,
                      'rms':            rms,
                      'ni':             len(self.hjd)
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
            mean_amp = np.mean(clean_amplitudes)
        except ValueError:
            print 'No Clean amplitudes for starid %s' % starid
            logger.error("No Clean amplitudes for starid %s" % starid)
            return
        
        clean_sigmas = lc.clean_amplitude/mean_amp
        
        if (issimilar(lc.clean_period, lc.clean_error, lc.pdm_period, lc.pdm_error) or \
             issimilar(lc.ls_period, lc.ls_error, lc.pdm_period, lc.pdm_error)):
            logstring = "clean: %.2f+-%.2f\tls: %.2f+-%.2f\tpdm: %.3f+-%.3f" % (lc.clean_period, 
                                                                 lc.clean_error, 
                                                                 lc.ls_period,
                                                                 lc.ls_error,
                                                                 lc.pdm_period,
                                                                 lc.pdm_error
                                                                 )
            print logstring
            logger.info(logstring)
            self.setperiod(starid, period)
        else:
            logstring = '%.1f < %.1f sigmas or theta %.2f > %.2f' % (clean_sigmas, self.sigma_limit, lc.pdm_theta, self.theta_limit)
            print logstring
            logger.info(logstring)
            
            self.setperiod(starid, np.nan)
            return
            
        #make the plots    
        star = {'tab':0, 'bv':bv}
        
        plt.subplot(511) ##################################################
        plt.title('%s (%d) V = %.2f B-V=%.2f S/N = %.1f' % (starid, star['tab'], vmag, bv,snr))
        self.plot_lightcurve()
        
        plt.subplot(512) ##################################################
        self.plot_clean(clean_periods, clean_amplitudes, clean_residuals, bv, 
                        lc.clean_period, lc.clean_error)
        
        plt.subplot(513) ##################################################
        self.plot_lombscargle(ls_periods, ls_powers, ls_window, bv,
                              lc.ls_period, lc.ls_error)

        plt.subplot(514) ##################################################
        plt.axvline(period, color='b')
        self.plot_pdm(pdm_periods, pdm_thetas, bv,
                      lc.pdm_period,lc.pdm_error)
        
        plt.subplot(515) ##################################################
        self.plot_phase(period)
        
        plt.savefig(config.plotpath+'%s(%d).pdf' % (starid,star['tab']))
        plt.close()

    def analysis(self):
        """
        perform various period and frequency analysis on each star with a 
        lightcurve
        """
        print 'Analysis'

        
        #set plot attributes
        plt.style.use('a4paper.mplstyle')
        
        for starid,vmag,bv in zip(self.stars, self.vmag, self.bv):
            self._analyzestar(starid, vmag, bv)

if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description='NGC 6633 analysis')
    parser.add_argument('-a', '--analysis', action='store_true', help='perform analysis')
    parser.add_argument('--make', action='store_true', help='make lightcurves')
#    parser.add_argument('filter', default='V', type=str, help='filter color to process')

    args = parser.parse_args()

    if args.make:
        analysis = Analysis(minperiod = 1.25, maxperiod = 20)
        analysis.load_candidates(config.datapath+'candidates.txt')
        analysis.make_lightcurves()
        analysis.save_candidates(config.datapath+'candidates_good.txt')
    
    if args.analysis: 
        analysis = Analysis()
        analysis.getstars()
        analysis.analysis()

        


    