#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 13, 2017

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import numpy as np
import matplotlib.pyplot as plt
from dbtable import DBTable
import config

N = 200

class IsoCrossCorr(object):
    '''
    This class makes a crosscorrelation of synthetical isochrones with the 
    actual color magnitude diagram.
    '''


    def __init__(self):
        '''
        Constructor
        '''
        from datasource import DataSource
        self.wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost)
    
    def loadisochrone(self):
        isofile = config.datapath+'yapsi_l_X0p703812_Z0p016188_600Myr.dat'
        Mvs, bvs = np.genfromtxt(isofile, usecols=[5,7], unpack=True)
        self.iso = np.zeros((N, N))
        for mvi, bvi in zip(Mvs, bvs):
            vmag = int((mvi+4.0)*(N/20.0))
            bv = int((bvi+0.2)*(N/2.0))
            if bv in np.arange(1, N-1) and vmag in np.arange(1, N-1):
                for dx in [-1, 0, 1]:
                    for dy in [-1, 0, 1]:
                        self.iso[vmag + dx, bv + dy] = 1.0
        
        plt.imshow(self.iso, extent=[0.0,2.0,20.0,0.0], interpolation='none', aspect='auto')   
        plt.show()
    
    def loaddata(self):
        
        alls = DBTable(self.wifsip, 'ngc6633', condition='NOT bv is NULL AND vmag<20 and (ni>0 or member or hiltner_member)')#
        self.stars = np.zeros((N, N))
        for star in alls:
            vmag = int((star['vmag']+4.0)*(N/20.0))
            bv = int((star['bv']+0.2)*(N/2.0))
            if bv<N and vmag<N:
                self.stars[vmag, bv] += 1.0
        
        plt.imshow(np.log(self.stars+1.0), extent=[0.0,2.0,20.0,0.0], interpolation='none', aspect='auto')   
        plt.show()
        
    def crosscorrelate(self):
        from scipy import signal
        #self.stars = np.log(self.stars+1.0)
        face  =  self.stars - self.stars.mean()
        template = self.iso - self.iso.mean()
        #template = np.fliplr(np.flipud(template))
        corr = signal.correlate2d(face, template, boundary='fill', mode='same')
        #corr = np.rot90(corr,k=2)
        corr = np.rot90(np.transpose(corr),k=2)
        y, x = np.unravel_index(np.argmax(corr), corr.shape)  # find the match
        bv, vmag = x/(N/2.0), y/(N/20.0)
        print bv, vmag
        plt.imshow(corr,extent=[0.0,2.0,20.0,0.0], interpolation='none', aspect='auto')
        plt.plot(bv, vmag, 'ro')
        plt.plot(0.165,7.71, 'go')
        plt.show()
        
if __name__ == '__main__':
    icc = IsoCrossCorr()
    icc.loaddata()
    icc.loadisochrone()
    icc.crosscorrelate()