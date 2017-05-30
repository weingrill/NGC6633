#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 4, 2017

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import numpy as np
import config
from datasource import DataSource
import pickle

from astropy.coordinates import SkyCoord,search_around_sky   # @UnresolvedImport
from astropy import units as u

from matplotlib import pyplot as plt

class CorrFrames(object):
    '''
    Applies the correction of frames based on the list of reference stars taken 
    from APASS9 and downloaded accorindgly into the table. See apass.py for DB 
    injection.
    '''


    def __init__(self, filtercol):
        if filtercol not in ['B','V']:
            raise(ValueError,'Wrong filter color')
        self.filtercol= filtercol
        self.wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost)

    @property
    def nobjids(self):
        return len(self.objids)
    
    @property
    def nrefstars(self):
        return len(self.refstarids)
    
    def getobjids(self):
        
        query = "SELECT objid FROM frames WHERE object LIKE 'NGC 6633 rot %%' AND filter='%s' ORDER BY objid" % self.filtercol
        result = self.wifsip.query(query)
        self.objids = [r[0] for r in result]  
        self.safepickle('corr_objids%s.pickle' % self.filtercol, self.objids)  
        print 'objids: ', self.nobjids
        
    def getreferencestars(self):
        query = 'SELECT starid, ra, dec, "B", "V" ' + \
                ' FROM referencestars ' + \
                ' WHERE referencestars.coord <@ box(point(276.44583,6.1416),point(277.17916,6.8749)) ' + \
                ' AND "V">0 and "B">0 ' + \
                ' ORDER BY starid;'
        stars = self.wifsip.query(query)
        
        self.refstarids = [r[0] for r in stars]      
        self.refra = [r[1] for r in stars]      
        self.refdec = [r[2] for r in stars] 
        if self.filtercol =='B':     
            self.refmag = np.array([r[3] for r in stars])      
        if self.filtercol =='V':     
            self.refmag = np.array([r[4] for r in stars])      
        self.safepickle('corr_starids%s.pickle' % self.filtercol, self.refstarids)  
        print 'reference stars: ', self.nrefstars

    def safepickle(self, filename, data):
        with open(config.datapath+filename, 'wb') as picklefile:
            pickle.dump(data, picklefile)
    
    def getstars(self, objid):
        query = """SELECT alphawin_j2000, deltawin_j2000, mag_isocor 
                    FROM phot
                    WHERE objid = '%s' AND flags=0
                    ORDER BY star;""" % objid
        result = self.wifsip.query(query)
        self.ra = np.array([r[0] for r in result])
        self.dec = np.array([r[1] for r in result])
        self.mag = np.array([r[2] for r in result])
        
    def setcorr(self, objid, corr):
        query = """UPDATE frames SET corr = %f WHERE objid = '%s';""" % (corr, objid)
        # this is a rude hack
        query = query.replace('nan', 'NULL')
        self.wifsip.execute(query)
        
    def matchframe(self, objid):
        # get imaging data
        #image_data = fetch_imaging_sample()
        c = SkyCoord(ra=self.ra*u.degree, dec=self.dec*u.degree)  
        catalog = SkyCoord(ra=self.refra*u.degree, dec=self.refdec*u.degree)  
        #idx, d2d, d3d = c.match_to_catalog_sky(catalog)  
        #matches = catalog[idx]
        #print len(matches)
        
        #idxc, idxcatalog, d2d, d3d = catalog.search_around_sky(c, 1.5*u.arcsec)
        idxc, idxcatalog, d2d, d3d = search_around_sky(c, catalog, 0.6*u.arcsec)
        if len(idxc) > 9:
            print len(idxc),
            oc = self.mag[idxc]-self.refmag[idxcatalog]
            meanoc = np.mean(oc)
            stdoc = np.std(oc)
            #plt.scatter(self.refmag[idxcatalog], oc)
            #plt.axhline(np.mean(oc), linestyle='-')
            #plt.axhline(np.median(oc), linestyle='-.')
            
            #plt.axhline(np.mean(oc)-np.std(oc), linestyle='--')
            #plt.axhline(np.mean(oc)+np.std(oc), linestyle='--')
            
            #perform sigma clipping
            oc = oc[abs(oc-meanoc) < 2*stdoc]
            
            
            meanoc = np.mean(oc)
            stdoc = np.std(oc)
            oc = oc[oc < meanoc+stdoc]
            #plt.axhline(np.mean(oc), linestyle='-', color='r')
            #plt.axhline(np.median(oc), linestyle='-.', color='r')
            
            #plt.axhline(np.mean(oc)-np.std(oc), linestyle='--', color='r')
            #plt.axhline(np.mean(oc)+np.std(oc), linestyle='--', color='r')
            
            #plt.grid()
            #plt.show()
    
            #plt.close()
            print '%.3f' % stdoc,
            print '%.3f' % np.mean(oc),
            print '%.3f' % np.std(oc)
            corr = np.mean(oc)
        else:
            print '!'
            corr = np.nan
        self.setcorr(objid, corr)
        
    
    def corrframe(self, objid):
        ''' match frame with reference '''
        self.getstars(objid)
        self.matchframe(objid)
    
    def corrframes(self):
        for objid in self.objids:
            print objid,
            self.corrframe(objid)
    
if __name__ == '__main__':
    cf = CorrFrames('V')
    cf.getobjids()
    cf.getreferencestars()
    cf.corrframes()