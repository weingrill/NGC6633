#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Oct 25, 2013

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import config
def filename2starid(filename):
    import os
    basename = os.path.basename(filename)
    return os.path.splitext(basename)[0].rstrip('(0)')

def set_null():
    query = """UPDATE ngc6633
               SET good=NULL;"""
    wifsip.execute(query)

def set_good(starid):
    query = """UPDATE ngc6633
               SET good=True
               WHERE starid like '%s';""" % starid
    wifsip.execute(query)

def set_bad(starid):
    query = """UPDATE ngc6633
               SET good=FALSE
               WHERE starid like '%s';""" % starid
    wifsip.execute(query)

def delete():
    """
    move already marked good files to respective folder
    """
    from glob import glob
    import os
    
    def _deletefile(filename):
        try:
            os.remove(filename)
        except OSError:
            print "Can't delete %s" % filename
            
    
    deletepdfs = [os.path.basename(pathname) for pathname in glob(config.plotpath+'delete/2014*.pdf')]
    #print goodpdfs
    for filename in deletepdfs:
        starid = filename.replace('(0).pdf', '')
        print filename, starid
        _deletefile(config.lightcurvespath+'pdm/%s.pdm' % starid)
        _deletefile(config.lightcurvespath+'ls/%s.ls' % starid)
        _deletefile(config.lightcurvespath+'clean/%s.clean' % starid)
        _deletefile(config.lightcurvespath+'%s.dat' % starid)
        _deletefile(config.plotpath+'delete/'+filename)


def _update(targetfolder, basefolder='', comparefolder=''):
    from glob import glob
    import os
    
    if comparefolder=='':
        comparefolder=targetfolder
    if not targetfolder in ['', 'good', 'bad', 'extra']:
        raise ValueError
    
    basepdfs = [os.path.basename(pathname) for pathname in glob(config.plotpath+'%s/2014*.pdf' % basefolder)]
    selectpdfs = [os.path.basename(pathname) for pathname in glob(config.plotpath+'%s/2014*.pdf' % comparefolder)]
    
    for basefile in basepdfs:
        
        if basefile in selectpdfs:
            print basefile
            os.rename(config.plotpath+basefile, config.plotpath+'%s/%s' % (targetfolder,basefile))

def update_good():
    """
    move already marked good files to respective folder
    """
    _update('good')

def update_bad():
    """
    move already marked bad files to respective folder
    """
    _update('bad')

def update_extra():
    """
    move already marked extra files to respective folder
    """
        
    _update('extra')

if __name__ == '__main__':
    from glob import glob
    from datasource import DataSource
    
    #delete()
    update_good()
    update_bad()
    update_extra()
    
    wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost)
    
    #set_null()
    
    goodpdfs = glob(config.plotpath+'good/*.pdf')
    goodids = [filename2starid(p) for p in goodpdfs]
    
    goodidstring = "('"+"', '".join(goodids)+"')"
    badpdfs = glob(config.plotpath+'bad/*.pdf')
    badids = [filename2starid(p) for p in badpdfs]
    badidstring = "('"+"', '".join(badids)+"')"
    
    params = {'good': goodidstring, 'bad': badidstring}
    query = "UPDATE ngc6633 SET good = NULL; " +\
            "UPDATE ngc6633 SET good = TRUE WHERE starid in %(good)s; " % params + \
            "UPDATE ngc6633 SET good = FALSE WHERE starid in %(bad)s;" % params
    wifsip.execute(query)
    

    wifsip.commit()
    wifsip.close()
    print '%d stars marked as good' % len(goodpdfs)
    print '%d stars marked as bad' % len(badpdfs)