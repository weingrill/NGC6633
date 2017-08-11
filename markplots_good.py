'''
Created on Oct 25, 2013

@author: jwe <jweingrill@aip.de>
'''
import config
def filename2starid(filename):
    import os
    basename = os.path.basename(filename)
    return os.path.splitext(basename)[0].rstrip('(0)')

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

def update_good():
    from glob import glob
    import os
    
    basepdfs = [os.path.basename(pathname) for pathname in glob(config.plotpath+'/2014*.pdf')]
    goodpdfs = [os.path.basename(pathname) for pathname in glob(config.plotpath+'good/2014*.pdf')]
    #print goodpdfs
    for base in basepdfs:
        
        if base in goodpdfs:
            print base
            os.rename(config.plotpath+base, config.plotpath+'good/'+base)

def update_bad():
    from glob import glob
    import os
    
    basepdfs = [os.path.basename(pathname) for pathname in glob(config.plotpath+'/2014*.pdf')]
    goodpdfs = [os.path.basename(pathname) for pathname in glob(config.plotpath+'bad/2014*.pdf')]
    #print goodpdfs
    for base in basepdfs:
        
        if base in goodpdfs:
            print base
            os.rename(config.plotpath+base, config.plotpath+'bad/'+base)

if __name__ == '__main__':
    from glob import glob
    from datasource import DataSource
    
    update_good()
    update_bad()

    
    wifsip = DataSource(database=config.dbname, user=config.dbuser, host=config.dbhost)
    goodpdfs = glob(config.plotpath+'good/*.pdf')
    for p in goodpdfs:
        #print p
        objid = filename2starid(p)
        print objid,'= good'
        set_good(objid)

    badpdfs = glob(config.plotpath+'bad/*.pdf')
    for p in badpdfs:
        #print p
        objid = filename2starid(p)
        print objid,'= bad'
        set_bad(objid)


    wifsip.commit()
    wifsip.close()
    print '%d stars marked as good' % len(goodpdfs)
    print '%d stars marked as bad' % len(badpdfs)