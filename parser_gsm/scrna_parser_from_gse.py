"""
author: Rongbin Zheng
resource: derive from Cistrome DB geo parser framework
purpose: parse single cell FNA-seq data from GEO database
"""
import os,sys
import json, re, time
import urllib.request, urllib.parse, urllib.error
import traceback
from datetime import datetime
import pickle
import subprocess
from operator import itemgetter
import random
import importlib
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

from xml.dom.minidom import parseString

#AUTO-load classifiers
#a trick to get the current module
_modname = globals()['__name__']
_this_mod = sys.modules[_modname]

_ppath = "/".join(_this_mod.__file__.split("/")[:-1])

# from django.utils.encoding import smart_str

# import scrna_parser_sampleDetail
# import pubmed
# ## build mysql environment
# import env
# models = env.models


""" ========== main script ========== """

def getSyncLog(infoStr):
    """ouput the record to DoneGsmXml.log file
    """
    os.system('echo "[%s] %s"' % (time.strftime('%H:%M:%S'), infoStr))


### GDS interface
def getGDSSamples(date_region=False):
    """Will run the predefined query and return a list of GDS ids
    NOTE: this returns ALL GDS samples which are of SRA type i.e.
    ALL CHIP-SEQ, RNA-SEQ, etc.
    """
    #expireDate = now - 30 days in seconds
    #ref: http://stackoverflow.com/questions/7430928/python-comparing-date-check-for-old-file
    # _expireDate = time.time() - 60 * 60 * 24 * 30

    ret = []
    #
    # #TRY: to read a file first -- IF IT IS NOT STALE
    path = os.path.join(_ppath, "gdsSamples.txt")
    
    #REAL URL
    URL = """http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=SRA[Sample%20Type]%20AND%20gse[Entry%20Type]%20AND%20(homo%20sapiens[Organism]%20OR%20mus%20musculus[Organism])&retmax=100000&usehistory=y"""
    if date_region:
        maxTime = date_region.split('-')[1]
        minTime = date_region.split('-')[0]
        print(date_region)
        URL = """http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=SRA[Sample%20Type]%20AND%20gse[Entry%20Type]%20AND%20(homo%20sapiens[Organism]%20OR%20mus%20musculus[Organism])&mindate=%s&maxdate=%s&datetype=pdat&retmax=100000&usehistory=y"""%(minTime, maxTime)
    try:
        getSyncLog("getGDSSample: %s" % URL) # output record
        f = urllib.request.urlopen(URL)
        root = ET.fromstring(f.read())
        f.close()

        #Get the IDList
        tmp = root.findall("IdList/Id")
        ret = [i.text for i in tmp]

        #write to disk
        f = open(path, "w")
        for l in ret:
            f.write("%s\n" % l)
        f.close()
        print("Refresh %s"%path)
    except:
        print("Exception in user code:")
        print('-' * 60)
        traceback.print_exc(file=sys.stdout)
        print('-' * 60)
    return ret
def gse_idToAcc(gdsId):
    """Given a GDS id, e.g. 300982523, tries to give a GDS accession, e.g.
    GSM982523

    NOTE: there is an algorithm: acc = "GSM"+gdsId[1:] (strip leading 0s)
    """
    #Cut = dropping of the "3" (which indicates sample) and removal of leading
    #leading 0s
    cut = gdsId[1:].lstrip("0")
    return "GSE%s" % cut

def proxyInstead(link, using=False):
    """using proxy to aviod forbidden
    """
    context = ''
    if using:
        #using proxy first
        try: # using proxy first, or using the read ip
            agent = [x.rstrip() for x in open('./pickle_file/proxy.txt')]
            proxy = {'http':'http://%s'%random.sample(agent, 1)[0]}
            urlf = urllib.request.urlopen(link, proxies = proxy)
            getSyncLog('.')
        except:
            urlf = urllib.request.urlopen(link)
            proxy = {'proxy':'local IP'}
            getSyncLog('.') # use for record, so that we can know what happened if error occured
        # check whether we get the correct inf
        context = urlf.read()
        urlf.close()
        if ('404 - File or directory not found' in context) or ('ERR_ACCESS_DENIED' in context) or (context.strip() == ''):
            urlf = urllib.request.urlopen(link)
            context = urlf.read()
            urlf.close()
            proxy = {'proxy':'local IP'}
            getSyncLog('.')
        getSyncLog('%s: %s'%(list(proxy.values())[0], link))
        context = context.decode(encoding='utf-8',errors='ignore')
        return context
    try:
        # time.sleep(0.3)
        urlf = urllib.request.urlopen(link)
        context = urlf.read()
        context = context.decode(encoding='utf-8',errors='ignore')
        urlf.close()
        getSyncLog('local IP: %s'%link)
        return context
    except:
        print('link problem: %s'%link)
    return None

def isXML(doc):
    """TEST if it is a valid geo XML record
    NOTE: first lines are-
    <?xml version="1.0" encoding="UTF-8" standalone="no"?>
    """
    f = doc.split("\n")
    return f[0].strip() == """<?xml version="1.0" encoding="UTF-8" standalone="no"?>"""

  
def getGeoXML(accession, path='geo_gse'):
    """HANDLES GEO XML records--i.e. our GEO XML librarian!
    Given a GEO ACCESSION ID, return the xml record for it
    (making the urllib call)"""
    #path pattern: EXAMPLE-GSE1126513 geo/GSE1126/GSE1126513
    #path = os.path.join(_ppath, ddir)
    if not os.path.exists(path):
        os.mkdir(path)
    subdir = os.path.join(path, accession[:7])
    if not os.path.exists(subdir):
        os.mkdir(subdir)
    path = os.path.join(subdir, "%s.xml" % accession)
    if os.path.exists(path):
        f = open(path)
        docString = f.read()
        f.close()
    else:
        #print accession
        URL = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s&view=quick&form=xml&targ=self" % accession
        try:
            #print "getGeoXML: %s" % URL
            #signal.alarm(180)
            docString = proxyInstead(link=URL)
            if not isXML(docString): # try again
                #signal.alarm(180)
                getSyncLog('.')
                docString = proxyInstead(link=URL)
            if isXML(docString):
                #write to file
                f = open(path, "w")
                f.write(docString)
                f.close()
                #getSyncLog(proxy.values()[0]+'\t'+accession + '\n')# output record
            else:
                print(accession)
                print("ERROR: accession is NOT xml. (The accession may be deleted from GEO repository)")
                f1 = open('gsm_notXML.txt', 'a')
                f1.write(accession + '\n')
                f1.close()
                return None
        except:
            print("Exception in user code:")
            print('-' * 60)
            traceback.print_exc(file=sys.stdout)
            print('-' * 60)
            docString = None
    return docString

def _sync_gse():
    ## get all GDS ids of GSE from API
    gdsSamples = getGDSSamples()
    getSyncLog('start: There are %s GDS Samples in sum'%(len(gdsSamples)))#
    local = [x.strip() for x in open('geo_gse_collection.txt', 'r')]
    out = open('gse_gdsSample_all.txt', 'w')
    out.write('\n'.join(gdsSamples))
    out.close()
    # convert to gds id to GSE and download XML file
    cnt = 0
    one_percent = len(gdsSamples)/100
    for gds in gdsSamples:
        cnt += 1
        if cnt % one_percent == 0:
            getSyncLog("%s%%"%(cnt/one_percent))
            time.sleep(3)
        gseid = gse_idToAcc(gds)
        if gseid in local:
            continue
        gseXML = getGeoXML(gseid)
        out = open('geo_gse_collection.txt', 'a')
        out.write(gseid+'\n')
        out.close() 
        time.sleep(0.03) # sleep to avoid IP blocking

    getSyncLog('done!')#

if __name__ == '__main__':
    _sync_gse()


