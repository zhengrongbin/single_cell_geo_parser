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

from django.utils.encoding import smart_str

import scrna_parser_sampleDetail
import pubmed
## build mysql environment
import env
models = env.models

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
    URL = """http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=SRA[Sample%20Type]%20AND%20gsm[Entry%20Type]%20AND%20(homo%20sapiens[Organism]%20OR%20mus%20musculus[Organism])&retmax=100000&usehistory=y"""
    if date_region:
        maxTime = date_region.split('-')[1]
        minTime = date_region.split('-')[0]
        print(date_region)
        URL = """http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=SRA[Sample%20Type]%20AND%20gsm[Entry%20Type]%20AND%20(homo%20sapiens[Organism]%20OR%20mus%20musculus[Organism])&mindate=%s&maxdate=%s&datetype=pdat&retmax=100000&usehistory=y"""%(minTime, maxTime)
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

def gsm_idToAcc(gdsId):
    """Given a GDS id, e.g. 300982523, tries to give a GDS accession, e.g.
    GSM982523

    NOTE: there is an algorithm: acc = "GSM"+gdsId[1:] (strip leading 0s)
    """
    #Cut = dropping of the "3" (which indicates sample) and removal of leading
    #leading 0s
    cut = gdsId[1:].lstrip("0")
    return "GSM%s" % cut
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

  
def getGeoXML(accession, path='geo'):
    """HANDLES GEO XML records--i.e. our GEO XML librarian!
    Given a GEO ACCESSION ID, return the xml record for it
    (making the urllib call)"""
    #path pattern: EXAMPLE-GSM1126513 geo/GSM1126/GSM1126513
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

def readGeoXML(path, docString=None):
    """
    Input: a file path or a string--default is to use the path

    Tries to read in the geo xml record,
    **KEY: REMOVES the xmlns line
    Returns the xml record text WITHOUT the xmlns line!!!
    """
    if docString:
        f = docString.split("\n")
    else:
        if path:
            f = open(path)
        else:
            f = ''
    tmp = []
    try:
        for l in f:
            if l.find("xmlns=\"http://www.ncbi.nlm.nih.gov/geo/info/MINiML\"") == -1:
                tmp.append(l)
    except: # sometimes, the xml is not saved as utf-8 in local, do getGeoXml agaim
        gsmid = os.path.basename(path).rstrip('.xml')
        ag = getGeoXML(accession=gsmid)
        f = open(os.path.join('./geo/'+gsmid[:7]+'/'+gsmid+'.xml'))
        tmp = []
        for l in f:
            if l.find("xmlns=\"http://www.ncbi.nlm.nih.gov/geo/info/MINiML\"") == -1:
                tmp.append(l)
    if not docString and not isinstance(f, str):
        f.close()
    return "".join(tmp)


def _getFieldXML(sample_path, fields = ['Sample/Library-Strategy', 'Sample/Description', 'Sample/Data-Processing',
    'Sample/Channel/Extract-Protocol', "Sample/Title", 'Sample/Channel/Source', 'Sample/Channel/Characteristics']):
    """
    get text of items in filds list from XML
    we need these info for matching key words, like single cell
    """
    text = readGeoXML(sample_path)
    ## parse XML content
    rt = ET.fromstring(text)
    info = {}
    for field in fields:
        tmp = rt.findall(field)
        if tmp:
            info[field] = tmp[0].text.strip()
        else:
            pass
    return info


def _matchKeyWord(xmlContent, key, fileds=False):
    ## match key words in a specific xml content
    ## return the match words
    res = {}
    if not fileds:
        fileds = list(xmlContent.keys()) # use all fields if not speficify
    for field in xmlContent.keys():
        if field in fileds:
            tmp = list(set(re.findall(r'%s'%key, xmlContent[field].replace('-', '').replace('_', ''), re.I)))
            if tmp:
                res[field]= tmp # a list in dict
    return res

def _match_scRNAseq(xmlContent):
    """
    match key words, like single cell, single cell RNA-seq, and sequencing platform,etc
    return a dict, contains geo items and matched keys
    """
    if not xmlContent:
        return {}
    #filter bulk RNA-seq in single cell RNA-seq GSE, uploder mentioned single cell but actually bulk RNA-seq 
    ## bulk RNA-seq apears in sample description, title, and source name
    for key in ['bulk rnaseq', 'bulk']:
        bulk1 = _matchKeyWord(xmlContent, key = 'bulk rnaseq', fileds = ['Sample/Description', 'Sample/Title', 'Sample/Channel/Source'])
        if bulk1:
            return {}
    ## or library stratey : bulk RNA-seq appears in any fields 
    bulk2 = _matchKeyWord(xmlContent, key = "library strategy: bulk rnaseq")
    if bulk2:
        return {} # return empty, if bulk RNA-seq appears in Title 
    match_res = {}
    ## 1. match with single cell words, remove special characters, like '-'
    for key1 in ['single cell', 'scrnaseq', 'singlecell rnaseq']:
        tmp = _matchKeyWord(xmlContent, key1)
        if tmp:
            for i in tmp.keys():
                match_res[i].extend(tmp[i]) if i in match_res.keys() else match_res.update({i:tmp[i]})
    # print(match_res)
    tmp = []
    ## 2. match with platform words
    keys = {}
    with open('platform.txt') as f:
        for line in f:
            line = line.rstrip().split('\t')
            keys[line[0]] = line[1] # platforms of sequecing, like 10X Genomics, Smart-seq2
    for key2 in keys.keys():
        tmp = _matchKeyWord(xmlContent, keys[key2]) # result of key word matching, a list in dict
        if tmp:
            for i in tmp.keys():
                # print(tmp)
                match_res[i].extend(tmp[i]) if i in match_res.keys() else match_res.update({i:tmp[i]})
    return match_res


def _checkType(acc, sample_path, type_need):
    """
    check whether it is single cell RNAseq or single cell ATAC-seq
    based on key words matching
    """ 
    ## read in XML to get sample description
    ret = {}
    if sample_path and os.path.isfile(sample_path):
        #NOTE: we need readGeoXML to process
        xmlContent = _getFieldXML(sample_path)
        # parse single cell
        if 'sc-rna-seq' in type_need:
            # single cell RNA-seq, Library-strategy nmush be RNA-Seq
            if ( 'Sample/Library-Strategy' in xmlContent.keys() ) and ( xmlContent['Sample/Library-Strategy'] == 'RNA-Seq' ):
                res = _match_scRNAseq(xmlContent)
                if res:
                    ret[acc] = res
                    return ret
    if sample_path and not os.path.isfile(sample_path):
        xml = getGeoXML(accession=acc)
        ret = _checkType(acc = acc, sample_path = sample_path, type_need = type_need)
        return ret
    return None
            

def getGeoSamples_byType(ddir="geo", ttype=["sc-rna-seq", "sc-atac-seq"], unique_ids=False, refresh=False):
    """A filter for our Geo model; searches our db for the specific sample
    type.
    NOTE: hones in on Library-Strategy tag

    Returns a list of samples fitting the specified

    NOTE: building this up takes time, around 10 secs almost 1 minute!
    TRY: caching the result, reading from a cached file takes only 1 minute
    Store them in files by the .ttype--in the local dir
    """
    ret = {}
    if not unique_ids: #qury all local samples
        #NEED to generate the file, and make the call recursively
        #actually, just one level of recursion b/c geo is pretty flat
        p = os.path.join(ddir)
        ls = os.listdir(p)
        ls = [x for x in ls if x.startswith('GSM') and x != 'GSM'] # for real gsm ids
        for df in ls:
            path = os.path.join(p, df)
            if os.path.isfile(path): #it's a file--check if it's ChIP-Seq
                print(path)
                typo = _checkType(acc=df.split(".")[0], sample_path=path, type_need=ttype) # check whether the seq type is chip-seq, atac, or dnase
                if typo:
                    ret.update(typo)
            else:
                #it's a dir recur
                newd = os.path.join(ddir, df)
                newdict = getGeoSamples_byType(ddir=newd, ttype=ttype, refresh=refresh)
                if newdict and (type(newdict) == type(ret)):
                    ret = dict(ret, **newdict)
    elif unique_ids: # query just for the specified gsm
        for gsmid in unique_ids:
            p = os.path.join(ddir+'/'+gsmid[:7]+'/'+gsmid+'.xml')
            typo = _checkType(acc=gsmid, sample_path=p, type_need=ttype)
            if typo:
                ret.update(typo)
    else:
        pass
    return ret


def getGeoSamples_byTypes(path, ddir="geo", datatype=False, gsmids=False, refresh=False): #ttypes = ["ATAC-Seq"]): #"ChIP-Seq", "DNase-Hypersensitivity"]):
    ret = []
    if not refresh and os.path.exists(path):
        ret = pickle.load(open(path))
        return ret
#    for t in ttypes:
    if datatype and gsmids:
        ret = getGeoSamples_byType(ddir=ddir, ttype=datatype, unique_ids=gsmids, refresh=refresh)
    elif datatype and not gsmids:
        ret = getGeoSamples_byType(ddir=ddir, ttype=datatype, refresh=refresh)
    elif gsmids and not datatype:
        ret = getGeoSamples_byType(ddir=ddir, unique_ids=gsmids, refresh=refresh)
    else:
        ret = getGeoSamples_byType(ddir=ddir, refresh=refresh)
    # pickle.dump(ret, open(path, "w"))
    return ret


def sync_samples(fsave, fill_or_not=False, DataType=False, dateRegion = False, refresh=False):
    """
    get IDs from GDS and convert to GSM id for given time stamp or most recent 100000
    DataType: this script are allowed to parse sc-rna-seq or sc-atac-seq
    """
    f = open(fsave, 'w')
    f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%('GSM', 
        'species', 'GSE', 'PMID', 'Paper', 'Title', 'CellType', 'Tissue', 'Disease',
        'Cell_Pop', 'Release_Date', 'fields_Annotation', 'fields_dataType', 'Key_Match')+'\n')
    f.close()

    if refresh:
        if DataType:
            datatype = ['sc-rna-seq', 'sc-atac-seq']
        else:
            datatype = False
        getSyncLog("# 1. resync the repository from Internet")# output record to log file
        local_repo_path = "repository_samples.pickle"
        local_db_samples = set(models.Samples.objects.values_list('unique_id', flat=True))
    
        gdsSamples = getGDSSamples(dateRegion) # get GDS ids
        begin = datetime.now()
        getSyncLog('start %s: There are %s GDS Samples in sum'%(dateRegion, len(gdsSamples)))# output record to log file

        gsm_collect = [] # collect GSM ids
        one_percent = len(gdsSamples)/100
        cnt = 0
        for gdsid in gdsSamples:
            cnt += 1
            if cnt % one_percent == 0:
                getSyncLog("%s%%"%(cnt/one_percent))
            gsm = gsm_idToAcc(gdsid) # convert GDS id to GSM id
            if gsm and (gsm not in local_db_samples) and (gsm not in gsm_collect):
                geoXML = getGeoXML(gsm) # download XML file
                
                """if get one xml of gsm, then parser the information immediately"""
                getType = getGeoSamples_byTypes(path=local_repo_path, datatype=datatype, gsmids=[gsm], refresh=refresh)
                if getType:
                    for i in getType.keys():
                        # parse sample annotation
                        list_sample = scrna_parser_sampleDetail.update_one_sample(gsmid=gsmid, ddir=xmlPath)
                        try:
                            list_sample.append(str(getType[i])) # add matched key words
                            f = open(fsave, 'a')
                            f.write('\t'.join(list_sample)+'\n')
                            f.close()
                        except:
                            getSyncLog("Error when writing in table: %s" % s)
                else:
                    out = open(fsave+'_others.txt', 'a')
                    out.write(iterm+'\n')
                    out.close()
                # out = open(fsave, 'a')
                # if getType:
                #     for i in getType.keys():
                #         out.write(str(i)+'\t'+str(getType[i])+'\n')
                # else:
                #     out.write(str(gsmid)+'\tNone'+'\n')
                # out.close()  
                    # gsm_collect.append(gsm) # collect gsm of just parsed samples
                    #getSyncLog('parse details: %s'%gsm)
                    # parse_detail_linker(fsave, fill_or_not, need_added_samples=getType, FactorType = DataType)

def sync_samples_from_gsm_factor(infile, gsm_col, factor_col, fsave, xmlPath='geo', fill_or_not = False, refresh = False):
    getSyncLog("try to add samples based on outside table which contain gsm ID and factor name")
    f = open(fsave, 'w')
    f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%('GSM', 
        'species', 'GSE', 'PMID', 'Paper', 'Title', 'CellType', 'Tissue', 'Disease',
        'Cell_Pop', 'Release_Date', 'fields_Annotation', 'fields_dataType', 'Key_Match')+'\n')
    f.close()
    need_added_samples = [x.rstrip().split('\t') for x in open(infile)]
    local_db_samples = set(models.Samples.objects.values_list('unique_id', flat=True))
    getSyncLog('totally %d samples need to be added'%len(need_added_samples))
    n = 1
    for iterm in need_added_samples:
        print(n) # let me know which the process 
        if (not refresh) and (iterm[int(gsm_col)] in local_db_samples):
            continue
        if factor_col:
            factor_know = iterm[int(factor_col)]
        else:
            factor_know = False
        # get seqtype
        gsmid = iterm[int(gsm_col)]
        print(gsmid)
        getType = getGeoSamples_byTypes(path = "repository_samples.pickle", datatype = ['sc-rna-seq', 'sc-atac-seq'],
         gsmids=[gsmid], refresh=refresh, ddir = xmlPath)
        if getType:
            for i in getType.keys():
                # parse sample annotation
                list_sample = scrna_parser_sampleDetail.update_one_sample(gsmid=gsmid, ddir=xmlPath)
                try:
                    list_sample.append(str(getType[i])) # add matched key words
                    f = open(fsave, 'a')
                    f.write('\t'.join(list_sample)+'\n')
                    f.close()
                except:
                    getSyncLog("Error when writing in table: %s" % s)
        else:
            out = open(fsave+'_others.txt', 'a')
            out.write('\t'.join(iterm)+'\n')
            out.close()


def getLocalGeo(fsave, fill_or_not=False, xmlPath="geo", DataType=False, refresh = False):
    """This function can be used if we have some xml file of GEO, then 
    recursion all local existed geo xml and get we wanted samples (ChIP, ATAC, DNase)
    """
    f = open(fsave, 'w')
    f.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%('GSM', 
        'species', 'GSE', 'PMID', 'Paper', 'Title', 'CellType', 'Tissue', 'Disease',
        'Cell_Pop', 'Release_Date', 'fields_Annotation', 'fields_dataType', 'Key_Match')+'\n')
    f.close()
    getSyncLog("# 1. go through all the xml and check the type in the path of %s"%xmlPath)
    local_repo_path = "repository_samples.pickle"
    if DataType:
        datatype = ['sc-rna-seq', 'sc-atac-seq']
    else:
        datatype = False
    local_repo_samples_dict = getGeoSamples_byTypes(path=local_repo_path, datatype=datatype, ddir=xmlPath, refresh=refresh)
    out = open(fsave+'_gsm.txt', 'w')
    out.write('\n'.join(list(set(local_repo_samples_dict.keys()))))
    out.close()
    local_repo_samples = set(local_repo_samples_dict.keys())

    getSyncLog("# 2. calculate new samples")

    # local_db_samples = set(models.Samples.objects.values_list('unique_id', flat=True))

    # getSyncLog("There are %d samples in local repo." % len(local_repo_samples))
    # getSyncLog("There are %d samples in local db." % len(local_db_samples))
    # need_added_samples = sorted(list(local_repo_samples - local_db_samples))
    for gsmid in local_repo_samples:
        print(gsmid)
        getType = getGeoSamples_byTypes(path = "repository_samples.pickle", datatype = ['sc-rna-seq', 'sc-atac-seq'], gsmids=[gsmid], refresh=refresh)
        if getType:
            for i in getType.keys():
                # parse sample annotation
                list_sample = scrna_parser_sampleDetail.update_one_sample(gsmid=gsmid, ddir=xmlPath)
                try:
                    list_sample.append(str(getType[i])) # add matched key words
                    f = open(fsave, 'a')
                    f.write('\t'.join(list_sample)+'\n')
                    f.close()
                except:
                    getSyncLog("Error when writing in table: %s" % s)
        # else:
        #     out = open(infile+'_others.txt', 'a')
        #     out.write('\t'.join(iterm)+'\n')
        #     out.close()


