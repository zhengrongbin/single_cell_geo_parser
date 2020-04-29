import json
import os
import re
import sys
import urllib.request, urllib.parse, urllib.error
import traceback
from datetime import datetime
import time
import pickle
import pickle
import subprocess
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

import env
models = env.models

import pubmed
import scrna_parser  


def _parse_a_field(description_dict, a_field, DCmodel, max_create_length=100, new=False):
    if  not description_dict.get(a_field, None):
        return None
    if len(description_dict.get(a_field, "")) > 0:
        if DCmodel in [models.CellTypes]: #remove '-' when match cell types, like T-cells > T cells
            description_dict_new = {}
            for v in list(description_dict.keys()):
                description_dict_new[v] = description_dict[v].replace('-', ' ')
            description_dict = description_dict_new  
        result_searched_by_name = sorted(DCmodel.objects.extra(
            where={"%s REGEXP CONCAT('([^a-zA-Z0-9]|^)', `name`, '([^a-rt-zA-RT-Z0-9]|$)')"},
                                  params=[description_dict[a_field]]),
            key=lambda o: len(o.name),
            reverse=True)

        if result_searched_by_name and len(result_searched_by_name[0].name.strip()) > 0:
            return result_searched_by_name[0]

        if DCmodel not in [models.Factors, models.Aliases] :
            result_searched_by_aliases = sorted(DCmodel.objects.exclude(aliases=None).exclude(aliases="").extra(
                where={"%s REGEXP CONCAT('([^a-zA-Z0-9]|^)', `aliases`, '([^a-rt-zA-RT-Z0-9]|$)')"},
                params=[description_dict[a_field]]),
                    key=lambda o: len(o.name),
                    reverse=True)
            if result_searched_by_aliases and len(result_searched_by_aliases[0].name.strip()) > 0:
                return result_searched_by_aliases[0]

    if new and (len(description_dict.get(a_field, "")) > 0) and (len(description_dict.get(a_field, "")) <= max_create_length):
        #return description_dict[a_field]
        if DCmodel in [models.CellLines, models.CellTypes, models.TissueTypes]: # these three will check wheather they appare in each other
            return description_dict[a_field]
        else:
            # ret, created = DCmodel.objects.get_or_create(name=description_dict[a_field])
            # if created:
            #     ret.status = 'new'
            return description_dict[a_field]


    return None

def _parse_fields(description_dict, strict_fields, greedy_fields, DCmodel, greedy_length=100):
    for sf in strict_fields:
        ret = _parse_a_field(description_dict, sf, DCmodel, greedy_length, new=False)
        if ret:
            return ret

            
    for gf in greedy_fields:
        ret = _parse_a_field(description_dict, gf, DCmodel, greedy_length, new=True)
        if ret:
            return ret
    if DCmodel in [models.Factors, models.Aliases]: #do not use all description information for factor part, since factor should be parse more serious.
        return None
    tmp = {}
    for x in list(description_dict.keys()): # using all description to parse information
        if x not in ['antibody', 'chip_antibody', 'last update date', 'release date']:
            tmp[x] = description_dict[x]
    characteristics = {'characteristics':' '.join(list(tmp.values()))}
    ret = _parse_a_field(characteristics, 'characteristics', DCmodel, greedy_length, new = False)
    if ret:
        return ret

    if DCmodel in [models.CellLines]: # search cell line by pubic database information
        ret = geo_parser_newVersion.search_cellline_from_out(characteristics, 'characteristics')
        if ret:
            if ret and (str(ret) not in ['OF', 'IP']):
                return ret
    return None

def parseCellType(description_dict):
    return _parse_fields(description_dict,
                         ['cell type', 'cell lineage', 'cell', 'cell line', 'source name', 'cell description', 'title', ],
                         ['cell type'],
                         models.CellTypes)


def parseCellPop(description_dict):
    return _parse_fields(description_dict, ['cell', 'source name', 'cell description', 'title', 'cell type', 'cell lineage'], [], models.CellPops)


def parseTissue(description_dict):
    return _parse_fields(description_dict,
                         ['tissue', 'tissue type', 'tissue depot', 'source name', 'cell description', 'title', 'cell type', 'cell lineage','cell', 'cell line'],
                         ['tissue', 'tissue type'],
                         models.TissueTypes)


def parseStrain(description_dict):
    return _parse_fields(description_dict,
                         ['strain', 'strain background', 'source name', 'cell description', 'title'],
                         ['strain'],
                         models.Strains)


def parseDisease(description_dict):
    return _parse_fields(description_dict,
                         ['disease', 'tumor stage', 'cell karotype', 'source name', 'title'],
                         ['disease'],
                         models.DiseaseStates)


def parseReleaseTime(description_dict):
    # a trick to convert time string into a datetime object
    time_struct = time.strptime(description_dict["release date"], "%Y-%m-%d")
    return datetime.fromtimestamp(time.mktime(time_struct))


def search_between_table(description_dict, strict_fields, serious_fields, DCmodel):
    for sf in strict_fields:
        sea = _parse_a_field(description_dict, sf, DCmodel, 100, new=False)
        if sea:
            return sea
    tmp = {}
    for k in list(description_dict.keys()):
        if k not in ['antibody', 'chip_antibody', 'last update date', 'release date']:
            tmp[k] = description_dict[k]
    characteristics = {'characteristics':' '.join(list(tmp.values()))}
    sea = _parse_a_field(characteristics, 'characteristics', DCmodel, 100, new = False)
    if sea:
        return sea
    if DCmodel in [models.CellTypes]:
        """try metamap if no cell type"""
        feature = characteristics['characteristics'].replace('(', ' ').replace(')', ' ')
        try:
            content = subprocess.getoutput('echo "%s" | metamap -y -I '%(feature))
            content = [x for x in content.split('\n') if x.endswith('[Cell]')]
            if content:
                proCell = content[0][content[0].index(':')+1:].rstrip('[Cell]').strip()
                if proCell:
                    Cell = proCell
                    if '(' in proCell:
                        Cell = proCell[:proCell.index('(')].strip('(').strip()
                    if Cell.lower() not in ['cell', 'cancer', 'tumor', 'clone', 'cell line', 'cells', 'human cell line', 'mouse cell line', 'cellline', 'celllines']:
                        if (not _parse_a_field({'cellType':Cell}, 'cellType', models.CellPops, 100, new = False)) and (not _parse_a_field({'cellType':Cell}, 'cellType', models.CellLines, 100, new = False)) : # in case metamap result in cell pop
                            print("metamap cell type")
                            return Cell
        except:
            pass
    if serious_fields:
        for sf in serious_fields:
            tmp = description_dict.get(sf, "")
            if tmp and (DCmodel == models.CellTypes) and (not _parse_a_field({'cell':tmp}, 'cell', models.CellLines, 100, new = False)) and (not _parse_a_field({'cell':tmp}, 'cell', models.TissueTypes, 100, new = False)) and (not _parse_a_field({'cell':tmp}, 'cell', models.CellPops, 100, new = False)):
                return tmp
            if tmp and (DCmodel == models.CellLines) and (not _parse_a_field({'cell':tmp}, 'cell', models.CellTypes, 100, new = False)):
                return tmp
            if tmp and (DCmodel == models.TissueTypes) and (not _parse_a_field({'cell':tmp}, 'cell', models.CellTypes, 100, new = False)) and (not _parse_a_field({'cell':tmp}, 'cell', models.CellLines, 100, new = False)):
                return tmp
            if tmp and (DCmodel not in [models.CellTypes, models.CellLines, models.TissueTypes]):
                return tmp
    return None

def parseAndsearch(description_dict, field):
    tmp_sea_cellLine = search_between_table(description_dict,
                         field, ['cell line'],
                         models.CellLines)
    tmp_sea_cellType = search_between_table(description_dict,
                         field, ['cell type'],
                         models.CellTypes)
    tmp_sea_tissueType = search_between_table(description_dict,
                         field, ['tissue', 'tissue type'],
                         models.TissueTypes)
    tmp_sea_cellpop = search_between_table(description_dict,
                         field, [],
                         models.CellPops)
    tmp_sea_disease = search_between_table(description_dict,
                         field, [],
                         models.DiseaseStates)
    if str(tmp_sea_tissueType).lower() in ['primary tumor', 'tumor']:
        tmp_sea_tissueType = None
    return {'cellType':tmp_sea_cellType, 'cellLine':tmp_sea_cellLine, 'tissueType':tmp_sea_tissueType, 'cellpop':tmp_sea_cellpop, 'disease':tmp_sea_disease}               
        

def parseGeoInfo(accession, ddir='geo'):
    """parse necessary information (detailed description from the Geo XML file"""
    xmlString = scrna_parser.readGeoXML(None, scrna_parser.getGeoXML(accession, path=ddir))
    tree = ET.fromstring(xmlString)
    ret = {}
    for node in tree.findall("Sample/Channel/Characteristics"):
        if node.get("tag"):
            ret[node.get("tag").replace("_", " ")] = node.text.strip()
    if not ret:
        ret["Characteristics"] = tree.findall("Sample/Channel/Characteristics")[0].text.strip()
    ret["source name"] = tree.find("Sample/Channel/Source").text.strip()
    ret["title"] = tree.find("Sample/Title").text.strip()
    ret['last update date'] = tree.find("Sample/Status/Last-Update-Date").text.strip()
    ret['release date'] = tree.find("Sample/Status/Release-Date").text.strip()
    return ret

def cleanCategory(s):
    """Given a string, replaces ' ' with '_'
    '/', '&', '.', '(', ')'with ''
    """
    tmp = s.replace(" ", "_")
    for bad in ['/', '&', '.', '(', ')', ',']:
        tmp = tmp.replace(bad, "")
    return tmp


def postProcessGeo(accession, ddir='geo', docString=None):
    """post processes the GEO record to feed into the classifiers
    uses docString IF it is available, otherwise will read the record
    """
    #ignore these tags
    ignore = ["Growth-Protocol", "Extract-Protocol", "Treatment-Protocol"]
    _thresh = 10 #max 10 words

    path = os.path.join(_ppath, "geoPost")
    if not os.path.exists(path):
        os.mkdir(path)
    subdir = os.path.join(path, accession[:7])
    if not os.path.exists(subdir):
        os.mkdir(subdir)
    path = os.path.join(subdir, "%s.txt" % accession)

    if os.path.exists(path):
        f = open(path)
        docString = f.read()
        f.close()
        return docString
    else:
        #need to build the doc
        if not docString:
            #read the document from geo
            docString = scrna_parser.getGeoXML(accession, ddir)
            #process the string

        if not docString:
            return None

        text = scrna_parser.readGeoXML(None, docString=docString)

        try:
            root = ET.fromstring(text)
        except:
            print("Could not parse: %s" % accession)
            print('-' * 60)
            traceback.print_exc(file=sys.stdout)
            print('-' * 60)
            return None

        #3. collect all of the information under Sample/Channel
        tmp = []

        #Things in sample to get:
        ls = ['Title', 'Source', 'Library-Strategy', 'Library-Source',
              'Library-Selection', 'Supplementary-Data']
        for t in ls:
            tag = root.findall("Sample/%s" % t)
            for tt in tag:
                tmp.append("%s(%s)" % (t.upper(), tt.text.strip().upper()))

        channels = root.findall("Sample/Channel")

        for c in channels:
            for child in c:
                category = ""

                if child.tag in ignore:
                    continue

                #Special case: Characteristic--take the tag attrib
                if child.tag == "Characteristics":
                    if "tag" in child.attrib and child.attrib["tag"]:
                        category = child.attrib["tag"].lstrip().rstrip()
                else:
                    category = child.tag.lstrip().rstrip()

                #convert categories like "cell line" to "CELL_LINE"
                #tmp.append("%s(%s)" % (category.replace(" ","_").upper(),
                #                       child.text.strip()))
                val = child.text.strip()
                #THRESHOLD: values can be at most 10 words
                if len(val.split()) <= _thresh:
                    tmp.append("%s(%s)" % (cleanCategory(category).upper(),
                                           val.upper()))
                else:
                    #take first 10 words
                    tmp.append("%s(%s)" % (cleanCategory(category).upper(),
                                           " ".join(val.split()[:_thresh]).upper()))


        #4. write the information to file
        f = open(path, "w")
        f.write("%s" % smart_str("\n".join(tmp)))
        f.close()
    
        return "\n".join(tmp)

### HELPER fns
def getFromPost(geoPost, cat):
    """tries to search for cat(VALUE) returns VALUE if found, otherwise ""
    NOTE: categories can be in UPPER or lowercase, eg. TITLE or title
    """
    m = re.search("%s\((.+)\)" % cat.upper(), geoPost)
    if m:
        return m.group(1)
    else:
        return ""
### Librarian fns
def gsmToGse(gsmid):
    """Given a gsmid, will try to get the geo series id (GSE) that the
    sample is associated with; if it is associated with several GSEs, then
    returns the first.

    STORES the GSE ID, e.g. 200030833 in the file
    NOTE: if we want GSEs then we need to translate IDs to GSEXXXXX just
          like we do for GSMs above

    uses this query:
    http://www.ncbi.nlm.nih.gov/gds/?term=gse%5BEntry+Type%5D+AND+GSM764990%5BGEO+Accession%5D&report=docsum&format=text
    """
    ret = None
    path = os.path.join(_ppath, "GSM_GSE")
    if not os.path.exists(path):
        os.mkdir(path)
    subdir = os.path.join(path, gsmid[:7])
    if not os.path.exists(subdir):
        os.mkdir(subdir)
    path = os.path.join(subdir, "%s.txt" % gsmid)
    if os.path.exists(path):
        f = open(path)
        ret = f.read().strip()
        f.close()
    else:
        #This URL is slow!
        #NOTE: for every ncbi query, try to find the eutils equivalent!
        #--it's faster
        #URL = "http://www.ncbi.nlm.nih.gov/gds/?term=gse[Entry Type] AND %s[GEO Accession]&report=docsum&format=text" % gsmid
        URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=gse[Entry%20Type]%20AND%20{acc}[GEO%20Accession]".format(acc=gsmid)

        try:
            #print "gsmToGse: %s" % URL
            #signal.alarm(180)
            urlf = scrna_parser.proxyInstead(link=URL)
            root = ET.fromstring(urlf)

            #Get the IDList
            tmp = root.findall("IdList/Id")
            if tmp:
                ret = tmp[0].text.strip()
                f = open(path, "w")
                f.write(ret)
                f.close()

        except:
            # try again:
            time.sleep(0.3)
            try:
                URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=gse[Entry%20Type]%20AND%20{acc}[GEO%20Accession]".format(acc=gsmid)
                urlf = scrna_parser.proxyInstead(link=URL)
                root = ET.fromstring(urlf)
                #Get the IDList
                tmp = root.findall("IdList/Id")
                if tmp:
                    ret = tmp[0].text.strip()
                    f = open(path, "w")
                    f.write(ret)
                    f.close()
            except:
                print("gsmToGse")
                print("URL is" + URL)
                print('-' * 60)
                traceback.print_exc(file=sys.stdout)
                print('-' * 60)

    return ret


def gseToPubmed(gseid):
    """Given a gseid, will try to get the pubmed id
    """
    ret = None
    path = os.path.join(_ppath, "GSE_PUB")
    if not os.path.exists(path):
        os.mkdir(path)
    subdir = os.path.join(path, gseid[:6])
    if not os.path.exists(subdir):
        os.mkdir(subdir)
    path = os.path.join(subdir, "%s.txt" % gseid)
    if os.path.exists(path):
        f = open(path)
        ret = f.read().strip()
        f.close()
    else:
        URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?id=%s&db=pubmed&dbfrom=gds" % gseid

        try:
            #print "gseToPubmed: %s" % URL
            #signal.alarm(180)
            urlf = scrna_parser.proxyInstead(link=URL)
            root = ET.fromstring(urlf)

            #Get the IDList
            tmp = root.findall("LinkSet/LinkSetDb/Link/Id")
            if tmp:
                ret = tmp[0].text.strip()
                f = open(path, "w")
                f.write(ret)
                f.close()
        except:
            print("gsmToGse")
            print('-' * 60)
            traceback.print_exc(file=sys.stdout)
            print('-' * 60)
    return ret

def gse_idToAcc(gdsId):
    """Given a GDS id, e.g. 200030833, tries to give a GDS accession, e.g.
    GSE30833

    NOTE: there is an algorithm: acc = "GSE"+gdsId[1:] (strip leading 0s)
    """
    #Cut = dropping of the "2" (which indicates series) and removal of leading
    #leading 0s
    cut = gdsId[1:].lstrip("0")
    return "GSE%s" % cut

def update_one_sample(gsmid, ddir='geo', parse_fields=['other_ids', 'paper', 'name', 'species', 'description', 'antibody', 'factor',
                                           'cell type', 'cell line', 'cell pop', 'tissue', 'strain', 'disease','update date','release date']):
    """Given a gsmid, tries to create a new sample--auto-filling in the
    meta fields


    If overwrite is True and there is the sample that has the same gsmid, this function will overwrite that sample

    NOTE: will try to save the sample!!

    Returns newly created sample
    """

    print('+++ description')
    description_dict = parseGeoInfo(gsmid, ddir)

    geoPost = postProcessGeo(gsmid, ddir=ddir)
    if not geoPost:
        return None

    if 'species' in parse_fields:
        print('+++ species')
        if getFromPost(geoPost, "organism") == "HOMO SAPIENS":
            species = models.Species.objects.get(pk=1)
        else:
            species = models.Species.objects.get(pk=2)

    if ('other_ids' in parse_fields) or ('paper' in parse_fields):
        print('+++ other IDs')
        gseId = gsmToGse(gsmid)
        if not gseId: # sometimes the gse parser failed, but can get it with some tries, strange!
            time.sleep(0.3)
            print('+++ again gse')
            gseId = gsmToGse(gsmid)
        pmid = gseToPubmed(gseId) if gseId else None

    if 'other_ids' in parse_fields:
        import json 
        idList = {'gse': gseId, 'pmid': pmid}
        print(idList)
        other_ids = json.dumps(idList)
        series_id = None
        try:
            if idList['gse']:
                series_id = gse_idToAcc(idList['gse'])#str(idList['gse'][4:])
        except:
            print('cannot find GSE_id')
    
    paper = None
    if 'paper' in parse_fields and pmid:
        print('+++ paper')
        try:
            paper = pubmed.getOrCreatePaper(pmid)
        except:
            paper = None

    if 'name' in parse_fields:
        print('+++ title')
        name = getFromPost(geoPost, "title")

    if 'species' in parse_fields:
        print('+++ species')
        if getFromPost(geoPost, "organism") == "HOMO SAPIENS":
            species = models.Species.objects.get(pk=1)
        else:
            species = models.Species.objects.get(pk=2)

    #HERE is where I need to create a classifier app/module
    #FACTOR, platform, species--HERE are the rest of them!

    if 'description' in parse_fields:
        print('+++ add description')
        description = json.dumps(description_dict)

    if 'cell type' in parse_fields:
        print('+++ cell type')
        tmp_celltype = None
        # get first parsed cell type information
        searchCellType = parseAndsearch(description_dict, ['cell type', 'cell lineage', 'cell', 'cell line', 'source name', 'cell description', 'title'])
        if searchCellType['cellType'] and (str(searchCellType['cellType']).upper() not in [str(searchCellType['cellLine']).upper(), str(searchCellType['tissueType']).upper()]):
            tmp_celltype = searchCellType['cellType'] # use the cell type if parsed information not in other tables. else use "None" defined before
    
    if 'tissue' in parse_fields:
        print('+++ tissue')
        tmp_tissue = None
        searchTisssue = parseAndsearch(description_dict, ['tissue', 'tissue type', 'tissue depot', 'source name', 'cell description', 'title', 'cell type', 'cell lineage','cell', 'cell line'])
        if searchCellType['tissueType'] and (str(searchCellType['tissueType']).upper() not in [str(searchCellType['cellLine']).upper(), str(searchCellType['cellType']).upper(), str(searchCellType['cellpop']).upper(), str(searchCellType['disease']).upper()]):
            tmp_tissue = searchTisssue['tissueType']
        else:
            if tmp_tissue:
                test_tissue = parseAndsearch({'cell type':str(tmp_tissue)}, ['cell type']) # test parsed tissue information whether in cell type table
            else:
                test_tissue = {'cellType':None}
            if test_tissue['cellType']: # means parsed tissue information in cell type table, then ignore tissue
                tmp_tissue = None

    if 'disease' in parse_fields:
        print('+++ disease')
        disease_state = parseDisease(description_dict)
        if searchCellType['disease']:
            disease_state = searchCellType['disease']

    if 'cell pop' in parse_fields:
        print('+++ cell pop')
        cell_pop = parseCellPop(description_dict)
        if searchCellType['cellpop']:
            cell_pop = searchCellType['cellpop']

    if 'release date' in parse_fields:
        geo_release_date = parseReleaseTime(description_dict)
    # fields for checking sample type
    sample_path = os.path.join(ddir+'/'+gsmid[:7]+'/'+gsmid+'.xml')
    os.system('echo %s'%sample_path)
    xmlContent = scrna_parser._getFieldXML(sample_path)

    res = [gsmid, str(species), str(series_id), str(pmid), str(paper), str(name), str(tmp_celltype),
        str(tmp_tissue), str(disease_state), str(cell_pop), str(geo_release_date), str(description_dict), str(xmlContent)]
    # time.sleep(0.3)
    return res
