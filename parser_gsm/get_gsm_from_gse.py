import os,sys
import urllib.request
import re

gseFile = sys.argv[1] # one-column file of GSE id

def _parse_gsm_from_gse(gse):
    url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s'%gse
    gse_handler = urllib.request.urlopen(url)
    gse_html = gse_handler.read()
    gse_html = gse_html.decode('utf-8', errors='ignore')
    gse_handler.close()
    regexp = re.compile('GSM[0-9]*')
    gsminfor = re.findall(r'GSM[0-9]*', gse_html)
    if gsminfor:
        gsm = list(set(gsminfor))
        return(gsm)
    return([])

gse_ids = [x.strip() for x in open(gseFile)]
for gse in gse_ids:
    try:
        gsms = _parse_gsm_from_gse(gse)
        if len(gsms)>0:
            print(gse, len(gsms))
            for gsm in gsms:
                out = open(gseFile+'_gsm.txt', 'a')
                out.write(gsm+'\t'+gse+'\n')
                out.close()
    except:
        traceback.print_exc(file=sys.stdout)
