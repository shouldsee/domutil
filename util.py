# from modeller import *
# from modeller.scripts import complete_pdb
# ### Enable this line to reduce verbositys
# # log.none()
import sys,os
import StringIO
import contextlib
import re 
import numpy as np

full = lambda p: os.path.expandvars(p)

p_nb = re.compile("Number of non-bonded pairs \(excluding 1-2,1-3,1-4\): *([0-9]*)")
p_energy=re.compile("Current energy *: *([0-9,\.,-]*)")
p_atomCount = re.compile("Number of all, selected real atoms *: ([0-9, ]{10})")
p_resCount = re.compile("Number of all residues in MODEL *: *([0-9]*)")
p_header = re.compile("NAME.*?\n")
p_hmmlen = re.compile('LENG  (\d+)\n')
p_cathdomain = re.compile("([0-9,a-z,A-Z]{7})")
p_cathFAheader=re.compile('\|([0-9,a-z,A-Z]*)\/')

# notlist = lambda lst: [not x for x in lst]
notlist = lambda lst: [not x for x in lst]
listNOT = lambda lst: [not x for x in lst]
listAND=lambda alst,blst: [ x & y for x,y in izip(alst,blst)]
listANDNOT=lambda alst,blst: [ x and not y for x,y in izip(alst,blst)]
list2dict = lambda lst: {v:i for i,v in enumerate(lst) }
get_vname = lambda var:  [ k for k,v in locals().items() if v is var][0]


a_href = lambda text,url: "<a href='%s'>%s</a>" % (url,text)


def reset_database_connection():
    from django import db
    db.close_old_connections()

import cPickle as pk
def pk_load(fname, cache_dir = 'data/'):
    fname = cache_dir + fname
    return pk.load(open(fname, 'rb'))
    
# def pk_dump(v, alias, cache_dir = 'data/'):
#     vname,var = v
#     # vname = get_vname( var )
#     cdir = cache_dir + vname + '/'
#     if not os.path.isdir(cdir):
#         os.makedirs(cdir)

#     fname = cache_dir + vname + '/' + alias
#     pk.dump( var, open(fname,'wb'))

def pk_dump(var, fname, cache_dir = 'data/'):
    fname = cache_dir + fname
    cdir = os.path.dirname(fname)
    if not os.path.isdir(cdir):
        os.makedirs(cdir)
    pk.dump( var, open(fname,'wb'))

levels=[ None,
'root',
'Class',
'arch',
'topo',
'homsf',
's35',
's60',
's95',
's100'];




forward_field = {
    'S':'cath_node__id',
    'H':'cath_node__parent__id',
    'T':'cath_node__parent__parent__id',
    'A':'cath_node__parent__parent__parent__id',
    'C':'cath_node__parent__parent__parent__parent__id',
}

reverse_field = {
    'S':'hmmprofile__hits',
    'H':'classification__hmmprofile__hits',
    'T':'classification__classification__hmmprofile__hits',
    'A':'classification__classification__classification__hmmprofile__hits',
    'C':'classification__classification__classification__classification__hmmprofile__hits'
}



@contextlib.contextmanager
def stdoutIO(stdout=None):
    oldout = sys.stdout
    olderr = sys.stderr
    if stdout is None:
        stdout = StringIO.StringIO()
    sys.stdout = stdout
    sys.stderr = stdout
    yield stdout
    sys.stdout = oldout
    sys.stderr = olderr


import urllib2

def download(url):
    request = urllib2.Request( url )
    response = urllib2.urlopen(request)
    s = response.read()
    print "finshed download %s" % url
    return s
def get_gzip(url = 'http://download.cathdb.info/cath/releases/daily-release/newest/cath-b-s35-newest.gz'):

# putative_s35_url = 'http://download.cathdb.info/cath/releases/daily-release/newest/cath-b-s35-newest.gz'
# url = 'http://download.cathdb.info/cath/releases/daily-release/newest/cath-b-s35-newest.gz'
    request = urllib2.Request( url)
    response = urllib2.urlopen(request)
    request.add_header('Accept-encoding', 'gzip')
    request.add_header('Accept-encoding', 'gz')
    response = urllib2.urlopen(request)
    if response.info().get('content-type') == 'application/x-gzip':
#     if 1:
        buf = StringIO.StringIO(response.read())
        f = gzip.GzipFile(fileobj=buf)
        data = f.read()# len(f)
    else:
        data = response.read()
    # type(f)
    # help(f)
    # lines = data.splitlines()
    mapdict = {}
    return data

import sys
import time
class counter():
    def __init__(self, lst,  per = 100, fper = 1, INF = False, stdout = sys.stdout,
    	ifprint = 1,
        threshold = 0.1,
        step = 1,
        behave = '',
        prefix = '' ):
        # import time
        # self.lst = list(lst);
        # self.imax= len(lst)
        if INF:
            self.imax = -1
        else:
            self.imax= sum( 1 for _ in lst)
        self.i   = 0
        self.f   = 0
        self.per = per
        self.fper = fper
        self.flst = []
        self.e = None
        self.stdout = stdout
        self.ifprint = ifprint
        self.prefix = prefix 
        self.t0 = time.time()
        self.threshold = threshold 
        self.step = step
        self.behave = behave


    def count(self, step = None, callback = None):
        if not self.i % self.per and self.ifprint:
            msg = '%d of %d'%(self.i,self.imax)
            msg = self.prefix + msg
            print >> sys.__stdout__, msg
            print >> self.stdout, msg
        # if not step:
        #     step = self.step
        self.i += (step or self.step)
        
    def fail(self, msg, ins = None, e = None):
     #    if not self.f % self.fper:
        if msg:
            msg = self.prefix + msg
            print >> sys.__stdout__, msg
            print >> self.stdout, msg
        self.f += 1
        self.flst += [ins]
        try:
            self.e = e
        except:
            pass
    def summary(self, behave = '', INPUT = ''):
        behave = behave or self.behave
        self.imax = self.i
        self.t1 = time.time()
        self.dur = self.t1 - self.t0
        try:
            self.frate = self.f / float( self.i )
        except:
            self.frate = 0.
        print >> self.stdout, '\n[SUMMARY]:%s' % self.prefix
        print >> self.stdout, '\t[Task]: finshed %s from %s ' % (behave, INPUT)
        print >> self.stdout, '\t[Failrate]%d instances of %d failed, ( %2.2f%% )' % (self.f, self.i, 100 * self.frate )
        print >> self.stdout, '\t[Duration]:Ended after %.4f sec' % ( self.dur )   # len(lst)d
        if self.frate > self.threshold:
            msg = 'fail rate is too high: expected < 10%%, actual: %2.2f%%' % (100 * self.frate)
            raise Exception( msg )

##### Parallel processing
import multiprocessing as mp
from multiprocessing.managers import BaseManager, SyncManager
class MyManager(SyncManager): pass
MyManager.register('counter',counter)


from tempfile import TemporaryFile , mkdtemp

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def split_file(fname, linecount = None, number = None):
    if linecount and number:
        raise Exception('You can only specify one of linecount or number')
    f_handles = []
    tempdir = mkdtemp(prefix = '/tmp/feng')
    
    with open(fname,'r') as f:
        lines = f.readlines()
    lcount = len(lines)
       
    if number:
#         with open(fname) as f:
#             lines = f.readlines()
        lcount = len(lines)
        nlcount = lcount // number + 1
        idxs = range(0, lcount, nlcount)
        for i,idx in enumerate(idxs):
#             f = TemporaryFile()
            temp_fname = tempdir + '/%d'%i
            f = open(temp_fname, 'w+')
            f.write( ''.join(lines[idx: idx+nlcount]) )
            f.seek(0)
            f_handles.append(temp_fname)
    if linecount:
        idx = 0
        i   = 0
        nlcount = linecount
        while 1:
            temp_fname = tempdir + '/%d'%i
            f = open(temp_fname, 'w+' )
            f.write( ''.join(lines[idx:idx + nlcount]) )
            f.seek(0)
            f_handles.append(temp_fname)
            if idx + nlcount+1 >= lcount: 
                break
            else:
                idx += nlcount
                i   += 1
    return f_handles
if __name__ == '__test__':
    ##### !!! "test.fa" is yet to be deposited !!! #####
    print "[TESTING]:split_file()"
    INPUT=full('$SEQlib/test.fa') 
    lc1 = len(open(INPUT,'r').read())
    f_handles = split_file( INPUT, linecount = 1000 )
    lc2 = len(''.join( open(f,'r').read() for f in f_handles))
    f_handles = split_file( INPUT, number = 10 )
    lc3 = len(''.join( open(f,'r').read() for f in f_handles))
    assert  lc3 == lc2 == lc1
    print "[PASSED]:split_file()"




import csv

def csv_listener( q, fname):
    '''listens for messages on the q, writes to file. '''
    ## Read "fname" from global. "fname" file must exists prior to call.

    f = open(fname, 'a') 
    c = csv.writer(f)
    while 1:
        m = q.get()
        if m == 'kill':
            # f.write('killed \n')
            break
        elif m == 'clear':
            f.truncate();    
        else:## Asumme a csv row
            row = m
            c.writerow( row )
        f.flush()
    f.close()



def worker( i, q, slist):
    # global wait, waitname
    # pdbfile = onlyfiles[i];
    (wait,waitname,reset,fname,env,args) = slist

    pdbfile = i;
    import os 

    pdbname = os.path.basename(pdbfile);
    if pdbname.split(".")[-1] in ["bak","csv"]:
        return 
    if wait:
        nDOPE = 0
        if pdbname == waitname:
            q.put('start');
            nDOPE = get_nDOPE( pdbfile, env = env, auto_complete = args.auto_complete)
        row = [pdbname, nDOPE ];

    else:    
        # print("\n\n//Testing structure from %s" % pdbfile)
        nDOPE = get_nDOPE( pdbfile, env = env, auto_complete = args.auto_complete)                
        row = [pdbname, nDOPE] ;
    q.put( row );
    
    return row



### Statistics helper functions

def MAD_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    if not isinstance(points,np.ndarray):
        points = np.array(points)
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh


def cov2corr(c):
    c = np.copy(c)
    try:
        d = np.diag(c)
    except ValueError:
        # scalar covariance
        # nan if incorrect value (nan, inf, 0), 1 otherwise
        return c / c
    stddev = np.sqrt(d.real)
    c /= stddev[:, None]
    c /= stddev[None, :]

    # Clip real and imaginary parts to [-1, 1].  This does not guarantee
    # abs(a[i,j]) <= 1 for complex arrays, but is the best we can do without
    # excessive work.
    np.clip(c.real, -1, 1, out=c.real)
    if np.iscomplexobj(c):
        np.clip(c.imag, -1, 1, out=c.imag)

    return c







" Here starts ISS utils"
from Bio import SearchIO
#### import hmmpf database
# print("finished")
hmmpf_file = full("$SEQlib/cath-S35-hmm3.lib")
# hmmpf_file = full("$SEQlib/tmp.hmmpf")

def search_idx(q_acc, fname, acc_func):
    # fname = f
    with open( full(fname) , 'r') as f:
        line = f.readline()
        idx = 0
        while line:
#             header = line.rstrip("\n")[6:]
#             acc = seqheader_parse_cath(header)["acc"]
            acc = acc_func(line)
            if acc == q_acc:
                break
            
            line = f.readline()
            idx += 1
        if not line:
            return None
        else:
            return idx
        
        

        
import subprocess
def index_hmmlib(hmmlib_file):
    idxfile = hmmlib_file + '.idx'
    if not os.path.isfile(idxfile):
        subprocess.call(['cat',hmmlib_file,'|','grep','NAME','>',idxfile])
    f = open(idxfile,'r')
    lines = f.readlines()
    f.close() 
    return lines


def search_lib( q_acc, hmmlib_file = hmmpf_file, acc_func = None):
    
    ##### check whether an index file exists *.idx
    idxfile = hmmlib_file + '.idx'
    if not os.path.isfile(idxfile):
        subprocess.call(['cat',hmmlib_file,'|','grep','NAME','>',idxfile])
    if not acc_func:
        acc_func = lambda line: seqheader_parse_cath(line.rstrip("\n")[6:])["acc"]
    
    idx = search_idx(q_acc ,idxfile, acc_func)
    print idx

    if not idx:
        return None
    else:
        with open( full(hmmlib_file) ,'r') as f:
            line = 1
            cidx = 0
#             ln = 0
            buf = ""
            while line:
                line = f.readline()                
                if cidx == idx:
                    buf += line
                if line.rstrip('\n')=='//':
                    cidx += 1
                    if buf:
                        return buf 



#### import sequence database
def seqheader_parse_cath(header):
    lst = header.split("|")
#     dbname = lst[0]
    dbname = "CATH"
    version = 'v' + lst[1]
    acc = lst[2].split("/")[0] 
    jdict = {
        "seqDB":{
            "name":   dbname,
            "version":version,        
            },
        "acc": acc,
        }
    return jdict

def parse_hmmlib(fname):
    with open(fname, "r") as f:
        buf = ''
        while 1:
            line = f.readline()
            buf += line
            if line == '//\n':
                yield buf
                buf = ''
            if not line:
                break

# import gc
def hsp2jdict(hsp,query = None, simple = False):
    jdict = hsp.__dict__
    if not simple:
	    jdict["query_id"] = query.id
	#     print hsp.hit_id
	    jdict["target_id"] = hsp.hit_id
    # jdict["target_id"] = sequence.objects.get( acc = hsp.hit_id).id
    it = jdict.pop('_items')
    jdict.pop('domain_index')
    jdict["start"]=jdict.pop("env_start")
    jdict["end"]=jdict.pop("env_end")
    jdict["logCevalue"] = max(-1000,np.log10(jdict.pop("evalue_cond")))
    jdict["logIevalue"] = max(-1000,np.log10(jdict.pop("evalue")))
    
    return jdict



def hmmsearch(hmm,seqDB_curr = None, seqDB_file = None, tmpdir = "/tmp", 
	tbl_format = 'hmmsearch3-domtab'):
#     if ["seqDB_curr"] in locals().keys():
    if not seqDB_curr:
        seqDB_curr = seqDB.objects.get(name = 'CATH')
#         seqDB_file = seqDB_curr.filepath
    if seqDB_file:
        seqDB_file = full(seqDB_file)
    

    tmphmm = "%s/hmm" % tmpdir
    with open(tmphmm,'w') as f:
        f.write(hmm.text)
    cmd = ["hmmsearch","--noali","-o",tmpdir+"/log","--domtblout",tmpdir + "/domtbl", 
                                     tmphmm, seqDB_file,]
    # print ' '.join(cmd)
    try:
        qtext = subprocess.check_output( cmd )
    except Exception as e:
    	print 'CMD is: %s' % (' '.join(cmd))
        raise e 
    # del qtext
    # gc.collect()

    parser = SearchIO.parse( tmpdir + "/domtbl", tbl_format )
    
    q_hits = next(parser,None)
    if q_hits:
	    q_hits.id = seqheader_parse_cath(q_hits.id)["acc"]
	    def acc_mapper(hit):
	#         hit.query
	        hit.id = seqheader_parse_cath(hit.id)["acc"]
	        return hit
	    q_hits = q_hits.hit_map(acc_mapper)
	    return q_hits
    
	# oldhits = hmm.hit4hmm2hsp_set.all()
	# for hit in q_hits:
	#     hsp = hit[0] ### Assume only one dom per hit
	#     jdict = hsp2jdict(hsp, query = hmm )
	#     jdict["target_id"] = sequence.objects.get( acc = jdict["target_id"]).id
	#     if not oldhits.filter(**jdict).exists():
	#         hit_db = hit4hmm2hsp(**jdict)
	#         hit_db.save()

#     dis_q = hmm.hit4hmm2hsp_set.values_list("target",flat = True).distinct()
    
#     q = hmm.hit4hmm2hsp_set.all()

#     if dis_q.count() < q.count():
#         q = q.exclude(id__in=list(dis_q.values_list("id",flat = True)) )
# #         q.delete()
# #         dis_q.update()
        
    # return hmm.hit4hmm2hsp_set

def batch_qs(qs, batch_size=1000):
    """
    Returns a (start, end, total, queryset) tuple for each batch in the given
    queryset.
    
    Usage:
        # Make sure to order your querset
        article_qs = Article.objects.order_by('id')
        for start, end, total, qs in batch_qs(article_qs):
            print "Now processing %s - %s of %s" % (start + 1, end, total)
            for article in qs:
                print article.body
    """
    total = sum( 1 for _ in qs)

    for start in range(0, total, batch_size):
        end = min(start + batch_size, total)
        # yield (start, end, total, qs[start:end])
        yield qs[start:end]
def acc_mapper_cath(hit):
#         hit.query
    hit.id = seqheader_parse_cath(hit.id)["acc"]
    return hit    



def ISS_normalise(hc1, hc2, hcboth):
    # if len(hc1) != 1 or isinstance(hc1,np.array):
    # else:
    log1 =  np.log10( hc1 + 1) 
    log2 =  np.log10( hc2 + 1) 
    log3 =  np.log10( hcboth + 1) 
    norm = (log1 + log2) / 2 - log3 
    return norm

def ISS_normalise_new(hc1, hc2, hcboth, geoavg = None):
    # if len(hc1) != 1 or isinstance(hc1,np.array):
    # else:
    log3 =  np.log10( hcboth + 1) 
    if geoavg:
        log4 = np.log10( geoavg + 1)
        norm = log3 / log4
    else:
        log1 =  np.log10( hc1 + 1) 
        log2 =  np.log10( hc2 + 1) 
        norm = 2 * log3 / (log1 + log2)  
    return norm


##### iterating over COO sparse matrix (10 times faster than DOK matrix)

import random
import itertools
from itertools import izip

# def using_nonzero(x):
#     rows,cols = x.nonzero()
#     for row,col in zip(rows,cols):
#         ((row,col), x[row,col])

# def using_coo(x):
#     cx = scipy.sparse.coo_matrix(x)    
#     for i,j,v in zip(cx.row, cx.col, cx.data):
#         (i,j,v)

def using_tocoo(x):
    cx = x.tocoo()    
    it = zip(cx.row, cx.col, cx.data)
    return it
    # for i,j,v in zip(cx.row, cx.col, cx.data):
    #     (i,j,v)

def using_tocoo_izip(x):
    cx = x.tocoo() 
    it = itertools.izip(cx.row, cx.col, cx.data)
    return it   
    # for i,j,v in itertools.izip(cx.row, cx.col, cx.data):
    #     yield (i,j,v)

def sort_coo(m,order = 1):
    tuples = izip(m.row, m.col, m.data)
    return sorted(tuples, key=lambda x: order * (x[2]) )



####################################################
####### Toy function to expand hmm into lists  #####
####################################################
hmms2hit_ids = lambda hmms: np.expand_dims(
    np.array(
    [list(
        hmm.hits.values_list("id", flat = True)
    )
     for hmm in hmms]
), axis = 1)

hmms2hit_ids_para = lambda hmms,pool: np.expand_dims(
    np.array(
    [list(
        hmm.hits.values_list("id", flat = True)
    )
     for hmm in hmms]
), axis = 1)


def wrap_mat(Da,Db):
    from scipy.sparse import dok_matrix
    l  = Da.shape[0]
    l2 = Db.shape[0]
    assert l == l2,'shape of two matrix don''t match, Da.shape: %s, Db.shape: %s ' % (l,l2)
    OUTPUT = dok_matrix( (l,l))
    it = using_tocoo_izip( Db )
    D_curr = Da.todok()

    d = dict()
    c = counter([],INF = 1,per = 10000)
    for x,y,v2 in it:
        c.count()
        v1 = D_curr[x,y]
        d[(x,y)] = [v1,v2]
    OUTPUT.update(d)
    c.summary("zipping sparse matrix")
    return OUTPUT




TEMPLATE_STRING_IF_INVALID = 'No attr:'

numeric_test = re.compile("^\d+$")
def getattribute_none(value, arg):
    """Gets an attribute of an object dynamically from a string name"""
    if hasattr(value, str(arg)):
        attr = getattr(value, str(arg));
        if callable(attr):
            return attr();
        else:
            return attr
    elif hasattr(value, 'has_key') and value.has_key(arg):
        return value[arg]
    elif numeric_test.match(str(arg)) and len(value) > int(arg):
        return value[int(arg)]
    else:
        return 'None'

def getattribute(value, arg, debug = False):
    """Gets an attribute of an object dynamically from a string name"""
    if arg.endswith("?"):
        arg = arg[:-1]
        callit = 1
    else:
        callit = 0

    if hasattr(value, str(arg)):
        attr = getattr(value, str(arg));
        # tp = type(attr) 
        # if callable(attr) and tp != "django.db.models.manager.Manager":
        if callit:
            return attr();
        else:
            return attr
    elif hasattr(value, 'has_key') and value.has_key(arg):
        return value[arg]
    elif numeric_test.match(str(arg)) and len(value) > int(arg):
        return value[int(arg)]
    else:
        return TEMPLATE_STRING_IF_INVALID + ' ' + arg

def getattribute_iter(value, args, debug = False,):
    args = args.split('__');
    # value 
    for arg in args:
        if arg:
            value = getattribute(value,arg, debug = debug);
        if debug:
            raise Exception(value)
        else:
            pass

    return value


def assert_expr(act, expr = '', behave = None, **kwargs):
    #### remember to pass your data_dict in kwargs, especially what got mentioned in "expr" ###########
    #######################################################################
    ####  example: expr = " == exp "  ==>  fullexpr = "act == exp"
    ####           "exp" needs to be passed as a keyword argument assert_expr(... , exp = "12345") 
    ###################################################################################
#     expr0 = expr[:]
    if not expr:
        expr = '== True '
        
    fullexpr = ' act ' + expr
    basemsg = "Expecting:%s, Actual:%s " % (expr,act)
    success = eval(fullexpr,locals(),kwargs)
    assert( success ),"[FAILED]:" + basemsg
    print "[PASSED]:" + basemsg
    
import itertools
def sortANDgroup(lst,key = None):
    lst = sorted(lst,key = keyF)
    return itertools.groupby(lst,key = keyF)



import collections
from scipy.sparse import *
def concat_dok( Ds, func = sum, args = [],):
    from scipy.sparse import dok_matrix
    if func == max:
        args = []
    elif func == sum:
        args = [()]
    
    same = 1
    s0 = Ds[0].shape
    keys = set()
    for D in Ds:        
        same *= D.shape == s0
        assert same,'size checking failed'
        keys = keys | set(D.keys())
    
        
    OUTPUT = dok_matrix( Ds[0].shape )
    it = keys

    d = dict()
    c = counter([],INF = 1,per = 10000)
    for k in it:
#         v = sum( (D.get(k,(0,)) for D in Ds),() )
        v = func( (D.get(k,(0,)) for D in Ds), *args )
        c.count()
#         v1 = D_curr[x,y]
        d[k] = v
    OUTPUT.update(d)
    c.summary("concatenating sparse matrix")
    return OUTPUT

def matrify( lst, flat = False, l = None):
    if not l:
        l = max(itertools.izip(*lst))
        if isinstance(l,tuple):
            l = l[0]
#     l = self.l
    OUTPUT = dok_matrix( (l + 1, l + 1), dtype = 'int')
    it = (tuple(x) for x in lst)
    count = collections.Counter(it)
    if not flat:
        for k,v in count.iteritems():
            count[k] = (v,)
    OUTPUT.update(count)
    return OUTPUT
