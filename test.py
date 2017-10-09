from ..utils import *






#############################
####   !!!!!!!!!!!!!   ######
####     DEPRECATED    ######
#######!!!!!!!!!!!!!!#######
#############################
def test__raw(D_raw, hmms):
    # next(counts.iteritems()[0]
    # it = counts.iteritems()
    it = using_tocoo_izip(D_raw)
    c = counter([],INF = 1)
    for x,y,v in it:
#         print x,y
        # hmm1 = hmms.get(id = x + 1)
        # hmm2 = hmms.get(id = y + 1)
        hmm1 = hmms.get(id = x )
        hmm2 = hmms.get(id = y )
        hmm1hits = hmm1.hits.values_list('id')
        hmm2hits = hmm2.hits.values_list('id')
        interhits = set(hmm1hits) & set(hmm2hits)
        intercount = len(interhits)
#         print v
#         print intercount
        msg = '[OK] %s against %s overlaps %d, with %d from '%(hmm1,hmm2, intercount, v )
        print msg
        assert v == len(interhits),'[ERROR] %s against %s overlaps %d, while D_raw returns %d '%(hmm1,hmm2, intercount, v )
        c.count()
        if c.i == 5:
            break
            
def test__norm(D_curr, hmms, norm_func = ISS_normalise):
    it = using_tocoo_izip(D_curr)
    c = counter([],INF = 1)
    for x,y,v_act in it:
        # hmm1 = hmms.get(id = x + 1)
        # hmm2 = hmms.get(id = y + 1)
        hmm1 = hmms.get(id = x )
        hmm2 = hmms.get(id = y )
        hmm1hits = hmm1.hits.values_list('id')
        hmm2hits = hmm2.hits.values_list('id')
        interhits = set(hmm1hits) & set(hmm2hits)
        intercount = len(interhits)
        
        v_exp = norm_func( len(hmm1hits), len(hmm2hits), intercount)
#         print v
#         print intercount
        msg = '[OK] %s against %s overlaps:: Expected:%s, Actual:%s'%(hmm1,hmm2, v_exp, v_act )
        print msg
        assert v_exp == v_act,'[ERROR] %s against %s overlaps:: Expected %s, Actual: %s '%(hmm1,hmm2, v_exp, v_act )
        c.count()
        if c.i == 5:
            break


def test__key(D_curr, s_exp = 1.0):
    it = using_tocoo_izip(D_curr)
    length = sum( 1 for _,_,_ in it)
    it = using_tocoo_izip(D_curr)
    sortcount = sum( x < y for x,y,_ in it  )
    s_act = float(sortcount)/length
#     s_exp = 1.0
    msg = 'Sorted rate: [ x < y for index (x,y) ) ] is %2.2f %%, Expected: %2.2f %%' % ( s_act*100, s_exp*100 )
    assert  s_exp == s_act, '[ERROR]' + msg
    print "[OK]" + msg





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






def test__raw(D_raw, reverse_dict, letter = 'S', seqDB_curr = None,depth = 5,debug = 0,**kwargs):
    # next(counts.iteritems()[0]
    # it = counts.iteritems()
    it = using_tocoo_izip(D_raw)
    c = counter([],INF = 1)#
    
    sids = set(seqDB_curr.sequence_set.values_list('id',flat = True))
    rv_field = reverse_field[letter]
    nodes = classification.objects.filter(level__letter = letter)
    
    for x,y,v in it:
#         print x,y
        ###### New routine
        hmm1 = nodes.filter(id = reverse_dict[x])
        hmm2 = nodes.filter(id = reverse_dict[y])
        hmm1hits = hmm1.values_list( rv_field, flat = True)
        hmm2hits = hmm2.values_list( rv_field, flat = True)
        ###### Old routine
#         hmm1 = hmms.get(id = x + 1)
#         hmm2 = hmms.get(id = y + 1)
#         hmm1hits = hmm1.hits.values_list('id')
#         hmm2hits = hmm2.hits.values_list('id')
        hmm1hits = set(hmm1hits) & sids
        hmm2hits = set(hmm2hits) & sids
        interhits = set(hmm1hits) & set(hmm2hits)
    
#         interhits = set(hmm1hits) & set(hmm2hits)
        exp = len(interhits)
        act = v
        
#         print v
#         print intercount
        msg = '%s against %s overlaps: Expected: %d, Actual: %d from '%(hmm1[0],hmm2[0], exp, act )
        success = ( v == len(interhits) )
        errmsg = "[ERROR] " + msg

        if not debug:
            assert success, errmsg
            if depth:
                print "[OK] " + msg
        else:
            if success:
                print "[OK] " + msg
            else:
                print "[ERROR] " + msg

        c.count()
        if c.i == depth:
            break
# if 
# test__raw( OUTPUT, reverse_dict = reverse_dict, letter = 'H', seqDB_curr = sDB)
def test__norm(D_curr, reverse_dict,  letter = 'S', seqDB_curr = None, dnum = 2,
               norm_func = ISS_normalise):
    it = using_tocoo_izip(D_curr)
    c = counter([],INF = 1)

    sids = set(seqDB_curr.sequence_set.values_list('id'))
    rv_field = reverse_field[letter]
    # nodes = classification.objects.filter(level__letter = letter)
    nodes = classification.objects
    
    for x,y,v in it:
        # hmm1 = hmms.get(id = x + 1)
        # hmm2 = hmms.get(id = y + 1)
        node1__id = reverse_dict[x]
        node2__id = reverse_dict[y]
        hmm1 = nodes.filter(id = node1__id)
        hmm2 = nodes.filter(id = node2__id)
        hmm1hits = hmm1.values_list( rv_field )
        hmm2hits = hmm2.values_list( rv_field )
        hmm1hits = set(hmm1hits) & sids
        hmm2hits = set(hmm2hits) & sids
        interhits = set(hmm1hits) & set(hmm2hits)
    
        intercount = len(interhits)
        
        exp = norm_func( len(hmm1hits), len(hmm2hits), intercount)
        act = v
#         print v
#         print intercount
        try:
            msg = '%s against %s overlaps: Expected: %f, Actual: %f from '%(hmm1[0],hmm2[0], exp, act )
        except:
            raise Exception("failed hmm1:%s , hmm2:%s" %(node1__id,node2__id))

#         print msg
        
        success = ( round(act, dnum) == round(exp, dnum) )
        errmsg = "[ERROR] " + msg
#         assert success, errmsg
#         print "[OK] " + msg
        if success:
            print "[OK] " + msg
        else:
            print  errmsg
        c.count()
        if c.i == 5:
            break