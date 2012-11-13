import urllib2
import math # voluntarily avoiding numpy...
import numpy as np # ... did not work in the end ;-)
import cPickle # to keep track of targets
import os

"""
example:

import simbad

# get the info from simbad:
s = simbad.query(['Polaris', 'delta Cep', 'eta Aql'])

# pretty print for org mode
simbad.prettyPrint(s)
"""

default_site = 'http://simbad.u-strasbg.fr/' # default site is Strasbourg
alternate_site = 'http://simbad.cfa.harvard.edu/' # alternate is Harvard
simbad_site = default_site
dbfile = 'simbad.dpy'

def query(identifiers, debug=False, closest=False, around=0):
    """
    identifiers is list of SIMBAD identifier

    if around is set, search around target for the specified value (in minutes)
    
    returns a list of dictionnaries with parameters (self
    explanatory). DIAM is the angular diameter in mas, estimated from
    photometry. TRANSIT MONTH is the moonth for which the objects
    transit at midnight.

    Note:
    - if magnitudes are not present, -99. is returned 
    - if VSINI (km/s) not found, set to -1
    - if PM(A,D) or PLX not found, set to 0
    - IRAS FLUXES are in Jy, for 12, 25, 60, 100
    """
    global simbad_site, dbfile
    
    if not isinstance(identifiers, list): # make a scalar into a list
        identifiers = [identifiers]
    if closest and len(identifiers)>1:
        print 'ID:', identifiers, len(identifiers)
        print 'closest=True only for single query...'
        return 
        
    ngroup = 40 # max number of objects to query at the same time
    if len(identifiers)>ngroup: # slice list
        # group by ngroup
        res = []
        while len(identifiers)>0:
            #print len(identifiers), 'objects to query'
            res.extend(query(identifiers[:ngroup], debug=debug, closest=closest,
                             around=around))
            identifiers = identifiers[ngroup:]
        return res
    
    ###########################################################
    ### from here, assumes identifiers is a list of strings ###    ###########################################################
     
    # check if it is in the DataBase
    if os.path.isfile(dbfile) and around==0:
        dbf = open(dbfile)
        db = cPickle.load(dbf)
        dbf.close()
        res = []
        for i in identifiers:
            if db.has_key(i):
                #print i, 'found in local DB'
                res.append(db[i])
            else:
                #print i, 'NOT found in local DB'
                res.append({})
    else:
        res = [{} for i in identifiers]
            
    # -- all target found in database
    if all([r.has_key('IDENTIFIER') for r in res]):
        return res
    
    rt_ = '%0D%0A' # cariage return
    plus_ = '%2B' # + in the URL
    separator = ';'
    format_ = "format+object+form1+\""+separator+"+%25IDLIST(1)+"+separator+"+%25COO(A+D)+"+separator+"+%25OTYPE+"+separator+"+%25SP+"+separator+"+%25PM(A+D)+"+separator+"+%25PLX(V+E)+"+separator+"+%25FLUXLIST(B)+"+separator+"+%25FLUXLIST(V)+"+separator+"+%25FLUXLIST(R)+"+separator+"+%25FLUXLIST(J)+"+separator+"+%25FLUXLIST(H)+"+separator+"+%25FLUXLIST(K)+"+separator+"+%25MEASLIST(rot;|F)+"+separator+"%25MEASLIST(iras;|F)"+separator+"%25MEASLIST(JP11;|F)"+"\""

    url = 'simbad/sim-script?submit=submit+script&script='+format_
    
    Nquery = 0
    IDquery = []
    for k,i in enumerate(identifiers):
        if not res[k].has_key('IDENTIFIER'):
            Nquery+=1
            IDquery.append(i)
            obj = i.replace('+', plus_)
            obj = obj.replace('_', ' ')
            obj = obj.replace(' ', '+')
            if ':' in i: # these must be coordinates!
                url = url+rt_+'query+coo+'+obj+'+radius%3D5s'
            elif around>0:
                url = url+rt_+'query+around+'+obj+'+radius%3D'+str(around)+'m'
            else:
                url = url+rt_+'query+id+'+obj

    if debug:
        print simbad_site+url
    try:
        lines = urllib2.urlopen(simbad_site+url, timeout=20).read()
    except:
        simbad_site = alternate_site
        print 'switching to alternate server...'
        try:
            lines = urllib2.urlopen(simbad_site+url, timeout=20).read()
        except:
            raise NameError('servers do not respond OR no internet connection')
               
    if debug:
        print lines
    lines = lines.split('\n')

    # go to data
    for k, l in enumerate(lines):        
        if ':error:' in l:
            #print '  ERROR:', lines[k+2]
            #print '------------------------------'
            #print lines
            return None
        if ':data:' in l:
            lines = lines[k+1:]
            break

    lines = filter(lambda x: len(x)>0, lines)

    if len(lines)!=Nquery and not closest and around==0:
        print '  ERROR: too many/few results!'
        return None
    
    if debug:
        print lines

    # read every line which is a different object 
    for k, l in enumerate(lines):
        obj = {}
        if around>0:
            obj['IDENTIFIER'] = 'around: '+identifiers[0]
        else:
            obj['IDENTIFIER'] = IDquery[k]
        
        obj['NAME'] = l.split(separator)[1].strip()        

        if '-' in l.split(separator)[2]:
            l_ra = l.split(separator)[2].split('-')[0]
            l_dec = '-'+l.split(separator)[2].split('-')[1]
        else:
            l_ra = l.split(separator)[2].split('+')[0]
            l_dec = '+'+l.split(separator)[2].split('+')[1]

        obj['RA'] = l_ra.strip()
        obj['DEC'] = l_dec.strip()
        
        if len(l_ra.split())==3:
            obj['RA.h'] = (float(l_ra.split()[0])+
                           float(l_ra.split()[1])/60.+
                           float(l_ra.split()[2])/3600.)
        elif len(l_ra.split())==2:
            obj['RA.h'] = (float(l_ra.split()[0])+
                           float(l_ra.split()[1])/60.)
        else:
            obj['RA.h'] = float(l_ra.split()[0])
                          
        obj['RA D'] = obj['RA.h']*15
        
        if len(l_dec.split())==3:
            obj['DEC.d'] = abs(float(l_dec.split()[0]))+\
                           float(l_dec.split()[1])/60.+\
                           float(l_dec.split()[2])/3600.
        elif len(l_dec.split())==2:
            obj['DEC.d'] = abs(float(l_dec.split()[0]))+\
                           float(l_dec.split()[1])/60.
        else:
            obj['DEC.d'] = abs(float(l_dec.split()[0]))
            
        obj['DEC.d'] = math.copysign(obj['DEC.d'],
                                     float(l_dec.split()[0]))
            
        # 15th Jan at midnight is ~ LST 6:00
        obj['TRANSIT MONTH'] = int(round((obj['RA.h']-6.00)/2.-1, 0))%12+1
        obj['TYPE'] = l.split(separator)[3].split('~')[0].strip()
        obj['SPTYPE'] = l.split(separator)[4].strip().split()[0]

        try:
            obj['PMA'] = float(l.split(separator)[5].split()[0])/1000.
            obj['PMD'] = float(l.split(separator)[5].split()[1])/1000.
        except:
            obj['PMA'] = 0.0
            obj['PMD'] = 0.0

        try:
            obj['PLX'] = float(l.split(separator)[6].split()[0])/1000.
            obj['EPLX'] = float(l.split(separator)[6].split()[1])/1000.
        except:
            obj['PLX'] = 0.0
            obj['EPLX'] = 0.0

        mags = ['B','V','R','J','H','K']
        for j, m in enumerate(mags):
            try:
                obj[m+'MAG'] = float(l.split(separator)[7+j].split()[1])
            except:
                try:
                    # take first number
                    tmp = l.split(separator)[7+j]
                    for i in range(len(tmp)):
                        if tmp[i].isdigit():                            
                            break
                    obj[m+'MAG'] = float(tmp[i:].split()[0])
                except:
                    obj[m+'MAG'] = np.nan
        try:
            obj['VSINI'] = float(l.split(separator)[13].split('|')[0].split()[0])
        except:
            obj['VSINI'] = -1 # failed
        iras_wl = ['12um', '24um', '60um', '100um']
        
        obj['IRAS'] = dict(zip(iras_wl, np.zeros(len(iras_wl)))) 
        for i,j in enumerate(iras_wl):
            try:
                obj['IRAS'][j] = float(l.split(separator)[14].split('|')[i].split()[0])
            except:
                obj['IRAS'][j] = np.nan
                
        JP11_wl = ['U', 'B', 'V', 'R', 'I', 'J', 'K', 'L', 'M', 'N', 'H']
        obj['JP11'] = dict(zip(JP11_wl, np.zeros(len(JP11_wl)))) 
        for i,j in enumerate(JP11_wl):
            try:
                obj['JP11'][j] = float(l.split(separator)[15].split('|')[i].split()[0])
            except:
                obj['JP11'][j] = np.nan
        if np.isnan(obj['KMAG']) and not np.isnan(obj['JP11']['K']):
            obj['KMAG']= obj['JP11']['K']
        
         
        res[identifiers.index(IDquery[k])] = obj
        if closest:
            break
        
    if around>0:
        for k in range(len(res)):
            res[k]['DIST D'] = math.sqrt( (res[0]['DEC.d']-res[k]['DEC.d'])**2+
                                          math.cos(res[0]['DEC.d']*3.1415/180)**2*
                                          (res[0]['RA D']-res[k]['RA D'])**2)
            res[k]['DIST S'] =  res[k]['DIST D']*3600              
    res = addApproxDiam(res, verbose=False)
    
    if around==0:
        try:
            if not isinstance(db, dict):
                db = {}
        except:
            db = {}
            
        for k,i in enumerate(IDquery):
            db[i]= res[k]
            
        dbf = open(dbfile, 'w')
        cPickle.dump(db, dbf, 2)
        dbf.close()
    return res

def prettyPrint(dics, keys=['NAME', 'RA', 'DEC', 'VMAG', 'KMAG']):
    """
    dics is the result of query. optimized for .org mode of emacs
    """
    title='|'
    for k in keys:
        title = title+k+'|'
    print title
    print '|-'
    for d in dics:
        line = '|'
        for k in keys:
            line = line+str(d[k])+'|'
        print line
    return

def addApproxDiam(dics, verbose=True):
    """
    add the approximated diameter estimated using V-K

    uses 'BMAG', 'VMAG', 'JMAG', 'HMAG', 'KMAG' keyword of the
    dictionnary
    """
    # surface brightness relations for dwarf stars
    # from Kervella et al. 2004
    k04 = {}
    #           coef0  coef1  error
    k04['BV']=[.9095, .4889, .0918]
    k04['BJ']=[.3029, .5216, .0307]
    k04['BH']=[.2630, .5134, .0189]
    k04['BK']=[.2538, .5158, .0100]
    k04['VJ']=[.3547, .5310, .0475]
    k04['VH']=[.2893, .5148, .0185]
    k04['VK']=[.2753, .5175, .0101]
    k04['JK']=[.5256, .5097, .0575]

    for k, d in enumerate(dics): # for each star
        diams = []
        errs = []
        for coul in k04.keys(): # for each color
            # check magnitudes are valid, compute diameter and error
            if d.has_key(coul[0]+'MAG') and d[coul[0]+'MAG']>-90 and\
            d.has_key(coul[1]+'MAG') and d[coul[1]+'MAG']>-90:
                diams.append(diamSurfBri(d[coul[0]+'MAG'], d[coul[1]+'MAG'],
                                         k04[coul]))
                errs.append(k04[coul][2]*diams[-1])
        if len(diams)>1:
            # weighted average\
            dics[k]['DIAM'] = reduce(lambda x,y: x+y, [diams[i]*errs[i]
                                                       for i in range(len(diams))])/\
                                                       reduce(lambda x,y: x+y, errs)
            dics[k]['DIAM'] = round(dics[k]['DIAM'],
                                    int(-math.log10(dics[k]['DIAM']) +3))
        elif len(diams)==1:
            dics[k]['DIAM'] = round(diams[0], int(-math.log10(diams[0])+3))
        else:
            dics[k]['DIAM'] = 0            
        if verbose: 
            print dics[k]['NAME'], '|', dics[k]['DIAM']
    return dics

def diamSurfBri(c0, c1, coef):
    """
    surface brightness formula for magnitudes c0 and c1 (color c1-c0)
    
    see Kervella & Fouque 2008
    """
    return 10**(coef[1] + coef[0]*(c0-c1) - 0.2*c0)

def createEdb(s):
    """
    convert simbad dic into a xEphem format (list of stringes)
    
    exple: Polaris,f|M|F7,2:31:48.704,89:15:50.72,2.02,2000
    """
    if isinstance(s, list):
        return [createEdb(x) for x in s]
    
    res = "%s,f|S|%s,%s,%s,%3.1f,2000" % (s['IDENTIFIER'],
                                          s['SPTYPE'][:4],
                                          s['RA'].strip().replace(" ", ":"),
                                          s['DEC'].strip().replace(" ", ":"),
                                          s['VMAG'])
    return res

def createRdb(s):
    """
    create RDB line for APES input catalog. Fields are tab separated, missing
    values are replaced by ~.
    
    stype:  'T' or 'R'
    system_id: same for 2 stars
    star_id
    alpha, delta in dd/hh:mm:ss.ss
    dalpha, ddelta: ?
    epoch: J2000
    equinox: 2000
    coord_syst: ircs
    mualpha, mudelta: propermotion, in mas/yr
    Tint_max: ?
    stdev_Dphi: ?
    parallax: in mas
    SP_type: spectral type
    Teff: in K
    lambda_eff: in microns?
    magV, magK, magH:
    MS_tar, r, MP_1, T0_1, period_1, ecc_1, a_1, inc_1, omega_1, OMEGA_1,
        MP_2, T0_2, period_2, ecc_2, a_2, inc_2, omega_2, OMEGA_2, MP_3, T0_3,
        period_3, ecc_3, a_3, inc_3, omega_3, OMEGA_3: parameters of astrometric
        signal
    
    """
    pass
    return
