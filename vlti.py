import numpy as np
import os
import simbad
from matplotlib import pyplot
import itertools

#      station   p        q        E        N       A0
layout_orientation = -18.984 # degrees
layout = {'A0':(-32.0010, -48.0130, -14.6416, -55.8116, 129.8495),
          'A1':(-32.0010, -64.0210, -9.4342, -70.9489, 150.8475),
          'B0':(-23.9910, -48.0190, -7.0653, -53.2116, 126.8355),
          'B1':(-23.9910, -64.0110, -1.8631, -68.3338, 142.8275),
          'B2':(-23.9910, -72.0110, 0.7394, -75.8987, 158.8455),
          'B3':(-23.9910, -80.0290, 3.3476, -83.4805, 166.8295),
          'B4':(-23.9910, -88.0130, 5.9449, -91.0303, 150.8275),
          'B5':(-23.9910, -96.0120, 8.5470, -98.5942, 174.8285),
          'C0':(-16.0020, -48.0130, 0.4872, -50.6071, 118.8405),
          'C1':(-16.0020, -64.0110, 5.6914, -65.7349, 142.8465),
          'C2':(-16.0020, -72.0190, 8.2964, -73.3074, 134.8385),
          'C3':(-16.0020, -80.0100, 10.8959, -80.8637, 150.8375),
          'D0':(0.0100, -48.0120, 15.6280, -45.3973, 97.8375),
          'D1':(0.0100, -80.0150, 26.0387, -75.6597, 134.8305),
          'D2':(0.0100, -96.0120, 31.2426, -90.7866, 150.8275),
          'E0':(16.0110, -48.0160, 30.7600, -40.1959, 81.8405),
          'G0':(32.0170, -48.0172, 45.8958, -34.9903, 65.8357),
          'G1':(32.0200, -112.0100, 66.7157, -95.5015, 129.8255),
          'G2':(31.9950, -24.0030, 38.0630, -12.2894, 73.0153),
          'H0':(64.0150, -48.0070, 76.1501, -24.5715, 58.1953),
          'I1':(72.0010, -87.9970, 96.7106, -59.7886, 111.1613),
          'J1':(88.0160, -71.9920, 106.6481, -39.4443, 111.1713),
          'J2':(88.0160, -96.0050, 114.4596, -62.1513, 135.1843),
          'J3':(88.0160, 7.9960, 80.6276, 36.1931, 124.4875),
          'J4':(88.0160, 23.9930, 75.4237, 51.3200, 140.4845),
          'J5':(88.0160, 47.9870, 67.6184, 74.0089, 164.4785),
          'J6':(88.0160, 71.9900, 59.8101, 96.7064, 188.4815),
          'K0':(96.0020, -48.0060, 106.3969, -14.1651, 90.1813),
          'L0':(104.0210, -47.9980, 113.9772, -11.5489, 103.1823),
          'M0':(112.0130, -48.0000, 121.5351, -8.9510, 111.1763),
          'U1':(-16.0000, -16.0000, -9.9249, -20.3346, 189.0572),
          'U2':(24.0000, 24.0000, 14.8873, 30.5019, 190.5572),
          'U3':(64.0000, 48.0000, 44.9044, 66.2087, 199.7447),
          'U4':(112.0000, 8.0000, 103.3058, 43.9989, 209.2302)}

# ---------- horizons ---------------
hfiles = os.listdir('VLTI_horizons')
hfiles = filter(lambda x: '.horizon' in x, hfiles)
horizon = {}
for h in hfiles:
    f = open(os.path.join('VLTI_horizons', h), 'r')
    lines = f.read().split('\n')[:-3]
    az = np.array([float(l.split()[0]) for l in lines])
    el = np.array([float(l.split()[1]) for l in lines])
    horizon[h[:2]] = (az, el)
    f.close()
        
vlti_longitude =  -70.40479659 # degrees
vlti_latitude  =  -24.62794830 # degrees

opdModels='./VLTI_opd'

def projBaseline(T1, T2, target, lst, ip1=1, ip2=3, min_altitude=20, max_OPD=95):
    """
    T1, T2: telescopes stations (e.g. 'A0', 'U1', etc.)
    target : [ra, dec] in decimal hour and decimal deg or name for SIMBAD
    lst   : decimal hour (scalar or array)
    ip1, ip2: input channels (optional)

    min_altitude (default 30 degrees) sets the minimum altitude for
    observation. Checks also for shadowinf from the UTs.

    max_OPD (default 95m) maximum OPD stroke of the delay lines.

    return a dictionnary. units are m and degrees 
    """
    if isinstance(target, str):
        s = simbad.query(target)[0]
        radec = [s['RA.h'], s['DEC.d']]
    else:
        radec = target
       
    # -- hour angle
    ha = (np.array(lst) - radec[0]) *360/24.0

    # -- alt-az
    dec = radec[1]
    d2r = lambda x: x*np.pi/180
    tmp1 =  np.sin(d2r(dec)) * np.sin(d2r(vlti_latitude)) +\
           np.cos(d2r(dec)) * np.cos(d2r(vlti_latitude)) *\
           np.cos(d2r(ha))
    alt = np.arcsin(tmp1)*180/np.pi
    
    tmp1 = np.cos(d2r(dec)) * np.sin(d2r(ha))
    tmp2 = -np.cos(d2r(vlti_latitude)) * np.sin(d2r(dec)) + \
           np.sin(d2r(vlti_latitude)) * np.cos(d2r(dec)) * \
           np.cos(d2r(ha));  
    az = (360-np.arctan2(tmp1, tmp2)*180/np.pi)%360

    b = [layout[T1][2]-layout[T2][2],
         layout[T1][3]-layout[T2][3],
         0.0] # assumes you *DO NOT* combine ATs and UTs
    
    # projected baseline
    ch_ = np.cos(ha*np.pi/180.)
    sh_ = np.sin(ha*np.pi/180.)
    cl_ = np.cos(vlti_latitude*np.pi/180.)
    sl_ = np.sin(vlti_latitude*np.pi/180.)
    cd_ = np.cos(radec[1]*np.pi/180.)
    sd_ = np.sin(radec[1]*np.pi/180.)

    # (u,v) coordinates in m
    u = ch_ * b[0] -  sl_*sh_ * b[1] + cl_*sh_ * b[2]
    v = sd_*sh_ * b[0] + (sl_*sd_*ch_+cl_*cd_) * b[1] -\
        (cl_*sd_*ch_ - sl_*cd_) * b[2] 
  
    # optical delay, in m
    d = -b[0]*cd_*sh_ -\
        b[1]*(sl_*cd_*ch_ - cl_*sd_) +\
        b[2]*(cl_*cd_*ch_ + sl_*sd_)
    d = np.array(d)
    
    opl1 = layout[T1][4] + 0.12*(ip1-1)
    opl2 = layout[T2][4] + 0.12*(ip2-1)
    opd = d+opl2-opl1

    observable = (alt>min_altitude)* \
                 (np.abs(opd) < max_OPD)*\
                 (alt>np.interp(az%360, horizon[T1][0], horizon[T1][1]))*\
                 (alt>np.interp(az%360, horizon[T2][0], horizon[T2][1]))

    if isinstance(alt, np.ndarray):
        observable = np.array(observable)

    airmass = alt*0+99
    czt = np.cos(np.pi*(90-alt)/180.)
    airmass = 1/czt*(1-.0012*(1/czt**2-1))

    if np.cos(d2r(dec))!=0:
        parang = 180*np.arcsin(np.sin(d2r(az))*
                               np.cos(d2r(vlti_latitude))/np.cos(d2r(dec)))/np.pi
    else:
        parang = 0.
    
    res = {'u':u, 'v':v, 'opd':opd, 'alt':alt, 'az':az,
           'observable':observable, 'parang':parang,
           'lst':lst, 'ra':radec[0],
           'dec':radec[1], 'B':np.sqrt(u**2+v**2),
           'PA':np.arctan2(u, v)*180/np.pi,
           'airmass':airmass, 'baseline':T1+T2,
           'opd':opd}
    
    return res

def nTelescopes(telescopes, target, lst, ip=None, plot=False,
                min_altitude=20, max_OPD=95):
    """
    target can be a string with the name of the target resolvable by
    simbad.
    """
    tmp = []
    if ip is None:
        ip = range(2*len(telescopes))[1::2]
    if isinstance(target, str):
        s = simbad.query(target)
        target_name=target
        target = [s[0]['RA.h'], s[0]['DEC.d']]
    else:
        # assume [ra dec]
        target_name = '%2d:%2d %3d:%2d' % (int(target[0]),
                                           int(60*(target[0]-int(target[0]))),
                                           int(target[1]),
                                           int(60*(abs(target[1])-int(abs(target[1]))))
                                           )
    res = {}
    for i in range(len(telescopes)+1):
        for j in range(len(telescopes))[i+1:]:
            tmp = projBaseline(telescopes[i],telescopes[j],
                               target, lst,ip1=ip[i], ip2=ip[j],
                               min_altitude=min_altitude, max_OPD=max_OPD)
            if not 'lst' in res.keys():
                # init
                for k in ['lst', 'airmass', 'observable',
                          'ra', 'dec', 'alt', 'az', 'parang']:
                    res[k] = tmp[k]
                res['baseline'] = [tmp['baseline']]
                for k in ['B', 'PA', 'u', 'v', 'opd']:
                    res[k] = {}
            else:
                # update
                res['observable'] = res['observable']*tmp['observable']
                res['baseline'].append(tmp['baseline'])
                
            for k in ['B', 'PA', 'u', 'v', 'opd']:
                res[k][tmp['baseline']] = tmp[k]
    
    if plot:
        where = lambda key: [res[key][k] for k in range(len(res[key])) 
                             if res['observable'][k]]  
        where2 = lambda key1, key2: [res[key1][key2][k] for k in range(len(res['lst']))
                                     if res['observable'][k]]
        pyplot.figure(0, figsize=(12,5))
        pyplot.clf()
        pyplot.subplot(121)
        pyplot.plot(res['lst'], res['alt'], '.', color='0.5')
        pyplot.plot(where('lst'), where('alt'), 'ko')
        pyplot.ylim(0,90)
        pyplot.xlabel('lst (h)')
        pyplot.ylabel('altitude (deg)')
        pyplot.title(target_name)
        ax = pyplot.subplot(122, polar=True)
        for b in res['baseline']:
            B = where2('B', b)
            B.extend(where2('B', b))
            PA = where2('PA', b)
            PA.extend([180+x for x in where2('PA', b)])
            pyplot.plot([x*3.1415/180 for x in PA],B, '.', label=b)
        pyplot.legend(ncol=1, prop={'size':8}, numpoints=1)
        ax.get_xaxis().set_visible(False)
    else:
        return res

def skyCoverage(T, uv_coverage=False, ha_coverage=True, Bmax=None, max_OPD=95,
                min_altitude=20):
    """
    plot the sky coverage of a VLTI configuration (2+ telescopes). By default,
    will display only the sky coverage (from 0 to 70 degrees in zenithal angle)
    and the drawing of the baselines on the platform.
    
    - T is a list of telescopes, such as ['A0', 'K0', 'G1']
    
    - uv_coverage= toggle the computation of the (u,v) coverage

    - Bmax= set the maximum baseline for the (u,v) plot. By default, it will be
        adjusted to the data, but in case you want to compare configurations, it
        is better to have the same range.
    
    - max_OPD=  maximum OPD correction for each baseline, in meters
        (default=False)
    
    LIMITATIONS:
    - there is no limitations in the number of telescopes, 2 is the minimum.
    - I did not implement the rules to avoid to combine telescopes in the same
        light duct (['A1', 'A0'] will work)
    - I did not implement VCM limitations
    - AT+UT is not accurate, even though it will plot something!
    - there is an option (max_OPD=) to give the maximum OPD compensation.
        Currently, I set it to 95m (it is really 100m, from 11m to 111m). You
        can play with the parameter if you want to see the effect on the double
        pass.
    """
    dec = np.linspace(-90, 80, 100)
    res_alt = []
    res_az = []
    for d in dec:
        lst = np.linspace(-12,12,1+101*np.cos(d*np.pi/180))
        x = nTelescopes(T, [0.0, d], lst, min_altitude = 20, max_OPD=max_OPD)
        w = np.where((x['alt']>20)*(1-np.float_(x['observable'])))
        res_alt.extend(x['alt'][w])
        res_az.extend(x['az'][w])
    res_alt = np.array(res_alt)
    res_az = np.array(res_az)

    if uv_coverage or ha_coverage:
        color_map = 'Spectral'        
        dec = np.linspace(-85,40,25)
        tracks = []
        for d in dec:
            lst = np.linspace(-12,12,97)[:-1] # one tick per 15 minutes
            x = nTelescopes(T, [0.0, d], lst,min_altitude =min_altitude, max_OPD=max_OPD)
            tracks.append(x)
 
    from matplotlib import pyplot
    
    if uv_coverage:
        pyplot.figure(1, figsize=(17,5))
    elif ha_coverage:
        pyplot.figure(1, figsize=(15,6))
    else:
        pyplot.figure(1, figsize=(9,6))    
    pyplot.clf()
    pyplot.subplots_adjust(wspace=0.05, top=0.9, left=0.02,right=0.98,bottom=0.1)

    # ------------------------------
    if uv_coverage:
        ax = pyplot.subplot(141,polar=1)
    elif ha_coverage:
        ax = pyplot.subplot(131,polar=1)
    else:
        ax = pyplot.subplot(121,polar=1)
        
    pyplot.title('Sky Coverage')
    pyplot.plot((res_az-90)*np.pi/180, 90-res_alt, 'xk', alpha=0.6)
    pyplot.text(np.pi/2, 60, 'N', size='x-large', color='k',
                horizontalalignment='center',
                verticalalignment='center',)
    pyplot.text(0.0, 60, 'E', size='x-large', color='k',
                horizontalalignment='center',
                verticalalignment='center')
    pyplot.text(-np.pi/2, 60, 'S', size='x-large', color='k',
                horizontalalignment='center',
                verticalalignment='center')
    pyplot.text(np.pi, 60, 'W', size='x-large', color='k',
                horizontalalignment='center',
                verticalalignment='center')
    if uv_coverage or ha_coverage:
        for tr in tracks:
            w = np.where(tr['observable'])
            pyplot.plot((tr['az'][w]-90)*np.pi/180, 90-tr['alt'][w],'o',
                color=pyplot.cm.get_cmap(color_map)((tr['dec']-min(dec))/np.ptp(dec)),
                linewidth=3, alpha=0.5)
    ax.get_xaxis().set_visible(False)
    # ------------------------------
    if uv_coverage:
        ax = pyplot.subplot(142)
    elif ha_coverage:
        ax = pyplot.subplot(132)
    else:
        ax = pyplot.subplot(122)

    ax.set_axis_off()
    pyplot.title(reduce(lambda x,y: x+'-'+y, T))
    ax.set_aspect('equal')
    # -- stations:
    #pyplot.plot([layout[k][2] for k in layout.keys() if not 'U' in k],
    #    [layout[k][3] for k in layout.keys() if not 'U' in k], 'w.',
    #    markersize=10)
    pyplot.plot([layout[k][2] for k in layout.keys() if 'U' in k],
        [layout[k][3] for k in layout.keys() if 'U' in k], 'o', color='y',
        markersize=20)
    # -- name of station
    for k in layout.keys():
        pyplot.text(layout[k][2], layout[k][3], k,
                    horizontalalignment='center',
                    verticalalignment='center',
                    color='k', fontsize=9, weight='black')
    # -- baselines:
    for i in range(len(T)+1):
        for j in range(len(T))[i+1:]:
            pyplot.plot([layout[T[i]][2], layout[T[j]][2]],
                [layout[T[i]][3], layout[T[j]][3]],'-', linewidth=3,
                alpha=0.3, color='k')
    pyplot.ylim(-120,120)
    pyplot.xlim(-40,160)
    
    # ------------------------------
    if not uv_coverage and not ha_coverage:
        return
    
    if ha_coverage:
        ax = pyplot.subplot(133)
        for tr in tracks:
            w = np.where(tr['observable'])
            pyplot.plot(tr['lst'][w], tr['dec']*np.ones(len(w[0])), 'o',
                        color=pyplot.cm.get_cmap(color_map)((tr['dec']-min(dec))/np.ptp(dec)))
        pyplot.xlim(-8, 8)
        pyplot.ylim(-90,45)
        pyplot.xlabel('Hour Angle (h)')
        pyplot.ylabel('declination ($^o$)')
        pyplot.grid(None, axis='both')
        return
    
    ax = pyplot.subplot(143,polar=1)
    pyplot.title('u,v coverage')
    for tr in tracks:
        w = np.where(tr['observable'])
        for b in tr['baseline']:
            pyplot.plot((tr['PA'][b][w]-90)*np.pi/180, tr['B'][b][w],'-',
                color=pyplot.cm.get_cmap(color_map)((tr['dec']-min(dec))/np.ptp(dec)),
                alpha=0.5, markersize=6)
            pyplot.plot((tr['PA'][b][w]+90)*np.pi/180, tr['B'][b][w],'-',
                color=pyplot.cm.get_cmap(color_map)((tr['dec']-min(dec))/np.ptp(dec)),
                alpha=0.5, markersize=6)
    if not Bmax is None:
        ax.set_rmax(Bmax)
    ax.get_xaxis().set_visible(False)
    pyplot.subplot(144)
    for tr in tracks:
        w = np.where(tr['observable'])
        pyplot.plot(tr['dec']*np.ones(2), [0,np.diff(tr['lst'][w]).sum()],
                    color=pyplot.cm.get_cmap(color_map)((tr['dec']-min(dec))/np.ptp(dec)),
                        linewidth=5)
    pyplot.xlim(-85, 45)
    pyplot.ylim(0,10)
    pyplot.xlabel('declination (deg)')
    pyplot.title('observability (hours)')
    return

def plotAllConfigs():
    quads = [['A1','G1','K0','J3'],
             ['D0','H0','G1','I1'],
             ['A1','B2','C1','D0'],
             ['U1','U2','U3','U4']]

    for q in quads:
        for t1 in q:
            for t2 in q:
                if t1<t2:
                    print t1, t2
                    skyCoverage([t1,t2])
                    pyplot.savefig(t1+'-'+t2+'.png')
                    for t3 in q:
                        if t2<t3:
                            print t1, t2, t3
                            skyCoverage([t1,t2,t3])
                            pyplot.savefig(t1+'-'+t2+'-'+t3+'.png')
                            for t4 in q:
                                if t3<t4:
                                    print t1, t2, t3, t4
                                    skyCoverage([t1,t2,t3,t4])
                                    pyplot.savefig(t1+'-'+t2+'-'+t3+'-'+t4+'.png')

def findOpdConfig(telescopes, IP=None, respectOrder=False, bestForFINITO=False, compact=False):
    """
    look in ./VLTI_opd (opdModels variable) for OPD models and check if there is
    a combination of models that can be used.
    
    telescopes is a list of stations, such as "A1", "K0" or "U2". 
    
    IP is an optional list of requested IPs for the telescopes. default IP
    values (used when IP==None) are [1,3], [1,3,5] or [1,3,5,7] for 2, 3 or 4
    telescopes configurations.
    
    Note that the order in which telescopes are given does matter, unless
    'respectOrder' is set to False. Then, any telescope can go to any IP.
    
    returns a list of tuples containing the len(telescopes)-1 models for a
    particular configuration. 
    """
    
    models = loadOpdModels()
    if IP is None:
        IP = range(2*len(telescopes))[1::2]
    if respectOrder:
        suitables = []
        for ref in range(len(telescopes)): # --reference telescope
            # -- all the other telescopes:
            idx = filter(lambda x: x!=ref, range(len(telescopes)))
            tmp = []
            for k in idx:
                tmp.extend(filter(lambda x: (x['S1']==telescopes[ref] and
                                             x['S2']==telescopes[k] and
                                             x['IP1']==IP[ref] and
                                             x['IP2']==IP[k]) , models))
            # -- if a model was found for each pair with the reference:
            if all([isinstance(x, dict) for x in tmp]) and \
                len(tmp)==len(telescopes)-1:
                suitables.append(tuple([x['name'] for x in tmp]))
    else:
        # -- loop on the permutations of IP and apply function
        # -- respecting the order this time
        suitables = []
        for ips in itertools.permutations(IP):
            suitables.extend(findOpdConfig(telescopes, IP=ips,
                                           respectOrder=True, compact=False))
    # -- now check that it forms an actual triangle
    if len(telescopes)==3:
        newSuitables = []
        for k,s in enumerate(suitables):
            if s[0].split('-')[0]==s[1].split('-')[0] or\
               s[0].split('-')[1].split('_')[0]==s[1].split('-')[0] or\
               s[0].split('-')[0]==s[1].split('-')[1].split('_')[0] or\
               s[0].split('-')[1].split('_')[0]==s[1].split('-')[1].split('_')[0]:
                newSuitables.append(s)
        suitables=newSuitables
    
    # -- in case of a triangle, optimize for FINITO (track on 2 shortest baselines)
    if len(telescopes)==3 and bestForFINITO and len(suitables)>1:
        #print 'choosing best configuration for FINITO out of', len(suitables),
        #print 'configurations'
        # -- conpute baselines length at zenith
        b = [projBaseline(telescopes[0], telescopes[1],[0,-24.0], [0.0])['B'][0],
             projBaseline(telescopes[1], telescopes[2],[0,-24.0], [0.0])['B'][0],
             projBaseline(telescopes[0], telescopes[2],[0,-24.0], [0.0])['B'][0]]
        #print 'baselines length at Zenith:'
        #print '  '*4, telescopes[0]+'-'+telescopes[1], round(b[0],1), 'm'
        #print '  '*4, telescopes[1]+'-'+telescopes[2], round(b[1],1), 'm'
        #print '  '*4, telescopes[0]+'-'+telescopes[2], round(b[2],1), 'm'
        
        if max(b)==b[0]:
            T0=telescopes[2]
        elif max(b)==b[1]:
            T0=telescopes[0]
        else:
            T0=telescopes[1]
        # -- keep models for which the 
        bestForF = filter(lambda x: T0+'IP3' in x[0] or T0+'IP3' in x[1],  suitables)
        if len(bestForF)>0:
            suitables = bestForF
        else:
            print 'failed to optimize for FINITO:', T0, 'should be on IP3'
            
    if compact:
        # -- return a more compact version of the result
        res = []
        for s in suitables:
            tmp = []
            for x in s:
                tmp.extend([x.split('-')[0],x.split('-')[1].split('_')[0]])
            tmp = list(set(tmp))
            tmp = np.array(tmp)[np.argsort([x[-3:] for x in tmp])]
            res.append('-'.join(list(tmp)))
        return res
    return suitables  
    
def loadOpdModels():
    global opdModels
    files = os.listdir(opdModels)
    files = filter(lambda x: 'model.dat' in x, files)
    res = []
    for f in files:
        T1 = f.split('-')[0]
        T2 = f.split('-')[1].split('_')[0]
        res.append({'name':f,
               # Delay Line, Telescope, Station and Input Channel
               'DL1':int(T1[2]),'T1':T1[3:6],'S1':T1[6:8],'IP1':int(T1[10]),
               'DL2':int(T2[2]),'T2':T2[3:6],'S2':T2[6:8],'IP2':int(T2[10])})
    return res
        
    