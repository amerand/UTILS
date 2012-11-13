import vlti
import ephem
import numpy as np
from matplotlib import pyplot as plt

def vltiAscii(T, targets, date=None, moonDist=20.0, useColor=True,
              windDir=None, lstStep=5, sortByRA=False, useUnicode=False,
              plotOpl=False):
    """
    T: list of telescopes, e.g. ['A0', 'K0', 'G1']

    targets: list of target Simbad can resolve.
    
    date=None -> Now (default)
    date='2012/03/20'
    date='2012/03/20 03:00:00' (UT time)
    
    moonDist: minimum distance to the Moon (degrees)
    
    windDir: direction of the wind in case of wind restriction. Can be
    in degrees (0 is North and 90 is East, same convention as Paranal
    ASM) or given litteraly: 'N', 'SW', 'NE' etc. None (default) means
    no restrictions.
    """
    vlt = ephem.Observer()
    vlt.lat = '-24:37:38'
    vlt.long = '-70:24:15'
    vlt.elevation = 2635.0
  
    """
    http://stackoverflow.com/questions/287871/print-in-terminal-with-colors-using-python
    """
    colors = {'GREEN':'\033[92m', 'GREENBG':'\033[42m',
              'BLUE':'\033[94m',  'BLUEBG':'\033[44m',
              'RED':'\033[91m',   'REDBG':'\033[41m',
              'YELLOW':'\033[93m','YELLOWBG':'\033[43m',
              'GRAY':'\033[90m',
              'CYAN':'\033[96m',
              'MAGENTA':'\033[95m',
              'NO_COL':'\033[0m' }
        
    if not useColor:
        for k in colors.keys():
            colors[k] = ''

    if isinstance(windDir, str):
        direction = {'':None, 'N':0, 'S':180, 'E':90, 'W':270,
                     'NE':45, 'NW':315, 'SW':225, 'SE':135}
        windDir = direction[windDir]

    windArrow = {'0':8659, '45':8665, '90':8656, '135':8662,
                 '180':8657, '225':8663, '270':8658, '315':8664,
                 '360':8659}
    
    if useUnicode and not windDir is None:
        windK = np.argmin([np.abs(float(k)-windDir) for k in windArrow.keys()])
        windChr = unichr(windArrow[windArrow.keys()[windK]])
    else:
        windChr='w'
    # ephemeris of the Sun
    sun = ephem.Sun()
    sun_eph={}

    # now
    if date is None:
        vlt.date = ephem.date(ephem.now())
    else:
        vlt.date= date
    sun_eph['NOW'] = vlt.date
        
    sun_eph['SUNSET'] = vlt.next_setting(sun) # sun sets
    vlt.horizon = '-18:0'
    sun_eph['TWI END']=vlt.next_setting(sun) # TWI ends
    vlt.horizon = '0'
    sun_eph['SUNRISE'] = vlt.next_rising(sun) # sun sets 
    vlt.horizon = '-18:0'
    sun_eph['TWI START']=vlt.next_rising(sun) # TWI starts
    vlt.horizon = '0:0'

    # -- check that the events are in the right order
    for k in sun_eph.keys():
        if float(sun_eph[k]) > float(vlt.date) + 0.5:
            sun_eph[k] = ephem.date(sun_eph[k]-1)
        if float(sun_eph[k]) < float(vlt.date) - 0.5:
            sun_eph[k] = ephem.date(sun_eph[k]+1)
    
    # offset between LST and UT
    lstOffset = float(vlt.sidereal_time())/(2*3.141592) - \
                float(vlt.date)%1 - 0.5
    
    # position of the Moon
    moon = ephem.Moon()
    moon.compute(vlt)
    # unicode sympol of the moon: 0 is full, 50 is new
    moonChr = {'0':9673, '25':9681, '40':9789, '50':9712,
               '60':9790, '75':9680, '100':9673}
    moonChr = {'00':9673, '25':5613, '40':9789, '50':9712,
               '60':9790, '75':5610, '100':9673}
    
    moonK = np.argmin([np.abs(float(k)-moon.phase) for k in moonChr.keys()])
    
    if useUnicode:
        # does not work :(
        #moonC = unichr(moonChr[moonChr.keys()[moonK]])
        moonC = unichr(9790); moonC=unichr(9686)
        upChr = unichr(10003); upChr = unichr(9624)
        downChr = unichr(8901)
        nowChr = unichr(9731)
        lineChr = unichr(8901)
        tickChr = [unichr(9312), unichr(9313), unichr(9314)]
        transitChr = unichr(9615)
        sunChr = [unichr(9728), unichr(9728)]
        twiChr = [unichr(9732), unichr(9732)]
        notobsChr = unichr(9889)
    else:
        moonC = 'm'
        upChr = '='
        downChr = '_'
        notobsChr = '-'
        nowChr = '@'
        lineChr = '-'
        tickChr = [':', ':', ':']
        transitChr = 'T'
        sunChr = ['<', '>']
        twiChr = ['(', ')']
        notobsChr = '^'
        
    sun_eph_lst = {}
    for k in sun_eph.keys():
        sun_eph_lst[k] = (sun_eph[k]+lstOffset-0.5)%1*24

    kz= ['SUNSET', 'TWI END', 'TWI START', 'SUNRISE']
    for k in range(3):
        #print kz[k], sun_eph_lst[kz[k]]
        if sun_eph_lst[kz[k+1]]<sun_eph_lst[kz[k]]:
            sun_eph_lst[kz[k+1]] += 24
        
    lst = np.arange(np.floor(sun_eph_lst['SUNSET'])-0,
                    np.ceil(sun_eph_lst['SUNRISE'])+1./lstStep,
                    1./lstStep)
    sun_eph_lst['1/2'] = 0.5*sun_eph_lst['TWI START'] +\
                         0.5*sun_eph_lst['TWI END']
    sun_eph_lst['1/4'] = 0.25*sun_eph_lst['TWI START'] +\
                         0.75*sun_eph_lst['TWI END']
    sun_eph_lst['3/4'] = 0.75*sun_eph_lst['TWI START'] +\
                         0.25*sun_eph_lst['TWI END']
    ## if useUnicode:
    ##     lstString1 = []
    ##     for l in lst:
    ##         if np.abs(int(l)-l)<1e-6:
    ##             print l,l%24.0, (l%24)>20, 
    ##             if (l%24)<21:
    ##                 lstString1.append(unichr(9311+int(l%24)))
    ##             else:
    ##                 lstString1.append(unichr(12860+int(l%24)))
    ##         else:
    ##             lstString1.append(' ')
    ##     lstString2 = ''
    ## else:
   
    lstString1 = [str(int((l%24+0.01)/10.))
                  if np.abs(round(l,0)-l)<0.01 else ' ' for l in lst]
    lstString2 = [str(int((l%24+0.01)%10.))
                  if np.abs(round(l,0)-l)<0.01 else ' ' for l in lst]
    
    lstString3 = [colors['GRAY']+lineChr+colors['NO_COL'] for l in lst]
    lstString3[np.abs(lst-sun_eph_lst['SUNSET']).argmin()] =\
                                colors['YELLOW']+sunChr[0]+colors['NO_COL']
    lstString3[np.abs(lst-sun_eph_lst['SUNRISE']).argmin()] =\
                                colors['YELLOW']+sunChr[1]+colors['NO_COL']
    lstString3[np.abs(lst-sun_eph_lst['TWI START']).argmin()] =\
                                colors['YELLOW']+twiChr[1]+colors['NO_COL']
    lstString3[np.abs(lst-sun_eph_lst['TWI END']).argmin()] =\
                                colors['YELLOW']+twiChr[0]+colors['NO_COL']
    lstString3[np.abs(lst-sun_eph_lst['1/4']).argmin()] =\
                                colors['YELLOW']+tickChr[0]+colors['NO_COL']
    lstString3[np.abs(lst-sun_eph_lst['1/2']).argmin()] =\
                                colors['YELLOW']+tickChr[1]+colors['NO_COL']
    lstString3[np.abs(lst-sun_eph_lst['3/4']).argmin()] =\
                                colors['YELLOW']+tickChr[2]+colors['NO_COL']
    lstString3[np.abs(lst-sun_eph_lst['NOW']).argmin()] =\
                                colors['GREEN']+nowChr+colors['NO_COL']
    
    obsStrings = []
    if isinstance(targets, str): # single target -> list of one element
        targets = [targets]
        
    RAs = []
    
    if plotOpl:
        plt.figure(0)
        plt.clf()
        kp=0
        
    for target in targets:
        min_altitude = 25.
        tmp = vlti.nTelescopes(T, target, lst, min_altitude=min_altitude)
        RAs.append(tmp['ra'])
        dra = (moon.ra*12/np.pi-tmp['ra'])*np.cos(tmp['dec']*np.pi/180)*15
        ddec = (moon.dec*180/np.pi-tmp['dec'])
        dist = np.sqrt(dra**2+ddec**2)
        
        obsString = []
        for k in range(len(tmp['alt'])):
            if  tmp['observable'][k]:
                if useUnicode:
                    obsString.append(unichr(8320+int(tmp['alt'][k]/10.)))
                    #obsString.append(unichr(10111+int(tmp['alt'][k]/10.)))
                else:
                    obsString.append(str(int(tmp['alt'][k]/10.)))
                    
            elif tmp['alt'][k]>=min_altitude:
                obsString.append(notobsChr)
            else:
                obsString.append(downChr)
        
        for k in range(len(obsString)-1)[::-1]:
            if obsString[k]==obsString[k+1]:
                if obsString[k]==downChr: # too low
                    obsString[k+1]=downChr
                elif obsString[k]==notobsChr: # not observable
                    obsString[k+1]=notobsChr
                elif dist < moonDist:
                    obsString[k+1]=moonC
                else:
                    obsString[k+1]=upChr
                        
        iT = tmp['alt'].argmax()
        if not iT==0 and not iT==len(tmp['alt']):
            obsString[iT] = transitChr

        # wind restriction:
        if not windDir is None:
            # AZ convetion are different between wind and Target!
            for k in range(len(obsString)):
                if obsString[k] == upChr and \
                       np.abs((-tmp['az'][k]-windDir)%360-180)<=90:
                    obsString[k] = windChr
            
        if useColor:
            for k in range(len(obsString)):
                if obsString[k]==downChr:
                    obsString[k]=colors['GRAY']+downChr+colors['NO_COL']
                elif obsString[k]==notobsChr:
                    obsString[k]=colors['RED']+notobsChr+colors['NO_COL']
                elif obsString[k]==moonC:
                    obsString[k]=colors['BLUE']+moonC+colors['NO_COL']
                elif obsString[k]==upChr:
                    obsString[k]=colors['GREEN']+upChr+colors['NO_COL']
                elif obsString[k]==transitChr:
                    obsString[k]=colors['YELLOW']+transitChr+colors['NO_COL']
                elif obsString[k]==windChr:
                    obsString[k]=colors['REDBG']+windChr+colors['NO_COL']
        obsString = reduce(lambda x,y: x+y, obsString)
        if isinstance(target, str):
            targ_name = target
        else:
            # assume [ra dec]
            targ_name = '%2d:%2d %3d:%2d' % (int(target[0]),
                                             int(60*(target[0]-int(target[0]))),
                                             int(target[1]),
                                             int(60*(abs(target[1])-
                                                     int(abs(target[1]))))
                                            )
        obsStrings.append(obsString+' '+colors['YELLOW']+
                          targ_name+colors['NO_COL'])
        if plotOpl:
            cols = ['r', 'g', 'b', 'c', 'm', 'y', '0.5']
            typ = ['-','--',':','-.']
            t = 0
            for k in tmp['opd'].keys():
                plt.plot(tmp['lst'], tmp['opd'][k], label=targ_name+' '+k,
                         color=cols[kp%len(cols)],
                         linestyle=typ[t%len(typ)])
                plt.plot(tmp['lst'][tmp['observable']],
                         tmp['opd'][k][tmp['observable']],
                         linestyle=typ[t%len(typ)],
                         color=cols[kp%len(cols)], linewidth=5, alpha=0.5)
                t+=1
            kp+=1
    
    if plotOpl:
        plt.legend(prop={'size':9})
        plt.xlabel('LST (h)')
        plt.ylabel('OPD (m)')
        plt.grid()
    if sortByRA:
        obsStrings = [obsStrings[k] for k in np.argsort(RAs)]
    
    lstString1 = reduce(lambda x,y: x+y, lstString1)
    if len(lstString2)>1:
        lstString2 = reduce(lambda x,y: x+y, lstString2)
    lstString3 = reduce(lambda x,y: x+y, lstString3)
    
    print colors['BLUE'], '('*3, reduce(lambda x,y:x+'-'+y, T), \
          sun_eph['NOW'], moonC+':Moon (%2.0f%%)<'%(moon.phase), moonDist,
    if not windDir is None:
        print windChr+':wind restr. (%3d)'%(windDir),
    print ')'*3+colors['NO_COL']
    print colors['GRAY']+lstString1+' LST'+colors['NO_COL']
    if lstString2!='':
        print colors['GRAY']+lstString2+colors['NO_COL']
    print lstString3

    if useUnicode:
        # zodiac strip
        lstString4 = [unichr(9800+int(((l-3)%24)/2.)) if 
                      np.abs(int(l)-l)<1e-6 and np.abs((l+1)%2. < 1e-6)
                      else ' ' for l in lst]
        lstString4 = reduce(lambda x,y: x+y, lstString4)
        print lstString4
        
    for k in range(len(targets)):
        print obsStrings[k]
    
