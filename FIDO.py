import numpy as np
import math
import sys
from scipy.special import jv
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from Tkinter import *
from pylab import setp
import random 

# Geometry programs
def SPH2CART(sph_in):
    r = sph_in[0]
    colat = (90. - sph_in[1]) * dtor
    lon = sph_in[2] * dtor
    x = r * np.sin(colat) * np.cos(lon)
    y = r * np.sin(colat) * np.sin(lon)
    z = r * np.cos(colat)
    return [x, y, z]

def CART2SPH(x_in):
# calcuate spherical coords from 3D cartesian
# output lat not colat
    r_out = np.sqrt(x_in[0]**2 + x_in[1]**2 + x_in[2]**2)
    colat = np.arccos(x_in[2] / r_out) * 57.29577951
    lon_out = np.arctan(x_in[1] / x_in[0]) * 57.29577951
    if lon_out < 0:
        if x_in[0] < 0:
            lon_out += 180.
        elif x_in[0] > 0:
            lon_out += 360. 
    elif lon_out > 0.:
	    if x_in[0] < 0:  lon_out += 180. 
    return [r_out, 90. - colat, lon_out]

def rotx(vec, ang):
# Rotate a 3D vector by ang (input in degrees) about the x-axis
    ang *= dtor
    yout = np.cos(ang) * vec[1] - np.sin(ang) * vec[2]
    zout = np.sin(ang) * vec[1] + np.cos(ang) * vec[2]
    return [vec[0], yout, zout]

def roty(vec, ang):
# Rotate a 3D vector by ang (input in degrees) about the y-axis
    ang *= dtor
    xout = np.cos(ang) * vec[0] + np.sin(ang) * vec[2]
    zout =-np.sin(ang) * vec[0] + np.cos(ang) * vec[2]
    return [xout, vec[1], zout]

def rotz(vec, ang):
# Rotate a 3D vector by ang (input in degrees) about the y-axis
	ang *= dtor
	xout = np.cos(ang) * vec[0] - np.sin(ang) * vec[1]
	yout = np.sin(ang) * vec[0] + np.cos(ang) * vec[1]
	return [xout, yout, vec[2]]



# New functions
def isinCME(vec_in, CME_shape):
    # Check and see if the requested point is actually in the CME and return
    # the cylindrical radial distance (from center of FR)
    # Function assumes vec_in is in CME Cartesian coord sys
    thetas = np.linspace(-math.pi/2, math.pi/2, 101)
    # determine the xz positions of the rope axis
    xFR = CME_shape[0] + CME_shape[1] * np.cos(thetas)
    zFR = CME_shape[3] * np.sin(thetas)
    dists2 = (vec_in[0] - xFR)**2 + vec_in[1]**2 + (vec_in[2] - zFR)**2
    myb2 = np.min(dists2)
    minidxs = np.where(dists2 == myb2)
    # unwrap
    minidx = minidxs[0]
    temp = thetas[np.where(dists2 == myb2)]
    mythetaT = temp[0]
    #testing 
    #if dists2[minidx[0]]/CME_shape[2]<1: 
    #    print dists2[minidx]/CME_shape[2], mythetaT
    # add a second iteration to refine B
    # see which side of minidx the actual minimum is on
    if len(minidx) == 1: # if perfectly symmetric can get two equi dists at back edge
        if minidx < len(dists2) - 1: # check to make sure not already at edge
            if dists2[minidx-1] < dists2[minidx+1]: startidx = minidx - 1
            else:  startidx = minidx + 1
    	# repeat the same procedure on the side with the acutal min
            if dists2[minidx-1] != dists2[minidx+1]:
                thetas2 = np.linspace(thetas[startidx], thetas[minidx], 101)
                xFR2 = CME_shape[0] + CME_shape[1] * np.cos(thetas2)
                zFR2 = CME_shape[3] * np.sin(thetas2)
                dists2 = (vec_in[0] - xFR2)**2 + vec_in[1]**2 + (vec_in[2] - zFR2)**2
                myb2 = np.min(dists2)
                minidxs = np.where(dists2 == myb2)
                minidx = minidxs[0]
                temp = thetas2[np.where(dists2 == myb2)]
                mythetaT = temp[0]
    myb = np.sqrt(myb2)
    CME_crossrad = CME_shape[2]
    if (myb > CME_crossrad):
        myb = -9999.
        return myb, -9999, -9999., -9999, -9999
    else:
	# thetaP should swing quickly at closest approach to FR
        mythetaP = np.arcsin(vec_in[1] / myb)
        origTP = mythetaP
	# doesn't run through second loop sometimes
        if ('xFR2' not in locals()) and ('xFR2' not in globals()):
            xFR2 = xFR
	# convert thetaP to complement
        # had a random case that blew up here -> added try/except
        try: 
            if vec_in[0] < (xFR2[minidx] ):	
                if vec_in[1] > 0: mythetaP = math.pi - mythetaP
                else: mythetaP = -math.pi - mythetaP
        except:
            pass
    return myb, mythetaT, mythetaP, 0, CME_crossrad

def getBvector(CME_shape, minb, thetaT, thetaP):
    tdir = np.array([-(CME_shape[1] + minb * np.cos(thetaP)) * np.sin(thetaT), 0., (CME_shape[3] + minb * np.cos(thetaP)) * np.cos(thetaT)])
    pdir = np.array([-minb * np.sin(thetaP) * np.cos(thetaT), minb * np.cos(thetaP), -minb * np.sin(thetaP) * np.sin(thetaT)])
    tmag = np.sqrt(np.sum(tdir**2))
    pmag = np.sqrt(np.sum(pdir**2))
    tdir = tdir / tmag
    pdir = pdir / pmag
    return tdir, pdir

def update_insitu():
    CME_shape = np.zeros(4)
    # set up CME shape as [d, a, b, c]
    shapeC = np.tan(CMEAW*dtor) / (1. + CMESRB + np.tan(CMEAW*dtor) * (CMESRA + CMESRB)) 
    global CMEtilt # being stupid without this for some reason
    if np.abs(CMEtilt) < 1e-3: 
        CMEtilt = 1e-3
    dtorang = (CMElon - FFlon0) / np.sin(CMEtilt * dtor) * dtor
    CMEnose = 215. / np.sqrt((1 - (CMESRA + CMESRB) * shapeC * (1. - np.cos(dtorang)))**2 + (1 + CMESRB)**2 * shapeC**2 * np.sin(dtorang)**2)  - 10.
    CME_shape[3] = CMEnose * shapeC
    CME_shape[1] = CME_shape[3] * CMESRA
    CME_shape[2] = CME_shape[3] * CMESRB 
    CME_shape[0] = CMEnose - CME_shape[1] - CME_shape[2]
    obsBx = []
    obsBy = []
    obsBz = []
    tARR = []
    t = 0.
    CMEB = CMEB0
    FFlon = FFlon0
    thetaPprev = -42
    switch = 0
    flagExp = False
    while t < tmax:  
        tARR.append(t/3600.)
	# convert flyer position to CME Cartesian coordinates
	# have to do in loop to account for Earth orbit
	# get Sun xyz position
        FF_sunxyz = SPH2CART([215, FFlat, FFlon])
	# rotate to CMEcoord system
        temp = rotz(FF_sunxyz, -CMElon)
        temp2 = roty(temp, CMElat)
        FF_CMExyz = rotx(temp2, CMEtilt)
	# determine CME expansion and propagation
	# calculate CME shape
        if flagExp == False:
            CME_shape[3] = CMEnose * np.tan(CMEAW*dtor) / (1. + CMESRB + np.tan(CMEAW*dtor) * (CMESRA + CMESRB))
            CME_shape[1] = CME_shape[3] * CMESRA
            CME_shape[2] = CME_shape[3] * CMESRB 
        CME_shape[0] = CMEnose - CME_shape[1] - CME_shape[2]
	# check to see if the flyer is currently in the CME
	# if so, get its position in CME torus coords 
        minb, thetaT, thetaP, flagit, CME_crossrad = isinCME(FF_CMExyz, CME_shape)
	# need to check that thetaP doesn't jump as it occasionally will
        if flagit != -9999:
	    # define it the first time its used
            if thetaPprev==-42: 
                tempa = CME_shape[1]+CME_shape[2]
                tempb = CME_shape[3]+CME_shape[2]
                thetaPprev = thetaP
                if expansion_model == 'None': flagExp = True
            delThetaP = np.abs(thetaP - thetaPprev)
            if (delThetaP > 0.5) and (np.abs(delThetaP < 3.1)): thetaP = thetaPprev
            thetaPprev = thetaP
	# get the toroidal and poloidal magnetic field
        Btor = CMEB * jv(0, 2.4 * minb / CME_crossrad)
        Bpol = CMEB * CMEH * jv(1, 2.4 * minb / CME_crossrad)
	# save the magnetic field if in CME
        if flagit != -9999.:
            # convert this into CME Cartesian coordinates
            tdir, pdir = getBvector(CME_shape, minb, thetaT, thetaP)
            Bt_vec = Btor * tdir
            Bp_vec = Bpol * pdir
            Btot = Bt_vec + Bp_vec
            # convert to spacecraft coord system
            temp = rotx(Btot, -CMEtilt)
            temp2 = roty(temp, CMElat - FFlat) 
            BSC = rotz(temp2, CMElon - FFlon)
            #print minb/CME_crossrad, thetaT, thetaP, Btor, Bpol, BSC
            obsBx.append(BSC[0])
            obsBy.append(BSC[1])
            obsBz.append(BSC[2])
        else:
            obsBx.append(0.)
            obsBy.append(0.)
            obsBz.append(0.)
        #print t/3600, minb, CME_shape[2], CMEnose, CME_shape[1]
        t += dt
	# CME nose moves to new position
        CMEnose += CMEvr * dt / 7e5 # change position in Rsun
	# update the total magnetic field
        if flagExp == False:
            CMEB *= ((CMEnose - CMEvr * dt / 7e5) / CMEnose)**2
	# determine new lon of observer
        FFlon += dt * rotspeed
    print FFlon, rotspeed
    return obsBx, obsBy, obsBz, tARR

def pullGUIvals():
    global CMElat, CMElon, CMEtilt, CMEAW, CMESRA, CMESRB, CMEB0, CMEH, CMEvr, tshift, CMEstart, CMEend, CMEmid, plotCME
    global expansion_model
    expansionToggle = expansionToggleVAR.get()
    expansion_model = 'None'
    if expansionToggle == 0: expansion_model = 'Self-Similar'
    global Autonormalize
    Autonormalize = False
    if autonormVAR.get() == 1: Autonormalize = True
    
    try:
        CMElat = float(e1.get())
    except:
        if canprint: print("Need CME latitude")
        plotCME = False
    try:
        CMElon = float(e2.get())  
    except:
        if canprint: print("Need CME longitude")
        plotCME = False
    try:  
        CMEtilt = float(e3.get())
    except:
        if canprint: print("Need CME tilt")
        plotCME = False
    try:
        CMEAW   = float(e3b.get())
    except:
        if canprint: print("Need CME angular width")
        plotCME = False
    try:
        CMESRA  = float(e4.get())
    except:
        if canprint: print("Need CME shape A")
        plotCME = False
    try:
        CMESRB  = float(e5.get())
    except:
        if canprint: print("Need CME shape B")
        plotCME = False
    try:
        CMEB0   = float(e6.get())
    except:
        if canprint: print("Need CME flux rope magnitude B0")
        plotCME = False
    try:
        CMEH	= float(e7.get())
    except:
        if canprint: print("Need handedness +1 or -1")
        plotCME = False
    try:
        CMEvr   = float(e8.get())
    except:
        if canprint: print("Need CME velocity")
        plotCME = False
    try:
        tshift  = float(e9.get())
    except:
        tshift = 0
    global FFlat, FFlon0, plotstart, plotend, pad
    try:    
        FFlat   = float(eR1.get())
    except:
        if canprint: print("No spacecraft latitude, assuming at 0")
        FFlat = 0
    try:
        FFlon0  = float(eR2.get()) 
    except:
        if canprint: print("Need spacecraft longitude, assuming at 0")
        FFlon0 = 0
    try:
        CMEstart = float(eS1.get())
        plotstart = CMEstart - pad/24.
        if ISfilename != False:
            if CMEstart < minISdate:
                print('CME_start before first in situ time')
                plotCME = False
                
            CMEend = float(eS2.get())
            plotend = CMEend + pad/24.
            if CMEend > maxISdate:
                print('CME_end before last in situ time')
                plotCME = False
    except:
        canScore = False
        # set equal to plot start if have data otherwise it should stay at 0
        if (range_flag == False) and (plotstart > 0):
            CMEstart = plotstart
            CMEend = plotend
        if canprint: print('No CME start/stop time -> no score')
        if canprint: print('Assuming CME starts at plot start')
    
    
def update_plot():
    global CMElat, CMElon, CMEtilt, CMEAW, CMESRA, CMESRB, CMEB0, CMEH, CMEvr, tshift, CMEstart, CMEend, CMEmid, avg_obs_B
    global range_flag, plotstart, plotend
    global scale_flag, canScore
    plotCME = True
    if Launch_GUI:
        pullGUIvals()
    # originally programmed tilt to be deg clockwise from N
    # but counterclockwise from W more common in literature
    # rather than alter code for shape, just convert it here
    global othertilt
    othertilt = CMEtilt
    if CMEtilt >=-90:
        CMEtilt = 90 - CMEtilt
    else:
        CMEtilt = -180 - CMEtilt
        
    if plotCME==True:
        # define as globals so can use to calculate score
        global obsBx, obsBy, obsBz, obsB, tARR, totalscore
        obsBx, obsBy, obsBz, tARR = update_insitu()
        obsBx, obsBy, obsBz = np.array(obsBx), np.array(obsBy), np.array(obsBz)
        tARR = np.array(tARR) # tarr in hours from 0
        zerotime = (CMEstart + tshift/24.) # in days
        tARR = tARR/24. # also in days
    	# find the beginning and end of the actual magnetic cloud portion of the data 
    	# check if beginning and end are equal to 0
        if np.sum(np.abs(obsBx)) > 0.:
            if obsBx[0] == 0:
                MCidxi = np.min(np.where(np.abs(obsBx) > 0))
    		# trim off excess zeroes since finding starting point of CME difficult
                idx1 = np.max([0, MCidxi])
                obsBx = obsBx[idx1:] 
                obsBy = obsBy[idx1:] 
                obsBz = obsBz[idx1:] 
                tARR  = tARR[idx1:]
                tARR -= tARR[0] # restart the t array at 0
            if obsBx[-1] == 0:
                MCidxf = np.max(np.where(np.abs(obsBx) > 0))
                idx2 = MCidxf
                if idx2 > len(obsBx) - 1: idx2 = len(obsBx) -1 
                obsBx = obsBx[:idx2] 
                obsBy = obsBy[:idx2] 
                obsBz = obsBz[:idx2] 
                tARR  = tARR[:idx2]
            obsB = np.sqrt(obsBx**2 + obsBy**2 + obsBz**2)
            # scale the CME to match at midpoint (ignoring B0 with this)
            tARR += zerotime
            CMEmid = np.mean(tARR)
            if ISfilename != False:
                avg_obs_B = np.mean(d_Btot[np.where(np.abs(d_tUN - CMEmid) < 2./24.)])
            scale = 1.   
            if (Autonormalize==True) and (ISfilename !=False): 
                    cent_idx = np.where(np.abs(tARR - CMEmid) < 2./24.)[0]  
                    if len(cent_idx) > 0: 
                        avg_FIDO_B = np.mean(obsB[cent_idx])
                        scale = avg_obs_B / avg_FIDO_B
                    else:
                    	print('CME too short to autonormalize, reverting to B0')
                    obsBx *= scale
                    obsBy *= scale
                    obsBz *= scale
                    obsB = np.sqrt(obsBx**2 + obsBy**2 + obsBz**2)

            if canScore:
                scoreBx, scoreBy, scoreBz = calc_score()
                if canprint: print('score '+ str(totalscore))
        else:
                plotCME = False
                print ('No impact expected')
                totalscore = 9999.
                scoreBx, scoreBy, scoreBz = 9999., 9999., 9999.
 
    else:
	    print('No impact expected')
	    plotCME = False
        
    if range_flag == True:
        if ISfilename != False:
            plotstart = d_tUN[0]
            plotend = d_tUN[-1]
        else:
            plotstart = tARR[0]-0.2
            plotend = tARR[-1]+0.2
        range_flag = False
        
    ax2.clear()
    ax3.clear()
    ax4.clear()
    ax5.clear()
    # temporarily set to nice things so can easily determine min/max ranges
            
    if ISfilename != False:    
        ax2.plot(d_tUN, d_Btot, 'k', linewidth=4)
        ax3.plot(d_tUN, d_Bx, 'k', linewidth=4)
        ax4.plot(d_tUN, d_By, 'k', linewidth=4)
        ax5.plot(d_tUN, d_Bz, 'k', linewidth=4)
        if plotCME != False:
            maxBtot = np.max([np.max(d_Btot), np.max(obsB)])
            minBx = np.abs(np.min([np.min(d_Bx), np.min(-obsBx)]))
            maxBx = np.max([np.max(d_Bx), np.max(-obsBx)])
            minBy = np.abs(np.min([np.min(d_By), np.min(-obsBy)]))
            maxBy = np.max([np.max(d_By), np.max(-obsBy)])
            minBz = np.abs(np.min([np.min(d_Bz), np.min(obsBz)]))
            maxBz = np.max([np.max(d_Bz), np.max(obsBz)])
            
        else:
            maxBtot = np.max(d_Btot)
            minBx = np.abs(np.min(d_Bx))
            maxBx = np.max(d_Bx)
            minBy = np.abs(np.min(d_By))
            maxBy = np.max(d_By)
            minBz = np.abs(np.min(d_Bz))
            maxBz = np.max(d_Bz)
    elif plotCME != False:
        maxBtot = np.max(obsB)
        minBx = np.abs(np.min(-obsBx))
        maxBx = np.abs(np.max(-obsBx))
        minBy = np.abs(np.min(-obsBy))
        maxBy = np.abs(np.max(-obsBy))
        minBz = np.abs(np.min(obsBz))
        maxBz = np.abs(np.max(obsBz))
        
    else:
        maxBtot = 1
        minBx, minBy, minBz = 1, 1, 1 # set as neg later
        maxBx, maxBy, maxBz = 1, 1, 1
    
    if canScore:
        ax2.plot([CMEstart, CMEstart], [0., 1.3*maxBtot], 'k--', linewidth=2)
        ax2.plot([CMEend, CMEend], [0., 1.3*maxBtot], 'k--', linewidth=2)
        ax3.plot([CMEstart, CMEstart], [-1.3*minBx, 1.3*maxBx], 'k--', linewidth=2)
        ax3.plot([CMEend, CMEend], [-1.3*minBx, 1.3*maxBx], 'k--', linewidth=2)
        ax4.plot([CMEstart, CMEstart], [-1.3*minBy, 1.3*maxBy], 'k--', linewidth=2)
        ax4.plot([CMEend, CMEend], [-1.3*minBy, 1.3*maxBy], 'k--', linewidth=2)
        ax5.plot([CMEstart, CMEstart], [-1.3*minBz, 1.3*maxBz], 'k--', linewidth=2)
        ax5.plot([CMEend, CMEend], [-1.3*minBz, 1.3*maxBz], 'k--', linewidth=2)

    ax2.set_ylabel('B (nT)')
    setp(ax2.get_xticklabels(), visible=False)    
    ax3.set_ylabel('B$_x$ (nT)')
    setp(ax3.get_xticklabels(), visible=False)
    ax4.set_ylabel('B$_y$ (nT)')
    setp(ax4.get_xticklabels(), visible=False)
    ax5.set_ylabel('B$_z$ (nT)')
    if plotCME == True: 
        ax2.plot(tARR, obsB, 'r', linewidth=4)
        ax3.plot(tARR, -obsBx, 'r', linewidth=4)
        ax4.plot(tARR, -obsBy, 'r', linewidth=4)
        ax5.plot(tARR, obsBz, 'r', linewidth=4)
        if (ISfilename != False) and canScore:
            ax2.annotate('%0.2f'%(totalscore), xy=(1, 0), color='r', xycoords='axes fraction', fontsize=16, horizontalalignment='right', verticalalignment='bottom')    
            ax3.annotate('%0.2f'%(scoreBx), xy=(1, 0), color='r', xycoords='axes fraction', fontsize=16, horizontalalignment='right', verticalalignment='bottom')
            ax4.annotate('%0.2f'%(scoreBy), xy=(1, 0), color='r', xycoords='axes fraction', fontsize=16, horizontalalignment='right', verticalalignment='bottom')    
            ax5.annotate('%0.2f'%(scoreBz), xy=(1, 0), color='r', xycoords='axes fraction', fontsize=16, horizontalalignment='right', verticalalignment='bottom')  
    
    # indices calculation - first figure out rotating to GSM
    if (plotCME == True) and (show_indices == True):
        fracyear = CMEstart / 365.
        rotang = 23.856 * np.sin(6.289 * fracyear + 0.181) + 8.848
        GSMBx = []
        GSMBy = []
        GSMBz = []
        #print('Rotating by '+str(rotang)+' to GSM')
        for i in range(len(obsB)):
            vec = [-obsBx[i], -obsBy[i], obsBz[i]]
            GSMvec = rotx(vec, -rotang)
            GSMBx.append(GSMvec[0])
            GSMBy.append(GSMvec[1])
            GSMBz.append(GSMvec[2])
        # calculate Kp 
        GSMBy = np.array(GSMBy)
        GSMBz = np.array(GSMBz)
        Bt = np.sqrt(GSMBy**2 + GSMBz**2)
        thetac = np.arctan2(np.abs(GSMBy), GSMBz)
        dphidt = np.power(CMEvr, 4/3.) * np.power(Bt, 2./3.) * np.power(np.sin(thetac/2),8/3.) 
        # Emmons 2013 expression for Kp
        #Kp = 0.0002947* dphidt + 1
        # Mays/Savani expression, better behaved for high Kp
        Kp = 9.5 - np.exp(2.17676 - 5.2001e-5*dphidt)
        ax3.plot(tARR, GSMBx, 'b--', linewidth=4, zorder=0)
        ax4.plot(tARR, GSMBy, 'b--', linewidth=4, zorder=0)
        ax5.plot(tARR, GSMBz, 'b--', linewidth=4, zorder=0)
        ax6.plot(tARR, Kp, 'b', linewidth=4)
        ax6.plot([CMEstart, CMEstart], [0, 2+np.max(Kp)], 'k--', linewidth=2)
        ax6.plot([CMEend, CMEend], [0, 2+np.max(Kp)], 'k--', linewidth=2)
        ax6.set_ylim([0, np.max(Kp)+2])
      
    
    scl = 1.25
    ax2.set_ylim([0,scl*maxBtot])
    ax3.set_ylim([-1.*scl*minBx, scl*maxBx])
    ax4.set_ylim([-1.*scl*minBy, scl*maxBy])
    ax5.set_ylim([-1.*scl*minBz, scl*maxBz])
    ax2.set_xlim([plotstart, plotend])
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.subplots_adjust(right=0.8, wspace=0.001, hspace=0.25)
    # fix the bottom plot if added the indices plot
    if show_indices: 
        plt.subplots_adjust(hspace=0.05)
        setp(ax5.get_xticklabels(), visible=False)
        ax6.set_xlabel('Day of Year')
        ax6.set_ylabel('Kp Index')
    else:
        ax5.set_xlabel('Day of Year')
    plt.tight_layout()
    plt.subplots_adjust(top=0.94)
    plt.gcf().canvas.draw()

def calc_score():
    maxt = int(24*(CMEend-CMEstart))# -1
    FIDO_hrt  = np.zeros(maxt)
    FIDO_hrBx = np.zeros(maxt)
    FIDO_hrBy = np.zeros(maxt)
    FIDO_hrBz = np.zeros(maxt)
    ACE_hrBx = np.zeros(maxt)
    ACE_hrBy = np.zeros(maxt)
    ACE_hrBz = np.zeros(maxt)
    
    ACE_t = (d_tUN - CMEstart) * 24
    FIDO_hr = (tARR - CMEstart) *24
    for i in range(maxt):
        ishere = len(FIDO_hr[abs(FIDO_hr - (i+.5))  <= 0.5]) > 0
        if ishere:
            FIDO_hrt[i]  =  np.mean(FIDO_hr[abs(FIDO_hr - (i+.5))  <= 0.5])
            FIDO_hrBx[i] =  np.mean(-obsBx[abs(FIDO_hr - (i+.5))  <= 0.5])
            FIDO_hrBy[i] =  np.mean(-obsBy[abs(FIDO_hr - (i+.5))  <= 0.5])
            FIDO_hrBz[i] =  np.mean(obsBz[abs(FIDO_hr - (i+.5))  <= 0.5])
        if len(d_Bx[abs(ACE_t - (i+.5))  <= 0.5]) > 0:
            ACE_hrBx[i]  =  np.mean(d_Bx[abs(ACE_t - (i+.5))  <= 0.5])
            ACE_hrBy[i]  =  np.mean(d_By[abs(ACE_t - (i+.5))  <= 0.5])
            ACE_hrBz[i]  =  np.mean(d_Bz[abs(ACE_t - (i+.5))  <= 0.5])

    # determine total avg B at hourly intervals
    ACE_hrB = np.sqrt(ACE_hrBx**2 + ACE_hrBy**2 + ACE_hrBz**2)
    errX = np.abs((FIDO_hrBx[np.where(ACE_hrBx!=0)] - ACE_hrBx[np.where(ACE_hrBx!=0)])) / np.mean(ACE_hrB[np.where(ACE_hrBx!=0)])
    errY = np.abs((FIDO_hrBy[np.where(ACE_hrBx!=0)] - ACE_hrBy[np.where(ACE_hrBx!=0)])) / np.mean(ACE_hrB[np.where(ACE_hrBx!=0)])
    errZ = np.abs((FIDO_hrBz[np.where(ACE_hrBx!=0)] - ACE_hrBz[np.where(ACE_hrBx!=0)])) / np.mean(ACE_hrB[np.where(ACE_hrBx!=0)])
    scoreBx = np.mean(errX)
    scoreBy = np.mean(errY)
    scoreBz = np.mean(errZ)   
    global totalscore
    totalscore = np.mean(np.sqrt(errX**2+errY**2+errZ**2))

    if tARR[-1] < CMEend - .5/24.: totalscore += 5.
    overamt = tARR[-1] - CMEend
    if overamt > 1/24.: totalscore += 0.1 * (24 * overamt - 1) 
    return scoreBx, scoreBy, scoreBz
        
def save_plot():
    if canprint: print('saving as '+my_name)
    plt.savefig(my_name+'.png')
    f1 = open(my_name+'.txt', 'w')
    if ISfilename != False: 
        f1.write('insitufile: '+ISfilename+' \n')
    f1.write('%-13s %8.2f \n' % ('CME_lat: ', CMElat))
    f1.write('%-13s %8.2f \n' % ('CME_lon: ', CMElon))
    f1.write('%-13s %8.2f \n' % ('CME_tilt: ', othertilt))
    f1.write('%-13s %8.2f \n' % ('CME_AW: ', CMEAW))
    f1.write('%-13s %8.2f \n' % ('CME_Ashape: ', CMESRA))
    f1.write('%-13s %8.2f \n' % ('CME_Bshape: ', CMESRB))
    f1.write('%-13s %8.2f \n' % ('CME_vr: ', CMEvr))
    f1.write('%-13s %8.2f \n' % ('CME_B0: ', CMEB0))
    f1.write('%-13s %8.2f \n' % ('CME_pol: ', CMEH))
    f1.write('%-13s %8.2f \n' % ('tshift: ', tshift))
    f1.write('%-13s %8.2f \n' % ('Sat_lat: ', FFlat))
    f1.write('%-13s %8.2f \n' % ('Sat_lon: ', FFlon0))
    f1.write('%-13s %8.2f \n' % ('CME_start: ', CMEstart))
    f1.write('%-13s %8.2f \n' % ('CME_stop: ', CMEend))
    f1.write('Launch_GUI: '+ str(Launch_GUI)+  '\n')
    f1.write('Autonormalize: '+ str(Autonormalize)+  '\n')
    f1.write('Save_Profile: '+ str(Save_Profile)+  '\n')
    f1.write('Expansion_Model: '+ expansion_model+  '\n')
    f1.write('Silent: '+ str(not canprint)+ '\n')
    f1.write('Indices: '+ str(show_indices))
    f1.close()
    print tshift/24.
    if Save_Profile == True:
            f1 = open(my_name+'.dat', 'w')
            if canprint: print('saving profile') 
            for i in range(len(obsBx)):
                f1.write('%10.5f %10.4f %10.4f %10.4f \n' % (tARR[i], obsBx[i], obsBy[i], obsBz[i]))
            f1.close()
    

def get_inputs(inputs):
    # take a file with unsorted input values and return a dictionary.
    # variable names have to match their names below
    possible_vars = ['insitufile', 'Sat_lat', 'Sat_lon', 'CME_lat', 'CME_lon', 'CME_tilt', 'CME_AW', 'CME_vr', 'tshift', 'CME_Ashape', 'CME_Bshape', 'CME_B0', 'CME_pol', 'CME_start', 'CME_stop', 'Autonormalize', 'Launch_GUI', 'Save_Profile', 'Expansion_Model', 'Silent', 'Indices']
    
    # if matches add to dictionary
    input_values = {}
    for i in range(len(inputs)):
        temp = inputs[i]
        if temp[0][:-1] in possible_vars:
            input_values[temp[0][:-1]] = temp[1]
    return input_values



# Parameters that define the simulation
global FFlat, FFlon0, CMElat, CMElon, CMEtilt, CMEAW, CMESRA, CMESRB, CMEvr, CMEB0, CMEH, tshift, CMEstart, CMEend

global tmax, dt
tmax = 80 * 3600. # maximum time of observations
dt = 1 * 60. # time between spacecraft obs

# useful global variables
global rsun, dtor, radeg, kmRs
rsun  =  7e10		 # convert to cm, 0.34 V374Peg
dtor  = 0.0174532925  # degrees to radians
radeg = 57.29577951    # radians to degrees
kmRs  = 1.0e5 / rsun # km (/s) divided by rsun (in cm)

global rotspeed
rotspeed = 1./ 3600. / 24 / 365. * 360 

# Get the CME number
global my_name

if len(sys.argv) < 2: 
    #sys.exit("Need an input file")
    print('No input file, running without in situ data and starting with defaults')
    input_values = []
    my_name = 'temp'
else:
    input_file = sys.argv[1]
    inputs = np.genfromtxt(input_file, dtype=str)
    my_name = input_file[:-4]
    input_values = get_inputs(inputs)

# Get in situ filename
global ISfilename, canScore
ISfilename = False
canScore = False
if 'insitufile' in input_values:
    if input_values['insitufile'] != 'NONE':
        ISfilename = input_values['insitufile']
        canScore = True

# Print things to command line?
global canprint 
canprint = True
if 'Silent' in input_values: 
    if input_values['Silent']=='True': canprint = False
if canprint: print('Files will be saved as '), my_name

# Pop up a GUI or just save a file
global Launch_GUI
Launch_GUI = False
if 'Launch_GUI' in input_values:
   if input_values['Launch_GUI'] == 'True': Launch_GUI = True
   
# Option to print a file with simulation results
global Save_Profile
Save_Profile = False
if 'Save_Profile' in input_values:
    if input_values['Save_Profile'] == 'True': Save_Profile = True
      
# Autonormalizing magnitude
global Autonormalize
Autonormalize = False
if 'Autonormalize' in input_values:
    if input_values['Autonormalize'] == 'True': Autonormalize = True
    
# Pick expansion mode
global expansion_model
expansion_model = 'None'
if 'Expansion_Model' in input_values:
    if input_values['Expansion_Model'] == 'None':
        expansion_model = 'None'
    elif input_values['Expansion_Model'] == 'Self-Similar':
        expansion_model = 'Self-Similar'        
    else:
        sys.exit('Expansion_Model should be None or Self-Similar')
        
# Determine if we are plotting kp index
global show_indices
show_indices = False
if 'Indices' in input_values:
    if input_values['Indices'] == 'True':
        show_indices = True

        
   
# Establish a GUI root if needed
if Launch_GUI == True:
    root = Tk()
    root.configure(background='gray75')
    
fig2 = plt.figure()
# set up the panels of the figure depending on whether we 
# want to show Kp or not
if show_indices == False:
    ax2  = fig2.add_subplot(411)
    ax3  = fig2.add_subplot(412, sharex=ax2)
    ax4  = fig2.add_subplot(413, sharex=ax2)
    ax5  = fig2.add_subplot(414, sharex=ax2)
else:
    ax2  = fig2.add_subplot(511)
    ax3  = fig2.add_subplot(512, sharex=ax2)
    ax4  = fig2.add_subplot(513, sharex=ax2)
    ax5  = fig2.add_subplot(514, sharex=ax2)
    ax6  = fig2.add_subplot(515, sharex=ax2)
    setp(ax5.get_xticklabels(), visible=False)
    
setp(ax2.get_xticklabels(), visible=False)
setp(ax3.get_xticklabels(), visible=False)
setp(ax4.get_xticklabels(), visible=False)
    
plt.subplots_adjust(right=0.8, wspace=0.001, hspace=0.0001)
plt.tight_layout()

# More GUI, complains if haven't def fig right above this so splitting root/canvas parts
if Launch_GUI:
    # Add fig to GUI canvas
    canvas = FigureCanvasTkAgg(fig2, master=root)
    canvas.get_tk_widget().grid(row=0, column=2, rowspan=30) #.grid(row=0,column=0)

    # CME parameters
    Label(root, text='CME Parameters', bg='gray75').grid(column=0,row=0, columnspan=2)
    Label(root, text='CME Lat (deg):', bg='gray75').grid(column=0, row=1)
    e1 = Entry(root, width=10)
    e1.grid(column=1,row=1)
    Label(root, text='CME Lon (deg):', bg='gray75').grid(column=0, row=2)
    e2 = Entry(root, width=10)
    e2.grid(column=1, row=2)
    Label(root, text='Tilt from W (deg):', bg='gray75').grid(column=0, row=3)
    e3 = Entry(root, width=10)
    e3.grid(column=1, row=3)
    Label(root, text='Angular Width (deg):', bg='gray75').grid(column=0, row=4)
    e3b = Entry(root, width=10)
    e3b.grid(column=1, row=4)
    Label(root, text='CME vr (km/s):', bg='gray75').grid(column=0, row=5)
    e8 = Entry(root, width=10)
    e8.grid(column=1, row=5)


    Label(root, text='Time Shift (hr):', bg='gray75').grid(column=0, row=16)
    e9 = Entry(root, width=10)
    e9.grid(column=1, row=16)


    # Torus Parameters
    Label(root, text="Torus Shape Parameters", bg='gray75').grid(column=0, row=6, columnspan=2)
    Label(root, text='A:', bg='gray75').grid(column=0, row=7)
    e4 = Entry(root, width=7)
    e4.grid(column=1, row=7)
    Label(root, text='B:', bg='gray75').grid(column=0, row=8)
    e5 = Entry(root, width=7)
    e5.grid(column=1, row=8)
    
    # Flux rope params
    Label(root, text='Force Free Parameters', bg='gray75').grid(column=0, row=11, columnspan=2)
    Label(root, text='B0:', bg='gray75').grid(column=0, row=12)
    e6 = Entry(root, width=10)
    e6.grid(column=1, row=12)
    Label(root, text='Pol. Direction:', bg='gray75').grid(column=0, row=13)
    e7 = Entry(root, width=10)
    e7.grid(column=1, row=13)
    

    # check button for autonormalizing magnitude
    global autonormVAR  
    autonormVAR = IntVar()
    autonormVAR.set(0)
    if Autonormalize: autonormVAR.set(1)
    Label(root, text='Autonormalize', bg='gray75').grid(column=3, row=0, columnspan=2)
    normCheck = Checkbutton(root, bg='gray75', var=autonormVAR).grid(column=3, row=1, columnspan=2)

    global expansionToggleVAR
    expansionToggleVAR = IntVar()
    Label(root, text='Expansion Profile:', bg='gray75').grid(row=2, column=3,columnspan=2)
    Label(root, text='Self-Similar', bg='gray75').grid(row=3, column=3,columnspan=1)
    Radiobutton(root, variable=expansionToggleVAR, value=0, bg='gray75').grid(column=3,row=4)
    Label(root, text='None', bg='gray75').grid(row=3, column=4,columnspan=1)
    Radiobutton(root, variable=expansionToggleVAR, value=1, bg='gray75').grid(column=4,row=4)
    expansionToggleVAR.set(1)
    if expansion_model == "Self-Similar": expansionToggleVAR.set(0)

    # CME start stop time
    Label(root, text='Observed CME Boundaries', bg='gray75').grid(column=3,row=5, columnspan=2)
    Label(root, text='Start (DOY)', bg='gray75').grid(column=3,row=6)
    eS1 = Entry(root, width=10)
    eS1.grid(column=4, row=6)
    Label(root, text='Stop (DOY)', bg='gray75').grid(column=3,row=7)
    eS2 = Entry(root, width=10)
    eS2.grid(column=4, row=7)


    # FIDO parameters
    Label(root, text='FIDO position', bg='gray75').grid(column=3,row=10, columnspan=2)
    Label(root, text='FIDO Lat:', bg='gray75').grid(column=3, row=11)
    eR1 = Entry(root, width=10)
    eR1.grid(column=4, row=11)
    Label(root, text='FIDO Lon:', bg='gray75').grid(column=3, row=12)
    eR2 = Entry(root, width=10)
    eR2.grid(column=4, row=12)

    Label(root, text='Update Plot', bg='gray75').grid(column=3,row=14, columnspan=2)
    draw_button = Button(root, command = update_plot, bg='gray75')
    draw_button.grid(row=15,column=3, columnspan=2)

    Label(root, text='Save Plot', bg='gray75').grid(column=3,row=16, columnspan=1)
    print_button = Button(root, bg='black', command = save_plot)
    print_button.grid(row=17,column=3, columnspan=1)

    Label(root, text='Quit', bg='gray75').grid(column=4,row=16, columnspan=1)
    quit_button = Button(root, command = root.quit)
    quit_button.grid(row=17, column=4, columnspan=2)

# Set sim params to default, replace with any given values
FFlat    = 0. 
FFlon0   = 0.
CMElat   = 0.
CMElon   = 0.
CMEtilt  = 0.
CMEAW    = 45. 
CMESRA   = 0.75
CMESRB   = 0.35
CMEvr    = 440.
CMEB0    = 25.
CMEH     = 1
tshift   = 0.
CMEstart = 0.
CMEend   = 0.

if 'Sat_lat' in input_values: FFlat = float(input_values['Sat_lat'])
if 'Sat_lon' in input_values: FFlon0 = float(input_values['Sat_lon'])
if 'CME_lat' in input_values: CMElat = float(input_values['CME_lat']) 
if 'CME_lon' in input_values: CMElon = float(input_values['CME_lon'])
if 'CME_tilt' in input_values: CMEtilt = float(input_values['CME_tilt']) 
if 'CME_AW' in input_values: CMEAW = float(input_values['CME_AW'])
if 'CME_Ashape' in input_values: CMESRA = float(input_values['CME_Ashape'])
if 'CME_Bshape' in input_values: CMESRB = float(input_values['CME_Bshape'])   
if 'CME_vr' in input_values:CMEvr = float(input_values['CME_vr'])
if 'tshift' in input_values: tshift = float(input_values['tshift'])
if 'CME_B0' in input_values: CMEB0 =float(input_values['CME_B0'])
if "CME_pol" in input_values: 
    if input_values['CME_pol'][0] == '-': CMEH = -1
if 'CME_start' in input_values: CMEstart =float(input_values['CME_start'])
if 'CME_stop' in input_values: CMEend = float(input_values['CME_stop'])

# see if we have IS and start/stop -> can calc score
# otherwise fix so will plot and flag not to score
global CMEmid, plotstart, plotend
global range_flag, scale_flag
pad = 3
scale_flag = False
range_flag = False
plotstart = CMEstart - pad/24.
if CMEend == 0:
    print('!!!Have CME start but not stop!!!')
    print('!!!Defaulting to duration of a day!!!')
    print('!!!May cause error if insufficient in situ data!!!')
    print('!!!Should not autonormalize or use score unless fix CME stop!!!')
    CMEend = CMEstart + 1
plotend   = CMEend + pad/24.
CMEmid    = 0.5 * (CMEstart + CMEend)

# insert values if using GUI, otherwise just set globals
if Launch_GUI:
    eR1.insert(0, FFlat)
    eR2.insert(0, FFlon0)
    e1.insert(0, CMElat) 
    e2.insert(0, CMElon)
    e3.insert(0, CMEtilt) 
    e3b.insert(0, CMEAW)
    e4.insert(0, CMESRA) 
    e5.insert(0, CMESRB)   
    e8.insert(0, CMEvr)
    e9.insert(0, tshift)
    e6.insert(0, CMEB0)
    e7.insert(0, CMEH)
    eS1.insert(0, CMEstart)
    eS2.insert(0, CMEend)




# read in data (coded to work with ACE)
global d_t, d_Btot, d_Bx, d_By, d_Bz, Wind_t, Wind_B, Wind_Bx, Wind_By, Wind_Bz

if (ISfilename != False):
    # Don't need skip_header for prepped daya
    data = np.genfromtxt(ISfilename, dtype=np.float, skip_header=44)
    if (range_flag==False):
        i_date = int(plotstart)
        i_hour = int(plotstart % 1 * 24)
        f_date = int(plotend) 
        f_hour = int(plotend % 1 % 1 * 24)
        # determine the initial and final index of the desired time period
        try:
            iidx = np.min(np.where(data[:,0] >= plotstart))
        except:
            sys.exit('CME_start outside in situ data range')
        try:
            fidx = np.max(np.where(data[:,0] <= plotend))
        except:
            sys.exit('CME_stop outside in situ data range')
    else:
        iidx = 0
        fidx = -2
    d_fracdays = data[iidx:fidx+1,0]
    d_Bx = data[iidx:fidx+1,1]
    d_By = data[iidx:fidx+1,2]
    d_Bz = data[iidx:fidx+1,3]
    d_Btot = np.sqrt(d_Bx**2 + d_By**2 + d_Bz**2)
    d_t = (d_fracdays - d_fracdays[0]) * 24
    d_tUN = d_fracdays
    global minISdate, maxISdate
    minISdate, maxISdate = np.min(d_tUN), np.max(d_tUN)
    # calculate the average magnetic field in the middle, used to normalize out B0
    global avg_obs_B
    if range_flag == False:
        avg_obs_B = np.mean(d_Btot[np.where(np.abs(d_tUN - CMEmid) < 2./24.)])
else:
    d_Btot, d_Bx, d_By, d_Bz = 0., 0., 0., 0.

# run the initial conditions
update_plot()

if Launch_GUI != False:
	root.mainloop()
	root.quit()
else:
	save_plot()

