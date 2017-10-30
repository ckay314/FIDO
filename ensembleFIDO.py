import numpy as np
import math
import sys
import os
import pickle
import dateit as DI
from scipy.special import jv
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from Tkinter import *
from pylab import setp
import random 
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
##from mpl_toolkits.axes_grid1 import make_axes_locatable
#from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

# These are things I am copying from ForeCAT and reusing--------------------------|
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

#---end of boring copied portion----------------------------------------------|


# New functions
def isinCME(vec_in, CME_shape):
	# Check and see if the requested point is actually in the CME and return
	# the cylindrical radial distance (from center of FR)
	# Function assumes vec_in is in CME Cartesian coord sys
    thetas = np.linspace(-math.pi/2, math.pi/2, 101)
    # determine the xz positions of the rope axis
    if CMEshapeToggle==0:
	# assume CME_shape = [rcenter, a, b, c]
        xFR = CME_shape[0] + CME_shape[1] * np.cos(thetas)
        zFR = CME_shape[3] * np.sin(thetas)
    if CMEshapeToggle==1:
    # assume shape is [h, b, rho, 0]
        bGCS   = CME_shape[1]
        rhoGCS = CME_shape[2]
        X0GCS  = (rhoGCS + bGCS * kappa**2*np.cos(thetas)) / (1-kappa**2)
        RGCS   = np.sqrt((bGCS**2*kappa**2 - rhoGCS**2)/(1-kappa**2) + X0GCS**2)
        #print np.mean(X0GCS**2), (bGCS**2*kappa**2 - rhoGCS**2)/(1-kappa**2)
        xFR    = bGCS + X0GCS*np.cos(thetas) 
        zFR    = X0GCS * np.sin(thetas)
    dists2 = (vec_in[0] - xFR)**2 + vec_in[1]**2 + (vec_in[2] - zFR)**2
    myb2 = np.min(dists2)
    minidxs = np.where(dists2 == myb2)
    minidx = minidxs[0]
    temp = thetas[np.where(dists2 == myb2)]
    mythetaT = temp[0]
	#for i in range(len(thetas)):
	#	print vec_in, xFR[i], zFR[i], thetas[i]*radeg, np.sqrt(dists2[i])
	# add a second iteration to refine B
	# see which side of minidx the actual minimum is on
    if minidx < len(dists2) - 1: # check to make sure not already at edge
        if dists2[minidx-1] < dists2[minidx+1]: startidx = minidx - 1
        else:  startidx = minidx + 1
		# repeat the same procedure on the side with the acutal min
        if dists2[minidx-1] != dists2[minidx+1]:
            thetas2 = np.linspace(thetas[startidx], thetas[minidx], 101)
            if CMEshapeToggle==0:
                xFR2 = CME_shape[0] + CME_shape[1] * np.cos(thetas2)
                zFR2 = CME_shape[3] * np.sin(thetas2)
            if CMEshapeToggle==1:
                X0GCS2  = (rhoGCS + bGCS * kappa**2*np.cos(thetas2)) / (1-kappa**2)
                RGCS2   = np.sqrt((bGCS**2*kappa**2 - rhoGCS**2)/(1-kappa**2) + X0GCS2**2)
                xFR2    = bGCS + X0GCS2*np.cos(thetas2) 
                zFR2    = X0GCS2 * np.sin(thetas2)

            dists2 = (vec_in[0] - xFR2)**2 + vec_in[1]**2 + (vec_in[2] - zFR2)**2
			#check if behind FR
			# need to do this
            myb2 = np.min(dists2)
            minidxs = np.where(dists2 == myb2)
            minidx = minidxs[0]
            temp = thetas2[np.where(dists2 == myb2)]
            mythetaT = temp[0]
			#print xFR2[minidx], zFR2[minidx]
    myb = np.sqrt(myb2)
    if CMEshapeToggle==0: CME_crossrad = CME_shape[2]
    if CMEshapeToggle==1: 
        if ('RGCS2' not in locals()) and ('RGCS2' not in globals()):
            RGCS2 = RGCS
        CME_crossrad = RGCS2[np.where(dists2 == myb2)]
    #print myb, CME_crossrad, CME_shape
    if (myb > CME_crossrad):
		#print 'Outside CME, distance is ', myb, CME_shape[2]
        myb = -9999.
        return myb, -9999, -9999., -9999, -9999
    else:
		# thetaP should swing quickly at closest approach to FR
        mythetaP = np.arcsin(vec_in[1] / myb)
        origTP = mythetaP
        if ('xFR2' not in locals()) and ('xFR2' not in globals()):
            xFR2 = xFR
        # had a random case that blew up here -> added try/except
        try: 
            if vec_in[0] < (xFR2[minidx] ):	
			    #print 'flipped it'
                if vec_in[1] > 0: mythetaP = math.pi - mythetaP
                else: mythetaP = -math.pi - mythetaP
        except:
            pass
	# return the distance, the toroidal, and the poloidal angle of the 
	#print myb, mythetaT*radeg, mythetaP * radeg, CME_shape[2] 
    return myb, mythetaT, mythetaP, 0, CME_crossrad

def getBvector(CME_shape, minb, thetaT, thetaP):
    if CMEshapeToggle==0:
	    tdir = np.array([-(CME_shape[1] + minb * np.cos(thetaP)) * np.sin(thetaT), 0., (CME_shape[3] + minb * np.cos(thetaP)) * np.cos(thetaT)])
	    pdir = np.array([-minb * np.sin(thetaP) * np.cos(thetaT), minb * np.cos(thetaP), -minb * np.sin(thetaP) * np.sin(thetaT)])
    if CMEshapeToggle==1:
        hGCS   = CME_shape[0]
        bGCS   = CME_shape[1]
        rhoGCS = CME_shape[2]
        sT = np.sin(thetaT)
        cT = np.cos(thetaT)
        sP = np.sin(thetaP)
        cP = np.cos(thetaP)
        X0GCS  = (rhoGCS + bGCS * kappa**2*cT) / (1-kappa**2)
        RGCS   = minb  #np.sqrt((bGCS**2*kappa**2 - rhoGCS**2)/(1-kappa**2) + X0GCS**2)
        dXdt   = - bGCS * kappa**2  * sT / (1 - kappa**2)
        dRdt   = X0GCS * dXdt / RGCS
        tdir = np.array([-sT*(X0GCS + RGCS * cP) + cT* (dXdt+ dRdt *cP), dRdt*sP, cT*(X0GCS + RGCS * cP) + sT* (dXdt+ dRdt *cP)])
        pdir = np.array([-RGCS * sP * cT, RGCS * cP, -RGCS * sP * sT])
    tmag = np.sqrt(np.sum(tdir**2))
    pmag = np.sqrt(np.sum(pdir**2))
    tdir = tdir / tmag
    pdir = pdir / pmag
    return tdir, pdir

def update_insitu():
    CME_shape = np.zeros(4)
    if CMEshapeToggle==0:
        # set up CME shape as [d, a, b, c]
        shapeC = np.tan(CMEAW*dtor) / (1. + CMESRB + np.tan(CMEAW*dtor) * (CMESRA + CMESRB)) 
        dtorang = (CMElon - FFlon0) / np.sin(CMEtilt * dtor) * dtor
        CMEnose = 215. / np.sqrt((1 - (CMESRA + CMESRB) * shapeC * (1. - np.cos(dtorang)))**2 + (1 + CMESRB)**2 * shapeC**2 * np.sin(dtorang)**2)  - 10.
        CME_shape[3] = CMEnose * shapeC
        CME_shape[1] = CME_shape[3] * CMESRA
        CME_shape[2] = CME_shape[3] * CMESRB 
        CME_shape[0] = CMEnose - CME_shape[1] - CME_shape[2]
    if CMEshapeToggle==1:
        # set up CME shape as [h, b, rho, <empty>]
        # alphaGCS is in rads, AW not
        # calculate smart starting point of nose - always seems to be near 215 so not really useful
        #dtorang = (CMElon - FFlon0) / np.sin(CMEtilt * dtor) * dtor
        #boverh = 1./ np.cos(alphaGCS)
        #Xoverh = 1./(1. - kappa**2) * (np.tan(alphaGCS) + kappa**2*np.cos(dtorang)/np.cos(alphaGCS))
        #Roverh = np.sqrt(1./(1-kappa**2)*(kappa**2/np.cos(alphaGCS)  - np.tan(alphaGCS)**2)+ Xoverh**2)
        #XplusR = Xoverh + Roverh
        #est_nose = 215. / np.sqrt(XplusR**2 + 1/np.cos(alphaGCS)**2 + 2*np.cos(dtorang)/np.cos(alphaGCS)*XplusR) 
        #print est_nose * (1. + np.sin(alphaGCS)) / ((1-kappa) * np.cos(alphaGCS))
        CMEnose=215 
        hGCS =  CMEnose* (1-kappa) * np.cos(alphaGCS) / (1. + np.sin(alphaGCS))
        bGCS = hGCS / np.cos(alphaGCS)
        rhoGCS = hGCS * np.tan(alphaGCS)  
        #print alphaGCS*radeg, hGCS, bGCS, rhoGCS
        # can't calc X0 or R here since depends on tor_angle
        CME_shape = [hGCS, bGCS, rhoGCS, 0.]
    obsBx = []
    obsBy = []
    obsBz = []
    tARR = []
    t = 0.
    CMEB = CMEB0
    FFlon = FFlon0
    thetaPprev = -42
    switch = 0
    #print 'in update insitu', CMElat, CMElon, CMEtilt, CMEAW, CMEB0, kappa
    while t < tmax:  
        tARR.append(t/3600.)
		# convert flyer position to CME Cartesian coordinates
		# have to do in loop to account for Earth orbit
		# get Sun xyz position
        FF_sunxyz = SPH2CART([215., FFlat, FFlon])
		# rotate to CMEcoord system
        temp = rotz(FF_sunxyz, -CMElon)
        temp2 = roty(temp, CMElat)
        FF_CMExyz = rotx(temp2, CMEtilt)
        #print FF_CMExyz
		# determine CME expansion and propagation
		# calculate CME shape
        if CMEshapeToggle==0:
            CME_shape[3] = CMEnose * np.tan(CMEAW*dtor) / (1. + CMESRB + np.tan(CMEAW*dtor) * (CMESRA + CMESRB))
            CME_shape[1] = CME_shape[3] * CMESRA
            CME_shape[2] = CME_shape[3] * CMESRB 
            CME_shape[0] = CMEnose - CME_shape[1] - CME_shape[2]
        if CMEshapeToggle==1: 
            CME_shape[0] = CMEnose* (1-kappa) * np.cos(alphaGCS) / (1. + np.sin(alphaGCS))
            CME_shape[1] = CME_shape[0] / np.cos(alphaGCS)
            CME_shape[2] = CME_shape[0] * np.tan(alphaGCS)  
		# check to see if the flyer is currently in the CME
		# if so, get its position in CME torus coords 
        minb, thetaT, thetaP, flagit, CME_crossrad = isinCME(FF_CMExyz, CME_shape)
		# need to check that thetaP doesn't jump as it occasionally will
        if flagit != -9999:
			# define it the first time its used
            if thetaPprev==-42: 
                thetaPprev = thetaP
                #print 'init cross section radius: ', CME_crossrad
                #print CME_shape
            delThetaP = np.abs(thetaP - thetaPprev)
 			#print thetaPprev, thetaP, delThetaP, delThetaP>0.5
            if delThetaP > 0.5: thetaP = thetaPprev
            thetaPprev = thetaP
			#print minb, CME_shape, FF_CMExyz, print temp2
		# get the toroidal and poloidal magnetic field
        if fluxRopeToggle==0:
            Btor = CMEB * jv(0, 2.4 * minb / CME_crossrad)
            Bpol = CMEB * CMEH * jv(1, 2.4 * minb / CME_crossrad)
        if fluxRopeToggle==1:
			Btor = CMEB * (cTau - (minb/CME_crossrad)**(1+cN))
			Bpol = -CMEB * (minb/CME_crossrad)**(1+cM) / cC1 * (-np.sign(CMEH)) * ((cN+1.)/(cM+2.))     
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
			#if switch==0:  # use this to check values of initial point in MC
			#	print CME_shape, minb, thetaT, thetaP
			#	switch = 1
        else:
            obsBx.append(0.)
            obsBy.append(0.)
            obsBz.append(0.)
        #print t/3600, minb, CME_shape[2], CMEnose, CME_shape[1]
        t += dt
		# CME nose moves to new position
        CMEnose += CMEvr * dt / 7e5 # change position in Rsun
		# update the total magnetic field
        CMEB *= ((CMEnose - CMEvr * dt / 7e5) / CMEnose)**2
		# determine new lon of observer
        FFlon += dt / 3600. / 24 / 365. * 360 
    #print 'through CME'
    return obsBx, obsBy, obsBz, tARR

def update_plot():
    global CMElat, CMElon, CMEtilt, CMEAW, CMESRA, CMESRB, CMEB0, CMEH, CMEvr, tshift
    # get the value of the radio buttons
    global fluxRopeToggle, CMEshapeToggle
    fluxRopeToggle = fluxRopeToggleVAR.get()
    CMEshapeToggle = CMEshapeToggleVAR.get()
    CMElat = float(e1.get())
    CMElon = float(e2.get())
    CMEtilt = float(e3.get())
    CMEAW   = float(e3b.get())
    if CMEshapeToggle==0:
        CMESRA  = float(e4.get())
        CMESRB  = float(e5.get())
    global alphaGCS, kappa
    if CMEshapeToggle==1:
        kappa   = float(eGCS.get())
        alphaGCS = CMEAW *dtor - np.arcsin(kappa)
        #b_unit = (1. - kappa)/(1+np.sin(alphaGCS))
        #p_unit = np.sin(alphaGCS) * (1. - kappa)/(1+np.sin(alphaGCS))
        #torAR = np.arctan(p_unit/(1-kappa**2) + kappa*b_unit/np.sqrt(1-kappa**2))
        #print alphaGCS*radeg, CMEAW, torAR * radeg
        if alphaGCS < 0:
            print 'alpha less than 0, increase AW or decrease kappa'
            CMElat = 99.
        #print alphaGCS
    if fluxRopeToggle==0:
        CMEB0   = float(e6.get())
        CMEH	= float(e7.get())
    CMEvr   = float(e8.get())
    tshift  = float(e9.get())
    global FFlat, FFlon0
    FFlat   = float(eR1.get())
    FFlon0  = float(eR2.get()) 
    global cTau, cC1, cM, cN
    cC1 = 0
    if fluxRopeToggle==1:
        CMEB0     = float(ec1.get()) 
        cM      = float(ec2.get())
        cN      = float(ec3.get())
        cTau    = float(ec4.get())
        cC1     = float(ec5.get())
    print CMElat, CMElon, CMEtilt, CMEAW, CMESRA, CMESRB, CMEB0, CMEH, CMEvr, tshift, cC1
    # define as globals so can use to calculate score
    global obsBx, obsBy, obsBz, obsB, tARR
    obsBx, obsBy, obsBz, tARR = update_insitu()
    obsBx, obsBy, obsBz = np.array(obsBx), np.array(obsBy), np.array(obsBz)
    tARR = np.array(tARR) 
    tshift = CMEstart + tshift/24.
    tARR = tARR/24.
    
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
			#print len(obsBx), idx2, MCidxf
            obsBx = obsBx[:idx2] 
            obsBy = obsBy[:idx2] 
            obsBz = obsBz[:idx2] 
            tARR  = tARR[:idx2]
    obsB = np.sqrt(obsBx**2 + obsBy**2 + obsBz**2)

    if autofitT == True:
        oldlength = np.max(tARR)
        tARR = obsCMEtime / oldlength * tARR

    # scale the CME to match at midpoint (ignoring B0 with this)
    avg_FIDO_B = np.mean(obsB[np.where(np.abs(tshift+tARR - CMEmid) < 2./24.)])
    scale = 1.
    if autonormVAR.get()==1: scale = avg_obs_B / avg_FIDO_B
    #print 'magnitude', scale*avg_obs_B, scale*CMEB0
    obsBx *= scale
    obsBy *= scale
    obsBz *= scale
    obsB = np.sqrt(obsBx**2 + obsBy**2 + obsBz**2)
    
    scoreBx, scoreBy, scoreBz = calc_score()
    print 'initial score', totalscore
    print 'Bx, By, Bz',  scoreBx, scoreBy, scoreBz
    print ''
      
    #ax2.clear()
    #ax3.clear()
    #ax4.clear()
    #ax5.clear()
    #ax2.plot(d_t, d_Btot, 'k', linewidth=3)
    maxBtot = np.max([np.max(d_Btot), np.max(obsB)])
    ax2.plot([RCstart, RCstart], [0., 1.3*maxBtot], 'k--', linewidth=2)
    ax2.plot([RCend, RCend], [0., 1.3*maxBtot], 'k--', linewidth=2)
    if flag_it == False:
        ax2.plot([TNCstart, TNCstart], [0., 1.3*np.max(maxBtot)], '--', color='dimgrey', linewidth=2)
        ax2.plot([TNCend, TNCend], [0., 1.3*np.max(maxBtot)], '--', color='dimgrey',linewidth=2)
        ax2.plot(Wind_t, Wind_B, color='dimgrey', linewidth=4, zorder=0)
    ax2.plot(d_tUN, d_Btot, 'k', linewidth=4, zorder=1)
    ax2.plot(tshift + tARR, obsB, color=CMEcolor, linewidth=4, zorder=CMEnumber)
    ax2.set_ylabel('B (nT)')
    setp(ax2.get_xticklabels(), visible=False)

    minBx = np.abs(np.min([np.min(d_Bx), np.min(-obsBx)]))
    maxBx = np.max([np.max(d_Bx), np.max(-obsBx)])
    ax3.plot([RCstart, RCstart], [-1.3*minBx, 1.3*maxBx], 'k--', linewidth=2)
    ax3.plot([RCend, RCend], [-1.3*minBx, 1.3*maxBx], 'k--', linewidth=2)
    if flag_it == False:
        ax3.plot([TNCstart, TNCstart], [-1.3*minBx, 1.3*maxBx], '--', color='dimgrey',linewidth=2)
        ax3.plot([TNCend, TNCend], [-1.3*minBx, 1.3*maxBx], '--', color='dimgrey', linewidth=2)
        ax3.plot(Wind_t, Wind_Bx, color='dimgrey', linewidth=4, zorder=0)
    ax3.plot(d_tUN, d_Bx, 'k', linewidth=4, zorder=1)
    ax3.plot(tshift + tARR, -obsBx, color=CMEcolor, linewidth=4, zorder=CMEnumber)
    #ax3.annotate('%0.2f'%(scoreBx), xy=(1, 0), color='r', xycoords='axes fraction', fontsize=16,
                #horizontalalignment='right', verticalalignment='bottom')    
    ax3.set_ylabel('B$_x$ (nT)')
    setp(ax3.get_xticklabels(), visible=False)

    minBy = np.abs(np.min([np.min(d_By), np.min(-obsBy)]))
    maxBy = np.max([np.max(d_By), np.max(-obsBy)])
    ax4.plot([RCstart, RCstart], [-1.3*minBy, 1.3*maxBy], 'k--', linewidth=2)
    ax4.plot([RCend, RCend], [-1.3*minBy, 1.3*maxBy], 'k--', linewidth=2)
    if flag_it == False:
        ax4.plot([TNCstart, TNCstart], [-1.3*minBy, 1.3*maxBy], '--', color='dimgrey', linewidth=2)
        ax4.plot([TNCend, TNCend], [-1.3*minBy, 1.3*maxBy], '--', color='dimgrey', linewidth=2)
        ax4.plot(Wind_t, Wind_By, color='dimgrey', linewidth=4, zorder=0)
    ax4.plot(d_tUN, d_By, 'k', linewidth=4, zorder=1)
    #ax4.annotate('%0.2f'%(scoreBy), xy=(1, 0), color='r', xycoords='axes fraction', fontsize=16,
                #horizontalalignment='right', verticalalignment='bottom')    
    ax4.set_ylabel('B$_y$ (nT)')
    ax4.plot(tshift + tARR, -obsBy, color=CMEcolor, linewidth=4, zorder=CMEnumber)
    setp(ax4.get_xticklabels(), visible=False)

    minBz = np.abs(np.min([np.min(d_Bz), np.min(obsBz)]))
    maxBz = np.max([np.max(d_Bz), np.max(obsBz)])
    ax5.plot([RCstart, RCstart], [-1.3*minBz, 1.3*maxBz], 'k--', linewidth=2)
    ax5.plot([RCend, RCend], [-1.3*minBz, 1.3*maxBz], 'k--', linewidth=2)
    if flag_it == False:
        ax5.plot([TNCstart, TNCstart], [-1.3*minBz, 1.3*maxBz], '--', color='dimgrey',linewidth=2)
        ax5.plot([TNCend, TNCend], [-1.3*minBz, 1.3*maxBz], '--', color='dimgrey',linewidth=2)
        ax5.plot(Wind_t, Wind_Bz, color='dimgrey', linewidth=4, zorder=0)
    ax5.plot(d_tUN, d_Bz, 'k', linewidth=4, zorder=1)
    #ax5.annotate('%0.2f'%(scoreBz), xy=(1, 0), color='r', xycoords='axes fraction', fontsize=16,
                #horizontalalignment='right', verticalalignment='bottom')    
    ax5.set_ylabel('B$_z$ (nT)')

    ax5.set_xlabel('Day of Year')
    ax5.plot(tshift + tARR, obsBz, color=CMEcolor, linewidth=4, zorder=CMEnumber)
    scl = 1.25
    ax2.set_ylim([0,scl*maxBtot])
    ax3.set_ylim([-1.*scl*minBx, scl*maxBx])
    ax4.set_ylim([-1.*scl*minBy, scl*maxBy])
    ax5.set_ylim([-1.*scl*minBz, scl*maxBz])
    ax2.set_xlim([plotstart, plotend])
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.subplots_adjust(right=0.8, wspace=0.001, hspace=0.25)
    fig2.suptitle(datestr, fontsize=16)
    plt.tight_layout()
    plt.subplots_adjust(top=0.94)
    plt.gcf().canvas.draw()

def calc_score():
    # calculate the score for Bx By and Bz with ACE
    #maxt      = int(np.max(tARR*24))-1
    maxt = int(24*(CMEend-CMEstart))# -1
    #print maxt
    FIDO_hrt  = np.zeros(maxt)
    FIDO_hrBx = np.zeros(maxt)
    FIDO_hrBy = np.zeros(maxt)
    FIDO_hrBz = np.zeros(maxt)
    ACE_hrBx = np.zeros(maxt)
    ACE_hrBy = np.zeros(maxt)
    ACE_hrBz = np.zeros(maxt)
    
    # set to start checking at beginning or end of obs CME
    #global CMEstart, CMEend
    
    #FIDO_hr  = (tARR + tshift) * 24 
    #print CMEstart, tshift, d_tUN[0]
    ACE_t = (d_tUN - CMEstart) * 24
    FIDO_hr = (tARR + tshift - CMEstart) *24
    #print CMEend, CMEstart, maxt
    #ACE_t    = d_tUN *24 - FIDO_hr[0]  # normalize so FIDO_hr[0] = 0
    #FIDO_hr -= FIDO_hr[0]
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
    global run_results
    run_results = [FIDO_hrt, FIDO_hrBx, FIDO_hrBy, FIDO_hrBz]

    # determine total avg B at hourly intervals
    # determine total avg B at hourly intervals
    ACE_hrB = np.sqrt(ACE_hrBx**2 + ACE_hrBy**2 + ACE_hrBz**2)
    errX = np.abs((FIDO_hrBx[np.where(ACE_hrBx!=0)] - ACE_hrBx[np.where(ACE_hrBx!=0)])) / np.mean(ACE_hrB[np.where(ACE_hrBx!=0)])
    errY = np.abs((FIDO_hrBy[np.where(ACE_hrBx!=0)] - ACE_hrBy[np.where(ACE_hrBx!=0)])) / np.mean(ACE_hrB[np.where(ACE_hrBx!=0)])
    errZ = np.abs((FIDO_hrBz[np.where(ACE_hrBx!=0)] - ACE_hrBz[np.where(ACE_hrBx!=0)])) / np.mean(ACE_hrB[np.where(ACE_hrBx!=0)])
    global scoreBx, scoreBy, scoreBz
    scoreBx = np.mean(errX)
    scoreBy = np.mean(errY)
    scoreBz = np.mean(errZ)   
    global totalscore
    totalscore = np.mean(np.sqrt(errX**2+errY**2+errZ**2))
    
    if np.mean(np.abs(FIDO_hrBx)) < 0.0001:
        print 'missed'
        scoreBx, scoreBy, scoreBz = 9.99, 9.99, 9.99
        totalscore = 99.
    # had a random case that blew up here -> added try/except
    try:
        # add to total score if CME is too short by more than 30 mins
        if tARR[-1] + tshift < CMEend - .5/24.: totalscore += 5.
        # penalize if over an hour too long
        #if tARR[-1] + tshift > CMEend + 1./24.: totalscore += 5.
        # penalize proportional to time over beyond an hour
        overamt = tARR[-1] + tshift - CMEend
        #print overamt
        if overamt > 1/24.: totalscore += 0.1 * (24 * overamt - 1) 
        
    except:
        totalscore += 10.  # return bad score if fails
    return scoreBx, scoreBy, scoreBz
    
def find_bf():
    prevCMElat = float(e1.get())
    prevCMElon = float(e2.get())
    prevCMEtilt = float(e3.get())
    prevCMEAW   = float(e3b.get())
    #prevCMEB0   = float(e6.get()) 
    prevCMESRA  = float(e4.get())
    prevCMESRB  = float(e5.get())
    #prevkappa   = float(eGCS.get())
    prevC1     = float(ec5.get())
    print prevCMElat, prevCMElon, prevCMEtilt, prevCMEAW, prevCMESRA, prevCMESRB, prevC1
    
    for bf_counter in range(2000):
        global totalscore
        oldscore = totalscore
    
        global CMElat, CMElon, CMEtilt, CMEAW, CMEB0, CMESRA, CMESRB, kappa, cC1      
        newCMElat = prevCMElat + 0.1 * (random.random() - 0.5)
        newCMElat = np.max([latrange[0], newCMElat])
        newCMElat = np.min([latrange[1], newCMElat])
        CMElat = newCMElat
        
        newCMElon = prevCMElon + 0.1 * (random.random() - 0.5)
        newCMElon = np.max([lonrange[0], newCMElon])
        newCMElon = np.min([lonrange[1], newCMElon])
        CMElon = newCMElon
        
        newCMEtilt = prevCMEtilt + 0.1 * (random.random() - 0.5)
        newCMEtilt = np.max([tiltrange[0], newCMEtilt])
        newCMEtilt = np.min([tiltrange[1], newCMEtilt])
        CMEtilt = newCMEtilt

        newCMEAW = prevCMEAW + 0.1 * (random.random() - 0.5)
        newCMEAW = np.max([AWrange[0], newCMEAW])
        newCMEAW = np.min([AWrange[1], newCMEAW])
        CMEAW = newCMEAW
 
        #newCMEB0 = prevCMEB0 + 0.25 * (random.random() - 0.5)
        #newCMEB0 = np.max([B0range[0], newCMEB0])
        #newCMEB0 = np.min([B0range[1], newCMEB0])
        #CMEB0 = newCMEB0

        newCMESRA = prevCMESRA + 0.005 * (random.random() - 0.5)
        newCMESRA = np.max([Arange[0], newCMESRA])
        newCMESRA = np.min([Arange[1], newCMESRA])
        CMESRA = newCMESRA

        newCMESRB = prevCMESRB + 0.001 * (random.random() - 0.5)
        newCMESRB = np.max([Brange[0], newCMESRB])
        newCMESRB = np.min([Brange[1], newCMESRB])
        CMESRB = newCMESRB
        
        #newkappica = prevkappa + 0.005 * (random.random() - 0.5)
        #newkappica = np.max([kapparange[0], newkappica])
        #newkappica = np.min([kapparange[1], newkappica])
        #kappa = newkappica
        #global alphaGCS
        #alphaGCS = CMEAW *dtor - np.arcsin(kappa)
        
        newC1 = prevC1 + 0.01 * (random.random() - 0.5)
        newC1 = np.max([C1range[0], newC1])
        newC1 = np.min([C1range[1], newC1])
        cC1 = newC1
        
        
     
        global obsBx, obsBy, obsBz, obsB, tARR
        obsBx, obsBy, obsBz, tARR = update_insitu()
        obsBx, obsBy, obsBz = np.array(obsBx), np.array(obsBy), np.array(obsBz)
        tARR = np.array(tARR) / 24.

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
			    #print len(obsBx), idx2, MCidxf
                obsBx = obsBx[:idx2] 
                obsBy = obsBy[:idx2] 
                obsBz = obsBz[:idx2] 
                tARR  = tARR[:idx2]
        obsB = np.sqrt(obsBx**2 + obsBy**2 + obsBz**2)

        # scale the CME to match at midpoint (ignoring B0 with this)
        avg_FIDO_B = np.mean(obsB[np.where(np.abs(tshift+tARR - CMEmid) < 2./24.)])
        scale = 1.
        if autonormVAR.get()==1: scale = avg_obs_B / avg_FIDO_B
        obsBx *= scale
        obsBy *= scale
        obsBz *= scale
        obsB = np.sqrt(obsBx**2 + obsBy**2 + obsBz**2)

        scoreBx, scoreBy, scoreBz = calc_score()
        
        print bf_counter, oldscore, totalscore, CMElat, CMElon, CMEtilt, CMEAW, CMESRA, CMESRB
        

        if totalscore > 1.001 * oldscore:
            print 'rejected'
            totalscore = oldscore
        elif np.isnan(totalscore):
            print 'rejected'
            totalscore = oldscore            
        else:
            prevCMElat  = CMElat
            prevCMElon  = CMElon
            prevCMEtilt = CMEtilt
            prevCMEAW   = CMEAW
           # prevCMEB0   = CMEB0
            prevCMESRA  = CMESRA
            prevCMESRB  = CMESRB
            #prevkappa   = kappa
            prevC1      = cC1
        #print ''
 
    e1.delete(0, END)
    e1.insert(0, prevCMElat)
    e2.delete(0, END)
    e2.insert(0, prevCMElon)
    e3.delete(0, END)
    e3.insert(0, prevCMEtilt)
    e3b.delete(0, END)
    e3b.insert(0, prevCMEAW)
    #e6.delete(0, END)
    #e6.insert(0, prevCMEB0)
    e4.delete(0, END)
    e4.insert(0, prevCMESRA)
    e5.delete(0, END)
    e5.insert(0, prevCMESRB)
    #eGCS.delete(0, END)
    #eGCS.insert(0, prevkappa)
    #ec1.delete(0, END)
    #ec1.insert(0, prevCMEB0)
    ec5.delete(0, END)
    ec5.insert(0, prevC1)

    update_plot()
    print prevCMElat, prevCMElon, prevCMEtilt, prevCMEAW, prevCMESRA, prevCMESRB, prevC1
    
def save_plot():
	mystrid = str(myid+1)
	print 'saving FIDOCME'+mystrid
	plt.savefig('FIDOCME'+mystrid+'.png')
	f1 = open('FIDOCME'+mystrid+'.txt', 'w')

	f1.write('%-13s %8.2f \n' % ('CME_lat', CMElat))
	f1.write('%-13s %8.2f \n' % ('CME_lon', CMElon))
	f1.write('%-13s %8.2f \n' % ('CME_tilt', CMEtilt))
	f1.write('%-13s %8.2f \n' % ('CME_AW', CMEAW))
	f1.write('%-13s %8.2f \n' % ('CME_SRA', CMESRA))
	f1.write('%-13s %8.2f \n' % ('CME_SRB', CMESRB))
	f1.write('%-13s %8.2f \n' % ('CME_vr', CMEvr))
	f1.write('%-13s %8.2f \n' % ('CME_Btor', CMEB0))
	f1.write('%-13s %8.2f \n' % ('CME_bpol', CMEH))
	tadjust = float(e9.get())
	f1.write('%-13s %8.2f \n' % ('t_shift', tadjust))

	f1.write('%-13s %8.2f \n' % ('FIDO_lat', FFlat))
	f1.write('%-13s %8.2f \n' % ('FIDO_lon', FFlon0))
	f1.close()

random.seed(42)

# Parameters to play with (like string for a cat)
global FFlat, FFlon0, CMElat, CMElon, CMEtilt, CMEAW, CMESRA, CMESRB, CMEvr, CMEB0, CMEH

global tmax, dt
tmax = 80 * 3600. # maximum time of observations
dt = 1 * 60. # time between spacecraft obs

# useful global variables
global rsun, dtor, radeg, kmRs
rsun  =  7e10		 # convert to cm, 0.34 V374Peg
dtor  = 0.0174532925  # degrees to radians
radeg = 57.29577951    # radians to degrees
kmRs  = 1.0e5 / rsun # km (/s) divided by rsun (in cm)


# Get the CME number
global myid
myid =int(sys.argv[1]) - 1

# open the corresponding CME file
mynumSTR = str(myid+1)
if len(mynumSTR) < 2: mynumSTR='0'+mynumSTR
myfile = 'CME'+mynumSTR+'.txt'

# set up the GUI
root = Tk()
root.configure(background='gray75')
fig2 = plt.figure()
canvas = FigureCanvasTkAgg(fig2, master=root)
canvas.get_tk_widget().grid(row=0, column=2, rowspan=30) #.grid(row=0,column=0)

# set up the panels of the figure
ax2  = fig2.add_subplot(411)
ax2.set_ylabel('B (G)')
setp(ax2.get_xticklabels(), visible=False)
ax3  = fig2.add_subplot(412, sharex=ax2)
ax3.set_ylabel('B$_x$ (G)')
setp(ax3.get_xticklabels(), visible=False)
ax4  = fig2.add_subplot(413, sharex=ax2)
ax4.set_ylabel('B$_y$ (G)')
setp(ax4.get_xticklabels(), visible=False)
ax5  = fig2.add_subplot(414, sharex=ax2)
ax5.set_ylabel('B$_z$ (G)')
ax5.set_xlabel('Time (hours)')
plt.subplots_adjust(right=0.8, wspace=0.001, hspace=0.001)
plt.tight_layout()


# CME parameters
Label(root, text='CME Parameters', bg='gray75').grid(column=0,row=0, columnspan=2)
Label(root, text='CME Lat:', bg='gray75').grid(column=0, row=1)
e1 = Entry(root, width=10)
e1.grid(column=1,row=1)
Label(root, text='CME Lon:', bg='gray75').grid(column=0, row=2)
e2 = Entry(root, width=10)
e2.grid(column=1, row=2)
Label(root, text='Tilt from N (Deg):', bg='gray75').grid(column=0, row=3)
e3 = Entry(root, width=10)
e3.grid(column=1, row=3)
Label(root, text='Angular Width (Deg):', bg='gray75').grid(column=0, row=4)
e3b = Entry(root, width=10)
e3b.grid(column=1, row=4)
Label(root, text='CME vr (km/s):', bg='gray75').grid(column=0, row=5)
e8 = Entry(root, width=10)
e8.grid(column=1, row=5)
Label(root, text='Time Shift:', bg='gray75').grid(column=0, row=6)
e9 = Entry(root, width=10)
e9.grid(column=1, row=6)


# radio button to choose between torus and GCS
global CMEshapeToggleVAR
CMEshapeToggleVAR = IntVar()
CMEshapeToggleVAR.set(0)
Label(root, text='Flux Rope Shape:', bg='gray75').grid(row=7, column=0,columnspan=2)
Radiobutton(root, text='Torus', variable=CMEshapeToggleVAR, value=0, bg='gray75').grid(column=0,row=8)
Radiobutton(root, text='GCS', variable=CMEshapeToggleVAR, value=1, bg='gray75').grid(column=1,row=8)

# Torus Parameters
Label(root, text="Torus Parameters", bg='gray75').grid(column=0, row=9, columnspan=2)
Label(root, text='A:', bg='gray75').grid(column=0, row=10)
e4 = Entry(root, width=10)
e4.grid(column=1, row=10)
Label(root, text='B:', bg='gray75').grid(column=0, row=11)
e5 = Entry(root, width=10)
e5.grid(column=1, row=11)

# GCS parameter
Label(root, text="GCS Parameter", bg='gray75').grid(column=0, row=12, columnspan=2)
Label(root, text='kappa:', bg='gray75').grid(column=0, row=13)
eGCS = Entry(root, width=10)
eGCS.grid(column=1, row=13)

# Flux rope parameters
# radio button to choose between force-free and circular
global fluxRopeToggleVAR
fluxRopeToggleVAR = IntVar()
fluxRopeToggleVAR.set(0)
Label(root, text='Flux Rope Model:', bg='gray75').grid(row=0, column=3,columnspan=2)
Radiobutton(root, text='Force Free', variable=fluxRopeToggleVAR, value=0, bg='gray75').grid(column=3,row=1)
Radiobutton(root, text='Circular', variable=fluxRopeToggleVAR, value=1, bg='gray75').grid(column=4,row=1)

# check button for autonormalizing magnitude
global autonormVAR
autonormVAR = IntVar()
autonormVAR.set(1)
normCheck = Checkbutton(root, text='Auto Normalize', bg='gray75', var=autonormVAR).grid(column=3, row=2)

Label(root, text='Force Free Parameters', bg='gray75').grid(column=3, row=3, columnspan=2)
Label(root, text='B0:', bg='gray75').grid(column=3, row=4)
e6 = Entry(root, width=10)
e6.grid(column=4, row=4)
Label(root, text='Pol. Direction:', bg='gray75').grid(column=3, row=5)
e7 = Entry(root, width=10)
e7.grid(column=4, row=5)
Label(root, text='Circular Parameters', bg='gray75').grid(column=3, row=6, columnspan=2)
Label(root, text='B0:', bg='gray75').grid(column=3, row=7)
ec1 = Entry(root, width=10)
ec1.grid(column=4, row=7)
Label(root, text='m:', bg='gray75').grid(column=3, row=8)
ec2 = Entry(root, width=10)
ec2.grid(column=4, row=8)
Label(root, text='n:', bg='gray75').grid(column=3, row=9)
ec3 = Entry(root, width=10)
ec3.grid(column=4, row=9)
Label(root, text='tau:', bg='gray75').grid(column=3, row=10)
ec4 = Entry(root, width=10)
ec4.grid(column=4, row=10)
Label(root, text='C1:', bg='gray75').grid(column=3, row=11)
ec5 = Entry(root, width=10)
ec5.grid(column=4, row=11)

# FIDO parameters
Label(root, text='FIDO position', bg='gray75').grid(column=3,row=12)
Label(root, text='FIDO Lat:', bg='gray75').grid(column=3, row=13)
eR1 = Entry(root, width=10)
eR1.grid(column=4, row=13)
Label(root, text='FIDO Lon:', bg='gray75').grid(column=3, row=14)
eR2 = Entry(root, width=10)
eR2.grid(column=4, row=14)

print_button = Button(root, text="Save Plot", command = save_plot)
print_button.grid(row=31,column=1)
draw_button = Button(root, text="Update Plot", command = update_plot)
draw_button.grid(row=31,column=2)
BF_button = Button(root, text="Find Best Fit", command = find_bf)
BF_button.grid(row=32,column=2)
quit_button = Button(root, text="Quit", command = root.quit, fg='red')
quit_button.grid(row=31, column=3)


# preload last variables if file exists
# Load in CME data
data = np.genfromtxt(myfile, dtype=None)
datestr = data[0][1]
datestr = datestr[:10] + ' ' + datestr[11:]
print 'datestr', datestr
# make a new datestr that actually has word months
month = int(datestr[5:7])
yr = datestr[0:4]
day = datestr[8:10]
month_dict = {1:' January ', 2:' February ', 3:' March ', 4:' April ', 5:' May ', 6:' June ', 7:' July ', 8:' August ', 9:' September ', 10:' October ', 11:' November ', 12: ' December '}
datestr = day+month_dict[month]+yr + ' ' + datestr[11:]


CMEtime = data[1][1] # get DOY from the year frac 
print 'CME time', CMEtime

# insert values 
e1.insert(0, data[52][1]) # latitude
e2.insert(0, data[53][1]) # longitude
e3.insert(0, data[54][1]) # tilt
e3b.insert(0, data[55][1]) # AW

e4.insert(0, data[56][1])    # A
e5.insert(0, data[57][1])    # B
# add alternative for kappa

e8.insert(0, data[49][1])   # vcme
e6.insert(0, data[62][1])   # B0
#\e6.insert(0, 0.005)   # B0
used_pol = data[8][1]
global CMEH
CMEH = 1
if used_pol[1] == '-': CMEH = -1
e7.insert(0, CMEH) # H
    
e9.insert(0, 0) # tshift
eR1.insert(0, data[2][1])  #FIDO lat
eR2.insert(0, data[50][1])  #FIDO lon

ec1.insert(0,data[62][1])  #B0
ec2.insert(0,0)         #m
ec3.insert(0,1)         #n
ec4.insert(0,1)         #tau
ec5.insert(0,1)         #C1



global latrange, lonrange, tiltrange, AWrange, Arange, Brange, B0range, tshiftrange, C1range

latrange = [float(data[52][1])-5, float(data[52][1])+5.]
lonrange = [float(data[53][1])-5., float(data[53][1])+5.]
tiltrange = [float(data[54][1])-10., float(data[54][1])+10.]
# might have dif AW for tor/GCS so pick whatever currently used
the_AW = float(e3b.get())
AWrange  = [the_AW-10., the_AW+10.]
Brange = [0.01, 0.99]
Arange = [0.4, 1.0]
#kapparange = [np.min([myrow2[9]-0.25, 0.01]), myrow2[9]+0.25]
# opposite sign of force free H
C1range = [0.5, 2.0]

# get the in situ information
global RCstart, RCend, TNCstart, TNCend, plotstart, plotend
RCstart = float(data[43][1])
RCend   = float(data[44][1])

TNCstart = float(data[45][1])
TNCend = float(data[46][1])

flag_it = False
if TNCstart < 0:
	TNCstart = RCstart
	TNCend   = RCend
	flag_it  = True
global CMEstart, CMEend, CMEmid
CMEstart  = float(data[47][1])
CMEend    = float(data[48][1])
CMEmid    = 0.5 * (CMEstart + CMEend)
print 'RC start, end: ',  RCstart, RCend
print 'TNC start, end: ',  TNCstart, TNCend
print "start, mid, end:  ", CMEstart, CMEmid, CMEend
pad = 3
plotstart = np.min([TNCstart, RCstart]) - pad/24.
plotend   = np.max([TNCend, RCend]) + pad/24.
print plotstart, plotend
global obsCMEtime
obsCMEtime = CMEend - CMEstart

# read in ACE data
global d_t, d_Btot, d_Bx, d_By, d_Bz, Wind_t, Wind_B, Wind_Bx, Wind_By, Wind_Bz
filename = 'ACE_CME'+ str(myid+1) +'.dat'
i_date = int(plotstart)
i_hour = int(plotstart % 1 * 24)
f_date = int(plotend) 
f_hour = int(plotend % 1 % 1 * 24)
print i_date, i_hour, f_date, f_hour
data = np.genfromtxt(filename, skip_header = 44, dtype=np.float)
d_yr = data[:,0]
yrlen = 365
if (d_yr[0] % 4 == 0): yrlen = 366
d_days = data[:,1]
d_hours  = data[:,2]
# determine the initial and final index of the desired time period
iidx = np.min(np.where((d_days == i_date) & (d_hours == i_hour)))
fidx = np.min(np.where((d_days == f_date) & (d_hours == f_hour)))
d_fracdays = data[iidx:fidx+1,5]
d_Btot = data[iidx:fidx+1,6]
d_Bx = data[iidx:fidx+1,7]
d_By = data[iidx:fidx+1,8]
d_Bz = data[iidx:fidx+1,9]
d_t = (d_fracdays - d_fracdays[0]) * 24
d_tUN = d_fracdays

# add in Earth motion in Car Lon if not already included (is for 45 CME cases)
#FFlon0  += (d_tUN[0] - CMEtime) * 360. / 365. #account for change in Earths position due to orbit
#eR2.insert(0, FFlon0)

# calculate the average magnetic field in the middle, used to normalize out B0
global avg_obs_B
avg_obs_B = np.mean(d_Btot[np.where(np.abs(d_tUN - CMEmid) < 2./24.)])

# read in Wind data if we have it
if flag_it==False:
	dataW = np.genfromtxt('Wind_CME'+ str(myid+1) +'.txt', skip_header=102, dtype=None)
	Wind_t  = []
	Wind_Bx = []
	Wind_By = []
	Wind_Bz = []
	for j in range(len(dataW)):
		Wrow = dataW[j]
		date = Wrow[0]
		time = Wrow[1]
		mydate = (DI.fracitup(int(date[6:]), int(date[3:5]), int(date[0:2]), int(time[0:2]), int(time[3:5]) + float(time[6:])/60.) - int(date[6:])) * yrlen
		if (mydate >= plotstart) & (mydate <=plotend):
			Wind_t.append(mydate)
			Wind_Bx.append(Wrow[2])
			Wind_By.append(Wrow[3])
			Wind_Bz.append(Wrow[4])	
	Wind_t  = np.array(Wind_t)
	Wind_Bx = np.array(Wind_Bx)
	Wind_By = np.array(Wind_By)
	Wind_Bz = np.array(Wind_Bz)
	# remove any bad data points
	Wind_t = Wind_t[np.where(Wind_Bx > -1e10)]
	Wind_Bz = Wind_Bz[np.where(Wind_Bx > -1e10)]
	Wind_By = Wind_By[np.where(Wind_Bx > -1e10)]
	Wind_Bx = Wind_Bx[np.where(Wind_Bx > -1e10)]
	Wind_B = np.sqrt(Wind_Bx**2 + Wind_By**2 + Wind_Bz**2)	

#autonormVAR.set(0)


f1 = open('fullSetdata.pkl', 'rb')
LASCO_date, Frac_date, Earth_lat, Earth_Clon, CR_number, AR_number, BS_POL, CS_POL, Used_POL, PIL_Lat, PIL_Lon, PIL_tiltX, COR1_time, COR1_lat, COR1_Clon, COR1_Slon, COR1_tiltX, COR1_dist, COR1_kap, COR1_AW, COR2_time, COR2_lat, COR2_Clon, COR2_Slon, COR2_tiltX, COR2_dist, COR2_kap, COR2_AW, FCAT_lat, FCAT_lon, FCAT_tilt, FCAT_A, FCAT_B, FCAT_initR, FCAT_r1, FCAT_r2, FCAT_vmin, FCAT_vmax, FCAT_AW0, FCAT_AWm, FCAT_AWmax, FCAT_AWscl, FCAT_mass, RC_start, RC_end, Wind_start, Wind_end, CME_start, CME_end, InSitu_v,newElon, CMEquality, FIDO_lat, FIDO_lon, FIDO_tilt, FIDO_AW, FIDO_A, FIDO_B, FIDO_xscr, FIDO_yscr, FIDO_zscr, FIDO_scor, FIDO_B0, FIDOB0_xscr, FIDOB0_yscr, FIDOB0_zscr, FIDOB0_scor, FIDO0_xscr, FIDO0_yscr, FIDO0_zscr, FIDO0_scor, FGCS_A, FGCS_B, FGCS_xscr, FGCS_yscr, FGCS_zscr, FGCS_scor, BFFT_lat, BFFT_lon, BFFT_tilt, BFFT_AW, BFFT_A, BFFT_B, BFFT_xscr, BFFT_yscr, BFFT_zscr, BFFT_scor, BFCT0_scor, BFCT_lat, BFCT_lon, BFCT_tilt, BFCT_AW, BFCT_A, BFCT_B, BFCT_C1, BFCT_tau, BFCT_xscr, BFCT_yscr, BFCT_zscr, BFCT_scor = pickle.load(f1)
f1.close

#e6.delete(0, END)
#e6.insert(0, -21) # latitude
#ec1.delete(0, END)
#ec1.insert(0, -21) # latitude

global autofitT
autofitT = True


print "normal"
# FF tor
global CMEcolor, CMEnumber
CMEcolor = '#332288'
CMEnumber = 10

update_plot()
# open up the ensemble pickle
FCAT_ensemble = pickle.load(open('ensemble_CME'+mynumSTR+'SCS.pkl','rb'))

print FIDO_lat[myid], FIDO_lon[myid], FIDO_tilt[myid]

CMEcolor = 'r'
CMEnumber = 2

number_misses = 0
one_hour = 1./24.
num_hours = int(obsCMEtime / one_hour)
ts = [one_hour * (i + 0.5) for i in range(num_hours+1) ]
print one_hour / 2.
print ts, obsCMEtime
global run_results
FIDO_ensemble = []
scoretot = []
scorex = []
scorey = []
scorez = []
for i in range(100):
    my_FCAT = FCAT_ensemble[i]
    #print my_FCAT[-1,:]
    # insert values 
    e1.delete(0, END)
    e1.insert(0, my_FCAT[-1,0]) # latitude
    e2.delete(0, END)
    e2.insert(0, my_FCAT[-1,1]) # longitude
    e3.delete(0, END)
    e3.insert(0, my_FCAT[-1,2]) # tilt
    print i+1
    update_plot()
    
    # check if miss
    if np.isnan(np.mean(np.abs(run_results[1]))): number_misses += 1
    # calculate hourly average values otherwise
    else:
        FIDO_ensemble.append(run_results)
        scorex.append(scoreBx)
        scorey.append(scoreBy)
        scorez.append(scoreBz)
        scoretot.append(totalscore)
        

print 'num misses', number_misses
print 'Bx ', np.mean(scorex), np.std(scorex) 
print 'By ', np.mean(scorey), np.std(scorey) 
print 'Bz ', np.mean(scorez), np.std(scorez) 
print 'B ', np.mean(scoretot), np.std(scoretot) 


pickle.dump(FIDO_ensemble, open('ensembleFIDO'+str(myid+1)+'SCS.pkl', 'wb'))
plt.savefig('ensembleFIDO'+str(myid+1)+'SCS.png')
#plt.savefig('Dec2008_plot.png')
