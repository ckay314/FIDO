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
    plotCME = True
    try:
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
        if fluxRopeToggle==1:
            CMEB0     = float(ec1.get()) 
            cM      = float(ec2.get())
            cN      = float(ec3.get())
            cTau    = float(ec4.get())
            cC1     = float(ec5.get())
    except:
        plotCME = False
        print 'Error in CME params'
        
    if plotCME==True:
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

        # scale the CME to match at midpoint (ignoring B0 with this)
        avg_FIDO_B = np.mean(obsB[np.where(np.abs(tshift+tARR - CMEmid) < 2./24.)])
        scale = 1.
        if autonormVAR.get()==1: scale = avg_obs_B / avg_FIDO_B
        obsBx *= scale
        obsBy *= scale
        obsBz *= scale
        obsB = np.sqrt(obsBx**2 + obsBy**2 + obsBz**2)

        scoreBx, scoreBy, scoreBz = calc_score()
        print 'initial score', totalscore
      
    ax2.clear()
    ax3.clear()
    ax4.clear()
    ax5.clear()
    #ax2.plot(d_t, d_Btot, 'k', linewidth=3)
    # reassign so won't complain when looking for best plot range
    if plotCME == False:
        obsB, obsBx, obsBy, obsBz = d_Btot, d_Bx, d_By, d_Bz
    maxBtot = np.max([np.max(d_Btot), np.max(obsB)])
    ax2.plot([CMEstart, CMEstart], [0., 1.3*maxBtot], 'k--', linewidth=2)
    ax2.plot([CMEend, CMEend], [0., 1.3*maxBtot], 'k--', linewidth=2)
    #if flag_it == False:
    #    ax2.plot([TNCstart, TNCstart], [0., 1.3*np.max(maxBtot)], 'b--', linewidth=2)
    #    ax2.plot([TNCend, TNCend], [0., 1.3*np.max(maxBtot)], 'b--', linewidth=2)
    #    ax2.plot(Wind_t, Wind_B, 'b', linewidth=4)
    ax2.plot(d_tUN, d_Btot, 'k', linewidth=4)
    #print d_tUN
    if plotCME == True: ax2.plot(tshift + tARR, obsB, 'r', linewidth=4)
    ax2.set_ylabel('B (nT)')
    setp(ax2.get_xticklabels(), visible=False)

    minBx = np.abs(np.min([np.min(d_Bx), np.min(-obsBx)]))
    maxBx = np.max([np.max(d_Bx), np.max(-obsBx)])
    ax3.plot([CMEstart, CMEstart], [-1.3*minBx, 1.3*maxBx], 'k--', linewidth=2)
    ax3.plot([CMEend, CMEend], [-1.3*minBx, 1.3*maxBx], 'k--', linewidth=2)
    #if flag_it == False:
    #    ax3.plot([TNCstart, TNCstart], [-1.3*minBx, 1.3*maxBx], 'b--', linewidth=2)
    #    ax3.plot([TNCend, TNCend], [-1.3*minBx, 1.3*maxBx], 'b--', linewidth=2)
    #    ax3.plot(Wind_t, Wind_Bx, 'b', linewidth=4)
    ax3.plot(d_tUN, d_Bx, 'k', linewidth=4)
    if plotCME == True: 
        ax3.plot(tshift + tARR, -obsBx, 'r', linewidth=4)
        ax3.annotate('%0.2f'%(scoreBx), xy=(1, 0), color='r', xycoords='axes fraction', fontsize=16,
                horizontalalignment='right', verticalalignment='bottom')    
    ax3.set_ylabel('B$_x$ (nT)')
    setp(ax3.get_xticklabels(), visible=False)

    minBy = np.abs(np.min([np.min(d_By), np.min(-obsBy)]))
    maxBy = np.max([np.max(d_By), np.max(-obsBy)])
    ax4.plot([CMEstart, CMEstart], [-1.3*minBy, 1.3*maxBy], 'k--', linewidth=2)
    ax4.plot([CMEend, CMEend], [-1.3*minBy, 1.3*maxBy], 'k--', linewidth=2)
    #if flag_it == False:
    #    ax4.plot([TNCstart, TNCstart], [-1.3*minBy, 1.3*maxBy], 'b--', linewidth=2)
     #   ax4.plot([TNCend, TNCend], [-1.3*minBy, 1.3*maxBy], 'b--', linewidth=2)
    #    ax4.plot(Wind_t, Wind_By, 'b', linewidth=4)
    ax4.plot(d_tUN, d_By, 'k', linewidth=4)
    if plotCME == True: 
            ax4.plot(tshift + tARR, -obsBy, 'r', linewidth=4)
            ax4.annotate('%0.2f'%(scoreBy), xy=(1, 0), color='r', xycoords='axes fraction', fontsize=16,
                horizontalalignment='right', verticalalignment='bottom')    
    ax4.set_ylabel('B$_y$ (nT)')
    setp(ax4.get_xticklabels(), visible=False)

    minBz = np.abs(np.min([np.min(d_Bz), np.min(obsBz)]))
    maxBz = np.max([np.max(d_Bz), np.max(obsBz)])
    ax5.plot([CMEstart, CMEstart], [-1.3*minBz, 1.3*maxBz], 'k--', linewidth=2)
    ax5.plot([CMEend, CMEend], [-1.3*minBz, 1.3*maxBz], 'k--', linewidth=2)
    #if flag_it == False:
    #    ax5.plot([TNCstart, TNCstart], [-1.3*minBz, 1.3*maxBz], 'b--', linewidth=2)
    #    ax5.plot([TNCend, TNCend], [-1.3*minBz, 1.3*maxBz], 'b--', linewidth=2)
    #    ax5.plot(Wind_t, Wind_Bz, 'b', linewidth=4)
    ax5.plot(d_tUN, d_Bz, 'k', linewidth=4)
    ax5.set_ylabel('B$_z$ (nT)')

    ax5.set_xlabel('Day of Year')
    if plotCME == True: 
        ax5.plot(tshift + tARR, obsBz, 'r', linewidth=4)
        ax5.annotate('%0.2f'%(scoreBz), xy=(1, 0), color='r', xycoords='axes fraction', fontsize=16,
                horizontalalignment='right', verticalalignment='bottom')    
    
    
    scl = 1.25
    ax2.set_ylim([0,scl*maxBtot])
    ax3.set_ylim([-1.*scl*minBx, scl*maxBx])
    ax4.set_ylim([-1.*scl*minBy, scl*maxBy])
    ax5.set_ylim([-1.*scl*minBz, scl*maxBz])
    ax2.set_xlim([plotstart, plotend])
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.subplots_adjust(right=0.8, wspace=0.001, hspace=0.25)
    #fig2.suptitle(datestr, fontsize=16)
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
 
    latrange = [prevCMElat-5, prevCMElat+5.]
    lonrange = [prevCMElon-5., prevCMElon+5.]
    tiltrange = [prevCMEtilt-10., prevCMEtilt+10.]
    AWrange  = [prevCMEAW-10., prevCMEAW+10.]
    Brange = [0.01, 0.99]
    Arange = [0.4, 1.0]
    #kapparange = [np.min([myrow2[9]-0.25, 0.01]), myrow2[9]+0.25]
    # opposite sign of force free H
    C1range = [0.5, 2.0]
    
    for bf_counter in range(4000):
        global totalscore
        oldscore = totalscore
    
        global CMElat, CMElon, CMEtilt, CMEAW, CMEB0, CMESRA, CMESRB, kappa, cC1      
        newCMElat = prevCMElat + 0.05 * (random.random() - 0.5)
        newCMElat = np.max([latrange[0], newCMElat])
        newCMElat = np.min([latrange[1], newCMElat])
        CMElat = newCMElat
        
        newCMElon = prevCMElon + 0.05 * (random.random() - 0.5)
        newCMElon = np.max([lonrange[0], newCMElon])
        newCMElon = np.min([lonrange[1], newCMElon])
        CMElon = newCMElon
        
        newCMEtilt = prevCMEtilt + 0.05 * (random.random() - 0.5)
        newCMEtilt = np.max([tiltrange[0], newCMEtilt])
        newCMEtilt = np.min([tiltrange[1], newCMEtilt])
        CMEtilt = newCMEtilt

        newCMEAW = prevCMEAW + 0.05 * (random.random() - 0.5)
        newCMEAW = np.max([AWrange[0], newCMEAW])
        newCMEAW = np.min([AWrange[1], newCMEAW])
        CMEAW = newCMEAW
 
        #newCMEB0 = prevCMEB0 + 0.25 * (random.random() - 0.5)
        #newCMEB0 = np.max([B0range[0], newCMEB0])
        #newCMEB0 = np.min([B0range[1], newCMEB0])
        #CMEB0 = newCMEB0

        newCMESRA = prevCMESRA + 0.002 * (random.random() - 0.5)
        newCMESRA = np.max([Arange[0], newCMESRA])
        newCMESRA = np.min([Arange[1], newCMESRA])
        CMESRA = newCMESRA

        newCMESRB = prevCMESRB + 0.0005 * (random.random() - 0.5)
        newCMESRB = np.max([Brange[0], newCMESRB])
        newCMESRB = np.min([Brange[1], newCMESRB])
        CMESRB = newCMESRB
        
        #newkappica = prevkappa + 0.005 * (random.random() - 0.5)
        #newkappica = np.max([kapparange[0], newkappica])
        #newkappica = np.min([kapparange[1], newkappica])
        #kappa = newkappica
        #global alphaGCS
        #alphaGCS = CMEAW *dtor - np.arcsin(kappa)
        
        newC1 = prevC1 + 0.005 * (random.random() - 0.5)
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
        

        if totalscore > 1.0001 * oldscore:
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

def get_inputs(inputs):
    # take a file with unsorted input values and return a dictionary.
    # variable names have to match their names below
    possible_vars = ['insitufile', 'Earth_lat', 'Earth_lon', 'CME_lat', 'CME_lon', 'CME_tilt', 'CME_AW', 'CME_vr', 'tshift', 'CME_Ashape', 'CME_Bshape', 'CME_B0', 'CME_pol', 'CME_start', 'CME_stop', 'CME_shape_model', 'Flux_rope_model', 'Autonormalize', 'Circ_m', 'Circ_n', 'Circ_C1', 'Circ_Tau']
    
    # if matches add to dictionary
    input_values = {}
    for i in range(len(inputs)):
        temp = inputs[i]
        if temp[0][:-1] in possible_vars:
            input_values[temp[0][:-1]] = temp[1]
    return input_values


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
if len(sys.argv) < 2: sys.exit("Need an input file")

input_file = sys.argv[1]
inputs = np.genfromtxt(input_file, dtype=None)

input_values = get_inputs(inputs)

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
if 'CME_shape_model' in input_values:
    if input_values['CME_shape_model'] == 'Torus': CMEshapeToggleVAR.set(0)
    elif input_values['CME_shape_model'] == 'Torus': CMEshapeToggleVAR.set(1)

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
if 'Flux_rope_model' in input_values:
    if input_values['Flux_rope_model'] == 'FF': fluxRopeToggleVAR.set(0)
    elif input_values['Flux_rope_model'] == 'GCS': fluxRopeToggleVAR.set(1)
Label(root, text='Flux Rope Model:', bg='gray75').grid(row=0, column=3,columnspan=2)
Radiobutton(root, text='Force Free', variable=fluxRopeToggleVAR, value=0, bg='gray75').grid(column=3,row=1)
Radiobutton(root, text='Circular', variable=fluxRopeToggleVAR, value=1, bg='gray75').grid(column=4,row=1)

# check button for autonormalizing magnitude
global autonormVAR
autonormVAR = IntVar()
if 'Autonormalize' in input_values:
    if input_values['Autonormalize'] == 'True': autonormVAR.set(1)
    elif input_values['Autonormalize'] == 'False': autonormVAR.set(0)

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
#data = np.genfromtxt(myfile, dtype=None)
#datestr = data[0][1]
#datestr = datestr[:10] + ' ' + datestr[11:]
#print datestr
#CMEtime = data[1][1] # get DOY from the year frac 
#print 'CME time', CMEtime

# insert values 
if 'Earth_lat' in input_values:
    eR1.insert(0, input_values['Earth_lat'])
if 'Earth_lon' in input_values:
    eR2.insert(0, input_values['Earth_lat'])
if 'CME_lat' in input_values:
    e1.insert(0, input_values['CME_lat']) 
if 'CME_lon' in input_values:
    e2.insert(0, input_values['CME_lon'])
if 'CME_tilt' in input_values:
    e3.insert(0, input_values['CME_tilt']) 
if 'CME_AW' in input_values:
    e3b.insert(0, input_values['CME_AW'])
if 'CME_Ashape' in input_values:
    e4.insert(0, input_values['CME_Ashape']) 
if 'CME_Bshape' in input_values:
    e5.insert(0, input_values['CME_Bshape'])   
if 'CME_vr' in input_values:
    e8.insert(0, input_values['CME_vr'])
if 'tshift' in input_values:
    e9.insert(0, input_values['tshift'])
if 'CME_B0' in input_values:
    # stick in both for now
    e6.insert(0, input_values['CME_B0'])
    ec1.insert(0,input_values['CME_B0']) 
if 'Circ_m' in input_values:
    ec2.insert(0,input_values['Circ_m']) 
if 'Circ_n' in input_values:
    ec3.insert(0,input_values['Circ_n'])        
if 'Circ_Tau' in input_values:
    ec4.insert(0,input_values['Circ_Tau'])     
if 'Circ_C1' in input_values:
    ec5.insert(0,input_values['Circ_C1'])  

global CMEH
CMEH = 1
if 'CME_pol' in input_values:    
    used_pol = input_values['CME_pol']
    if used_pol[1] == '-': CMEH = -1
    e7.insert(0, CMEH) # H
    


# get the in situ information
#global RCstart, RCend, TNCstart, TNCend, plotstart, plotend
global CMEstart, CMEend, CMEmid, plotstart, plotend

if 'CME_start' in input_values:  CMEstart = float(input_values['CME_start'])
else: sys.exit('Add CME_start to input file')
if 'CME_stop' in input_values:  CMEend = float(input_values['CME_stop'])
else: sys.exit('Add CME_stop to input file')

pad = 3
plotstart = CMEstart - pad/24.
plotend   = CMEend + pad/24.
CMEmid    = 0.5 * (CMEstart + CMEend)
print "start, mid, end:  ", CMEstart, CMEmid, CMEend



# read in ACE data
global d_t, d_Btot, d_Bx, d_By, d_Bz, Wind_t, Wind_B, Wind_Bx, Wind_By, Wind_Bz

if 'insitufile' in input_values:
    filename = input_values['insitufile']
else: sys.exit('Add insitufile to input file')

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
#FFlon0  += (d_tUN[0] - CMEtime) * 360. / 365. #account for change in Earths position due to orbit
#eR2.insert(0, FFlon0)

# calculate the average magnetic field in the middle, used to normalize out B0
global avg_obs_B
avg_obs_B = np.mean(d_Btot[np.where(np.abs(d_tUN - CMEmid) < 2./24.)])

flag_it = False # turned off for now
# read in Wind data if we have it
'''if flag_it==False:
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
	Wind_B = np.sqrt(Wind_Bx**2 + Wind_By**2 + Wind_Bz**2)	'''

update_plot()
root.mainloop()
root.quit()

