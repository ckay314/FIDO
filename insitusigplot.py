from pylab import *
import matplotlib.pyplot  as pyplot
import pickle
import numpy as np
from matplotlib import gridspec

matplotlib.rcParams.update({'font.size':12})

idstr='43'
data = pickle.load(open('ensembleFIDO'+idstr+'.pkl', 'rb'))
data2 = pickle.load(open('ensembleFIDO'+idstr+'synch2.pkl', 'rb'))

nCMEs = len(data)
nCMEs2 = len(data2)
print nCMEs, nCMEs2
nT = len(data[0][0])
ts = data[0][0] # should have ts the same for all cases

means = np.zeros([nT, 4])
stds  = np.zeros([nT, 4])
mins  = np.zeros([nT, 4])
maxs  = np.zeros([nT, 4])
meds  = np.zeros([nT, 4])
means2 = np.zeros([nT, 4])
stds2  = np.zeros([nT, 4])
mins2  = np.zeros([nT, 4])
maxs2  = np.zeros([nT, 4])
meds2  = np.zeros([nT, 4])

totB = np.zeros([nCMEs, nT])
totB2 = np.zeros([nCMEs2, nT])

for i in range(nT):
    for j in range(nCMEs):
        totB[j,i] = np.sqrt(data[j][1][i]**2 + data[j][2][i]**2 + data[j][3][i]**2)
    for j in range(nCMEs2):
        totB2[j,i] = np.sqrt(data2[j][1][i]**2 + data2[j][2][i]**2 + data2[j][3][i]**2)

for i in range(nT):
    for j in range(3):
        # need to loop through first idx (CMEnum) since doesn't like slicing array that way
        indiv_vals = []
        for k in range(nCMEs):
            indiv_vals.append(data[k][j+1][i])
        indiv_vals2 = []
        for k in range(nCMEs2):
            indiv_vals2.append(data2[k][j+1][i])
        means[i,j] = np.mean(indiv_vals)
        stds[i,j]  = np.std(indiv_vals)
        mins[i,j]  = np.min(indiv_vals)
        maxs[i,j]  = np.max(indiv_vals)
        meds[i,j]  = np.median(indiv_vals)
        means2[i,j] = np.mean(indiv_vals2)
        stds2[i,j]  = np.std(indiv_vals2)
        mins2[i,j]  = np.min(indiv_vals2)
        maxs2[i,j]  = np.max(indiv_vals2)
        meds2[i,j]  = np.median(indiv_vals2)
    means[i,3] = np.mean(totB[:,i])
    stds[i,3]  = np.std(totB[:,i])
    mins[i,3]  = np.min(totB[:,i])
    maxs[i,3]  = np.max(totB[:,i])
    meds[i,3]  = np.median(totB[:,i])
    means2[i,3] = np.mean(totB2[:,i])
    stds2[i,3]  = np.std(totB2[:,i])
    mins2[i,3]  = np.min(totB2[:,i])
    maxs2[i,3]  = np.max(totB2[:,i])
    meds2[i,3]  = np.median(totB2[:,i])
        
    
    

# read in in situ data
data3  = np.genfromtxt('CME'+idstr+'.txt', dtype=None)
RCstart = float(data3[43][1])
RCend   = float(data3[44][1])
pad = 3
plotstart = RCstart - pad/24.
plotend   = RCend + pad/24.

filename = 'ACE_CME'+ idstr +'.dat'
i_date = int(plotstart)
i_hour = int(plotstart % 1 * 24)
f_date = int(plotend) 
f_hour = int(plotend % 1 % 1 * 24)
print i_date, i_hour, f_date, f_hour
data4 = np.genfromtxt(filename, skip_header = 44, dtype=np.float)
d_yr = data4[:,0]
yrlen = 365
if (d_yr[0] % 4 == 0): yrlen = 366
d_days = data4[:,1]
d_hours  = data4[:,2]
# determine the initial and final index of the desired time period
iidx = np.min(np.where((d_days == i_date) & (d_hours == i_hour)))
fidx = np.min(np.where((d_days == f_date) & (d_hours == f_hour)))
d_fracdays = data4[iidx:fidx+1,5]
d_Btot = data4[iidx:fidx+1,6]
d_Bx = data4[iidx:fidx+1,7]
d_By = data4[iidx:fidx+1,8]
d_Bz = data4[iidx:fidx+1,9]
d_t = (d_fracdays - d_fracdays[0]) * 24 - 3

# get ACE hourly avgs for determining percent time within +- 1 sig
ACE_hrB = np.zeros([nT,4])
for i in range(nT):
        if len(d_Bx[abs(d_t - (i+.5))  <= 0.5]) > 0:
            ACE_hrB[i,0]  =  np.mean(d_Bx[abs(d_t - (i+.5))  <= 0.5])
            ACE_hrB[i,1]  =  np.mean(d_By[abs(d_t - (i+.5))  <= 0.5])
            ACE_hrB[i,2]  =  np.mean(d_Bz[abs(d_t - (i+.5))  <= 0.5])
            ACE_hrB[i,3] = np.sqrt(ACE_hrB[i,0]**2 + ACE_hrB[i,1]**2 + ACE_hrB[i,2]**2)

percX = float(len(np.where(np.abs(ACE_hrB[:,0] - means[:,0]) < stds[:,0])[0])) / float(len(stds[:,0]))
percY = float(len(np.where(np.abs(ACE_hrB[:,1] - means[:,1]) < stds[:,1])[0])) / float(len(stds[:,1]))
percZ = float(len(np.where(np.abs(ACE_hrB[:,2] - means[:,2]) < stds[:,2])[0])) / float(len(stds[:,2]))
print percX, percY, percZ
percX2 = float(len(np.where(np.abs(ACE_hrB[:,0] - means2[:,0]) < stds[:,0])[0])) / float(len(stds2[:,0]))
percY2 = float(len(np.where(np.abs(ACE_hrB[:,1] - means2[:,1]) < stds[:,1])[0])) / float(len(stds2[:,1]))
percZ2 = float(len(np.where(np.abs(ACE_hrB[:,2] - means2[:,2]) < stds[:,2])[0])) / float(len(stds2[:,2]))
print percX2, percY2, percZ2

    
    
    
fig = figure()
gs = gridspec.GridSpec(4,2)
ax0a = fig.add_subplot(gs[0])
ax1a = fig.add_subplot(gs[2], sharex=ax0a)
ax2a = fig.add_subplot(gs[4], sharex=ax0a)
ax3a = fig.add_subplot(gs[6], sharex=ax0a)
ax0b = fig.add_subplot(gs[1], sharey=ax0a)
ax1b = fig.add_subplot(gs[3], sharex=ax0b, sharey=ax1a)
ax2b = fig.add_subplot(gs[5], sharex=ax0b, sharey=ax2a)
ax3b = fig.add_subplot(gs[7], sharex=ax0b, sharey=ax3a)

setp(ax0b.get_yticklabels(), visible=False)
setp(ax1b.get_yticklabels(), visible=False)
setp(ax2b.get_yticklabels(), visible=False)
setp(ax3b.get_yticklabels(), visible=False)
setp(ax0a.get_xticklabels(), visible=False)
setp(ax0b.get_xticklabels(), visible=False)
setp(ax1a.get_xticklabels(), visible=False)
setp(ax1b.get_xticklabels(), visible=False)
setp(ax2a.get_xticklabels(), visible=False)
setp(ax2b.get_xticklabels(), visible=False)

degree = unichr(176)
ax0a.set_title('Synoptic')
ax0b.set_title('Synchronic')
ax0a.set_ylabel('B (G)')
ax1a.set_ylabel('B$_x$ (G)')
ax2a.set_ylabel('B$_y$ (G)')
ax3a.set_ylabel('B$_z$ (G)')
ax3a.set_xlabel('Time (hours)')
ax3b.set_xlabel('Time (hours)')

ax0a.fill_between(ts, mins[:,3], maxs[:,3], color='Silver')
ax1a.fill_between(ts, mins[:,0], maxs[:,0], color='Silver')
ax2a.fill_between(ts, mins[:,1], maxs[:,1], color='Silver')
ax3a.fill_between(ts, mins[:,2], maxs[:,2], color='Silver')
ax0b.fill_between(ts, mins2[:,3], maxs2[:,3], color='Silver')
ax1b.fill_between(ts, mins2[:,0], maxs2[:,0], color='Silver')
ax2b.fill_between(ts, mins2[:,1], maxs2[:,1], color='Silver')
ax3b.fill_between(ts, mins2[:,2], maxs2[:,2], color='Silver')

ax0a.fill_between(ts, means[:,3]-stds[:,3], means[:,3]+stds[:,3], color='Gray')
ax1a.fill_between(ts, means[:,0]-stds[:,0], means[:,0]+stds[:,0], color='Gray')
ax2a.fill_between(ts, means[:,1]-stds[:,1], means[:,1]+stds[:,1], color='Gray')
ax3a.fill_between(ts, means[:,2]-stds[:,2], means[:,2]+stds[:,2], color='Gray')
ax0b.fill_between(ts, means2[:,3]-stds2[:,3], means2[:,3]+stds2[:,3], color='Gray')
ax1b.fill_between(ts, means2[:,0]-stds2[:,0], means2[:,0]+stds2[:,0], color='Gray')
ax2b.fill_between(ts, means2[:,1]-stds2[:,1], means2[:,1]+stds2[:,1], color='Gray')
ax3b.fill_between(ts, means2[:,2]-stds2[:,2], means2[:,2]+stds2[:,2], color='Gray')

ax0a.plot(d_t, d_Btot, 'g', linewidth=2.5)
ax1a.plot(d_t, d_Bx, 'g', linewidth=2.5)
ax2a.plot(d_t, d_By, 'g', linewidth=2.5)
ax3a.plot(d_t, d_Bz, 'g', linewidth=2.5)
ax0b.plot(d_t, d_Btot, 'g', linewidth=2.5)
ax1b.plot(d_t, d_Bx, 'g', linewidth=2.5)
ax2b.plot(d_t, d_By, 'g', linewidth=2.5)
ax3b.plot(d_t, d_Bz, 'g', linewidth=2.5)

ax0a.plot(ts, means[:,3], 'k--', linewidth=2.5)
ax1a.plot(ts, means[:,0], 'k--', linewidth=2.5)
ax2a.plot(ts, means[:,1], 'k--', linewidth=2.5)
ax3a.plot(ts, means[:,2], 'k--', linewidth=2.5)
ax0b.plot(ts, means2[:,3], 'k--', linewidth=2.5)
ax1b.plot(ts, means2[:,0], 'k--', linewidth=2.5)
ax2b.plot(ts, means2[:,1], 'k--', linewidth=2.5)
ax3b.plot(ts, means2[:,2], 'k--', linewidth=2.5)

ax0a.plot(ts, meds[:,3], 'k', linewidth=2.5)
ax1a.plot(ts, meds[:,0], 'k', linewidth=2.5)
ax2a.plot(ts, meds[:,1], 'k', linewidth=2.5)
ax3a.plot(ts, meds[:,2], 'k', linewidth=2.5)
ax0b.plot(ts, meds2[:,3], 'k', linewidth=2.5)
ax1b.plot(ts, meds2[:,0], 'k', linewidth=2.5)
ax2b.plot(ts, meds2[:,1], 'k', linewidth=2.5)
ax3b.plot(ts, meds2[:,2], 'k', linewidth=2.5)

ax0a.plot(ts, totB[0,:], 'b', linewidth=2.5)
ax1a.plot(ts, data[0][1], 'b', linewidth=2.5)
ax2a.plot(ts, data[0][2], 'b', linewidth=2.5)
ax3a.plot(ts, data[0][3], 'b', linewidth=2.5)
ax0b.plot(ts, totB[0,:], 'b', linewidth=2.5)
ax1b.plot(ts, data[0][1], 'b', linewidth=2.5)
ax2b.plot(ts, data[0][2], 'b', linewidth=2.5)
ax3b.plot(ts, data[0][3], 'b', linewidth=2.5)
ax0b.plot(ts, totB2[0,:], 'r', linewidth=2.5)
ax1b.plot(ts, data2[0][1], 'r', linewidth=2.5)
ax2b.plot(ts, data2[0][2], 'r', linewidth=2.5)
ax3b.plot(ts, data2[0][3], 'r', linewidth=2.5)

#ax1a.plot(ts, ACE_hrB[:,0], 'go', )
#ax2a.plot(ts, ACE_hrB[:,1], 'go', )
#ax3a.plot(ts, ACE_hrB[:,2], 'go', )

ax0a.set_xlim([-2, np.max(ts)+2])
ax0b.set_xlim([-2, np.max(ts)+2])

setp(xticks)

subplots_adjust(hspace=0.001, wspace=0.001)
savefig('sig_ensemble_CME'+idstr+'.png')
plt.show()
