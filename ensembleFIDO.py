import FIDO
import numpy as np

# inps order is FIDO lat [0], FIDO lon0 [1], CMElat [2], CMElon [3], 
# CMEtilt [4], CMEAW [5], CMESRA [6], CMESRB [7], CMEvr [8], CMEB0 [9], 
# CMEH [10], tshift [11],  CME start [12],  CME end [13]

# set up params to vary and their uncertainties
# typically vary [2,3,4,5,6,7,8,9] if going for all CME observables
uncIDs = [2,3,4,5,6,7,8,9]
# have uncs for every options, use above ids to pick relevant ones
# default uncs = [0, 0, 5, 10, 10, 10, 0.2, 0.1, 50, 5, 0, 0, 0, 0]
# but can vary as desired
uncs = [0, 0, 5, 10, 10, 10, 0.2, 0.1, 50, 5, 0, 0, 0, 0]

def getnewinput(origvec, uncs, uncIDs):
    newinput = np.zeros(len(origvec))
    # make a copy so doesn't change the orig
    newinput[:] = origvec[:]
    for idx in uncIDs:
        newinput[idx] = np.random.normal(origvec[idx],0.5*uncs[idx],1)[0] 
    return newinput

# set up FIDO to run
# read in the text file and grab labeled inputs
input_values = FIDO.readinputfile()

# make sure these are set right, overwrite input file values
input_values['Launch_Gui'] = 'False'
input_values['No_Plot'] = 'True'
input_values['Silent'] = 'True'
# may want to make sure not autonorming
input_values['Autonormalize'] = 'False'

# set up the general properties (how will FIDO be ran)
FIDO.setupOptions(input_values)    
# set up the actual simulation input params
inps = FIDO.getInps(input_values)

# need these if potentially including observed data
# as can still use observed data to scale magnitude 
# though this has not been thoroughly tested
# otherwise no reason to provide in situ data for comp
# as not currently outputing any fit scores
if FIDO.ISfilename!=False: 
    FIDO.checkStartStop(inps)
    FIDO.setupObsData(inps)  



# Ensemble loop
misses = 0
maxiters = 10
Bxs = [[] for i in range(maxiters)]
Bys = [[] for i in range(maxiters)]
Bzs = [[] for i in range(maxiters)]
Kps = [[] for i in range(maxiters)]
i = 0
while i < maxiters:
    newinps = getnewinput(inps,uncs,uncIDs)
    Bout, tARR = FIDO.run_case(newinps)
    if len(Bout[0])>0 :
        Kp, BoutGSM = FIDO.calc_indices(Bout, inps[12], inps[8])
        Kp = FIDO.hourify(tARR-inps[12], Kp)
        Bx = FIDO.hourify(tARR-inps[12], Bout[0])
        By = FIDO.hourify(tARR-inps[12], Bout[1])
        Bz = FIDO.hourify(tARR-inps[12], Bout[1])
        Bxs[i] = Bx
        Bys[i] = By
        Bzs[i] = Bz
        Kps[i] = Kp
        print np.mean(Kp), np.max(Kp)
        i += 1
    else:
        print 'miss'
        misses += 1
print 'Ran '+str(maxiters)+' impacts, had '+str(misses)+ ' misses'



# function to take a data vector and print to file in format 
# that can be read by np.genfromtxt
def saveoutput(vecin,lens,maxlen,fname):
    f1 = open(fname,'w')
    for i in range(maxiters):
        outvec = np.zeros(maxlen)
        outvec[:lens[i]] = vecin[i]
        outprint = ''
        for item in outvec: outprint = outprint +'{:4.2f}'.format(item)+' '
        f1.write(outprint+'\n')
    f1.close()
    


# Saving output
lens = [len(Bxs[i]) for i in range(maxiters)]
maxlen = np.max(lens)
saveoutput(Kps, lens, maxlen, 'ensembleresultsKp.dat')

    