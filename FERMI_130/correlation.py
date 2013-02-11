#===============================================================================
# Compute 2 point correlation information on a dataset
# Author: Eric Carlson 
# Date: Oct 9th, 2012
#===============================================================================
import matplotlib.pyplot as plt
import scipy.spatial.distance as dist
import math, pickle
import numpy as np




def RandomCorrelations(nPoints = 5000,width = 20.0):
    X1 = GenRandomSet(nPoints,width) # random points
    ND1,DDcount1,DDr1 = GenAutoHist(X1)  # Autocorrelate and bin
    pickle.dump((ND1,DDcount1,DDr1), open('./randAuto.pickle','wb'))  # Write to file
    
    # Same for cross correlations
    X2 = GenRandomSet(nPoints,width) # random points
    ND1a, ND1b,DDcount1,DDr1 = GenCrossHist(X1,X2)  # Autocorrelate and bin
    pickle.dump((ND1a, ND1b,DDcount1,DDr1), open('./randCross.pickle','wb'))  # Write to file
    return
                
                
def Correlation(X1,X2,X3, pulsar = ''):
    
    #=======================================================
    # Autocorrelations
    
    # Generate histograms
    ND1,DDcount1,DDr1 = GenAutoHist(X1)
    ND2,DDcount2,DDr2 = GenAutoHist(X2)
    ND3,DDcount3,DDr3 = GenAutoHist(X3)
    
    if pulsar != '':
        NDP,DDcountP,DDrP = GenAutoHist(pulsar)
        
    
    # Load random histogram 
    NR,RRcount,RRr = pickle.load(open('./randAuto.pickle','rb'))
    # Compute xi function 
    xi1 = xiAppend(NR,ND1,DDcount1,RRcount)
    xi2 = xiAppend(NR,ND2,DDcount2,RRcount)
    xi3 = xiAppend(NR,ND3,DDcount3,RRcount)
    
    if pulsar != '':
        xiP = xiAppend(NR,NDP,DDcountP,RRcount)
    
    #=======================================================
    # Plot
    fig = plt.figure(1, (8,10))
    fig.add_subplot(211)
    plt.subplot(211)
    plt.step(DDr1, xi1,label = '70-110 GeV',ls='-.')
    plt.step(DDr2, xi2,label = '120-140 GeV',ls='-.')
    plt.step(DDr3, xi3,label = '150-300 GeV',ls='-.')
    if pulsar != '':
        plt.step(DDrP, xiP,label = '1 Pulsar')
    
    
    CorrelateSimulations2()
    
    
    
    plt.xlabel(r'$r[^\circ]$')
    plt.ylabel(r'$\xi$')   
    #plt.ylim(-1,1)
    #plt.xlim(10**-1.2,10)
    plt.xscale('log')
    # Legend
    from matplotlib.font_manager import FontProperties
    fontP = FontProperties()
    fontP.set_size('small')
    plt.legend(loc=1, ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.,labelspacing = .2)
    plt.legend(('70-110 GeV','120-140 GeV','150-300 GeV','NFW','NFWDECAY','Galprop BG','1 Pulsar','2 Pulsar','3 Pulsar'),loc=1, ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.,labelspacing = .2)
    plt.ylim(-1,20)
    
    #=======================================================
    # X-Correlations
    NDa1,NDb1,DDcount1,DDr1 = GenCrossHist(X1,X2)
    ND2a1,ND2b1,DD2count1,DD2r1 = GenCrossHist(X2,X1)
    
    
    NDa2,NDb2,DDcount2,DDr2 = GenCrossHist(X1,X3)
    ND2a2,ND2b2,DD2count2,DD2r2 = GenCrossHist(X3,X1)
    
    NDa3,NDb3,DDcount3,DDr3 = GenCrossHist(X2,X3)
    ND2a3,ND2b3,DD2count3,DD2r3 = GenCrossHist(X3,X2)
    
    NRa, NRb,RRcount,RRr = pickle.load(open('./randCross.pickle','rb'))
    xi1 = xiAppendCross(NRa,NRb,NDa1,NDb1,DDcount1,RRcount)
    xi21 = xiAppendCross(NRa,NRb,ND2a1,ND2b1,DD2count1,RRcount)
    
    xi2 = xiAppendCross(NRa,NRb,NDa2,NDb2,DDcount2,RRcount)
    
    
    xi3 = xiAppendCross(NRa,NRb,NDa3,NDb3,DDcount3,RRcount)
    
    #=======================================================
    # Plot
    fig.add_subplot(212)
    plt.subplot(212)
    plt.step(DDr1, xi1,label = '110x130GeV',ls='--')
    #plt.step(DD2r1, xi21,label = '130x110GeV')
    
    plt.step(DDr2, xi2,label = '110x150 GeV',ls='--')
    
    plt.step(DDr3, xi3,label = '130x150 GeV',ls='--')
    plt.xlabel(r'$r[^\circ]$')
    plt.ylabel(r'$\xi$')   
    plt.ylim(-1,1)
    plt.xlim(10**-1.2,10)
    plt.xscale('log')
    # Legend
    plt.legend(loc=1, ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.,labelspacing = .2)
    
    plt.show()
    
    
def GetXi(X):
    
    #=======================================================
    # Autocorrelations
    
    # Generate histograms
    ND1,DDcount1,DDr1 = GenAutoHist(X)    
    # Load random histogram 
    NR,RRcount,RRr = pickle.load(open('./randAuto.pickle','rb'))
    # Compute xi function 
    xi = xiAppend(NR,ND1,DDcount1,RRcount)
    return xi
    


def CorrelateSimulations(numSims = 10,numPhotons=500):
    import MC
    
    pulsar1,pulsar2,pulsar3,NFW,NFWDECAY = [],[],[],[],[]
    
    # Pulsars
    profile = ('PULSAR',)
    fileOut = 'PULSARRateMap.pickle'
    mcSims = MC.RUN_PULSAR(numSims, fileOut, numPhotons=numPhotons,numPulsars=1,angularSize=10.0, outputSize=300, mcList='test' + '.pickle')
    for j in range(len(mcSims)):
        pulsar1.append([])
        for i in range(len(mcSims[j][0])):
            pulsar1[j].append((mcSims[j][0][i],mcSims[j][1][i]))
            
              
    profile = ('PULSAR',)
    fileOut = 'PULSARRateMap.pickle'
    mcSims = MC.RUN_PULSAR(numSims, fileOut, numPhotons=numPhotons,numPulsars=2,angularSize=10.0, outputSize=300, mcList='test' + '.pickle')
    for j in range(len(mcSims)):
        pulsar2.append([])
        for i in range(len(mcSims[j][0])):
            pulsar2[j].append((mcSims[j][0][i],mcSims[j][1][i]))
    
    
    profile = ('PULSAR',)
    fileOut = 'PULSARRateMap.pickle'
    mcSims = MC.RUN_PULSAR(numSims, fileOut, numPhotons=numPhotons,numPulsars=3,angularSize=10.0, outputSize=300, mcList='test' + '.pickle')
    for j in range(len(mcSims)):
        pulsar3.append([])
        for i in range(len(mcSims[j][0])):
            pulsar3[j].append((mcSims[j][0][i],mcSims[j][1][i]))
    
    
    profile = ('NFW',23.5,1.0)
    fileOut = 'NFWRateMap.pickle'
    mcSims = MC.RUN(numSims, fileOut, numPhotons=numPhotons,angularSize=10.0, outputSize=300, mcList='test' + '.pickle')
    for j in range(len(mcSims)):
        NFW.append([])
        for i in range(len(mcSims[j][0])):
            NFW[j].append((mcSims[j][0][i],mcSims[j][1][i]))
    
    
    profile = ('NFWDECAY',23.5,1.0)
    fileOut = 'NFWDECAYRateMap.pickle'
    mcSims = MC.RUN(numSims, fileOut, numPhotons=numPhotons,angularSize=10.0, outputSize=300, mcList='test' + '.pickle')
    for j in range(len(mcSims)):
        NFWDECAY.append([])
        for i in range(len(mcSims[j][0])):
            NFWDECAY[j].append((mcSims[j][0][i],mcSims[j][1][i]))

    

    for simSet in (pulsar1,pulsar2,pulsar3,NFW,NFWDECAY):
        xiTotal = []
        for i in range(10):
            xiTotal.append([])
        for sim in simSet:
            xi = GetXi(sim)
            for i in range(len(xi)):
                xiTotal[i].append(xi[i])
        xiFinal, xiFinalErr = [], []
        for i in range(len(xiTotal)):
            xiFinal.append(np.mean(xiTotal[i]))
            xiFinalErr.append(np.std(xiTotal[i]))
        plt.errorbar(bins, xiFinal, xiFinalErr,marker='o',ls = '')
    
    from matplotlib.font_manager import FontProperties
    fontP = FontProperties()
    fontP.set_size('small')
    plt.legend(('1 Pulsar','2 Pulsar', '3 Pulsar','NFW', 'NFW Decay'),loc=1, ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.,labelspacing = .2)
    plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel(r'$r[^\circ]$')
    plt.ylabel(r'$\xi$')
    #plt.savefig('correlation_puls_nfw_nfwdecay_500.pdf')
    plt.show()

#bins = np.logspace(-1.5,.69897,10)
bins = [0,.2,.4,1,2,3,4,5,6,7]
def CorrelateSimulations2(numSims = 10,numPhotons=50):
    
    import MC
    pulsar1,pulsar2,pulsar3,NFW,NFWDECAY,BGONLY = [],[],[],[],[],[]
    
    profile = ('NFW',23.5,1.0)
    fileOut = 'NFWRateMap.pickle'
    mcSims = MC.RUN(numSims, fileOut, numPhotons=numPhotons,angularSize=10.0, outputSize=300, mcList='test' + '.pickle')
    for j in range(len(mcSims)):
        NFW.append([])
        for i in range(len(mcSims[j][0])):
            NFW[j].append((mcSims[j][0][i],mcSims[j][1][i]))
    
    
    profile = ('NFWDECAY',23.5,1.0)
    fileOut = 'NFWDECAYRateMap.pickle'
    mcSims = MC.RUN(numSims, fileOut, numPhotons=numPhotons,angularSize=10.0, outputSize=300, mcList='test' + '.pickle')
    for j in range(len(mcSims)):
        NFWDECAY.append([])
        for i in range(len(mcSims[j][0])):
            NFWDECAY[j].append((mcSims[j][0][i],mcSims[j][1][i]))
            
            
    profile = ('NFWDECAY',23.5,1.0)
    fileOut = 'NFWDECAYRateMap.pickle'
    mcSims = MC.RUN_BG_ONLY(numSims, fileOut, numPhotons=numPhotons,angularSize=10.0, outputSize=300, mcList='test' + '.pickle')
    for j in range(len(mcSims)):
        BGONLY.append([])
        for i in range(len(mcSims[j][0])):
            BGONLY[j].append((mcSims[j][0][i],mcSims[j][1][i]))        
            
    
    # Pulsars
    profile = ('PULSAR',)
    fileOut = 'PULSARRateMap.pickle'
    mcSims = MC.RUN_PULSAR(numSims, fileOut, numPhotons=numPhotons,numPulsars=1,angularSize=10.0, outputSize=300, mcList='test' + '.pickle')
    for j in range(len(mcSims)):
        pulsar1.append([])
        for i in range(len(mcSims[j][0])):
            pulsar1[j].append((mcSims[j][0][i],mcSims[j][1][i]))
            
              
    profile = ('PULSAR',)
    fileOut = 'PULSARRateMap.pickle'
    mcSims = MC.RUN_PULSAR(numSims, fileOut, numPhotons=numPhotons,numPulsars=2,angularSize=10.0, outputSize=300, mcList='test' + '.pickle')
    for j in range(len(mcSims)):
        pulsar2.append([])
        for i in range(len(mcSims[j][0])):
            pulsar2[j].append((mcSims[j][0][i],mcSims[j][1][i]))
    
    
    profile = ('PULSAR',)
    fileOut = 'PULSARRateMap.pickle'
    mcSims = MC.RUN_PULSAR(numSims, fileOut, numPhotons=numPhotons,numPulsars=3,angularSize=10.0, outputSize=300, mcList='test' + '.pickle')
    for j in range(len(mcSims)):
        pulsar3.append([])
        for i in range(len(mcSims[j][0])):
            pulsar3[j].append((mcSims[j][0][i],mcSims[j][1][i]))

    

    for simSet in (NFW,NFWDECAY,BGONLY,pulsar1,pulsar2,pulsar3):
        xiTotal = []
        for i in range(10):
            xiTotal.append([])
        for sim in simSet:
            xi = GetXi(sim)
            for i in range(len(xi)):
                xiTotal[i].append(xi[i])
        xiFinal, xiFinalErr = [], []
        for i in range(len(xiTotal)):
            xiFinal.append(np.mean(xiTotal[i]))
            xiFinalErr.append(np.std(xiTotal[i]))
        plt.errorbar(bins, xiFinal, xiFinalErr,marker='o',ls = '')
    
#    from matplotlib.font_manager import FontProperties
#    fontP = FontProperties()
#    fontP.set_size('small')
#    plt.xscale('log')
#    plt.xlabel(r'$r[^\circ]$')
#    plt.ylabel(r'$\xi$')
    #plt.savefig('correlation_puls_nfw_nfwdecay_500.pdf')
    #plt.show()    
    
def GenAutoHist(X1):
    """Internal.  Generates a histogram of pairwise distances for autocorrelation"""
    
    #bins = np.logspace(-1.5,.69897,10)
    bins = [0,.2,.4,1,2,3,4,5,6,7]
    
    r1= []
    ND1 = 0 # num points in box
    for i in range(len(X1)):
        if -5<X1[i][0] <5 and -5<X1[i][1]<5: # if point is in box
            ND1+=1 # increment counter
            for j in range(len(X1)):         # loop over all points
                if -5<X1[j][0] <5 and -5<X1[j][1]<5 and j>i: #if 2nd point in box and this pair hasn't been computed
                    r1.append(math.sqrt( (X1[j][0]-X1[i][0])**2.0 + (X1[j][1]-X1[i][1])**2.0)) # append pairwise distance
                elif  (X1[i][0]>5 or X1[i][0]<-5) and (5<X1[i][1] or X1[i][1]<-5): # if outside box then count it.
                    r1.append(math.sqrt( (X1[j][0]-X1[i][0])**2.0 + (X1[j][1]-X1[i][1])**2.0))
    
    # now bin these
    DDhist1 = np.histogram(r1, bins = bins)
    bins = []
    for i in range(len(DDhist1[1])-1):
        bins.append((DDhist1[1][i]+DDhist1[1][i+1])/2.)
    DDcount1, DDr1 = DDhist1[0][:], bins # using outer boundaries
    
    return ND1,DDcount1,DDr1
    
def GenCrossHist(X1,X2):
    """Internal.  Generates a histogram of pairwise distances for cross-correlation"""
    r1= []
    ND1, ND2 = 0,0 # num points in box
    for i in range(len(X1)):
        if -5<X1[i][0] <5 and -5<X1[i][1]<5: # if point is in box
            ND1+=1 # increment counter
            for j in range(len(X2)):         # loop over all points in X2
                r1.append(math.sqrt( (X2[j][0]-X1[i][0])**2.0 + (X2[j][1]-X1[i][1])**2.0)) # append pairwise distance
                
                
    for i in range(len(X2)):
        if -5<X2[i][0] <5 and -5<X2[i][1]<5: # if point in x2 is in box
            ND2+=1 # increment counter
            for j in range(len(X1)):  # loop over points in x1, but avoid any points in box
                    if  (X1[j][0]>5 or X1[j][0]<-5) and (5<X1[j][1] or X1[j][1]<-5): # if outside box then count
                        r1.append(math.sqrt((X2[i][0]-X1[j][0])**2.0 + (X2[i][1]-X1[j][1])**2.0))
    # now bin these
    
    DDhist1 = np.histogram(r1, bins = bins)
    DDcount1, DDr1 = DDhist1[0][:], DDhist1[1][1:] # using outer boundaries
    
    return ND1,ND2,DDcount1,DDr1
    
    
def xiAppend(NR,ND,DDcount,RRcount):
    xi = []
    for i in range(len(DDcount)):
        if float(RRcount[i])!=0:
            xi.append(NR*(NR-1.0)/(ND*(ND-1))*float(DDcount[i])/float(RRcount[i])-1) 
        else:
            xi.append(-1.0)
    return xi

def xiAppendCross(NR1,NR2,ND1,ND2,DDcount,RRcount):
    xi = []
    for i in range(len(DDcount)):
        if float(RRcount[i])!=0:
            xi.append(NR1*NR2/(ND1*ND2)*float(DDcount[i])/float(RRcount[i])-1) 
        else:
            xi.append(-1.0)
    return xi

    
def Cross_Dist(X1,X2):
    '''
    Computes a pairwise distance list for two sets of coordinate pairs 
    '''
    d = []
    for i in X1:
        for j in X2:
            d.append(math.sqrt((i[0]-j[0])**2 +(i[1]-j[1])**2)) 
    return d
    
def GenRandomSet(length, width = 10.0):
    """
    Generates a random set of 'length' coordinate pairs over a square interval.
    """
    X = []
    for i in range(length):
        X.append((np.random.rand(2)-1.)*width/2.0)
    return X
    
    


import ParseFermi
events110 = ParseFermi.Import_File('photons.txt', energyRange = (70000,110000),lonRange=(-10,10),latRange = (-10,10))
events130 = ParseFermi.Import_File('photons.txt', energyRange = (120000,140000),lonRange=(-10,10),latRange = (-10,10))
events150 = ParseFermi.Import_File('photons.txt', energyRange = (150000,300000),lonRange=(-10,10),latRange = (-10,10))

X1,X2,X3 = [],[],[]
for i in events110:
    X1.append((i[1],i[2]))
for i in events130:
    X2.append((i[1],i[2]))
for i in events150:
    X3.append((i[1],i[2]))    





#print len(X1),len(X2),len(X3)
#CorrelateSimulations()
#CorrelateSimulations2() # No Pulsars


#Correlation(pulsar,NFW,NFWDECAY)
 

#X1,X2,X3 = GenRandomSet(100,width = 20.),GenRandomSet(250,width = 20.),GenRandomSet(500,width = 20.)
#RandomCorrelations(nPoints = 5000)
Correlation(X1,X2,X3)




