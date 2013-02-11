import pickle
import MC,MCSTATS,DBSCAN
import numpy as np
import sys, math

# Global parameters
angularSize = 10.0 # box surrounding galactic center
size = 300  # size of output ratemap
numSims = 10000 # Number of simulations to use.
numSims2 = 1000
numPhotons = 48
numPhotons2 = 500
outputSize = 300
#outputSize = 1000

def run_MC_NFW():
    """Generate Rate Map and Run Monte-Carlo NFW""" 
    #NFW Rate Map
    profile = ('NFW',23.5,1.0)
    fileOut = 'NFWRateMap.pickle'
    #map = MC.Gen_Annihilation_Map(angularSize, size, profile,fileOut)
    #MCSTATS.Plot_Rate_Map('NFWRateMap.pickle', angularSize, fileOut)
    # MC
    mcSims = MC.RUN(numSims, fileOut, numPhotons=numPhotons,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_NFW_48.pickle',flatLevel = 0.25)
    mcSims = MC.RUN(numSims2, fileOut, numPhotons=numPhotons2,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_NFW_500.pickle',flatLevel = 0.25)

def run_MC_EIN():
    """Generate Rate Map and Run Monte-Carlo Einasto"""
    # Einasto rate map
    profile = ('EIN',20,0.17)
    fileOut = 'EINRateMap.pickle'
    #map = MC.Gen_Annihilation_Map(angularSize, size, profile,fileOut)
    #MCSTATS.Plot_Rate_Map('EINRateMap.pickle', angularSize, fileOut)
    # MC
    mcSims = MC.RUN(numSims, fileOut, numPhotons=numPhotons,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_EIN_48.pickle',flatLevel = 0.25)
    mcSims = MC.RUN(numSims2, fileOut, numPhotons=numPhotons2,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_EIN_500.pickle',flatLevel = 0.25)
    
def run_MC_FLAT():
    """Generate Rate Map and Run Monte-Carlo Flat"""
    # Einasto rate map
    profile = ('FLAT',0,0)
    fileOut = 'FLATRateMap.pickle'
    #map = MC.Gen_Annihilation_Map(angularSize, size, profile,fileOut)
    #MCSTATS.Plot_Rate_Map('FLATRateMap.pickle', angularSize, fileOut)
    # MC
    mcSims = MC.RUN(numSims, fileOut, numPhotons=numPhotons,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_FLAT_48.pickle',flatLevel = 0.25)
    mcSims = MC.RUN(numSims2, fileOut, numPhotons=numPhotons2,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_FLAT_500.pickle',flatLevel = 0.25)

def run_MC_NFWDECAY():
    """Generate Rate Map and Run Monte-Carlo NFW but not squared"""
    # Einasto rate map
    profile = ('NFWDECAY',23.5,1.0)
    fileOut = 'NFWDECAYRateMap.pickle'
    #map = MC.Gen_Annihilation_Map(angularSize, size, profile,fileOut)
    #MCSTATS.Plot_Rate_Map('NFWDECAYRateMap.pickle', angularSize, fileOut)
    # MC
    mcSims = MC.RUN(numSims, fileOut, numPhotons=numPhotons,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_NFWDECAY_48.pickle',flatLevel = 0.25)
    mcSims = MC.RUN(numSims2, fileOut, numPhotons=numPhotons2,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_NFWDECAY_500.pickle',flatLevel = 0.25)

def run_MC_PULSAR(numPulsars=6):
    """Generate Rate Map and Run Monte-Carlo NFW but not squared"""
    # Pulsars
    profile = ('PULSAR',)
    fileOut = 'PULSARRateMap.pickle'
    #map = MC.Gen_Annihilation_Map(angularSize, size, profile,fileOut)
    #MCSTATS.Plot_Rate_Map('PULSARRateMap.pickle', angularSize, fileOut)
    # MC
    mcSims = MC.RUN_PULSAR(numSims, fileOut, numPhotons=numPhotons,numPulsars=numPulsars,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_PULSAR_48_' + str(numPulsars) + '.pickle',flatLevel = 0.25)
    mcSims = MC.RUN_PULSAR(numSims2, fileOut, numPhotons=numPhotons2,numPulsars=numPulsars,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_PULSAR_500_' + str(numPulsars) + '.pickle',flatLevel = 0.25)

def run_CLUSTER_ANALYSIS(mcFile,clusterFile, numAnalyze = 0):
    """Run cluster analysis on mcFile output from MC.RUN."""
    print 'Loading MC File...'
    mcSims = pickle.load(open(mcFile, "r" ))
    print 'Running Cluster Analysis...'
    meanHeight = []
    coph = []
    sigma= []
    num = []
    cent_X,cent_Y = [],[]
    
    if ((numAnalyze == 0) or (numAnalyze > len(mcSims))):
        numAnalyze =len(mcSims)
    
    print 'Analyzing ' + str(numAnalyze) + ' simulations...'
    
    for i in range(numAnalyze):
        #print 'Cluster ' ,i
        xVec = MCSTATS.Pixel_to_Degrees(mcSims[i][0], outputSize, angularSize)
        yVec = MCSTATS.Pixel_to_Degrees(mcSims[i][1], outputSize, angularSize)
        clusterResults = MCSTATS.Compute_Clusters(xVec, yVec,method = 'average',searchMethod='incon')
        meanHeight.append(clusterResults[0])
        coph.append(clusterResults[1])
        sigma.append(clusterResults[2])
        num.append(clusterResults[3])
        cent_X.append(np.average(xVec)-angularSize/2.)
        cent_Y.append(np.average(yVec)-angularSize/2.)
        if (i%10==0):
            sys.stdout.write('\r' + str(i+1)+'/'+str(len(mcSims)))
            sys.stdout.flush()
    
    pickle.dump((meanHeight,coph,sigma,num),open(clusterFile,'wb'))
    print" File Output to " + clusterFile
    
    
    
    
    
    
    
    
def run_DBSCAN_ANALYSIS(mcFile,resultsFile, numAnalyze = 0, profile):
    """Run cluster analysis on mcFile output from MC.RUN."""
    print 'Loading MC File...'
    mcSims = pickle.load(open(mcFile, "r" ))
    print 'Running Cluster Analysis...'
    
    
    clusterResults = []
    for i in range(numAnalyze):
        #print 'Cluster ' ,i
        xVec = MCSTATS.Pixel_to_Degrees(mcSims[i][0], outputSize, angularSize)
        yVec = MCSTATS.Pixel_to_Degrees(mcSims[i][1], outputSize, angularSize)
        
    MCSTATS.DBSCAN_Compute_Clusters(mcSims, profile, S_cut=1.0, flatLevel=0.25, numAnalyze = numAnalyze)
        
        if (i%5==0):
            sys.stdout.write('\r' + str(i+1)+'/'+str(len(mcSims)))
            sys.stdout.flush()
    
    pickle.dump(clusterResults,open(resultsFile,'wb')
    print "File Output to " + clusterFile
    return
    
    
    
    
    
    
    
    
    
    
    
def Gen_Plots(clusterFile):
    # Generate Plots
    (meanHeight,coph,sigma,num) = pickle.load(open(clusterFile,'r'))
    MCSTATS.Height_Histogram((meanHeight,coph,sigma,num),'Analysis_' + clusterFile ,bins=100)
    #MCSTATS.Plot_MC_Centroids(cent_X, cent_Y, fileOut='centroids')

def Gen_Plots2():
    """48 photon plots"""
    inputList = (('NFW','Clusters_MCOut_NFW_48.pickle'),
                 ('Flat','Clusters_MCOut_FLAT_48.pickle'), 
                 ('Ein','Clusters_MCOut_EIN_48.pickle'),
                 ('NFW Decay','Clusters_MCOut_NFW_48.pickle'),
                 ('Pulsar 3 ','Clusters_MCOut_PULSAR_48_3.pickle'),
                 ('Pulsar 6','Clusters_MCOut_PULSAR_48_6.pickle'),
                 ('Pulsar 10','Clusters_MCOut_PULSAR_48_10.pickle'),
                 ('Pulsar 20','Clusters_MCOut_PULSAR_48_20.pickle'),
                 ('Pulsar 30','Clusters_MCOut_PULSAR_48_30.pickle')
                 )
    clusterResults = [[],[]]
    for i in inputList:
        (meanHeight,coph,sigma,num) = pickle.load(open(i[1],'r'))
        clusterResults[0].append(i[0])
        clusterResults[1].append(meanHeight)
    MCSTATS.Height_Histogram2(clusterResults, 'CompareClustering',bins =100)
    
def Gen_Plots3():
    """4800 photon plots"""
    inputList = (('NFW','Clusters_MCOut_NFW_4800.pickle'),
                 ('Flat','Clusters_MCOut_FLAT_4800.pickle'), 
                 ('Ein','Clusters_MCOut_EIN_4800.pickle'),
                 ('NFW Decay','Clusters_MCOut_NFW_4800.pickle'),
                 ('Pulsar 3 ','Clusters_MCOut_PULSAR_4800_3.pickle'),
                 ('Pulsar 6','Clusters_MCOut_PULSAR_4800_6.pickle'),
                 ('Pulsar 10','Clusters_MCOut_PULSAR_4800_10.pickle'),
                 ('Pulsar 20','Clusters_MCOut_PULSAR_4800_20.pickle'),
                 ('Pulsar 30','Clusters_MCOut_PULSAR_4800_30.pickle')
                 )
    clusterResults = [[],[]]
    for i in inputList:
        (meanHeight,coph,sigma,num) = pickle.load(open(i[1],'r'))
        clusterResults[0].append(i[0])
        clusterResults[1].append(meanHeight)
    MCSTATS.Height_Histogram2(clusterResults, 'CompareClustering_4800',bins =15)
    


#################
# Pulsar Testing
#################
#profile = ('PULSAR',)
#fileOut = 'PULSARRateMap.pickle'
##map = MC.Gen_Annihilation_Map(angularSize, size, profile,fileOut)
##MCSTATS.Plot_Rate_Map('PULSARRateMap.pickle', angularSize, fileOut)
## MC
#mcSims = MC.RUN_PULSAR(2, fileOut, numPhotons=numPhotons,numPulsars=6,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_PULSAR_48.pickle')
#MCSTATS.Plot_MC_Positions('MCOut_PULSAR_48.pickle', 'pulsarTest')


#############################
# Cluster 'incon' Testing
#############################
#profile = ('NFW',23.5,1.0)
#fileOut = 'NFWRateMap.pickle'
#mcSims = MC.RUN(1, fileOut, numPhotons=50,angularSize=angularSize, outputSize=outputSize, mcList='MC_TEST_NFW_4800.pickle')


#profile = ('PULSAR',)
#fileOut = 'PULSARRateMap.pickle'
#mcSims = MC.RUN_PULSAR(1, fileOut, numPhotons=100, numPulsars = 10,angularSize=angularSize, outputSize=outputSize, mcList='MC_TEST_NFW_4800.pickle')
#
#print 'Running Cluster Analysis...'
#meanHeight = []
#coph = []
#sigma= []
#num = []
#cent_X,cent_Y = [],[]
#for i in range(len(mcSims)):
#    #print 'Cluster ' ,i
#    xVec = MCSTATS.Pixel_to_Degrees(mcSims[i][0], outputSize, angularSize)
#    yVec = MCSTATS.Pixel_to_Degrees(mcSims[i][1], outputSize, angularSize)
#    clusterResults = MCSTATS.Compute_Clusters(xVec, yVec,method = 'average',searchMethod='incon')
#    meanHeight.append(clusterResults[0])
#    coph.append(clusterResults[1])
#    sigma.append(clusterResults[2])
#    num.append(clusterResults[3])
#    cent_X.append(np.average(xVec)-angularSize/2.)
#    cent_Y.append(np.average(yVec)-angularSize/2.)
#    print clusterResults



#Gen_Plots2()
#Gen_Plots3()




# Pulsars




def DBSCAN_PULSAR(numPulsars = 6, numPhotons =48,eps = 0.201719,n_cluster = 3,nCore = 2, S_cut=1.0, numSims = 100,flatLevel = .25):
    profile = ('PULSAR',)
    fileOut = 'PULSARRateMap.pickle'
    
    BGTemplate = pickle.load(open('BGRateMap.pickle','r'))
    
    mcSims = MC.RUN_PULSAR(numSims, fileOut, numPhotons=numPhotons,numPulsars=numPulsars,angularSize=angularSize, outputSize=outputSize, mcList='MCOut_PULSAR_48_' + str(numPulsars) + '.pickle',flatLevel = 0.25)
    
    
    cluster_Count = []   # Mean Number of clusters found
    cluster_Scale = []   # Mean Cluster Scale weighted by significance
    cluster_S = []       # Mean Significance weighted by number of cluster members
    cluster_Members = [] # Mean number of cluster Members
    
    for sim in mcSims:
        xVec = MCSTATS.Pixel_to_Degrees(sim[0], outputSize, angularSize)
        yVec = MCSTATS.Pixel_to_Degrees(sim[1], outputSize, angularSize)
        X = []
        for i in range(len(xVec)):
            X.append((xVec[i],yVec[i]))
        
        (clusters, clusterCount, noiseCount) = DBSCAN.RunDBScan(X, eps, n_cluster,nCore = nCore, plot=False)
        
    
        #===========================================================================
        # Compute Cluster Properties
        #===========================================================================
        clusterSigs = []
        clusterDist = []
        clusterSigmas = []
        clusterMembers = []
        for cluster in clusters:
            S =  DBSCAN.Compute_Cluster_Significance(cluster, BGTemplate, len(xVec),flatLevel = flatLevel)
            if S>S_cut:
                d, sigma = DBSCAN.Compute_Cluster_Scale(cluster)  # compute the cluster Scale
                clusterSigs.append(S)
                clusterDist.append(d)
                clusterSigmas.append(sigma)
                clusterMembers.append(len(cluster))
                
                #print S, d, len(cluster)
        # If no clusters found
        if len(clusterSigs)==0:
            continue
        #===========================================================================
        # Compute Weighted Mean and Weighted Stdev 
        #===========================================================================
        clusterScale = np.average(clusterDist, weights = clusterSigs)
        weights = np.array(clusterSigs)/np.sum(clusterSigs)
        stdev = 0
        for i in range(len(clusterSigmas)):
            stdev += weights[i]**2 * clusterSigmas[i]**2 
        stdev = math.sqrt(stdev)
        #print clusterScale, stdev
        # Append to Master List
        cluster_Scale.append(clusterScale)
        cluster_S.append(np.average(clusterSigs,weights = clusterMembers))
        cluster_Count.append(len(clusters))
        cluster_Members.append(np.mean(clusterMembers))
        
    return cluster_Scale, cluster_S, cluster_Count, cluster_Members
    
    
def DBSCAN_FERMI(energyRange = (120000,140000),eps=0.201719,n_cluster = 3,nCore = 2, S_cut = 1.0, flatLevel = .25):
    """Import Fermi Data and run DBSCAN"""
    import ParseFermi
    
    # Load fermi events 
    events = ParseFermi.Import_File('photons.txt', energyRange = energyRange,lonRange=(-5,5),latRange = (-5,5)) #@UndefinedVariable
    # Get into correct form
    X = []
    for i in events:
        X.append((i[1],i[2]))
        
    # Run DBSCAN
    (clusters, clusterCount, noiseCount) = DBSCAN.RunDBScan(X, eps, n_cluster,nCore = 2,plot=False)
    # Load BG Template
    BGTemplate = pickle.load(open('BGRateMap.pickle','r'))
    
    #===========================================================================
    # Compute Cluster Properties
    #===========================================================================
    clusterSigs = []
    clusterDist = []
    clusterSigmas = []
    clusterMembers = []
    for cluster in clusters:
        S =  DBSCAN.Compute_Cluster_Significance(cluster, BGTemplate, len(events),flatLevel =flatLevel)
        if S>S_cut:
            d, sigma = DBSCAN.Compute_Cluster_Scale(cluster)  # compute the cluster Scale
            clusterSigs.append(S)
            clusterDist.append(d)
            clusterSigmas.append(sigma)
            clusterMembers.append(len(cluster))
    #===========================================================================
    # Compute Weighted Mean and Weighted Stdev 
    #===========================================================================
    clusterScale = np.average(clusterDist, weights = clusterSigs)
    weights = np.array(clusterSigs)/np.sum(clusterSigs)
    stdev = 0
    for i in range(len(clusterSigmas)):
        stdev += weights[i]**2 * clusterSigmas[i]**2 
    stdev = math.sqrt(stdev)
        
    return clusterScale, np.average(clusterSigs,weights = clusterMembers),len(clusters),np.mean(clusterMembers)


    


#eps=0.0790533330292  # 68% containment
eps=0.201719         # 95% containment    
#eps=0.376319987699  # 99% containment

##===============================================================================
## Testing Flatness
##===============================================================================
#import numpy as np
#import matplotlib.pyplot as plt
#
#flatLevel = np.linspace(0, 1, 30)
#sig1 = []
#sig2 = []
#sig3 = []
#for i in flatLevel:
#   sig1.append(DBSCAN_FERMI( energyRange = (100000,120000),eps = eps, n_cluster = 3, nCore = 2,flatLevel = i)[1])
#   sig2.append(DBSCAN_FERMI( energyRange = (120000,140000),eps = eps, n_cluster = 3, nCore = 2,flatLevel = i)[1])
#   sig3.append(DBSCAN_FERMI( energyRange = (140000,180000),eps = eps, n_cluster = 3, nCore = 2,flatLevel = i)[1])
#plt.scatter(flatLevel, sig1,label = '110 GeV',c = 'b')
#plt.scatter(flatLevel, sig2,label = '130 GeV',c = 'r')
#plt.scatter(flatLevel, sig3,label = '150 GeV',c = 'g')
#plt.ylabel(r'Significance')
#plt.xlabel('Flatness Level')
#plt.legend()
#plt.show()
    
    


#DBSCAN_FERMI( energyRange = (120000,140000),eps = eps, n_cluster = 3, nCore = 2)

#clusterData = DBSCAN_PULSAR(numPulsars = 3, numPhotons =48,eps = eps,n_cluster = 3, nCore = 2,numSims = 1000, S_cut = 1.0)+ ('3 Pulsars',)
#MCSTATS.DBSCAN_STATS((clusterData,), fileOut='')




print 'runMC.py beginning in mode ' + sys.argv[1]
if (sys.argv[1] == '0'):
    run_MC_NFW()
if (sys.argv[1] == '1'):
    run_MC_EIN()
if (sys.argv[1] == '2'):
    run_MC_FLAT()
if (sys.argv[1] == '3'):
    run_MC_NFWDECAY()
if (sys.argv[1] == '4'):
    run_MC_PULSAR(int(sys.argv[2]))
    
if (sys.argv[1] == '10'):
    run_CLUSTER_ANALYSIS(sys.argv[2],'Clusters_'+sys.argv[2], numAnalyze = 100)
    Gen_Plots('Clusters_'+sys.argv[2])
    

