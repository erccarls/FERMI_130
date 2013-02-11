#######################################################
# This is intended to extrapolate clustering scaling relations with photon count.
#######################################################

import pickle
import MC,MCSTATS
import numpy as np
import sys
import matplotlib.pyplot as plt  #@UnresolvedImport
import matplotlib.image as mpimg #@UnresolvedImport
import matplotlib.cm as cm #@UnresolvedImport
import matplotlib, scipy #@UnresolvedImport


## Global parameters
#angularSize = 10.0 # box surrounding galactic center
#size = 300  # size of output ratemap
#numSims = 100 # Number of simulations to use.
#outputSize = 300
#
#
#
#
#runList = (('NFWRateMap.pickle',('NFW',23.5,1.0)),
#           ('EINRateMap.pickle',('EIN',20,0.17)),
#           ('FLATRateMap.pickle',('FLAT',0,0)),
#           ('NFWDECAYRateMap.pickle',('NFWDECAY',23.5,1.0)),
#           ('PULSARRateMap.pickle',('PULSAR',6)),
#           ('PULSARRateMap.pickle',('PULSAR',10)),
#           ('PULSARRateMap.pickle',('PULSAR',20)),
#           ('PULSARRateMap.pickle',('PULSAR',40)))
#
#ResultsList = []
#for run in runList:
#    fileOut = run[0]
#    profile = run[1]
#    
#
#    scan_event  = []
#    scan_height = []
#    scan_sigma  = []
#    scan_num    = []
#    for x in range(4,40,2):
#        eventCount = int(x*x)
#        print 'Setting Photon Count to ' + str(eventCount)
#        
#        mcSims = []
#        if profile[0] !='PULSAR':
#            mcSims = MC.RUN(numSims, fileOut, numPhotons=eventCount,angularSize=angularSize, outputSize=outputSize, mcList='clustering_' + profile[0] + '_' + str(eventCount) + '.pickle')
#        else:
#            mcSims = MC.RUN_PULSAR(numSims, fileOut, numPhotons=eventCount,numPulsars=profile[1],angularSize=angularSize, outputSize=outputSize, mcList='clustering_' + profile[0] + str(profile[1]) + '_' + str(eventCount) + '.pickle')
#        
#        
#        print 'Running Cluster Analysis...'
#        
#        mcSims = []
#        if profile[0] !='PULSAR':
#            mcSims = pickle.load(open('clustering_' + profile[0] + '_' + str(eventCount) + '.pickle', 'rb')) 
#        else:
#            mcSims = pickle.load(open('clustering_' + profile[0] + str(profile[1]) + '_' + str(eventCount) + '.pickle', 'rb'))
#            
#        
#        meanHeight = []
#        coph = []
#        sigma= []
#        num = []
#        cent_X,cent_Y = [],[]
#        
#        for i in range(len(mcSims)):
#            #print 'Cluster ' ,i
#            xVec = MCSTATS.Pixel_to_Degrees(mcSims[i][0], outputSize, angularSize)
#            yVec = MCSTATS.Pixel_to_Degrees(mcSims[i][1], outputSize, angularSize)
#            clusterResults = MCSTATS.Compute_Clusters(xVec, yVec,method = 'average',searchMethod='incon')
#            meanHeight.append(clusterResults[0])
#            coph.append(clusterResults[1])
#            sigma.append(clusterResults[2])
#            num.append(clusterResults[3])
#            cent_X.append(np.average(xVec)-angularSize/2.)
#            cent_Y.append(np.average(yVec)-angularSize/2.)
#            #print clusterResults   
#            if (i%5==0):
#                sys.stdout.write('\r' + str(i+1)+'/'+str(len(mcSims)))
#                sys.stdout.flush()
#        print 
#        scan_event.append(eventCount)
#        scan_height.append(np.average(meanHeight))
#        scan_sigma.append(np.std(meanHeight))
#        scan_num.append(np.average(num))
#    if profile[0] != 'PULSAR':    
#        ResultsList.append((profile[0], scan_event,scan_height,scan_sigma,scan_num))
#    else:
#        ResultsList.append((profile[0] + str(profile[1]), scan_event,scan_height,scan_sigma,scan_num))
#pickle.dump(ResultsList, open('clusterScaling.pickle','wb'))


#######################################################
# Plotting
#######################################################
from matplotlib.backends.backend_pdf import PdfPages
ResultsList = pickle.load(open('clusterScaling.pickle','r'))



# Plot Clustering Scale vs Photon Count
for scan in ResultsList:
    plt.errorbar(scan[1], scan[2], scan[3],fmt='o',label=scan[0])

#plt.xlim(1e1,1e4)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Photon Count')
plt.ylabel(r'Clustering Scale $[^\circ]$')
plt.legend(ncol=3,prop={'size':6})


pp = PdfPages('clusterScaling.pdf')
plt.savefig(pp, format='pdf',dpi=300)
print "Figures saved to clusterScaling.pdf\n"
pp.close()


# Plot Clustering Scale vs Photon Count
plt.clf()
for scan in ResultsList:
    plt.plot(scan[1], scan[4],label=scan[0])

#plt.xlim(1e1,1e4)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Photon Count')
plt.ylabel(r'$N_{Clusters}$')
plt.legend(ncol=3,loc=4,prop={'size':6})

pp = PdfPages('clusterScaling2.pdf')
plt.savefig(pp, format='pdf',dpi=300)
print "Figures saved to clusterScaling2.pdf\n"
pp.close()






