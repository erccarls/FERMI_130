'''
I've attached a file containing every photon above 100 GeV observed within 10 degrees 
of Sgr A*. Each photon is listed on one row, and the columns are as follows:

0. photon energy, 
1. galactic longitude (degrees), 
2. galactic latitude (degrees), 
3. theta value (degrees, instrument coordinates), 
4. phi angle(degrees, instrument coordinates), 
5. zenith angle (degrees),
6. time of event (seconds of MET time), 
7. event class of photon, 
8. conversion type (front or back converting event).
'''

import csv
import numpy as np
import MCSTATS
import matplotlib.pyplot as plt  #@UnresolvedImport
import matplotlib.image as mpimg #@UnresolvedImport
import matplotlib.cm as cm
import matplotlib, scipy
import pickle,sys
from matplotlib.backends.backend_pdf import PdfPages

def Import_File(CSVFile,energyRange = (-100000,10000000), lonRange=(-1000,1000),latRange = (-1000,1000)):
    '''
    returns list of events with the following fields:
    0. Energy
    1. Galactic longitude
    2. Galactic latitude
    3. Time of event
    '''
    values= []
    csvread = csv.reader(open(CSVFile, 'rb'),delimiter = ' ' )
    for i in csvread:
        if (float(i[0])>energyRange[0] and float(i[0])<energyRange[1] and float(i[2])>latRange[0] and float(i[2])<latRange[1] and (float(i[1])>360.+lonRange[0] or float(i[1])<lonRange[1])):
            lon = 0 
            if float(i[1]) > 180: 
                lon = float(i[1])-360.
            else: 
                lon = float(i[1])  
            values.append((float(i[0]),lon,float(i[2]),float(i[5])))
    return values
    
    
def Run_Scan2(NumEvents = 40,fileOut = 'scan_48.pickle',eventWidth = 10.):
    start, stop, inc = 30000,200000, 10000
    X,Y,XERR,YERR = [],[],[],[]
    window = 10      # Starting energy window width
    windowInc = 1.06
    count = 0
    numTrials = len(range(start,stop,inc))
    outputEvents = [] # (mean,width)
    for i in range(start,stop,inc):
        j = 0
        events = []
        xerr = 0
        while (len(events) < NumEvents and j < 5000):
            events = Import_File('photons.txt', energyRange = (i-window*i/1.0e4-windowInc**j,i+window*i/1.0e4+windowInc**j),lonRange=(-eventWidth/2.,eventWidth/2.),latRange = (-eventWidth/2.,eventWidth/2.))
            j+=1
            #print (i-(window-windowInc**j),i+(window+windowInc**j))
            xerr = window*i/1.0e4+windowInc**j
            if (j==5000):
                print 'Not enough events in range.'      
        outputEvents.append(events)
        X.append(i/1.0e3) # append mean energy
        XERR.append(xerr/1.0e3) # append bin widtherror
        sys.stdout.write('\nEnergy:' + str(count +1 )+'/'+str(numTrials) + '  Event Count ' + str(len(events))) 
        sys.stdout.flush()
        count+=1
    pickle.dump((X,outputEvents,XERR), open(fileOut,'wb'))

def Plot_Scan(fileIn = 'scan_48.pickle' ,fileOut='ClusteringScan',method = 'incon', eventWidth = 10, clusterOut= ''):
    scan = pickle.load(open(fileIn, "r" ))
    X,outputEvents,XERR = scan[0],scan[1],scan[2]
    
    
    #===========================================================================
    # Cluster Analysis
    #===========================================================================
    scales = []
    sigs   = []
    count  = []
    members= []
    sigmas = []
    
    for events in outputEvents:
        x,y = [],[]
        for q in events:
            x.append(q[1])
            y.append(q[2])
        # Run cluster Analysis
        

        (cluster_Scale, cluster_S, cluster_Count, cluster_Members,cluster_stdevs) = MCSTATS.DBSCAN_Compute_Clusters(((x,y),))
        # Append to Lists
        if cluster_Count == []:
            scales.append(float('NaN'))
            sigs.append(float('NaN'))
            count.append(float('NaN'))
            members.append(float('NaN'))
            sigmas.append(float('NaN'))
        else:
            scales.append(cluster_Scale[0])
            sigs.append(cluster_S[0])
            count.append(cluster_Count[0])
            members.append(cluster_Members[0])
            sigmas.append(cluster_stdevs[0])
    
    #===========================================================================
    # Plotting
    #===========================================================================
    fig = plt.figure(figsize=(8,12))
    #ax1 = fig.add_subplot(2,1,1)

    fig.add_subplot(4,1,2)    
    plt.errorbar(X, scales, sigmas, XERR, fmt='o')
    plt.xlabel('Photon Energy [GeV]')
    plt.ylabel(r'$r_{cluster} [^\circ$]')

    plt.axvline(x=130,color='r', label = '130 GeV')
    plt.axvline(x=111,color='r',label = '111 GeV')
#    plt.text(132, .6, '130 GeV', rotation=90)
#    plt.text(113, .6, '111 GeV', rotation=90)
#    plt.text(220, 2.8, method)
#    plt.text(160, 2.8, 'Width ' + str(eventWidth) + r'$^\circ$')

    plt.ylim((0,.02))
    plt.xlim((0,200))
    
    fig.add_subplot(4,1,1)
    plt.scatter(X,sigs,marker = 'o')
    plt.xlabel('Photon Energy [GeV]')
    plt.ylabel(r'$S$')
    plt.axvline(x=130,color='r', label = '130 GeV')
    plt.axvline(x=111,color='r',label = '111 GeV')
    
    fig.add_subplot(4,1,3)
    plt.scatter(X,count,marker = 'o')
    plt.xlabel('Photon Energy [GeV]')
    plt.ylabel(r'$N_{clusters}$')
    plt.axvline(x=130,color='r', label = '130 GeV')
    plt.axvline(x=111,color='r',label = '111 GeV')
    
    fig.add_subplot(4,1,4)
    plt.scatter(X,members,marker = 'o')
    plt.xlabel('Photon Energy [GeV]')
    plt.ylabel(r'$N_{members}$')
    plt.axvline(x=130,color='r', label = '130 GeV')
    plt.axvline(x=111,color='r',label = '111 GeV')

    
    pp = PdfPages(fileOut + '.pdf')
    plt.savefig(pp, format='pdf')
    print "Figures saved to ", str(fileOut)+ '.pdf\n'
    pp.close()
    plt.savefig(fileOut + '.png', format='png')
    #plt.show()

#Run_Scan2(NumEvents=48,fileOut = 'scan_48.pickle')
Plot_Scan(fileIn = 'scan_40.pickle')

def plotAngScales(fileOut = 'clusters'):
    import os
    files = [file for file in os.listdir('./') if (file.lower().endswith('.ang_scan'))]
    WIDTH = []
    SCALE = [[],[],[],[]]
    for file in files:
        width = float(file[0:2])
        WIDTH.append(width)
        scan = pickle.load(open(file, "r" ))
        X,Y,XERR,YERR = scan[0],scan[1],scan[2],scan[3]
        b = [[], [] ,[] ,[]]
        for i in range(len(X)):
            if (45 < X[i]< 85):
                b[0].append(Y[i])
            elif (85 < X[i]< 125):
                b[1].append(Y[i])
            elif (125 < X[i]< 165):
                b[2].append(Y[i])
            elif (165 < X[i]< 205):
                b[3].append(Y[i])
        for j in range(4):
            SCALE[j].append(np.mean(b[j]))
    
    plt.figure(figsize=(8,4))
    c = ['r','b','g','y']
    for i in range(4):
        plt.scatter(WIDTH,SCALE[i],c=c[i])
    plt.legend(['45-85 GeV','85-125 GeV', '125-165 GeV','165-205 GeV'])
    plt.xlabel(r'Box Width [$^\circ$]')
    plt.ylabel(r'Clustering Scale [$^\circ$]')
    pp = PdfPages(fileOut + '.pdf')
    plt.savefig(pp, format='pdf')
    print "Figures saved to ", str(fileOut)+ '.pdf\n'
    pp.close()
    plt.savefig(fileOut + '.png', format='png')
    plt.show()

#plotAngScales(fileOut = 'width_vs_clusterScale')


#for i in np.arange(5,42,1):
#for i in np.arange(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])):
#    evtCount = int(30.*i*i/(100.))
    #Run_Scan2(NumEvents=evtCount,fileOut = 'scan_'+ str(evtCount) + '_' + str(i) + 'x' + str(i) + '.pickle', eventWidth=i)
#    Plot_Scan(fileIn = 'scan_'+ str(evtCount) + '_' + str(i) + 'x' + str(i) + '.pickle', fileOut = 'RADIUS_ClusteringScan_'+ str(i),method='bin',eventWidth = i, clusterOut=str(i) + '.ang_scan' )

#Run_Scan2(NumEvents=30,fileOut = 'scan_30.pickle')
#Run_Scan2(NumEvents=20,fileOut = 'scan_20.pickle')
#Run_Scan2(NumEvents=40,fileOut = 'scan_40.pickle')
#Run_Scan2(NumEvents=50,fileOut = 'scan_50.pickle')
#Plot_Scan(fileIn = 'scan_40.pickle', fileOut = 'ClusteringScan_40_events_incon',method='incon')
#Plot_Scan(fileIn = 'scan_30.pickle', fileOut = 'ClusteringScan_30_events_median',method='median')
#Plot_Scan(fileIn = 'scan_30.pickle', fileOut = 'ClusteringScan_30_events_window',method='window')
#Plot_Scan(fileIn = 'scan_30.pickle', fileOut = 'ClusteringScan_30_events_bin',method='bin')

#for i in (20,30,40,50):
#    Run_Scan(NumEvents=i,fileOut = 'scan_' + str(i) + '.pickle')
#    Plot_Scan(fileIn = 'scan_' + str(i) + '.pickle', fileOut = 'ClusteringScan_' + str(i) + '_events.pdf')




