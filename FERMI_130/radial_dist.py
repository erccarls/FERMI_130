#===============================================================================
# This is intended to look at the radial distribution of an input data set
# Author: Eric Carlson
# Date: Oct 3, 2012
#===============================================================================


import numpy as np
import math
import matplotlib.pyplot as plt  #@UnresolvedImport
import matplotlib.image as mpimg #@UnresolvedImport
import matplotlib.cm as cm #@UnresolvedImport
import matplotlib, scipy #@UnresolvedImport


import csv
import numpy as np
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

def GetDensityProfile(X,Y,width = 1.0, startStop = (0,10.), numSteps = 100):
    '''Gets the photon density at a given radius for an input set of event coordinates'''
    
    R = []
    # Calc photon radius 
    for i in range(len(X)):
        R.append(math.sqrt(X[i]**2 + Y[i]**2))
    
    # Calc Density profile based on sliding window
    radius = np.linspace(startStop[0],startStop[1],numSteps)
    density = []
    counts = []
    
    bins = []
    for i in range(len(radius)-1):
        count = 0
        r_low = radius[i]-width/2.0
        if r_low < 0: 
            r_low = 0
        r_high = radius[i]+width/2.0
        for j in R:
            if (r_low <j < r_high):
                count+=1.0
        area = math.pi*(r_high**2-r_low**2)
        bins.append((r_low+r_high)/2.0)
        density.append(float(count) /area)
        counts.append(count)
        
    return (bins, density, counts)
        

        
def PlotDensityProfiles(list):
    """Takes a list of input tuples (radius, density, name) and plots them"""
    
    from matplotlib.font_manager import FontProperties
    fontP = FontProperties()
    fontP.set_size('small')
    
    plt.subplot(2,1,1)    
    for i in list:
        if 'Fermi' in i[3]:
            plt.semilogy(i[0],i[1],label = i[3])
        else:
            plt.semilogy(i[0],i[1],'--',label = i[3])

    plt.xlabel(r'$r_{GC} [^\circ]$')
    plt.ylabel(r'$\rho$')
    plt.legend(loc=1, ncol=1, fancybox=True, shadow=False,prop=fontP,borderaxespad=0.,labelspacing = .2)
    
    
    plt.subplot(2,1,2)
    for i in list:
        if 'Fermi' in i[3]:
            plt.semilogy(i[0],i[2],label = i[3])
        #else:
            #plt.semilogy(i[0],i[2],'--',label = i[3])
    plt.xlabel(r'$r_{GC} [^\circ]$')
    plt.ylabel(r'$N_\gamma$')

    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages('ClusterVsRadius.pdf')
    plt.savefig(pp, format='pdf')
    print "Figures saved to ", 'ClusterVsRadius'+ '.pdf\n',
    pp.close()

    plt.show()
        
    

def Run():
    import MC
    plotList = []

    Fermilow = Import_File('photons.txt', energyRange = (30000,100000),lonRange=(-15,15),latRange = (-15,15))    
    Fermi110 = Import_File('photons.txt', energyRange = (100000,120000),lonRange=(-15,15),latRange = (-15,15))
    Fermi130 = Import_File('photons.txt', energyRange = (120000,140000),lonRange=(-15,15),latRange = (-15,15))
    Fermi150 = Import_File('photons.txt', energyRange = (140000,180000),lonRange=(-15,15),latRange = (-15,15))
    
    x = []
    y = []
    for i in Fermilow:
        x.append(i[1])
        y.append(i[2])
    (r, rho, counts) = GetDensityProfile(x,y,width = .8, startStop = (0,10.), numSteps = 100)
    rho = np.array(rho)
    rho = np.array(rho)/np.average(rho)
    plotList.append((r,rho,counts, 'Fermi 30-100'))
        
        
    x = []
    y = []
    for i in Fermi110:
        x.append(i[1])
        y.append(i[2])
    (r, rho, counts) = GetDensityProfile(x,y,width = .8, startStop = (0,10.), numSteps = 100)
    rho = np.array(rho)
    rho = np.array(rho)/np.average(rho)
    plotList.append((r,rho,counts, 'Fermi 110'))
    
    x = []
    y = []
    for i in Fermi130:
        x.append(i[1])
        y.append(i[2])
    (r, rho, counts) = GetDensityProfile(x,y,width = .8, startStop = (0,10.), numSteps = 100)
    rho = np.array(rho)/np.average(rho)
    plotList.append((r,rho,counts, 'Fermi 130'))
    
    x = []
    y = []
    for i in Fermi150:
        x.append(i[1])
        y.append(i[2])
    (r, rho, counts) = GetDensityProfile(x,y,width = .8, startStop = (0,10.), numSteps = 100)
    rho = np.array(rho)/np.average(rho)
    plotList.append((r,rho,counts, 'Fermi 150'))

    counts = 0
    NFW  = MC.Gen_Annihilation_Profile(10, 200,('NFW',23.5,1.0))
    x = np.linspace(0, 10, 200)
    y = []
    for i in NFW:
        y.append(i)
    y = np.array(y)/np.average(y)
    plotList.append((x,y,counts, 'NFW'))
    
    NFW  = MC.Gen_Annihilation_Profile(10, 200,('NFWDECAY',23.5,1.0))
    x = np.linspace(0, 10, 200)
    y = []
    for i in NFW:
        y.append(i)
    y = np.array(y)/np.average(y)
    plotList.append((x,y,counts, 'NFWDECAY'))
    
    NFW  = MC.Gen_Annihilation_Profile(10, 200,('EIN',20,0.17))
    x = np.linspace(0, 10, 200)
    y = []
    for i in NFW:
        y.append(i)
    y = np.array(y)/np.average(y)
    plotList.append((x,y,counts, 'Einasto'))
    
    NFW  = MC.Gen_Annihilation_Profile(10, 200,('PULSAR',))
    x = np.linspace(0, 10, 200)
    y = []
    for i in NFW:
        y.append(i)
    y = np.array(y)/np.average(y)
    plotList.append((x,y,counts, 'PULSAR'))
    
    
    
    PlotDensityProfiles(plotList)
    
    
Run()



                
                 
    
    
    
    
    
    