import subprocess as sp
import os
#import sys
#
#print 'Syntax: min max inc numProcesses'
#
#min, max, inc, numProc = float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4])
#
#numTrials = int((max-min)/inc)
#numPerProc = int(numTrials/numProc)
#
#for i in range(0,numProc):
#    sp.Popen('python ScanFermi.py ' + str(i*numTrials))

sp.Popen('python ScanFermi.py 5 10 1')
sp.Popen('python ScanFermi.py 10 15 1')
sp.Popen('python ScanFermi.py 15 20 1')
sp.Popen('python ScanFermi.py 20 26 1')
sp.Popen('python ScanFermi.py 26 32 1')
sp.Popen('python ScanFermi.py 32 37 1')
sp.Popen('python ScanFermi.py 37 42 1')

