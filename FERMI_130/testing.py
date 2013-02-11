#import numpy as np
#import matplotlib.pyplot as plt  #@UnresolvedImport
#
#data = np.random.poisson(9,1000)
#values, bins = np.histogram(data, 10)
#
#print len(bins[:-1]) 
#print len(values)
#
#fig = plt.figure(1)
#fig.add_subplot(121)
#
#plt.bar(bins[:-1], values, fill = False)
#
#fig.add_subplot(122)
#
#plt.bar(bins[:-1], values/2., fill = True)
#
#fig.add_subplot(121)
#plt.bar(bins[:-1], values/2., fill = True)
#
#
#plt.show()


import multiprocessing as mp
from multiprocessing import pool

import numpy as np

def f(x):
    np.random.seed()
    for i in range(0,1000):
        np.random.poisson(lam = 1.)
    return np.random.ranf()
    
if __name__ == '__main__':
    
    #x = [1,2,3,4]
    p = pool.Pool(6)
    
    
    
    
    x = np.random.ranf(200)
    
    test = map(f,x)
    #print test
    
    print test.count(test[0])
    
    test = p.map(f, x)
    
    print test.count(test[5])
    
    


