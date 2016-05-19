"""
To generate weighted partition functions to be used in 
fitting nearest neighbor models. These can be read into 
Ising.py to be used in minimization routines.

Will change methods to support additional models.
Develop base class and 
"""

import sympy as sp
import numpy as np
import math
import csv
import sys
csv.field_size_limit(sys.maxsize)


class PartitionFunction(object):
    def __init__(self, R='kcal/mol', T='K', **kwargs):
        #initialize first and last matrices for all calculations
        self.begin_matrix = sp.Martix([[0,1]])
        self.end_matrix = sp.Matrix([[1],[1]])
        
        #have a dict with key:values for args - support R, T (K) **is this even important here? no, i think not
        #RT = (float(0.001987)*float(298.15)) #Use R for kcal/mol and 25deg C in K
        #check for requirements and set to defaults
        if R not in kwargs:
            #default to kcal/mol
            kwargs[R] = float(0.001987)
        if T not in kwargs:
            kwargs[T] = float(298.15)        
                    
#have create_function() method, support single interface, single intrinsic value

RT = sp.Symbol('RT')
#set up definitions for intrinsic terms
#dGNa = sp.Symbol('dGNa')
#dGNb = sp.Symbol('dGNb')
dGA = sp.Symbol('dGA')
mi = sp.Symbol('mi')
mii = sp.Symbol('mii')
x = sp.Symbol('x')
#dGB = sp.Symbol('dGB')
#dGCb = sp.Symbol('dGCb')
#dGS = sp.Symbol('dGS')
#Kna = sp.Symbol('Kna')
#Knb = sp.Symbol('Knb')
Ka = sp.Symbol('Ka')
#Kb = sp.Symbol('Kb')
#Kcb = sp.Symbol('Kcb')
#Ks = sp.Symbol('Ks')

#set up definitions for coupling terms
#dGAB = sp.Symbol('dGAB')
dGBA = sp.Symbol('dGBA')
#dGBNa = Symbol('dGBNa')
#dGANb = Symbol('dGANb')
#Wab = sp.Symbol('Wab')
Wba = sp.Symbol('Wba')
#Wbna = sp.Symbol('Wbna')
#Wanb = sp.Symbol('Wanb')
np.exp = sp.Symbol('np.exp')

#define matricies to be used to calculate partition functions
#Na = sp.Matrix([[(Kna*Wbna),1],[Kna,1]])
#Nb = sp.Matrix([[(Knb*Wanb),1],[Knb,1]])
A = sp.Matrix([[(Ka*Wba),1],[Ka,1]])
#B = sp.Matrix([[(Kb*Wab),1],[Kb,1]])
#S = sp.Matrix([[(Ks*Wba),1],[Ks,1]])
#Cb = sp.Matrix([[(Kcb*Wab),1],[Kcb,1]])

#set up partition functions for each construct
#in total, there are a lot of proteins:
#[['BAB2A'], ['BAB3A'], ['BAB4A'], ['AB2ACb'], ['AB3ACb'], ['AB4ACb'], ['BAB2ACb'], ['BAB3ACb'], ['BAB4ACb'], ['NbAB2A'], ['NbAB3A'], ['NbAB4A'], ['NaBAB2'], ['NaBAB3'], ['NaBAB4'], ['NaBABACb'], ['NaBAB2ACb'], ['NaBAB3ACb'], ['NaBAB4ACb'], ['NaBAB2A'], ['NaBAB3A'], ['NaBAB4A']]

#B(AB)nA type
z_BAB2A = begin*A*A*A*A*A*A*end
z_BAB3A = begin*A*A*A*A*A*A*A*A*end
z_BAB4A = begin*A*A*A*A*A*A*A*A*A*A*end
#(AB)nA type
z_AB2A = begin*A*A*A*A*A*end
z_AB3A = begin*A*A*A*A*A*A*A*end
z_AB4A = begin*A*A*A*A*A*A*A*A*A*end
#(AB)nS type
z_ABS = begin*A*A*A*end
z_AB2S = begin*A*A*A*A*A*end
z_AB3S = begin*A*A*A*A*A*A*A*end
z_AB4S = begin*A*A*A*A*A*A*A*A*A*end
#B(AB)n type
z_BAB = begin*A*A*A*end
z_BAB2 = begin*A*A*A*A*A*end
z_BAB3 = begin*A*A*A*A*A*A*A*end
z_BAB4 = begin*A*A*A*A*A*A*A*A*A*end
#B(AB)nS type
z_BABS = begin*A*A*A*A*end
z_BAB2S = begin*A*A*A*A*A*A*end
z_BAB3S = begin*A*A*A*A*A*A*A*A*end
z_BAB4S = begin*A*A*A*A*A*A*A*A*A*A*end

#Here, the partitionlist will contain the individual partition functions for each construct, followed by the total number of helices,
#and finally, the capping type.  All of this information will be utilized in the following loop where the fitting functions for each protein will be developed
#partitionlist = [[z_AB4ACb, 10, 'Cb', 'AB4ACb']]   
partitionlist = [[z_BAB2A, 6, 'AB', 'BAB2A'], [z_BAB3A, 8, 'AB', 'BAB3A'], [z_BAB4A, 10, 'AB', 'BAB4A'], [z_AB2A, 5, 'AB', 'AB2A'], [z_AB3A, 7, 'AB', 'AB3A'], [z_AB4A, 9, 'AB', 'AB4A'], [z_BAB, 3, 'AB', 'BAB'], [z_BAB2, 5, 'AB', 'BAB2'], [ z_BAB3, 7, 'AB', 'BAB3'], [ z_BAB4, 9, 'AB', 'BAB4'], [z_ABS, 3, 'ABS', 'ABS'], [z_AB2S, 5, 'ABS', 'AB2S'], [z_AB3S, 7, 'ABS', 'AB3S'], [z_AB4S, 9, 'ABS', 'AB4S'], [z_BABS, 4, 'ABS', 'BABS'], [z_BAB2S, 6, 'ABS', 'BAB2S'], [z_BAB3S, 8, 'ABS', 'BAB3S'], [z_BAB4S, 10, 'ABS', 'BAB4S']]   

thetalist = []

for i in range(0,len(partitionlist)):
    y = np.array(partitionlist[i][0])
    pf = y.item(0)

    #define partial derivatives for development of fraction folded functions
    pdKa = sp.diff(y.item(0), Ka)
    #pdKb = sp.diff(y.item(0), Kb)
    #pdKs = sp.diff(y.item(0), Ks)
    #pdWab = sp.diff(y.item(0), Wab)
    pdWba = sp.diff(y.item(0), Wba)
    
    protein = partitionlist[i][3]
    
    #define partition function for construct and make necessary substitutions to put into format for fitting function
    q = pf.subs({Ka:(np.exp(-((dGA - (mi*x))/RT))), Wba:(np.exp(-((dGBA - (mii*x))/RT)))})
    
    #make substitutions in partial derivative terms necessary for fitting function
    pKa = pdKa.subs({Ka:(np.exp(-((dGA - (mi*x))/RT))), Wba:(np.exp(-((dGBA - (mii*x))/RT)))})
    #pKb = pdKb.subs({Ka:(np.exp(-((dGA - (mi*x))/RT))), Kb:(np.exp(-((dGB - (mi*x))/RT))), Ks:(np.exp(-((dGS - (mi*x))/RT))), Wab:(np.exp(-((dGAB - (mii*x))/RT))), Wba:(np.exp(-((dGBA - (mii*x))/RT)))})
    #pKs = pdKs.subs({Ka:(np.exp(-((dGA - (mi*x))/RT))), Kb:(np.exp(-((dGB - (mi*x))/RT))), Ks:(np.exp(-((dGS - (mi*x))/RT))), Wab:(np.exp(-((dGAB - (mii*x))/RT))), Wba:(np.exp(-((dGBA - (mii*x))/RT)))})
    #pWab = pdWab.subs({Ka:(np.exp(-((dGA - (mi*x))/RT))), Kb:(np.exp(-((dGB - (mi*x))/RT))), Ks:(np.exp(-((dGS - (mi*x))/RT))), Wab:(np.exp(-((dGAB - (mii*x))/RT))), Wba:(np.exp(-((dGBA - (mii*x))/RT)))})
    pWba = pdWba.subs({Ka:(np.exp(-((dGA - (mi*x))/RT))), Wba:(np.exp(-((dGBA - (mii*x))/RT)))})
    
    if partitionlist[i][2] == "ABS":
        nhelix = partitionlist[i][1]
        nint = partitionlist[i][1] - 1
        theta = ((1/((nhelix)*q))*(((np.exp(-((dGA - (mi*x))/RT)))*pKa)))
          
    elif partitionlist[i][2] == "AB":
        nhelix = partitionlist[i][1]
        nint = partitionlist[i][1] - 1
        theta = ((1/((nhelix)*q))*(((np.exp(-((dGA - (mi*x))/RT)))*pKa)))
    print "simplyfying..."
    thetalist.append([protein, sp.simplify(theta)])

tlist = [item for item in thetalist]

with open("PartitionFunctions_c34PR_mi_miip1_Regan.csv", "wb") as n:
    writer = csv.writer(n, delimiter=',')
    writer.writerows(tlist)    
            
           
        
