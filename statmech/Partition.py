"""
To generate weighted partition functions to be used in 
fitting nearest neighbor models. These can be read into 
Ising.py to be used in minimization routines.

Will change methods to support additional models.
Develop base class and 
"""
from os import path
import sys
import errno
import sympy as sp
import numpy as np
import glob
import json
import multiprocessing as mp

from settings import * # will import the base directories of package

if CONSTANTS_DIR not in sys.path:
    sys.path.insert(0, CONSTANTS_DIR)
    
from Constants import gas_constants

class IsingModel():
    '''
    Class for generating Ising models. 
    '''
    
    def __init__(self, num_cpus=2):
        self.intrinsic_energies = {}
        self.extended_names = []
        self.model = ''
        self.num_cpus = num_cpus
        
    def ising1D(self, capping_scheme=['R', 'NR', 'NRC', 'RC'], intrinsic_term='R', max_length=6, min_length=1):
        self.model = '1D Ising'        
        # generate extended names and create symbolic intrinsic energy terms 
        for item in capping_scheme:     
            # define energy terms
            terms = list(item)
            for t in terms:
                if t not in self.intrinsic_energies:
                    self.intrinsic_energies[t] = sp.Symbol(t)
            
            # define extended names
            if item==intrinsic_term:
                self.extended_names.extend([item.replace(intrinsic_term, intrinsic_term*l) for l in np.arange(min_length,max_length+1) ])
            else:
                self.extended_names.extend([item.replace(intrinsic_term, intrinsic_term*l) for l in np.arange(min_length,max_length)])
        
        # now shuffle so that the length of calculations is distributed more evenly when parallelized (i.e. longer partition functions are not all on one core)
        np.random.shuffle(self.extended_names)
        self.model_pfs = [IsingPartitionFunction(item, intrinsic_energies=self.intrinsic_energies) for item in self.extended_names]
        pool = mp.Pool(self.num_cpus)
        self.partition_functions = pool.map(evaluate_pf, self.model_pfs)
        self.model_pfs = []
        for pf in self.partition_functions:
            c, theta = pf
            c.simplified_theta = theta
            self.model_pfs.append(c) 
       
    def ising2D(self, capping_scheme=['AB', 'ABC', 'BAB', 'BABC'], intrinsic_term='AB', max_length=6, min_length=1):
        return None
   
    def save_model(self, filename):
        filepath = path.join(MODELS_DIR, filename)
        
        if filename in self.existing_models:
            decision = raw_input('Caution! There is an existing model with this name. Would you like to overwrite? [y/n] : ')
            decision = decision.lower()
            if decision == 'n':
                filename = raw_input('Please supply an alternative filename : ')
        
        model = {filename: {'type': self.model_type, 'parameters': self.parameters, 'partition_functions': self.pfs}}
        
        with open(filepath, 'wb') as f:
            json.dump(unicode(model), f)
            
class IsingPartitionFunction(): 
    '''
    Class for generating partition functions for Ising models
    
    Will accept a string name, and compute the partition function 
    '''
    
    def __init__(self, name, intrinsic_energies=None, R='J/mol', T=293.15, model='homopolymer', dgdx='intrinsic', intrinsic_term='R'):
        
        if not intrinsic_energies:
            raise ValueError("intrinsic_energies is None. Please supply a symbolic dictionary of intrinsic energies")
        
        # initialize first and last matrices for all calculations, and all variables
        self.begin_matrix = sp.Matrix([[0,1]])
        self.end_matrix = sp.Matrix([[1],[1]])
        self.T = T
        self.R = gas_constants[R]['value']
        self.RT = self.R*self.T
        self.name = name
        self.parameters = ['dG'+i for i in intrinsic_energies.keys()]
        self.parameters.append('x')
        self.intrinsic_matricies = {}
        self.intrinsic_energies = intrinsic_energies
        x = sp.Symbol('x')
        
        if dgdx=='intrinsic':
            mi = sp.Symbol('mi')
            self.parameters.append('mi')
            self.deltaGsubs = {self.intrinsic_energies[i]: sp.exp( (sp.Symbol('dG'+i) - (mi*x)) / self.RT) for i in self.intrinsic_energies.keys()}
        
        elif dgdx=='interfacial':
            mii = sp.Symbol('mii')
            self.parameters.append('mii')
            print("Currently, only intrinsic energy differentiation is supported.")
            
        if model=='homopolymer':
            # create one coupling term 
            W = sp.Symbol('W')
            self.deltaGsubs[W] = sp.exp( -( sp.Symbol('dGii') / self.RT) )
            self.intrinsic_matrices = {i:sp.Matrix([[(self.intrinsic_energies[i]*W),1],[self.intrinsic_energies[i],1]]) for i in self.intrinsic_energies.keys()} # here matrices are mapped to each intrinsic term. ex. R : <R_matrix for multiplication>
            self.parameters.append('dGii')
                           
        elif model=='heteropolymer':
            # create two coupling terms
            try:
                Wab = sp.Symbol('Wab')
                Wba = sp.Symbol('Wba')
                self.deltaGsubs[Wab] = sp.exp(-(sp.Symbol('dGAB')/ self.RT))
                self.deltaGsubs[Wba] = sp.exp(-(sp.Symbol('dGBA')/ self.RT))
                self.parameters.append('dGAB')
                self.parameters.append('dGBA')
                
                self.intrinsic_matrices = {'A': sp.Matrix([[(self.intrinsic_energies['A'] * Wba),1],[self.intrinsic_energies['A'],1]]),
                                           'B': sp.Matrix([[(self.intrinsic_energies['B'] * Wab), 1], [self.intrinsic_energies['B'],1]])}
                
                for rep in self.intrinsic_energies.keys():
                    if rep not in self.intrinsic_matrices:
                        if rep=='N':
                            self.intrinsic_matrices['N'] = sp.Matrix([[(self.intrinsic_energies['N'] * Wab), 1], [self.intrinsic_energies['N'],1]])
                        if rep=='C':
                            self.intrinsic_matrices['C'] = sp.Matrix([[(self.intrinsic_energies['C'] * Wba),1],[self.intrinsic_energies['C'],1]])

            except:
                raise ValueError("Could not recognize intrinsic terms. Please check to ensure your capping scheme includes N or C as the definitions for N and C terminal caps, and the internal units are labeled as A and B")      
        
        else:
            raise ValueError("Could not recognize model. Currently, 1D 'homopolymer' and 'heteropolymer' Ising models are are supported")
        
    
    def evaluate(self):
        '''
        Evaluates self.name 
        
        Accepts: None
        Returns: a partition function representation for self.name
        '''
        nrep = len(self.name)
        unique_reps = list(set(self.name))
        
        pf = self.begin_matrix # begin pf
        for r in self.name:
            pf = pf*self.intrinsic_matrices[r] # multiply through
        pf = (pf*self.end_matrix)[0] # finish pf
        
        partials = {'pd_'+l: sp.diff(pf, self.intrinsic_energies[l]).subs(self.deltaGsubs) for l in unique_reps} # calculate partial derivatives with respect to each energy term             
        q = pf.subs(self.deltaGsubs) # now substitute energy terms with dGs
        add_term = np.sum([self.deltaGsubs[self.intrinsic_energies[r]]*partials['pd_'+r] for r in unique_reps]) #calculate additional_term based on repeats in name
        self.simplified_theta = sp.simplify(((1/(float(nrep)*q))*(add_term))) # calculate theta and simplify term for faster arithmetic 
        return self.simplified_theta
        
    def simplify(self):
        return None
    
    
    def name2pf(self, name):
        nrep = len(name)
        unique_reps = list(set(name))
        pf = self.begin_matrix #n begin pf
        for char in name:
            pf = pf*self.intrinsic_matrices[char] # multiply through
        pf = (pf*self.end_matrix)[0] # finish pf
        partial_derivs = {'pd_'+letter: sp.diff(pf, self.intrinsic_energies[letter]).subs(self.deltaGsubs) for letter in unique_reps}              
        q = pf.subs(self.deltaGsubs) # now substitute energy terms with dGs
        add_term = np.sum([self.deltaGsubs[self.intrinsic_energies[rep]]*partial_derivs['pd_'+rep] for rep in unique_reps]) #calculate additional_term based on repeats in name
        simplified_theta = sp.simplify(((1/(float(nrep)*q))*(add_term))) # calculate theta and simplify term for faster arithmetic 
        return (name, simplified_theta)
    
        print("Calculating partition %s functions.... \n\nNOTE: This is a lengthy operation and should be rewritten using dill to allow serialization for multiprocessing" % (len(self.extended_names)))   
        self.pfs = dict(map(name2pf, self.extended_names)) 
        print("Created %s partition functions using a %s model." % (str(len(self.extended_names)), model))
        self.model_type = model
    
    
#############################################################################
'''                    GLOBAL NAMESPACE FUNCS                             '''
#############################################################################        
  
def evaluate_pf(pf):
    '''
    Function to be called when parallelizing partition functions
    
    Accepts: an instance of the IsingPartitionFunction class
    Returns: an evaluated partition function
    '''
    theta = pf.evaluate()
    
    return pf, theta




          
'''                  
#have create_function() method, support single interface, single intrinsic value

RT = sp.Symbol('RT')
dGA = sp.Symbol('dGA')
mi = sp.Symbol('mi')
mii = sp.Symbol('mii')
x = sp.Symbol('x')
Ka = sp.Symbol('Ka')


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
            
'''           
if __name__ == '__main__':
    cpus = mp.cpu_count()
    model = IsingModel(num_cpus=cpus)
    
    model.ising1D()
    #pfg.save_model('homopolymer.json')     
