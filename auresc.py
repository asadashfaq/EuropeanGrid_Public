#! /usr/bin/env python
from numpy import *
from pylab import is_numlike
from pylab import is_string_like
import sys
import numpy as np

alphas=[0.67431640625 , 0.7955078125 , 0.71640625 , 0.68291015625 , 0.755322265625 , 0.8490234375 , 0.70126953125 , 0.7873046875 , 0.739599609375 , 0.659423828125 , 0.64169921875 , 0.660888671875 , 0.690698242188 , 0.639697265625 , 0.688647460938 , 0.713232421875 , 0.6625 , 0.69501953125 , 0.71611328125 , 0.7541015625 , 0.817724609375 , 0.731689453125 , 0.64677734375 , 0.6509765625 , 0.696557617188 , 0.70712890625 , 0.7068359375]

class node:
    """ The node class loads data from the '.npz' countryfiles and stores it in these
        predefined structures. 

        NOTE: the unnormalized load signal can be queried by typing
        n.load*n.mean, but we have functions for querying wind and solar.
        Type n.get_wind() or n.get_solar.

        Individual values for gamma and alpha can also be set by n.set_alpha().
        DON'T MODIFY THEM DIRECTLY IN n.gamma OR n.alpha AS THEY WONT BE 
        RECALCLUATED UNLESS YOU TYPE self.__update"""
    def __init__(self,path,fileName,ID):
        self.id = ID
        data = load(path + fileName)
        self.gamma = 1.0#float(setup[ID][0])
        self.alpha = 0.7#float(setup[ID][1]) #Alpha should be expanded to a vector. completalpha() can be applied in update()
        self.load = 1000*array(map(double,data['L']))
        self.nhours = len(self.load)
        self.normwind = array(map(double,data['Gw']))
        self.normsolar = array(map(double,data['Gs']))
        self.mean = mean(self.load)
        self.balancing = np.zeros(self.nhours)
        self.curtailment = np.zeros(self.nhours)
        self.label = data['datalabel']
        self.mismatch = None
        self.colored_import = None #Set using self.set_colored_i_import()
        data.close()
        self.ESS = ESS(0.001,0.00001,0.5,self.nhours)
        self._update_()
        self.gen=np.zeros(self.nhours)


    def _update_(self):
        self.mismatch=(self.get_wind()+self.get_solar())-self.load
    
    def get_import(self):
        """Returns import power time series in units of MW."""
        return get_positive(get_positive(-self.mismatch) - self.balancing) #Balancing is exported if it exceeds the local residual load.
        
    def get_export(self):
        """Returns export power time series in units of MW."""
        return get_positive(self.mismatch) - self.curtailment #+ get_positive(self.balancing - get_positive(-self.mismatch))

    def get_localRES(self):
        """Returns the local use of RES power time series in units of MW."""
        return self.get_wind() + self.get_solar() - self.curtailment - self.get_export()

    def get_localBalancing(self):
        """Returns the local use of balancing power time series in units of MW."""
        return get_positive(-self.mismatch) - self.get_import()

    def get_wind(self):
        """Returns wind power time series in units of MW."""
        return self.mean*self.gamma*self.alpha*self.normwind

    def get_solar(self):
        """Returns solar power time series in units of MW."""
        return self.mean*self.gamma*(1.-self.alpha)*self.normsolar

    def set_gamma(self,gamma,operation='='):
        if operation == '=':
            self.gamma = gamma
        else:
            self.gamma *= gamma
        self._update_()
    
    def set_alpha(self,alpha):
        self.alpha=alpha
        self._update_()


class ESS:
    def __init__(self,pmax,emax,soco,nhours,inflow=0):
        self.SOC=np.zeros(nhours+1)
        self.power=np.zeros(nhours)
        self.SOC[0]=soco
        self.pmax=pmax
        self.pmin=-pmax
        self.emax=emax

    def update(self,t):
        self.SOC[t+1] = self.SOC[t] - self.power[t]/self.emax
        ## but then command is the same... SOC -= delta/self.emax


class Nodes:
    """ Contains a collection of Node() objects and tools for managing them.

        When called as Europe=Nodes() (or just, N=Nodes()) it builds nodes 
        from the country data files. If called as N=Nodes(filename) the it looks 
        for a pre-solved filename with balancing, curtailment, etc.

        For the time being, Flows are stored separately in npy files, which are 
        simply opened with np.load. Nodes can be saved with the inbuilt function save_nodes()"""


    def __init__(self, admat='./settings/admat.txt', path='./data/', files=['ISET_country_AT.npz', 'ISET_country_FI.npz', 'ISET_country_NL.npz', 'ISET_country_BA.npz', 'ISET_country_FR.npz', 'ISET_country_NO.npz', 'ISET_country_BE.npz', 'ISET_country_GB.npz', 'ISET_country_PL.npz', 'ISET_country_BG.npz', 'ISET_country_GR.npz', 'ISET_country_PT.npz', 'ISET_country_CH.npz', 'ISET_country_HR.npz', 'ISET_country_RO.npz', 'ISET_country_CZ.npz', 'ISET_country_HU.npz', 'ISET_country_RS.npz', 'ISET_country_DE.npz', 'ISET_country_IE.npz', 'ISET_country_SE.npz', 'ISET_country_DK.npz', 'ISET_country_IT.npz', 'ISET_country_SI.npz', 'ISET_country_ES.npz', 'ISET_country_LU.npz', 'ISET_country_SK.npz'], load_filename=None, fivenodes=0, full_load=False):
        self.cache=[]
        self.pathadmat=admat
        for i in range(len(files)):
            n=node(path,files[i],i)
            self.cache=append(self.cache,n)
        F=np.zeros((size(files),self.cache[0].nhours))
        
        for i in self.cache:
            i.set_alpha(alphas[i.id])

        if load_filename != None:
            self._load_nodes_(load_filename, full_load, path='./results/')

    def __getitem__(self,x):
        return self.cache[x]
        
    def __len__(self):
        return len(self.cache)

    def set_gammas(self,value):
        # to change a single node's gamma, just write XX.set_gamma(yy)
        # 'value' can be a single number or a vector
        if np.size(value)==1:
            for i in self.cache: i.set_gamma(value)
        elif np.size(value)!=np.size(self.cache):
            print "Wrong gamma vector size. ", np.size(value,0)," were received, ",np.size(self.cache)," were expected."
        else:
            for i in self.cache:
                i.set_gamma(value[i.id])

    def set_alphas(self,value):
        # to change a single node's alpha, just write XX.set_gamma(yy)
        # 'value' can be a single number or a vector
        if size(value)==1:
            for i in self.cache: i.set_alpha(value)
        elif size(value)!=size(self.cache):
            print "Wrong gamma vector size. ", size(value,0)," were received, ",size(self.cache)," were expected."
        else:
            for i in self.cache:
                i.set_alpha(value[i.id])

    def get_Ppos(self,t):
        ppos=np.zeros(size(self.cache),dtype='d')
        for i in self.cache:
            ppos[i.id]= float(max(min(i.ESS.pmax,floor(i.ESS.SOC[t]*i.ESS.emax)),0.0)  )
        return double(ppos)

    def get_Pneg(self,t):
        pneg=np.zeros(size(self.cache),dtype='d')
        for i in self.cache:
            pneg[i.id]=float(-max(min(-i.ESS.pmin,floor((1-i.ESS.SOC[t])*i.ESS.emax)),0.0) )
        return double(pneg)

    def set_ESS(self,powerfrac,hours):
        for i in self.cache:
            p=i.mean
            i.ESS.pmax=p*powerfrac
            i.ESS.pmin=-p*powerfrac
            i.ESS.emax=p*hours

    def up_ESS(self,t):
        for i in self.cache:
            i.ESS.update(t)

    def save_nodes(self,filename,path='./results/'):
        """Saves the contents of a Nodes instance to a npz file."""
        
        attribute = dir(self[0])
        save_str = []
        #Determine which attributes to be saved
        for attribute in dir(self[0]):
            if attribute[0]=='_':
                continue
            elif is_numlike(getattr(self[0],attribute)) or is_string_like(getattr(self[0],attribute)):
                save_str.append(attribute + '=' + 'array([self[i].'+attribute+' for i in arange(len(self))])')

        #Write save file
        eval('np.savez(path+filename,'+','.join(save_str)+')')

        print 'Saved nodes to file: ', path+filename
        sys.stdout.flush()
        
    def _load_nodes_(self,load_filename, full_load, path='./results/'):
        """Loads a Nodes instance from an npz file."""

        npzobj = np.load(path+load_filename)
        if full_load == False:
            for attribute in npzobj.files: ## this is real
                if (attribute == 'balancing') or (attribute =='curtailment') or (attribute == 'ESS'): ## this is my bullshit
                    for i in arange(len(self)):
                        setattr(self.cache[i],attribute,npzobj[attribute][i])
        if full_load == True:
            for attribute in npzobj.files:            
                for i in arange(len(self)):
                    setattr(self.cache[i],attribute,npzobj[attribute][i])
        
        for n in self.cache:
            n._update_()
        npzobj.close()
        print 'Loaded nodes from file: ', path+load_filename
        sys.stdout.flush()


