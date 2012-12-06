#! /usr/bin/env python
from pylab import *
from scipy import *
import numpy as np
from time import time
import sys, os
import gurobipy as gb
from auresc import *

def AtoKh(N):
    """ For a bettter understandig, go open admat.txt file in the settings folder. This function 
        transforms a transmission table A into the incidence matrix K of size N x L with K[n,l] 
        if link 'l' starts at node and -1 if it ends there. Also returns the h vector, which 
        holds the actual transmission limits.

        In this version, K is returned as row and column indices and values, and used to build 
        the problem in CPLEX. The final entry returns a list of links with names, which is very 
        useful when you get lost in the numbers."""

    Ad=np.genfromtxt(N.pathadmat,dtype='d')
    L=0
    listFlows=[]
    for j in range(len(Ad)):
        for i in range(len(Ad)):
            if i>j:
                if Ad[i,j] > 0:
                    L+=1
    K_values=[]
    K_column_indices=[]
    K_row_indices=[]
    h=np.zeros(L*2)
#    h=np.append(h,np.zeros(3*len(Ad)))
    L=0
    for j in range(len(Ad)):
        for i in range(len(Ad)):
            if i>j:
                if Ad[i,j] > 0:
                    K_values.extend([1,-1])
                    K_column_indices.extend([L,L])
                    K_row_indices.extend([j,i])
                    h[2*L]=Ad[i,j]
                    h[2*L+1]=Ad[j,i]
                    listFlows.append([str(N[j].label)+" to " +str(N[i].label), L])
                    L+=1
    #K=np.spmatrix(K_values,K_row_indices,K_column_indices)
    #K=coo_matrix((K_values,(K_row_indices,K_column_indices)))    
    return K_row_indices,K_column_indices,K_values,h, listFlows

def get_quant(quant=0.99,filename='results/copper_flows.npy'):
    f=np.load(filename)
    flows=[]
    for i in f:
        flows.append(i*(i>0.))
        flows.append(-i*(i<0.))
    a=np.zeros(len(flows))
    b=np.zeros(len(flows))
    hs=np.zeros(len(flows))
    for i in range(len(flows)):
        a=hist(flows[i],cumulative=True,bins=100,normed=True)
        for j in range(len(a[0])):
            if (a[0][j]>=quant):
                hs[i]=a[1][j]
                break
            if j == len(a[0]) - 1:
                hs[i.id]=max(a[1][:])
    return hs

###############################################################
###############################################################

def sdcpf(N,path='./settings/',copper=0,lapse=None,b=1.0,h0=None):


###############################################################
################ Loading Default values #######################
###############################################################

    if lapse == None:
        lapse=N[0].nhours
    K_row,K_col,K_val,H,LF=AtoKh(N) # dummy node has been deleted from admat.txt!!!

    Nlinks=len(K_row)/2
    Nnodes=len(N)
    B_0=1e6


    if (h0 != None):
        H=h0
    h_neg=b*-H[1:Nlinks*2:2]
    h_pos=b*H[0:Nlinks*2:2]

    if (copper == 1):
        h_neg=-1.e6*np.ones(Nlinks)
        h_pos=1.e6*np.ones(Nlinks)



###############################################################
###############################################################


###############################################################
################ Preparing vars for solver    #################
###############################################################

    flw=np.zeros(Nlinks)
    sto=np.zeros(Nnodes)
    S=np.zeros((Nnodes,lapse))
    F=np.zeros((Nlinks,lapse))

###############################################################
###############################################################


###############################################################
################# Setting up the model ########################
###############################################################

    network = gb.Model()    

    efes=['f'+str(i+1) for i in range(Nlinks)]
    eses=['s'+str(i+1) for i in range(Nnodes)]
    exes=['x'+str(i+1) for i in range(Nnodes)]
    Names = efes + eses + exes

    #upper and lower bounds
    h_lower=h_neg
    h_upper=h_pos
    s_upper=np.ones(Nnodes)*1e6
    s_lower=np.ones(Nnodes)
    x_upper=np.ones(Nnodes)*1e6
    x_lower=np.zeros(Nnodes)
    upper_b=np.concatenate((h_upper,s_upper,x_upper))
    lower_b=np.concatenate((h_lower,s_lower,x_lower))
    for i in range(len(Names)):
        network.addVar(lb=lower_b[i],ub=upper_b[i],name=Names[i])

    network.update()
    #adding default linear constraints
    cNames=[]
    for i in range(2*Nnodes): cNames.append(str(i))

    rol=[1.0,-1.0]
    for r in rol:
        for n in range(Nnodes):
            ind=[]
            val=[]
            for i in range(len(K_row)):
                if K_row[i]==n:
                    ind.append('f'+str(K_col[i]+1))
                    val.append(-1*r*K_val[i])
            ind.append('s'+str(n+1))
            ind.append('x'+str(n+1))
            val.append(r)
            val.append(-1.0)
            var=[]
            for i in ind:
                var.append(network.getVarByName(i))
            if r==1:
                network.addConstr(lhs=gb.LinExpr(val,var),sense='<',rhs=1e6,name=cNames[n])
            if r==-1:
                network.addConstr(lhs=gb.LinExpr(val,var),sense='<',rhs=1e6,name=cNames[n+27])
      
    network.setParam("OutputFlag",0)
    network.setParam("FeasibilityTol",1e-2)          
    network.update()

###############################################################
################# Run the time series  ########################
###############################################################


    start=time()
    for t in range(lapse):
        solution, B_opt=solve_flows(network,N,t,h_pos,h_neg,K_row,K_col,K_val,B_0)
        flw=solution[0:Nlinks]
        sto=solution[Nlinks:Nlinks+Nnodes]

        for i in range(Nlinks):
            F[i,t]=flw[i]
        for i in range(Nnodes):
            S[i,t]=sto[i]
        for i in N:
            i.ESS.power[t]=S[i.id,t]
        N.up_ESS(t)
        B_0=B_opt
        #print time()-start
        if (np.mod(t,2073)==0) and t>0:
                print '.',
                sys.stdout.flush()
    end=time()
    print "\nCalculation took %3.1f seconds." % (end-start)
    #print end-start
    sys.stdout.flush()




###############################################################
############### Assignment to Objects #########################
###############################################################

    start2=time()
    D=zeros(Nnodes)
    M=zeros((Nnodes,Nlinks))
    for i in range(len(K_row)):
        M[K_row[i]][K_col[i]]=K_val[i]
    for t in range(lapse):
        for i in N:
            D[i.id]=i.mismatch[t]
        tmp = dot(M,F[:,t])
        mis= D - tmp + [x for x in S[:,t]]
        balancing=[-mis[i] if mis[i]<0 else 0. for i in range(len(mis))]
        curtailment=[mis[i] if mis[i]>0 else 0. for i in range(len(mis))]
        for i in N:
            i.balancing[t] = balancing[i.id]
            i.curtailment[t] = curtailment[i.id]
            i.ESS.power[t]=S[i.id,t]
    end=time()
    print "Assigning balancing and curtailment took %3.1f seconds." % (end-start2)
    print "Complete calculation took %3.1f seconds." % (end-start)
    return N,F

###############################################################
###############################################################





def solve_flows(network,N,t,h_pos,h_neg,K_row,K_col,K_val,BC_0):
    relaxation = 0.1
    start=time()
    Nlinks=len(K_row)/2
    Nnodes=len(N)
    efes=['f'+str(i+1) for i in range(Nlinks)]
    eses=['s'+str(i+1) for i in range(Nnodes)]
    exes=['x'+str(i+1) for i in range(Nnodes)]
    Names=efes + eses + exes

    Delta=[i.mismatch[t] for i in N]
    s_upper=np.array(N.get_Ppos(t))
    s_lower=np.array(N.get_Pneg(t))
    rhs=[-i for i in Delta] + Delta
    #print len(rhs)

    for s in range(len(eses)):
        network.getVarByName(eses[s]).ub=float(s_upper[s])
        network.getVarByName(eses[s]).lb=float(s_lower[s])
    network.update()
    for r in range(len(rhs)):
        network.getConstrByName(str(r)).setAttr("rhs",float(rhs[r]))
    equis=[]
    for x in exes:
        equis.append(network.getVarByName(x))
    unos=[]
    for i in range(Nnodes):
        unos.append(1)
    network.setObjective(expr=gb.LinExpr(unos,equis),sense=1)

    network.update()
    network.optimize()
    BC_opt=network.objVal + 1.0
    #print network.objVal
    
    # add new constraint
    nunos=[]
    for o in range(Nnodes):
        nunos.append(1)
    network.addConstr(lhs=gb.LinExpr(nunos,equis),sense='<',rhs=BC_opt,name="newcon")

    # set new objective
    unos=[]
    for i in range(Nlinks):
        unos.append(1)
    Fs=[]
    for f in efes:
        Fs.append(network.getVarByName(f))
    Qobj=gb.QuadExpr()
    Qobj.addTerms(unos,Fs,Fs)
    network.setObjective(expr=Qobj,sense=1)
    #solve
    network.update()
    network.optimize()

    # get solution
    v=[]
    for i in network.getVars():
        v.append(i.x)
    
    # remove contraint
    network.remove(network.getConstrs()[-1])


    return v, BC_opt






