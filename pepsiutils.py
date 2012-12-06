    #! /usr/bin/env python
from aures import *

cdict = {'red':  ((0.0, .281250,.281250),
                  (0.5, 138/256., 138/256.),
                  (1.0, 0.91015, 0.91015)),
          'green': ((0.0, 0.03906, 0.03906),
                    (0.5, 155/256., 155/256.),
                    (1.0, 0.496093, 0.496093)),
          'blue': ((0.0, 0.23828, 0.23828),
                   (0.5, 15/256., 15/256.),
                   (1.0, 0.0078125, 0.0078125))}
my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

def set_right_alphas(N):
    alphas=[0.67431640625 , 0.7955078125 , 0.71640625 , 0.68291015625 , 0.755322265625 , 0.8490234375 , 0.70126953125 , 0.7873046875 , 0.739599609375 , 0.659423828125 , 0.64169921875 , 0.660888671875 , 0.690698242188 , 0.639697265625 , 0.688647460938 , 0.713232421875 , 0.6625 , 0.69501953125 , 0.71611328125 , 0.7541015625 , 0.817724609375 , 0.731689453125 , 0.64677734375 , 0.6509765625 , 0.696557617188 , 0.70712890625 , 0.7068359375]
    for n in N:
        n.set_alpha(alphas[n.id])


def biggestpair(H):
    H0=np.zeros((len(H))/2)
    for i in range(len(H0)):
        H0[i]=max(H[2*i],H[2*i+1])
    return H0

def get_optimal_mix_balancing(L, GW, GS, gamma=1., returnall=True, normalized=False):
    L, GW, GS = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2) 
    weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))

    l = weighed_sum(L)
    Gw = weighed_sum(GW)	
    Gs = weighed_sum(GS)

    mismatch = lambda alpha_w: gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - l
    res_load = lambda alpha_w: sum(get_positive(-mismatch(alpha_w)))

    alpha_w_opt = fmin(res_load,0.5)
    res_load_1p = lambda alpha_w: res_load(alpha_w)-(res_load(alpha_w_opt)+.01*sum(l))

    alpha_w_opt_1p_interval = array([brentq(res_load_1p, 0, alpha_w_opt),brentq(res_load_1p, alpha_w_opt, 1)])

    if normalized:
        mismatch_opt = mismatch(alpha_w_opt)
    else:
        mismatch_opt = mismatch(alpha_w_opt)*mean(sum(L,axis=0))
    res_load_opt = sum(get_positive(-mismatch_opt))

    if returnall:
    #Returns: alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt
        return alpha_w_opt, alpha_w_opt_1p_interval, res_load_opt, mismatch_opt
    else:
        return alpha_w_opt

def histmany(N,Ns):
    for n in N:
        a=-n.load
        b=n.mismatch
        d=n.curtailment-n.balancing
        c=[]
        for x in d:
            if x>=1 or x<=-1:
                c.append(x)
        f=Ns[n.id].curtailment-Ns[n.id].balancing
        e=[]
        for x in f:
            if x>=1 or x<=-1:
                e.append(x)
        bins=arange(min(a),max(b),(max(b)-min(a))/250)
        p2=hist(b,histtype='step',bins=bins,align='mid',label="Mismatch - after RE")
        p1=hist(a,histtype='step',bins=bins,align='mid',label="Load")
        p3=hist(c,histtype='step',bins=bins,align='mid',label="Residual Load/Spillage - after KF")
        p4=hist(e,histtype='step',bins=bins,align='mid',label="Residual Load/Spillage - after KF & S")
        ax=subplot(1,1,1)    
        plt.ylabel('Incidence (h)')
        plt.xlabel('Power (MW)')
        plt.title(r'Histogram of Power Flows : '+str(n.label))
        gcf().set_size_inches([9*1.5,3*1.75])
        plt.ylim(0,max(p1[0])*1.075)
        handles,labels=ax.get_legend_handles_labels()
        ax.legend(handles,labels)
        plt.grid(True)
        savefig('./figures/6h_' + str(n.label) +'.pdf', dpi=300)
        plt.close()


def histone(N):
    alphas=np.arange(0,1.1,0.1)
    ###### Setting colours
    jet=plt.get_cmap("dark2")
    cnorm=matplotlib.colors.Normalize(vmin=0,vmax=1)
    colormap=matplotlib.cm.ScalarMappable(norm=cnorm,cmap=my_cmap)
    ######
    N.set_alphas(0)
    for n in N:
        b=arange(1.25*min(n.mismatch),1.25*max(n.mismatch),(1.25*max(n.mismatch)-1.25*min(n.mismatch))/250)
        for a in alphas:
            n.set_alpha(a)
            hist(n.mismatch,histtype='step',bins=b,align='mid',color=colormap.to_rgba(a),label="Mismatch, Alpha = "+str(a))
        ax=subplot(1,1,1)    
        plt.ylabel('Incidence (h)')
        plt.xlabel('Power (MW)')
        plt.title(r'Multialpha Histogram : '+str(n.label))
        gcf().set_size_inches([9*1.5,3*1.75])
        #plt.ylim(0,max(b[0])*1.075)
        handles,labels=ax.get_legend_handles_labels()
        ax.legend(handles,labels)
        plt.grid(True)
        savefig('./figures/Multialpha_' + str(n.label) +'.pdf', dpi=300)
        plt.close()

def Case_A(betas=[0.0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,1.0,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.0,2.25,2.50,2.75,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.5,15.0,17.5,20.0,25.0,30.0]):
    N=Nodes()
    a,b,c,hI,e=AtoKh(N)
    hQ=get_quant(0.99)
    h=np.zeros(88)    
    for b in betas:
        for i in range(88):
            h[i]=min(b*hI[i],hQ[i])
        N,F=sdcpf(N,h0=h)
        N.save_nodes('Case_A_Beta_'+str(b))
        save('./results/'+'Flows_Case_A_Beta_'+str(b),F)


def Case_B(h0=None,betas=[0.0,0.05,0.1,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,1.0,1.05,1.1,1.2,1.3,1.4,1.5,2.0,2.5,3.0,4.0,5.0]):
    if h0 == None: h0=get_quant(.99)
    N=Nodes()
    for b in betas:
        N,F=sdcpf(N,b=b,h0=h0)
        N.save_nodes('Case_B_Beta_'+str(b))
        save('./results/'+'Flows_Case_B_Beta_'+str(b),F)


#scenario C
def Case_C(quants=[0.75,0.8,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.9920,.9940,.9960,.9980,.9999]):
    N=Nodes()
    for q in quants:
        h=get_quant(q)
#        if sum(h) >=1:
#            for i in range(len(h)):
#                if h[i]==0: h[i]=50
#            print h
        N,F=sdcpf(N,h0=h)
        N.save_nodes('Case_C_Quant_'+str(q))
        save('./results/'+'Flows_Case_C_Quant_'+str(q),F)

#betas=[0.0,0.25,0.5,0.75,1.0,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.0,2.25,2.50,2.75,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.5]
    
def Plot_A():
    betas=[0.0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,1.0,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.0,2.25,2.50,2.75,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.5]
    PlotA=np.zeros((len(betas),2))
    N=Nodes()
    h=np.zeros(44)

    K1,k2,k3,Hi,lF=AtoKh(N)
    HI=biggestpair(Hi)

    Hop=get_quant(.99)
    Hopt=biggestpair(Hop)

    j=0
    for b in betas:
        for i in range(len(Hopt)):
            h[i]=min(b*HI[i],Hopt[i])
        PlotA[j,0]=b*sum(h)/sum(HI)
        N=Nodes(load_filename='Case_A_Beta_'+str(b)+'.npz')
        a=0
        d=0
        for i in N:
            a+=np.sum(i.balancing)
            d+=i.mean*i.nhours
        c=a/d
        PlotA[j,1]=c
        j+=1
        del(N)
    save('./results/PlotA',PlotA)
    return PlotA


def Plot_B():
    betas=[0.0,0.05,0.1,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,1.0,1.05,1.1,1.2,1.3,1.4,1.5,2.0,2.5,3.0,4.0,5.0]
    N=Nodes()
    K1,k2,k3,Hac,lF=AtoKh(N)
    Hact=biggestpair(Hac)
    Hop=get_quant(.99)
    Hopt=biggestpair(Hop)
    PlotB=np.zeros((len(betas),2))
    j=0
    for b in betas:
        PlotB[j,0]=b*sum(Hopt)/sum(Hact)
        N=Nodes(load_filename='Case_B_Beta_'+str(b)+'.npz')
        a=0
        d=0
        for i in N:
            a+=sum(i.balancing)
            d+=i.mean*i.nhours
        c=a/d
        PlotB[j,1]=c
        j+=1
        del(N)
    save('./results/PlotB',PlotB)
    return PlotB

def Plot_C():
    quants=[0.5,0.75,0.8,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.9920,.9940,.9960,.9980,.9999]
    N=Nodes()
    k1,k2,k3,Hac,lF=AtoKh(N)
    Hact=biggestpair(Hac)
    PlotC=np.zeros((len(quants),2))
    j=0
    for q in quants:
        Hop=get_quant(q)
        Hopt=biggestpair(Hop)
        PlotC[j,0]=sum(Hopt)/sum(Hact)
        N=Nodes(load_filename='Case_C_Quant_'+str(q)+'.npz')
        a=0
        d=0
        for i in N:
            a+=sum(i.balancing)
            d+=i.mean*i.nhours
        c=a/d
        PlotC[j,1]=c
        j+=1
        del(N)
    save('./results/PlotC',PlotC)
    return PlotC

#N=Nodes()
#N,F=sdcpf(N,copper=1)
#N.save_nodes("copper")
#save("./results/copper_flows",F)

#N,F=sdcpf(N)
#N.save_nodes("TodayCap")
#save("./results/TodayCap_flows",F)

#N,F=sdcpf(N,h0=get_quant(0.99))
#N.save_nodes("99PCap")
#save("./results/99PCap_flows",F)

#Case_A()
#Case_B()
#Case_C(quants=[0.5])

Case_A(betas=[0.0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85])
Plot_A()

