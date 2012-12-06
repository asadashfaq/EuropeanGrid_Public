#! /usr/bin/env python
from auresc import *
import networkx as nx
import csv
import sys
from pylab import *
from scipy.optimize import brentq
from scipy.optimize import fmin
#import datetime as dt
from localutils import get_quant, ISO2LONG, biggestpair

colwidth = (3.425)
dcolwidth = (2*3.425+0.236) 
cdict = {'red':  ((0.0, 1.0, 1.0),
                  (0.5, 0.0, 0.50),
                  (1.0, 0.0, 0.0)),
          'green': ((0.0, 0.0, 0.0),
                    (0.5, 1.0, 0.0),
                    (1.0, 0.0, 0.0)),
          'blue': ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 0.0),
                   (1.0, 1.0, 1.0))}
my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
colors_countries = ['#00A0B0','#6A4A3C','#CC333F','#EB6841','#EDC951'] #Ocean Five from COLOURlovers.
#sugar_colours = ['#F8CA00','#BD1550','#E97F02','#490A3D','#8A9B0F','#134B7C']
sugar_colours = ['#F8CA00','#AD0550','#FF9F22','#490A3D','#9AAB0F','#003B6C']
order = [24,11,10,22,3,9,17,13,4,23,12,0,14,16,26,25,15,6,18,2,8,7,19,21,20,5,1]
rolor= ( (16+3)/255. , (4*16+11)/255. , (7*16+12)/255.)
def get_positive(x):
    """ This function returns the positive side of variable. """
    return x*(x>0.)  #Possibly it has to be x>1e-10.

def get_optimal_mix_balancing(L, GW, GS, gamma=1., returnall=True, normalized=False):
    L, GW, GS = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2) 
    weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))

    l = weighed_sum(L)
    Gw = weighed_sum(GW)	
    Gs = weighed_sum(GS)

    mismatch = lambda alpha_w: gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - l
    res_load = lambda alpha_w: sum(get_positive(-mismatch(alpha_w)))

    alpha_w_opt = fmin(res_load,0.5,disp=False)
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

def get_q(tipo,s=None, N=None,quant=0.99):

    if s<>None:
        hs=0
        a=hist(s,cumulative=True,bins=100,normed=True)
        for j in range(len(a[0])):
            if (a[0][j]>=quant):
                hs=a[1][j]
                break
            if j == len(a[0]) - 1:
                hs=max(a[1][:])
        plt.close()
        return hs

    if N==None:
        N=Nodes()
    hs=np.zeros(len(N))
    for i in N:
        if tipo=="balancing":
            a=hist(i.mismatch,cumulative=True,bins=100,normed=True)
        if tipo=="load":
            a=hist(i.load,cumulative=True,bins=100,normed=True)
        for j in range(len(a[0])):
            if (a[0][j]>=quant):
                hs[i.id]=a[1][j]
                break
            if j == len(a[0]) - 1:
                hs[i.id]=max(a[1][:])
    plt.close()
    return hs

def plotbars1(): # Figures 2a and 2b
    alphas=[]
    resopts=[]
    realres=[]
    alphas_inter=[]
    N=Nodes()
    names=['EU','']
    solar=[1.0,0.0]

    for o in order:
        names.append(str(N[o].label))
        solar.append(1.0)

    #gcf().set_size_inches([15,5])
    width=0.95
    EU_load = sum(n.load for n in N)
    EU_mean = sum(n.mean for n in N)
    EU_load/= EU_mean
    EU_wind = sum(n.normwind for n in N)/27.0
    EU_solar= sum(n.normsolar for n in N)/27.0
    alpha,alpha_inter,res_load,wtf2=get_optimal_mix_balancing(EU_load,EU_wind,EU_solar)
    alphas.append(alpha)
    alphas.append(0.000001)
    alphas_inter.append(alpha_inter)
    alphas_inter.append(array([0.0000001,0.00001]))
    resopts.append(res_load/70128)
    resopts.append(0.0)
    for o in order:
        alpha,alpha_inter,res_load,wtf2=get_optimal_mix_balancing(N[o].load/N[o].mean,N[o].normwind,N[o].normsolar)
        alphas.append(alpha)
        alphas_inter.append(alpha_inter)
        resopts.append(res_load/70128)
        realres.append(N[o].mean*res_load)
    ind=np.arange(29)

    ax = subplot(1,1,1)
    gcf().set_size_inches([ 1.3*dcolwidth , 1.3*dcolwidth*0.4])#[9*1.5,3*1.75])
    gcf().set_dpi(400)

    #rects0=bar(ind*1.25,solar,width=width,align='center',color=(.9,.9,.0))
    rects1=bar(ind*1.25,alphas,width=width,align='center',color=rolor)
    axis(ymin=0,ymax=1.0,xmin=amin(ind)-.875,xmax=amax(ind*1.25)+.875)
    p1=ax.plot([-2,36],[mean(alphas[2:]),mean(alphas[2:])],linestyle='dashed',linewidth=1.5,color='k')
    errorbar(0,alphas[0],yerr=[ [float(alphas[0]-alphas_inter[0][0]) , float(alphas_inter[0][1]-alphas[0]) ]],color='k',linewidth=2.0)
    for i in range(2,29):
        errorbar(i*1.25,alphas[i],yerr=float(alphas[i]-alphas_inter[i][0]),color='k',linewidth=1.5)
        #errorbar(i*1.25,alphas[i],yerr=[[ float(alphas[i]-alphas_inter[i][0]) , float(alphas_inter[i][1]-alphas[i]) ]],color='k')

    yticks([0.0,0.25,0.5,0.75,1.0])
    xticks(ind*1.25+.35,names,rotation=60,ha='right',va='top')
    ylabel(r'Wind fraction ($\alpha_W$)')
    plt.text(-0.7,0.925, "(a)" , size='x-large')
    plt.tight_layout()
    savefig('./figures/optbar.pdf', dpi=300)
    plt.close()

    ax = subplot(1,1,1)
    gcf().set_size_inches([ 1.3*dcolwidth , 1.3*dcolwidth*0.40])#[9*1.5,3*1.75])
    gcf().set_dpi(400)
    rects2=bar(ind*1.25,resopts,width=width,align='center',color=rolor)
    m=sum(realres)/(EU_mean*70128)
    p0=ax.plot([-0.875,29*1.25+0.875],[m,m],linestyle='dashed',linewidth=2.0,color='k')
    axis(xmin=amin(ind)-.875,xmax=amax(ind*1.25)+.875)
    yticks([0.00,0.10,0.20,0.30])
    ax.set_yticklabels(["0.00","0.10",'0.20','0.30'])
    xticks(ind*1.25+.35,names,rotation=60,ha='right',va='top')
    ylabel(r'Residual load [normalised]')
    plt.text(-0.7,0.35*0.925, "(b)" , size='x-large')
    plt.tight_layout()
    savefig('./figures/minres.pdf', dpi=300)
    plt.close()


def plotbars2(): #figure 2c
    N=Nodes()
    names=['EU','']
    bkg=[]
    one=[]
    ten=[]
    nin=[]
    ninin=[]

    for o in order:
        names.append(str(N[o].label))

    #gcf().set_size_inches([15,5])
    width=0.95
    EU_load = sum(n.load for n in N)
    EU_mean = sum(EU_load)/70128
    EU_wind = sum(n.get_wind() for n in N)
    EU_solar= sum(n.get_solar() for n in N)
    EU_mismatch=(EU_wind+EU_solar)-EU_load
    bkg.append(get_q("none",s=EU_load,quant=0.99)/EU_mean)
    bkg.append(0.0)
    one.append(-get_q("none",s=EU_mismatch,quant=0.01)/EU_mean)
    one.append(0.0)
    ten.append(-get_q("none",s=EU_mismatch,quant=0.10)/EU_mean)
    ten.append(0.0)
    nin.append(get_q("none",s=EU_mismatch,quant=0.90)/EU_mean)
    nin.append(0.0)
    ninin.append(get_q("none",s=EU_mismatch,quant=0.99)/EU_mean)
    ninin.append(0.0)

    for o in order:    
        bkg.append(get_q("none",s=N[o].load,quant=0.99)/N[o].mean)
        one.append(-get_q("none",s=N[o].mismatch,quant=0.01)/N[o].mean)
        ten.append(-get_q("none",s=N[o].mismatch,quant=0.10)/N[o].mean)
        nin.append(get_q("none",s=N[o].mismatch,quant=0.90)/N[o].mean)
        ninin.append(get_q("none",s=N[o].mismatch,quant=0.99)/N[o].mean)

    ax = subplot(1,1,1)
    gcf().set_size_inches([ 1.3*dcolwidth , 1.3*dcolwidth*0.4])#[9*1.5,3*1.75])
    gcf().set_dpi(400)

    ind=np.arange(29)

    rects0=bar(ind*1.25,bkg,width=1.0,align='center',alpha=1.0,color=(.80,.80,.80))
    rects1=bar(ind*1.25-0.25,one,width=0.5,align='center',color="r")
    rects2=bar(ind*1.25-0.25,ten,width=0.5,align='center',color="r", hatch="\\\\")
    rects3=bar(ind*1.25+0.25,ninin,width=0.5,align='center',color=sugar_colours[4])
    rects4=bar(ind*1.25+0.25,nin,width=0.5,align='center',color=sugar_colours[4], hatch="//")

    pp = (rects0, rects1, rects2, rects3, rects4)#,rectsQ)
    pp_txtlabels = ("99% Q Load","99% Q Deficit","90% Q Deficit","99% Q Excess", "90% Q Excess")
    leg = legend(pp,pp_txtlabels,loc='upper right',ncol=len(pp),fancybox=True)
    ltext  = leg.get_texts();
    setp(ltext, fontsize=7)    # the legend text fontsize

    axis(xmin=amin(ind)-.875,xmax=amax(ind*1.25)+.875)
    yticks([0.00,0.5,1.0,1.5,2.0,2.5,3.0])
    ax.set_yticklabels(['0.00','0.50','1.00','1.50','2.00','2.50','3.00'])

    #yticks([0.0,0.25,0.5,0.75,1.0])
    xticks(ind*1.25+.35,names,rotation=60,ha='right',va='top')
    ylabel(r'Mismatch quantiles [normalised]')
    plt.text(-0.7,3*0.925, "(c)" , size='x-large')
    plt.tight_layout()
    savefig('./figures/manybar.pdf', dpi=300)
    plt.close()

def smoothhist(N): ## The mismatch histograms with different alphas
    for n in N:
        lines=[]
        lo=-n.load/n.mean
        x0=-2
        x1=3.001
        alphas=[0.0,0.7,1.0]

    ######
        b=np.arange(x0,x1,(x1-x0)/250.)
        lines.append(hist(lo,bins=b)[0]/(70128.))
        for a in alphas:
            n.set_alpha(a)
            lines.append(hist(n.mismatch/n.mean,bins=b)[0]/(70128.))
        plt.close()
        ax=subplot(1,1,1)
        plot(b[0:-1],lines[0],label="Load", linewidth=2.0, color ='k')
        plot(b[0:-1],lines[1],label=r"VRES mix $\alpha_W$ = 0.0", linewidth=1.5, color=sugar_colours[2]) 
        plot(b[0:-1],lines[2],label=r"VRES mix $\alpha_W$ = 0.7", linewidth=1.5, color=sugar_colours[4])
        plot(b[0:-1],lines[3],label=r"VRES mix $\alpha_W$ = 1.0", linewidth=1.5, color=sugar_colours[5])
        plt.ylabel('P($\Delta$)')
        plt.xlabel(r'Mismatch power [normalised]')
        plt.text(-1.9,max(lines[0]*1.00)*0.95, "(a)" , size='x-large')
        #plt.title(r'M : '+str(n.label))
        gcf().set_size_inches([ dcolwidth , dcolwidth*0.40])#[9*1.5,3*1.75])
        gcf().set_dpi(400)
        plt.ylim(0,max(lines[0])*1.075)
        plt.yticks(np.arange(0,max(lines[0])*1.075,0.01))
        plt.xlim(-2,3.001)
        plt.xticks([-2,-1,0,1,2,3])
        handles,labels=ax.get_legend_handles_labels()
        leg=ax.legend(handles,labels,title=ISO2LONG(n.label)+" ("+str(n.label)+")")
        setp(leg.get_texts(),fontsize="small")
        plt.grid(which="major",axis='x')
        plt.tight_layout()
        #show()
        savefig('./figures/SMultialpha_' + str(n.label) +'.pdf', dpi=400)
        plt.close()

    #### Now the same for aggregated EU
    lines=[]

    EU_load = sum(n.load for n in N)
    EU_mean = mean(EU_load)
    EU_load/= EU_mean
    EU_delta = sum(n.mismatch for n in N)
    x0=-2
    x1=3
    b=arange(x0,x1,(x1-x0)/250.)
    lines.append(hist(-EU_load,bins=b)[0]/70128.)
    for a in alphas:
        N.set_alphas(a)
        EU_delta = sum(n.mismatch for n in N)
        lines.append(hist(EU_delta/EU_mean,bins=b)[0]/70128.)
    plt.close()
    ax=subplot(1,1,1)
    plot(b[0:-1],lines[0],label="Load", linewidth=2.0, color='k')
    plot(b[0:-1],lines[1],label=r"VRES mix $\alpha_W$ = 0.0", linewidth=1.5, color=sugar_colours[2]) 
    plot(b[0:-1],lines[2],label=r"VRES mix $\alpha_W$ = 0.7", linewidth=1.5, color=sugar_colours[4])
    plot(b[0:-1],lines[3],label=r"VRES mix $\alpha_W$ = 1.0", linewidth=1.5, color=sugar_colours[5])
    plt.ylabel('P($\Delta$)')
    plt.xlabel(r'Mismatch power [normalised]')
    plt.text(-1.9,max(lines[0]*1.00)*0.95, "(b)" , size='x-large')
    #plt.title(r'M : '+str(n.label))
    gcf().set_size_inches([ dcolwidth , dcolwidth*0.40])#[9*1.5,3*1.75])
    gcf().set_dpi(400)
    plt.ylim(0,max(lines[0])*1.075)
    plt.xlim(-2,3.001)
    plt.xticks([-2,-1,0,1,2,3])
    plt.yticks(np.arange(0,max(lines[0])*1.075,0.01))
    handles,labels=ax.get_legend_handles_labels()
    leg=ax.legend(handles,labels, title="Europe (EU)")
    setp(leg.get_texts(),fontsize="small")
    plt.grid(which="major",axis='x')
    plt.tight_layout()
    #show()
    savefig('./figures/SMultialpha_EU.pdf', dpi=300)
    plt.close()

def smoothhist2(M,N,O): #M=today, N=intermediate, O=99P
    for n in N:
        lines=[]
    ######
        lo=-n.load/n.mean
        #x0=1.750*min(n.mismatch/n.mean)
        #x1=1.750*max(n.mismatch/n.mean)
        x0=-2
        x1=3.001
        b=arange(x0,x1,(x1-x0)/250.)
        lines.append(hist(lo,bins=b)[0]/(70128.))
        lines.append(hist(n.mismatch/n.mean,bins=b)[0]/(70128.))

        u=M[n.id].curtailment-M[n.id].balancing
        v=[]
        for w in u:
            if w>=1 or w<=-1:
                v.append(w/n.mean)
        lines.append(hist(v,bins=b)[0]/(70128.))

        u=n.curtailment-n.balancing
        v=[]
        for w in u:
            if w>=1 or w<=-1:
                v.append(w/n.mean)
        lines.append(hist(v,bins=b)[0]/(70128.))

        u=O[n.id].curtailment-O[n.id].balancing
        v=[]
        for w in u:
            if w>=1 or w<=-1:
                v.append(w/n.mean)
        lines.append(hist(v,bins=b)[0]/(70128.))

        plt.close()
        ax=subplot(1,1,1)
        plot(b[0:-1],lines[0],label="Load", linewidth=1.5, color="r")
        plot(b[0:-1],lines[1],label="No transmission", linewidth=1.5, color=sugar_colours[2]) 
        plot(b[0:-1],lines[2],label="Present layout", linewidth=1.5, color=sugar_colours[4])
        plot(b[0:-1],lines[3],label="Intermediate layout", linewidth=1.5, color=sugar_colours[1])
        plot(b[0:-1],lines[4],label="99% Quantile layout", linewidth=1.5, color=sugar_colours[5])
        #plot(b[0:-1],lines[4],label="$\Delta$, alpha = 0.8", linewidth=1.5, color=sugar_colours[3])
        #plot(b[0:-1],lines[5],label="$\Delta$, alpha = 1.0", linewidth=1.5, color=sugar_colours[4])
        plt.ylabel('P($\Delta$)')
        plt.xlabel(r'Mismatch Power [normalised]')
        plt.tight_layout()
        if n.label=="ES":
            plt.text(-1.9,max(lines[0]*1.00)*0.95, "(b)" , size='x-large')
        if n.label=="DK":
            plt.text(-1.9,max(lines[0]*1.00)*0.95, "(a)" , size='x-large')
        #plt.text(-0.5,max(lines[0]*1.00), "Mismatch: "  + ISO2LONG(n.label), size='x-large')
        #plt.title(r'M : '+str(n.label))
        gcf().set_size_inches([ dcolwidth , dcolwidth*0.40])#[9*1.5,3*1.75])
        gcf().set_dpi(400)
        plt.ylim(0,max(lines[0])*1.075)
        plt.xlim(-2,3)
        plt.xticks([-2,-1,0,1,2,3])
        plt.yticks(np.arange(0,max(lines[0])*1.075,0.01))
        handles,labels=ax.get_legend_handles_labels()
        leg=ax.legend(handles,labels,title=ISO2LONG(n.label)+" ("+str(n.label)+")")
        setp(leg.get_texts(),fontsize="small")
        plt.grid(which="major",axis='x')
        plt.tight_layout()
        #show()
        savefig('./figures/MutlitLayout_' + str(n.label) +'.pdf', dpi=300)
        plt.close()


    #show()



def AtoKh_old(N,G=None,h0=None,pathadmat='./settings/admat.txt'):
    Ad=np.genfromtxt(pathadmat,dtype='d')
    L=0
    if G==None:
        G=nx.Graph()
    listFlows=[]
    for j in range(len(Ad)):
        for i in range(len(Ad)):
            if i>j:
                if Ad[i,j] > 0:
                    L+=1
    K=np.zeros((len(Ad),L))
    h=np.zeros(L*2)
    h=np.append(h,np.zeros(3*len(Ad)))
    L=0
    j=0
    i=0
    for j in range(len(Ad)):
        for i in range(len(Ad)):
            if i>j:
                if Ad[i,j] > 0:
                    K[j,L]=1
                    K[i,L]=-1
                    h[2*L]=Ad[i,j]
                    h[2*L+1]=Ad[j,i]
                    listFlows.append([ str(N[j].label) , str(N[i].label) , L ])
                    L+=1  
    if h0 == None:
        h0=h

#    for l in range(L):
#        w=max(h0[l*2],h0[l*2+1])
#        if l==41:
#            w=h0[l*2]
#        G.add_edge(listFlows[l][0], listFlows[l][1] , weight= w)

    return K,h0,listFlows


def draw_static_network(N=None,F=None,tit="1"): ## All the network figures
    close()
    if N==None:
        N=Nodes()
    G=nx.Graph()
    nodelist=[]

    for n in N:
        G.add_node(str(n.label))
        nodelist.append(str(n.label))

    K,h,ListF=AtoKh_old(N)
    
    if F<>None:
        h=F
    for l in ListF:
        w=max(h[l[2]*2],h[l[2]*2+1])
        if l[2]==41:
            w=h[l[2]*2]
        G.add_edge(l[0], l[1] , weight= w)

  
    pos=nx.spring_layout(G)


    pos['AT']=[0.55,0.45]
    pos['FI']=[1.0,1.0]
    pos['NL']=[0.40,0.85]
    pos['BA']=[0.65,0.15]
    pos['FR']=[0.15,0.60]
    pos['NO']=[0.5,1.05]
    pos['BE']=[0.275,0.775]
    pos['GB']=[0.15,0.85]
    pos['PL']=[0.75,0.8]
    pos['BG']=[0.9,0.0]
    pos['GR']=[0.7,0.0]
    pos['PT']=[0.0,0.15]
    pos['CH']=[0.4,0.45]
    pos['HR']=[0.75,0.3]
    pos['RO']=[1.0,0.15]
    pos['CZ']=[0.75,0.60]
    pos['HU']=[1.0,0.45]
    pos['RS']=[0.85,0.15]
    pos['DE']=[0.45,0.7]
    pos['IE']=[0.0,0.95]
    pos['SE']=[0.75,1.0]
    pos['DK']=[0.5,0.875]
    pos['IT']=[0.4,0.2]
    pos['SI']=[0.55,0.3]
    pos['ES']=[0.15,0.35]
    pos['LU']=[0.325,0.575]
    pos['SK']=[0.90,0.55]

    fig = figure(dpi=400,figsize=(colwidth,colwidth*0.65))

    ax2= fig.add_axes([.775,0.05,.95,.9]) #For displaying graph    
    ax2.plot([0.0,.08],[0.15,0.15],linewidth=0.75,color='k')
    ax2.plot([0.0,.08],[0.25,0.25],linewidth=1.0,color='k')
    ax2.plot([0.0,.08],[0.35,0.35],linewidth=1.25,color='k')
    ax2.plot([0.0,.08],[0.45,0.45],linewidth=1.5,color='k')
    ax2.plot([0.0,.08],[0.55,0.55],linewidth=1.75,color='k')
    ax2.plot([0.0,.08],[0.65,0.65],linewidth=2.0,color='k')
    ax2.plot([0.0,.08],[0.75,0.75],linewidth=2.5,color='k')
    ax2.plot([0.0,.08],[0.85,0.85],linewidth=3.0,color='k')
    ax2.plot([0.0,.08],[0.95,0.95],linewidth=3.5,color='k')
    ax2.plot([0.0,.08],[1.05,1.05],linewidth=4.0,color='k')
    ax2.text(0.09,0.175,"500 MW",fontsize=7)
    ax2.text(0.09,0.275,"1 GW",fontsize=7)
    ax2.text(0.09,0.375,"2 GW",fontsize=7)
    ax2.text(0.09,0.475,"3 GW",fontsize=7)
    ax2.text(0.09,0.575,"4.5 GW",fontsize=7)
    ax2.text(0.09,0.675,"7.5 GW",fontsize=7)
    ax2.text(0.09,0.775,"10 GW",fontsize=7)
    ax2.text(0.09,0.875,"15 GW",fontsize=7)
    ax2.text(0.09,0.975,"20 GW",fontsize=7)
    ax2.text(0.09,1.075,"> 20 GW",fontsize=7)
    ax2.axis([0.0,1.0,0.0,1.2])
   
    ax2.axis('off')


    ax1= fig.add_axes([-0.075,-0.05,.925,1.1]) #For displaying graph    
    nx.draw_networkx_nodes(G,pos,node_size=200,nodelist=nodelist,node_color="b",facecolor=(1,1,1))
    e0=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']<=500]
    e1=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>500 and d['weight']<=1000]
    e2=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1000 and d['weight']<=2000]
    e3=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>2000 and d['weight']<=3000]
    e4=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>3000 and d['weight']<=4500]
    e5=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>4500 and d['weight']<=7500]
    e6=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>7500 and d['weight']<=10000]
    e7=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>10000 and d['weight']<=15000]
    e8=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>15000 and d['weight']<=20000]
    e9=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>20000]
    
    nx.draw_networkx_edges(G,pos,edgelist=e0,width=.75,edgecolor=(0.4,0.4,0.4))
    nx.draw_networkx_edges(G,pos,edgelist=e1,width=1.0,edgecolor=(0.375,0.375,0.375))
    nx.draw_networkx_edges(G,pos,edgelist=e2,width=1.25,edgecolor=(0.35,0.35,0.35))
    nx.draw_networkx_edges(G,pos,edgelist=e3,width=1.5,edgecolor=(0.325,0.325,0.325))
    nx.draw_networkx_edges(G,pos,edgelist=e4,width=1.75,edgecolor=(0.3,0.3,0.3))
    nx.draw_networkx_edges(G,pos,edgelist=e5,width=2.0,edgecolor=(0.3,0.3,0.3))
    nx.draw_networkx_edges(G,pos,edgelist=e6,width=2.5,edgecolor=(0.3,0.3,0.3))
    nx.draw_networkx_edges(G,pos,edgelist=e7,width=3.0,edgecolor=(0.25,0.25,0.25))
    nx.draw_networkx_edges(G,pos,edgelist=e8,width=3.50,edgecolor=(0.2,0.2,0.2))
    nx.draw_networkx_edges(G,pos,edgelist=e9,width=4.0,edgecolor=(0.15,0.15,0.15))
    nx.draw_networkx_labels(G,pos,font_size=8,font_color='w',font_family='sans-serif')

    ax1.axis('off')
    #plt.tight_layout()
    savefig("./figures/network"+tit+".eps")
    #show() # display

def show_hist(link,tit,mean=None,filename='results/copper_flows.npy',b=250):
    plt.close()
    if mean==None:
        mean=1.0    
    flows=np.load(filename)   
    f=[]
#    ax=subplot(1,1,1)
    for i in flows[link]:
        if i>1 or i<-1:
            f.append(i/mean)
    zzone=[-get_quant(0.9999)[link*2+1],get_quant(0.9999)[link*2]]
    zone=[-get_quant(0.999)[link*2+1],get_quant(0.999)[link*2]]
    one=[-get_quant(0.99)[link*2+1],get_quant(0.99)[link*2]]
    five=[-get_quant(0.95)[link*2+1],get_quant(0.95)[link*2]]
    plt.close()
    ax=subplot(1,1,1)
    a=hist(f,bins=(max(f)-min(f))/(0.01*0.9/2.5),normed=0,histtype='stepfilled',color=rolor,weights=np.zeros((len(f)))+ 1./70128)
    plt.ylabel('P($F_l$)',size="large")
    plt.xlabel(r'Directed power flow [normalised]',size="large")
    plt.text(-1.46,.009,"France to Spain",size="large")
    #plt.title(r'Histogram of Power Flows : '+str(tit))
    vlines(zzone[0]/mean,0,0.0015,color='r',linewidth=1.2)
    plt.text(zzone[0]/mean*1.05,.0018,'0.01% Q',size=11)
    vlines(zone[0]/mean,0,0.003,color='r',linewidth=1.2)
    plt.text(zone[0]/mean*1.05,.0032,'0.1% Q',size=11)
    vlines(one[0]/mean,0,0.004,color='r',linewidth=1.2)
    plt.text(one[0]/mean*1.2,.0042,'1.0% Q',size=11)
    vlines(five[0]/mean,0,.0055,color='r',linewidth=1.2)
    plt.text(five[0]/mean*1.25,.0057,'5.0% Q',size=11)
    vlines(five[1]/mean,0,.0055,color='r',linewidth=1.2)
    plt.text(five[1]/mean*0.9,.0057,'95% Q',size=11)
    vlines(one[1]/mean,0,0.0040,color='r',linewidth=1.2)
    plt.text(one[1]/mean*0.8,0.0042,'99% Q',size=11)
    vlines(zone[1]/mean,0,0.003,color='r',linewidth=1.2)
    plt.text(zone[1]/mean*0.8,0.0032,'99.9% Q',size=11)
    vlines(zzone[1]/mean,0,0.0015,color='r',linewidth=1.2)
    plt.text(zzone[1]/mean*0.99,0.0017,'99.99% Q',size=11)
    plt.setp(ax.get_xticklabels(), fontsize=12) # rotation='vertical',
    plt.setp(ax.get_yticklabels(), fontsize=12)
    gcf().set_size_inches([1.3*dcolwidth,1.3*dcolwidth*0.5])
    plt.ylim(0,.01001)
    plt.xlim(-1.5,1.0)
    plt.grid(which="major",axis='x')
#    plt.axis([-5.35,5.65,0,0.5])#plt.axis([0.5, 1, 0, 38000])
    plt.tight_layout()
    savefig('./figures/' + str(link) +'.eps', dpi=400)
    plt.close()
    #show()
    filename='results/99PCap_flows.npy'
    plt.close()
    if mean==None:
        mean=1.0    
    flows=np.load(filename)   
    f=[]
#    ax=subplot(1,1,1)
    for i in flows[link]:
        if i>1 or i<-1:
            f.append(i/mean)
    print 
    plt.close()
    ax=subplot(1,1,1)
    a=hist(f,bins=(max(f)-min(f))/(0.01*0.9/2.5),normed=0,histtype='stepfilled',color=rolor,weights=np.zeros((len(f)))+ 1./70128)
    plt.ylabel('P($F_l$)',size="large")
    plt.xlabel(r'Transmission Magnitude ($\langle L_{DE} \rangle$)',size="large")
    plt.text(-0.48,2.35,"Germany to Denmark",size="large")
    #plt.title(r'Histogram of Power Flows : '+str(tit))
    plt.setp(ax.get_xticklabels(), fontsize=12) # rotation='vertical',
    plt.setp(ax.get_yticklabels(), fontsize=12)
    gcf().set_size_inches([1.3*dcolwidth,1.3*dcolwidth*0.5])
    plt.axis([-0.5,0.4,0.0,2.5])#plt.axis([0.5, 1, 0, 38000])
    plt.grid(True)
    plt.tight_layout()
    savefig('./figures/' + str(tit) +'_constr.eps', dpi=400)
    plt.close()
    #show()

def plot_allcases():
    plt.close()
    plota=np.load('./results/PlotA.npy')
    plotb=np.load('./results/PlotB.npy')
    plotc=np.load('./results/PlotC.npy')

    plota[:,0]*=74.83
    plotb[:,0]*=74.83
    plotc[:,0]*=74.83

#    plota[:,1]*=100
#    plotb[:,1]*=100
#    plotc[:,1]*=100


    ax=subplot(1,1,1)

    plt.xlabel('Total installed transmission capacity [GW]')
    plt.ylabel('Balancing energy [normalised]')


    m=np.min(plotc[:,1])
    p0=ax.plot([0,1500],[m,m],linestyle='dashed',color='k',linewidth=2.5)
    ax.vlines(74.83,0.0,.35,linestyle='dashed',color='k',linewidth=2.5)
    plotcc=np.append(plotc,[[900,m]],0)
    p1,=ax.plot(plota[:,0],plota[:,1],label='Interpolation A',linewidth=4.5, color=sugar_colours[1])
    p2,=ax.plot(plotb[:,0],plotb[:,1],label='Interpolation B',linewidth=4.25,color=sugar_colours[5], linestyle='dashed')
    p3,=ax.plot(plotcc[:,0],plotcc[:,1],label='Interpolation C',linewidth=4.5,color=sugar_colours[4])

    plt.axis([0, 900.1, .125, .27])
    plt.yticks([.15,.17,.19,.21,.23,.25,.27])
    plt.xticks([0,100,200,300,400,500,600,700,800,900])
    pp = (p1,p2,p3)
    pp_txtlabels = (r'Interpolation A',r'Interpolation B',r'Interpolation C')
    leg = legend(pp,pp_txtlabels,loc='upper right',ncol=1);
    ltext  = leg.get_texts();
    setp(ltext, fontsize=18)    # the legend text fontsize
    plt.text(74.83 + 10 ,.257,'Installed',size=16)
    plt.text(74.83 + 10 ,.250,'capacity',size=16)
    plt.text(74.83 + 10 ,.243,'2010-2011',size=16)
    plt.text(200 + 10 ,.144,'Minimum',size=16)
    plt.text(200 + 10 ,.137,'balancing',size=16)
    plt.text(200 + 10 ,.130,'energy',size=16)
    plt.tick_params(axis='both',which='both', labelsize=18)
    plt.xlabel('Total installed transmisison capacity [GW]',size=18)
    ax.yaxis.label.set_size(18)
    gcf().set_size_inches([1.2*dcolwidth,1.2*dcolwidth*0.65])




    plt.tight_layout()
    savefig('./figures/allcases.eps', dpi=300)
    plt.close()

def balbars(): 
    order = [24,11,10,22,3,9,17,13,4,23,12,0,14,16,26,25,15,6,18,2,8,7,19,21,20,5,1]

    load_noc=[]
    load_tod=[]
    load_int=[]
    load_99p=[]
    load_cop=[]

    M=Nodes(load_filename="Case_A_Beta_0.0.npz")
    N=Nodes(load_filename="TodayCap.npz")
    O=Nodes(load_filename="Case_B_Beta_0.4.npz")
    P=Nodes(load_filename="99PCap.npz")
#    Q=Nodes(load_filename="copper.npz")

    names=['EU','']
    for o in order:
        names.append(str(N[o].label))

    #gcf().set_size_inches([15,5])
    width=0.95

    EU_delta_noc=sum(n.balancing for n in M)
    EU_delta_tod=sum(m.balancing for m in N)
    EU_delta_int=sum(o.balancing for o in O)    
    EU_delta_99p=sum(n.balancing for n in P)
#    EU_delta_cop=sum(m.balancing for m in Q)


    EU_load=sum(n.nhours*n.mean for n in N)
    #load_today.append(sum(get_positive(EU_delta_TOD))/EU_load)
    #load_99Q.append(sum(get_positive(EU_delta_99Q))/EU_load)
    #load_NoC.append(sum(get_positive(EU_delta_NoC))/EU_load)
    load_noc.append(sum(EU_delta_noc)/EU_load)
    load_tod.append(sum(EU_delta_tod)/EU_load)
    load_int.append(sum(EU_delta_int)/EU_load)
    load_99p.append(sum(EU_delta_99p)/EU_load)
#    load_cop.append(sum(EU_delta_cop)/EU_load)

    load_noc.append(0)
    load_tod.append(0)
    load_int.append(0)
    load_99p.append(0)
#    load_cop.append(0)

    for o in order:
        load_noc.append(sum(get_positive(M[o].balancing))/(N[o].nhours*N[o].mean))
        load_tod.append(sum(get_positive(N[o].balancing))/(N[o].nhours*N[o].mean))
        load_int.append(sum(get_positive(O[o].balancing))/(N[o].nhours*N[o].mean))
        load_99p.append(sum(get_positive(P[o].balancing))/(N[o].nhours*N[o].mean))
#        load_cop.append(sum(get_positive(Q[o].balancing))/(N[o].nhours*N[o].mean))
    load_tod[17]=load_int[17]
    ind=np.arange(29)
    ax = subplot(1,1,1)

    rectsM=bar(ind*1.25,load_noc,width=width,align='center',color=(1,1,1),label="No transmission")
    rectsN=bar(ind*1.25,load_tod,width=width,align='center',color=rolor,alpha=0.5,label="Present layout")
    rectsO=bar(ind*1.25,load_int,width=width,align='center',color=rolor,alpha=1.0,label="Intermediate layout")
    rectsP=bar(ind*1.25,load_99p,width=width,align='center',color=rolor,alpha=1.0,label=r"$99^{\rm th}$ Q layout", hatch='//')
    #rectsQ=bar(ind*1.25,load_noc,width=width,align='center',color=(0,0,0),label="Unconstrained Flows")

    axis(ymin=0,ymax=.40001,xmin=amin(ind)-.875,xmax=amax(ind*1.25)+.875)
    gcf().set_size_inches([ 3*colwidth , 3*colwidth*0.50])#[9*1.5,3*1.75])
    gcf().set_dpi(400)
    yticks([0,.1,.2,.3,.4])
    xticks(ind*1.25+.35,names,rotation=60,ha='right',va='top')
    ylabel(r'Residual load [normalised]')
    pp = (rectsM,rectsN,rectsO,rectsP)#,rectsQ)
    pp_txtlabels = ("No transmission","Present layout","Intermediate layout",r"99% Q layout")#, "Unconstrained Flows")
    leg = legend(pp,pp_txtlabels,loc='upper right',ncol=len(pp));
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    plt.tight_layout()
    savefig('./figures/balbar.pdf', dpi=300)
    plt.close()

def findneighbours(i,matr=np.genfromtxt("./settings/admat.txt")):
    """ matr should be = np.genfromtxt("./settings/admat.txt")"""
    neighbours=[]
    for j in range(len(matr[i])):
        if matr[i][j]>0 and i<>j:
            neighbours.append(j)
    return neighbours

def scatterneigh(m,N,matr=np.genfromtxt("./settings/admat.txt")):
    close()
    neigh=findneighbours(m,matr)
    ymism=sum(N[i].mismatch*1 for i in neigh)
    ymism/=max(ymism)
    xmism=N[m].mismatch*1
    xmism/=max(xmism)
    xbin=int((max(xmism)-min(xmism))/0.05)
    ybin=int((max(ymism)-min(ymism))/0.05)
    #scatter(xmism,ymism)
    ax = subplot(1,1,1)
    hexbin(xmism,ymism,gridsize=(xbin,ybin),mincnt=1,cmap=matplotlib.cm.gist_heat_r,linewidths=1)
    xlabel(r'Mismatch: '+str(N[m].label))
    ylabel(r'Mismatch: Neighbours')
    ax.xaxis.set_ticks([0])
    ax.yaxis.set_ticks([0])
    ax.grid(1,linewidth=1.5, linestyle='solid')
    gcf().set_size_inches([ 1.5*colwidth , 1.5*colwidth*0.7])#[9*1.5,3*1.75])
    plt.tight_layout()
    savefig("./figures/scatter_type_1_"+str(N[m].label)+".eps")

def linkstable():
    N=Nodes()
    a,h,F=AtoKh_old(N)
    H=get_quant(0.99)
    CP=get_quant(1.0)
    D=H*0.4
    #for f in F:
    #    print f[0],"to",f[1]," & ",round(h[2*f[2]]/1000,2)," & ",round(D[2*f[2]]/1000,2)," & ",round(H[2*f[2]]/1000,2)," &",round(CP[2*f[2]]/1000,2)," & ",f[1],"to",f[0]," & ",round(h[1+2*f[2]]/1000,2)," & ",round(D[1+2*f[2]]/1000,2)," & ",round(H[1+2*f[2]]/1000,2)," &",round(CP[1+2*f[2]]/1000,2),"\\\\"
    print "Total of today's is ", sum(biggestpair(h))
    print "Total of Interme's is ", sum(biggestpair(D))
    print "Total of 99p's is ", sum(biggestpair(H))
    print "Total of copper plate's is ", sum(biggestpair(CP))

def draw_static_networks(N,H,F,tit="1"): ## All the network figures
    close()
    if N==None:
        N=Nodes()
    G=nx.Graph()
    J=nx.Graph()
    X=nx.Graph()
    nodelist=[]

    for n in N:
        G.add_node(str(n.label))
        J.add_node(str(n.label))
        X.add_node(str(n.label))
        nodelist.append(str(n.label))

    K,h,ListF=AtoKh_old(N)
    for l in ListF:
        w1=max(h[l[2]*2],h[l[2]*2+1])
        w2=max(H[l[2]*2],H[l[2]*2+1])
        w3=max(F[l[2]*2],F[l[2]*2+1])        
        if l[2]==41:
            w1=h[l[2]*2]
        G.add_edge(l[0], l[1] , weight= w1)
        J.add_edge(l[0], l[1] , weight= w2)
        X.add_edge(l[0], l[1] , weight= w3)

  
    pos=nx.spring_layout(G)
    pos['AT']=[0.55,0.45]
    pos['FI']=[1.0,1.0]
    pos['NL']=[0.40,0.85]
    pos['BA']=[0.65,0.15]
    pos['FR']=[0.15,0.60]
    pos['NO']=[0.5,1.05]
    pos['BE']=[0.275,0.775]
    pos['GB']=[0.15,0.85]
    pos['PL']=[0.75,0.8]
    pos['BG']=[0.9,0.0]
    pos['GR']=[0.7,0.0]
    pos['PT']=[0.0,0.15]
    pos['CH']=[0.4,0.45]
    pos['HR']=[0.75,0.3]
    pos['RO']=[1.0,0.15]
    pos['CZ']=[0.75,0.60]
    pos['HU']=[1.0,0.45]
    pos['RS']=[0.85,0.15]
    pos['DE']=[0.45,0.7]
    pos['IE']=[0.0,0.95]
    pos['SE']=[0.75,1.0]
    pos['DK']=[0.5,0.875]
    pos['IT']=[0.4,0.2]
    pos['SI']=[0.55,0.3]
    pos['ES']=[0.15,0.35]
    pos['LU']=[0.325,0.575]
    pos['SK']=[0.90,0.55]

    fig = figure(dpi=400,figsize=(colwidth,colwidth*2.4))

    rolor= ( (16+3)/255. , (4*16+11)/255. , (7*16+12)/255.)
    ax1= fig.add_axes([-0.125,.68,1.25,.31]) #For displaying graph    
    nx.draw_networkx_nodes(G,pos,node_size=200,nodelist=nodelist,node_color=rolor,facecolor=(1,1,1))
    e0=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']<=700]
    e1=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>700 and d['weight']<=1200]
    e2=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1200 and d['weight']<=1800]
    e3=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1800 and d['weight']<=2400]
    e4=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>2400 and d['weight']<=3300]
    e5=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>3300 and d['weight']<=4000]
    e6=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>4000 and d['weight']<=5500]
    e7=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>5500 and d['weight']<=8000]
    e8=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>8000 and d['weight']<=12000]
    e9=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>12000]
    ax1.text(-0.05,1.05,"(a)",fontsize=12)
    nx.draw_networkx_edges(G,pos,edgelist=e0,width=.5,edge_color='k',style='dotted')
    nx.draw_networkx_edges(G,pos,edgelist=e1,width=1.0,edge_color='k',alpha=0.4,style='dashed')
    nx.draw_networkx_edges(G,pos,edgelist=e2,width=1.5,edge_color='k',alpha=.7,style='dashed')
    nx.draw_networkx_edges(G,pos,edgelist=e3,width=1.00,edge_color='k',alpha=.2)
    nx.draw_networkx_edges(G,pos,edgelist=e4,width=2.00,edge_color='k',alpha=.2)
    nx.draw_networkx_edges(G,pos,edgelist=e5,width=3.00,edge_color='k',alpha=.2)
    nx.draw_networkx_edges(G,pos,edgelist=e6,width=4.00,edge_color='k',alpha=.2)
    nx.draw_networkx_edges(G,pos,edgelist=e7,width=5.00,edge_color='k',alpha=.2)
    nx.draw_networkx_edges(G,pos,edgelist=e8,width=5.00,edge_color='k',alpha=0.5)
    nx.draw_networkx_edges(G,pos,edgelist=e9,width=5.00,edge_color='k',alpha=1.0)
    nx.draw_networkx_labels(G,pos,font_size=8,font_color='w',font_family='sans-serif')
    ax1.axis('off')

    ax2= fig.add_axes([-0.125,.40,1.25,.31]) #For displaying graph    
    nx.draw_networkx_nodes(J,pos,node_size=200,nodelist=nodelist,node_color=rolor,facecolor=(1,1,1))
    e0=[(u,v) for (u,v,d) in J.edges(data=True) if d['weight']<=700]
    e1=[(u,v) for (u,v,d) in J.edges(data=True) if d['weight']>700 and d['weight']<=1200]
    e2=[(u,v) for (u,v,d) in J.edges(data=True) if d['weight']>1200 and d['weight']<=1800]
    e3=[(u,v) for (u,v,d) in J.edges(data=True) if d['weight']>1800 and d['weight']<=2400]
    e4=[(u,v) for (u,v,d) in J.edges(data=True) if d['weight']>2400 and d['weight']<=3300]
    e5=[(u,v) for (u,v,d) in J.edges(data=True) if d['weight']>3300 and d['weight']<=4000]
    e6=[(u,v) for (u,v,d) in J.edges(data=True) if d['weight']>4000 and d['weight']<=5500]
    e7=[(u,v) for (u,v,d) in J.edges(data=True) if d['weight']>5500 and d['weight']<=8000]
    e8=[(u,v) for (u,v,d) in J.edges(data=True) if d['weight']>8000 and d['weight']<=12000]
    e9=[(u,v) for (u,v,d) in J.edges(data=True) if d['weight']>12000]
    ax2.text(-.05,1.05,"(b)",fontsize=13)
    nx.draw_networkx_edges(G,pos,edgelist=e0,width=.5,edge_color='k',style='dotted')
    nx.draw_networkx_edges(G,pos,edgelist=e1,width=1.0,edge_color='k',alpha=0.4,style='dashed')
    nx.draw_networkx_edges(G,pos,edgelist=e2,width=1.5,edge_color='k',alpha=.7,style='dashed')
    nx.draw_networkx_edges(G,pos,edgelist=e3,width=1.00,edge_color='k',alpha=.2)
    nx.draw_networkx_edges(G,pos,edgelist=e4,width=2.00,edge_color='k',alpha=.2)
    nx.draw_networkx_edges(G,pos,edgelist=e5,width=3.00,edge_color='k',alpha=.2)
    nx.draw_networkx_edges(G,pos,edgelist=e6,width=4.00,edge_color='k',alpha=.2)
    nx.draw_networkx_edges(G,pos,edgelist=e7,width=5.00,edge_color='k',alpha=.2)
    nx.draw_networkx_edges(G,pos,edgelist=e8,width=5.00,edge_color='k',alpha=0.5)
    nx.draw_networkx_edges(G,pos,edgelist=e9,width=5.00,edge_color='k',alpha=1.0)
    nx.draw_networkx_labels(J,pos,font_size=8,font_color='w',font_family='sans-serif')
    ax2.axis('off')

    ax3= fig.add_axes([-0.125,0.12,1.25,.31]) #For displaying graph    
    nx.draw_networkx_nodes(X,pos,node_size=200,nodelist=nodelist,node_color=rolor,facecolor=(1,1,1))
    e0=[(u,v) for (u,v,d) in X.edges(data=True) if d['weight']<=700]
    e1=[(u,v) for (u,v,d) in X.edges(data=True) if d['weight']>700 and d['weight']<=1200]
    e2=[(u,v) for (u,v,d) in X.edges(data=True) if d['weight']>1200 and d['weight']<=1800]
    e3=[(u,v) for (u,v,d) in X.edges(data=True) if d['weight']>1800 and d['weight']<=2400]
    e4=[(u,v) for (u,v,d) in X.edges(data=True) if d['weight']>2400 and d['weight']<=3300]
    e5=[(u,v) for (u,v,d) in X.edges(data=True) if d['weight']>3300 and d['weight']<=4000]
    e6=[(u,v) for (u,v,d) in X.edges(data=True) if d['weight']>4000 and d['weight']<=5500]
    e7=[(u,v) for (u,v,d) in X.edges(data=True) if d['weight']>5500 and d['weight']<=8000]
    e8=[(u,v) for (u,v,d) in X.edges(data=True) if d['weight']>8000 and d['weight']<=12000]
    e9=[(u,v) for (u,v,d) in X.edges(data=True) if d['weight']>12000]
    ax3.text(-0.05,1.05,"(c)",fontsize=14)
    nx.draw_networkx_edges(G,pos,edgelist=e0,width=.5,edge_color='k',style='dotted')
    nx.draw_networkx_edges(G,pos,edgelist=e1,width=1.0,edge_color='k',alpha=0.4,style='dashed')
    nx.draw_networkx_edges(G,pos,edgelist=e2,width=1.5,edge_color='k',alpha=.7,style='dashed')
    nx.draw_networkx_edges(G,pos,edgelist=e3,width=1.00,edge_color='k',alpha=.2)
    nx.draw_networkx_edges(G,pos,edgelist=e4,width=2.00,edge_color='k',alpha=.2)
    nx.draw_networkx_edges(G,pos,edgelist=e5,width=3.00,edge_color='k',alpha=.2)
    nx.draw_networkx_edges(G,pos,edgelist=e6,width=4.00,edge_color='k',alpha=.2)
    nx.draw_networkx_edges(G,pos,edgelist=e7,width=5.00,edge_color='k',alpha=.2)
    nx.draw_networkx_edges(G,pos,edgelist=e8,width=5.00,edge_color='k',alpha=0.5)
    nx.draw_networkx_edges(G,pos,edgelist=e9,width=5.00,edge_color='k',alpha=1.0)
    nx.draw_networkx_labels(X,pos,font_size=8,font_color='w',font_family='sans-serif')
    ax3.axis('off')
    #plt.tight_layout()

    ax4= fig.add_axes([-0.075,0.0,1.5,.15]) #For displaying graph
    #ax4.plot([0,.65],[0.55,.55],color='k') 
#    ax4.vlines(0.06*1.05+0.02,0.55,0.65,color='k')   
#    ax4.vlines(0.12*1.05+0.02,0.55,0.65,color='k') 
#    ax4.vlines(0.18*1.05+0.02,0.55,0.65,color='k') 
#    ax4.vlines(0.24*1.05+0.02,0.55,0.65,color='k') 
#    ax4.vlines(0.30*1.05+0.02,0.55,0.65,color='k') 
#    ax4.vlines(0.36*1.05+0.02,0.55,0.65,color='k') 
#    ax4.vlines(0.42*1.05+0.02,0.55,0.65,color='k') 
#    ax4.vlines(0.48*1.05+0.02,0.55,0.65,color='k') 
#    ax4.vlines(0.54*1.05+0.02,0.55,0.65,color='k') 
#    ax4.vlines(0.60*1.05+0.02,0.55,0.65,color='k') 
    ax4.vlines(0.06*1.05+0.025,0.6,1.0,linewidth=0.5,color='k',linestyles='dotted')
    ax4.vlines(0.12*1.05+0.025,0.6,1.0,linewidth=1.0,color='k',alpha=0.4,linestyle='dashed')
    ax4.vlines(0.18*1.05+0.025,0.6,1.0,linewidth=1.5,color='k',alpha=0.7,linestyle='dashed')
    ax4.vlines(0.24*1.05+0.025,0.6,1.0,linewidth=1.0,color='k',alpha=0.2)
    ax4.vlines(0.30*1.05+0.025,0.6,1.0,linewidth=2.,color='k',alpha=0.2)
    ax4.vlines(0.36*1.05+0.025,0.6,1.0,linewidth=3.,color='k',alpha=0.2)
    ax4.vlines(0.42*1.05+0.025,0.6,1.0,linewidth=4.,color='k',alpha=0.2)
    ax4.vlines(0.48*1.05+0.025,0.6,1.0,linewidth=5.0,color='k',alpha=0.2)
    ax4.vlines(0.54*1.05+0.025,0.6,1.0,linewidth=5.0,color='k',alpha=0.5)
    ax4.vlines(0.60*1.05+0.025,0.6,1.0,linewidth=5.0,color='k',alpha=1.0)
    ax4.text(0.06*1.05+0.01,0.5,"$\leq$ 700 MW",fontsize=9,rotation=-60)
    ax4.text(0.12*1.05+0.01,0.5,"$\leq$ 1.2 GW",fontsize=9,rotation=-60)
    ax4.text(0.18*1.05+0.01,0.5,"$\leq$ 1.8 GW",fontsize=9,rotation=-60)
    ax4.text(0.24*1.05+0.01,0.5,"$\leq$ 2.4 GW",fontsize=9,rotation=-60)
    ax4.text(0.30*1.05+0.01,0.5,"$\leq$ 3.3 GW",fontsize=9,rotation=-60)
    ax4.text(0.36*1.05+0.01,0.5,"$\leq$ 4.0 GW",fontsize=9,rotation=-60)
    ax4.text(0.42*1.05+0.01,0.5,"$\leq$ 5.5 GW",fontsize=9,rotation=-60)
    ax4.text(0.48*1.05+0.01,0.5,"$\leq$ 8.0 GW",fontsize=9,rotation=-60)
    ax4.text(0.54*1.05+0.01,0.5,"$\leq$ 12 GW",fontsize=9,rotation=-60)
    ax4.text(0.60*1.05+0.01,0.5,"$\leq$ 30 GW",fontsize=9,rotation=-60)
    ax4.axis([0.0,1.0,0.0,1.2])
   
    ax4.axis('off')
    savefig("./figures/network"+tit+".pdf")
    #show() # display

def gormhex(N, gamma=1., path='./figures/'):

    matplotlib.rcParams['font.size'] = 10
    majorFormatter = FormatStrFormatter('%.1f')
    EU_mismatch=sum(n.mismatch for n in N)
    bin_size = .1
    for n in N:
        EU_no_n=EU_mismatch - n.mismatch
        print "now doing "+ISO2LONG(n.label)
        nm=n.mismatch/max(n.mismatch)
        EUm=EU_no_n/max(EU_no_n)
        figure(1); clf()
   
        gcf().set_dpi(300)
        gcf().set_size_inches([colwidth,2.8])
   
        #Calculate number of bins
        x_bins = int((max(nm)-min(nm))/bin_size)
        y_bins = int((max(EUm)-min(EUm))/bin_size)
        print "     now doing hexbining"   
        H = hexbin(nm,EUm,gridsize=(x_bins,y_bins),mincnt=1,cmap=matplotlib.cm.gist_heat_r,linewidths=0)
        print "     now done hexbining"        
        values = H.get_array()
        conf_lvl = [.9,.5,.1]
       
        levels=zeros(len(conf_lvl))
        sort_values = sort(values)
        cumsum_values = cumsum(sort_values)/sum(values)
        for j in arange(len(conf_lvl)):
            levels[j] = sort_values[argmax(cumsum_values>=(1-conf_lvl[j]))]
       
        #levels = mquantiles(values,quantiles)
        print levels
       
        paths = H.get_paths()
        xc, yc = zeros(len(values)), zeros(len(values))
        for k in arange(len(values)):
            xc[k] = mean(paths[k].vertices[0:6].transpose()[0])
            yc[k] = mean(paths[k].vertices[0:6].transpose()[1])
       
        X, Y = meshgrid(linspace(amin(xc),amax(xc),2*x_bins),linspace(amin(yc),amax(yc),2*y_bins))
       
        V = griddata(xc,yc,values,X, Y)
        V[1:-1,1:-1] = (V[1:-1,1:-1] + V[0:-2,1:-1] + V[2:,1:-1] + V[1:-1,0:-2] + V[1:-1,2:])/5.
        V[1:-1,1:-1] = (V[1:-1,1:-1] + V[0:-2,1:-1] + V[2:,1:-1] + V[1:-1,0:-2] + V[1:-1,2:])/5.
       
        #pcolor(X,Y,V)
       
        CS = contour(X,Y,V,levels,colors='k')
       
        ctxt = dict([(levels[k],'{0:.0f}%'.format(100*conf_lvl[k])) for k in arange(len(levels))])
        clabel(CS,fmt=ctxt, fontsize=8, inline=1, inline_spacing=10)
       
       
        axis(xmin=-1.5,xmax=2.3,ymin=-1.2,ymax=1.2)
        xticks(arange(-1.,2.51,1.))
        yticks(arange(-1.,1.2,.5))
        #colorbar()
        clim(0,500)
        #plot(mismatch_X[i],mismatch_EU_no_X,'k.')
   
        #Plot cross
        axhline(0,color='k',ls='-',lw=1.)
        axvline(0,color='k',ls='-',lw=1.)
   
       
   
        xlabel(r'$\Delta^{0:}$ [normalised]'.format('{'+str(n.label)+'}'))
        ylabel(r'$\Delta^{EU\backslash{} '+r'{0:}'.format(str(n.label))+r'}$ [normalised]')
        #ylabel(r'$\Delta^{EU\backslash X}$ [av.l.h.]')
   
        #[left, bottom, width, height]
        gca().set_position([.16, .135, .835, .86], which='both')
        #Text labels
        bb = gca().get_position()
       
       
        ### Text labels:
        shortage_perc = sum((n.mismatch<=0)*(EU_no_n<=0))*100/len(n.mismatch)
        import_perc = sum((n.mismatch<0)*(EU_no_n>0))*100/len(n.mismatch)
        export_perc = sum((n.mismatch>0)*(EU_no_n<0))*100/len(n.mismatch)
        excess_perc = sum((n.mismatch>=0)*(EU_no_n>=0))*100/len(n.mismatch)
       
        text(bb.xmin-.14, bb.ymin-.05,'Shortage ({0:.0f}%)'.format(shortage_perc),horizontalalignment='left',verticalalignment='top',transform = gca().transAxes,weight='bold',size='small')
        text(bb.xmin-.14, bb.ymax-.1,'Import ({0:.0f}%)'.format(import_perc),horizontalalignment='left',verticalalignment='bottom',transform = gca().transAxes,weight='bold',size='small')
        text(bb.xmax-.33, bb.ymin-.05,'Export ({0:.0f}%)'.format(export_perc),horizontalalignment='left',verticalalignment='top',transform = gca().transAxes,weight='bold',size='small')
        text(bb.xmax-.33, bb.ymax-.1,'Excess ({0:.0f}%)'.format(excess_perc),horizontalalignment='left',verticalalignment='bottom',transform = gca().transAxes,weight='bold',size='small')
       
        gca().xaxis.set_major_formatter(majorFormatter)
        gca().yaxis.set_major_formatter(majorFormatter)
   
        savefig('./figures/hexbin_gorm_'+str(n.label)+'.pdf',dpi=400)

#L=Nodes()

#M=Nodes(load_filename="Case_A_Beta_1.0.npz")

#N=Nodes(load_filename="Case_B_Beta_0.4.npz")

#O=Nodes(load_filename="99PCap.npz")

#smoothhist(M)
#smoothhist2(M,N,O)
##plotbars1()
plotbars2()
#gormhex(N)
#F=get_quant()
#H=0.4*F
#draw_static_networks(L,H,F)
#balbars()
#plot_allcases()
#show_hist(16,"FRES",L[4].mean)

