#! /usr/bin/env python
from auresc import *
#import networkx as nx
import csv
import sys
from pylab import *
from scipy.optimize import brentq
from scipy.optimize import fmin
#import datetime as dt

colwidth = (3.425)
double_column_width = (2*3.425+0.236)*1.5 #extrawide!
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
sugar_colours = ['#F8CA00','#AD0550','#FF9F22','#490A3D','#9AAB0F','#003B6C']

def get_positive(x):
    """ This function returns the positive side of variable. """
    return x*(x>0.)  #Possibly it has to be x>1e-10.

def ISO2LONG(x):
    table= np.load('./settings/ISET2ISO_country_codes.npy')
    for i in arange(len(table)):
        if table[i]["ISO"]==x: return table[i]['name']

def smoothhist(N):
    for n in N:
        lines=[]
        alphas=[0.0,0.4,0.6,0.8,1.0]
        cnorm=matplotlib.colors.Normalize(vmin=0,vmax=1)
        colormap=matplotlib.cm.ScalarMappable(norm=cnorm,cmap=my_cmap)
    ######
        N.set_alphas(0)
        lo=-n.load/n.mean
        x0=1.10*min(n.mismatch/n.mean)
        x1=1.10*max(n.mismatch/n.mean)
        b=arange(x0,x1,(x1-x0)/250)
        lines.append(100*hist(lo,bins=b)[0]/70128.)
        for a in alphas:
            n.set_alpha(a)
            lines.append(100*hist(n.mismatch/n.mean,bins=b)[0]/70128.)
        plt.close()
        ax=subplot(1,1,1)
        plot(b[0:-1],lines[0],label="Load", linewidth=1.5)
        plot(b[0:-1],lines[1],label="$\Delta$, alpha = 0.0", linewidth=1.5, color=sugar_colours[0]) 
        plot(b[0:-1],lines[2],label="$\Delta$, alpha = 0.4", linewidth=1.5, color=sugar_colours[1])
        plot(b[0:-1],lines[3],label="$\Delta$, alpha = 0.6", linewidth=1.5, color=sugar_colours[2])
        plot(b[0:-1],lines[4],label="$\Delta$, alpha = 0.8", linewidth=1.5, color=sugar_colours[3])
        plot(b[0:-1],lines[5],label="$\Delta$, alpha = 1.0", linewidth=1.5, color=sugar_colours[4])
        plt.ylabel('P($\Delta$) (%)')
        plt.xlabel(r'Mismatch Power ($\langle L \rangle $)')
        plt.tight_layout()
        plt.text(0.075,max(lines[0]*1.00), "Mismatch: "  + ISO2LONG(n.label), size='x-large')
        #plt.title(r'M : '+str(n.label))
        gcf().set_size_inches([ 3*colwidth , 3*colwidth*0.50])#[9*1.5,3*1.75])
        gcf().set_dpi(400)
        plt.ylim(0,max(lines[0])*1.075)
        plt.xlim(-2,4)
        handles,labels=ax.get_legend_handles_labels()
        ax.legend(handles,labels)
        plt.grid(True)
        #show()
        savefig('./figures/SMultialpha_' + str(n.label) +'.pdf', dpi=300)
        plt.close()
    #### Now the same for aggregated EU
    lines=[]
    N.set_alphas(0)
    EU_load = sum(n.load for n in N)
    EU_mean = mean(EU_load)
    EU_load/= EU_mean
    EU_delta = sum(n.mismatch for n in N)
    x0=1.10*min(EU_delta/EU_mean)
    x1=1.10*max(EU_delta/EU_mean)
    b=arange(x0,x1,(x1-x0)/250)
    lines.append(100*hist(-EU_load,bins=b)[0]/70128.)
    for a in alphas:
        N.set_alphas(a)
        EU_delta = sum(n.mismatch for n in N)
        lines.append(100*hist(EU_delta/EU_mean,bins=b)[0]/70128.)
    plt.close()
    ax=subplot(1,1,1)
    plot(b[0:-1],lines[0],label="Load", linewidth=1.5)
    plot(b[0:-1],lines[1],label="$\Delta$, alpha = 0.0", linewidth=1.5, color=sugar_colours[0]) 
    plot(b[0:-1],lines[2],label="$\Delta$, alpha = 0.4", linewidth=1.5, color=sugar_colours[1])
    plot(b[0:-1],lines[3],label="$\Delta$, alpha = 0.6", linewidth=1.5, color=sugar_colours[2])
    plot(b[0:-1],lines[4],label="$\Delta$, alpha = 0.8", linewidth=1.5, color=sugar_colours[3])
    plot(b[0:-1],lines[5],label="$\Delta$, alpha = 1.0", linewidth=1.5, color=sugar_colours[4])
    plt.ylabel('P($\Delta$) (%)')
    plt.xlabel(r'Mismatch Power ($\langle L \rangle $)')
    plt.text(0.075,max(lines[0]*1.00), "Mismatch: EU Aggregated", size='x-large')
    plt.tight_layout()
    gcf().set_size_inches([ 3*colwidth , 3*colwidth*0.5])
    gcf().set_dpi(400)
    plt.ylim(0,max(lines[0])*1.075)
    plt.xlim(-2,4)
    handles,labels=ax.get_legend_handles_labels()
    ax.legend(handles,labels)
    plt.grid(True)
    #show()
    savefig('./figures/SMultialpha_EU.pdf', dpi=300)
    plt.close()

def smoothhist2(N,M):
    for n in N:
        lines=[]
    ######
        lo=-n.load/n.mean
        x0=1.10*min(n.mismatch/n.mean)
        x1=1.10*max(n.mismatch/n.mean)
        b=arange(x0,x1,(x1-x0)/250)
        lines.append(100*hist(lo,bins=b)[0]/70128.)
        lines.append(100*hist(n.mismatch/n.mean,bins=b)[0]/70128.)
        u=n.curtailment-n.balancing
        v=[]
        for w in u:
            if w>=1 or w<=-1:
                v.append(w/n.mean)
        lines.append(100*hist(v,bins=b)[0]/70128.)
        u=M[n.id].curtailment-M[n.id].balancing
        v=[]
        for w in u:
            if w>=1 or w<=-1:
                v.append(w/n.mean)
        lines.append(100*hist(v,bins=b)[0]/70128.)
        plt.close()
        ax=subplot(1,1,1)
        plot(b[0:-1],lines[0],label="Load", linewidth=1.5)
        plot(b[0:-1],lines[1],label="$\Delta$, no transmission", linewidth=1.5, color=sugar_colours[3]) 
        plot(b[0:-1],lines[2],label="$\Delta$, actual caps.", linewidth=1.5, color=sugar_colours[1])
        plot(b[0:-1],lines[3],label="$\Delta$, $99^{\mathrm{th}}$ Q caps.", linewidth=1.5, color=sugar_colours[4])
        #plot(b[0:-1],lines[4],label="$\Delta$, alpha = 0.8", linewidth=1.5, color=sugar_colours[3])
        #plot(b[0:-1],lines[5],label="$\Delta$, alpha = 1.0", linewidth=1.5, color=sugar_colours[4])
        plt.ylabel('P($\Delta$) (%)')
        plt.xlabel(r'Mismatch Power ($\langle L \rangle $)')
        plt.tight_layout()
        plt.text(-0.5,max(lines[0]*1.00), "Mismatch: "  + ISO2LONG(n.label), size='x-large')
        #plt.title(r'M : '+str(n.label))
        gcf().set_size_inches([ 3*colwidth , 3*colwidth*0.50])#[9*1.5,3*1.75])
        gcf().set_dpi(400)
        plt.ylim(0,max(lines[0])*1.075)
        plt.xlim(-2,2.5)
        handles,labels=ax.get_legend_handles_labels()
        ax.legend(handles,labels)
        plt.grid(True)
        #show()
        savefig('./figures/AAA_' + str(n.label) +'.pdf', dpi=300)
        plt.close()
    #### Now the same for aggregated EU
    lines=[]
    EU_load = sum(n.load for n in N)
    EU_mean = mean(EU_load)
    EU_load/= EU_mean
    EU_delta = sum(n.mismatch for n in N)
    x0=1.10*min(EU_delta/EU_mean)
    x1=1.10*max(EU_delta/EU_mean)
    b=arange(x0,x1,(x1-x0)/250)
    lines.append(100*hist(-EU_load,bins=b)[0]/70128.)
    EU_delta = sum(n.mismatch for n in N)
    lines.append(100*hist(EU_delta/EU_mean,bins=b)[0]/70128.)
    u=sum(n.curtailment-n.balancing for n in N)
    v=[]
    for w in u:
        if w>=1 or w<=-1:
            v.append(w/EU_mean)
    lines.append(100*hist(v,bins=b)[0]/70128.)
    u=sum(m.curtailment-m.balancing for m in M)
    v=[]
    for w in u:
        if w>=1 or w<=-1:
            v.append(w/EU_load)
    lines.append(100*hist(v,bins=b)[0]/70128.)
    plt.close()
    ax=subplot(1,1,1)
    plot(b[0:-1],lines[0],label="Load", linewidth=1.5)
    plot(b[0:-1],lines[1],label="$\Delta$, no transmission", linewidth=1.5, color=sugar_colours[3]) 
    plot(b[0:-1],lines[2],label="$\Delta$, actual caps.", linewidth=1.5, color=sugar_colours[1])
    plot(b[0:-1],lines[3],label="$\Delta$, $99^{\mathrm{th}}$ Q caps.", linewidth=1.5, color=sugar_colours[4])
    #plot(b[0:-1],lines[4],label="$\Delta$, alpha = 0.8", linewidth=1.5, color=sugar_colours[3])
    #plot(b[0:-1],lines[5],label="$\Delta$, alpha = 1.0", linewidth=1.5, color=sugar_colours[4])
    plt.ylabel('P($\Delta$) (%)')
    plt.xlabel(r'Mismatch Power ($\langle L \rangle $)')
    plt.tight_layout()
    plt.text(-0.5,max(lines[0]*1.00), "Mismatch: Aggregated EU", size='x-large')
    #plt.title(r'M : '+str(n.label))
    gcf().set_size_inches([ 3*colwidth , 3*colwidth*0.50])#[9*1.5,3*1.75])
    gcf().set_dpi(400)
    plt.ylim(0,max(lines[0])*1.075)
    plt.xlim(-2,2.5)
    handles,labels=ax.get_legend_handles_labels()
    ax.legend(handles,labels)
    plt.grid(True)
    #show()
    savefig('./figures/AAA_EU.pdf', dpi=300)
    plt.close()

def esshist(N,Ns): ## The mismatch histograms with different alphas
    for n in N:
        lines=[]
        lo=-n.load/n.mean
        x0=-2
        x1=3.001
    ######
        nostor=[]
        sistor=[]
        for x in (n.curtailment-n.balancing):
            if x >1 or x<-1: nostor.append(x)
        for x in (Ns[n.id].curtailment-Ns[n.id].balancing):
            if x >1 or x<-1: sistor.append(x)
        b=np.arange(x0,x1,(x1-x0)/250.)
        lines.append(hist(lo,bins=b)[0]/(70128.))
        lines.append(hist(n.mismatch/n.mean,bins=b)[0]/(70128.))
        lines.append(hist(nostor/n.mean,bins=b)[0]/(70128.))
        lines.append(hist(sistor/n.mean,bins=b)[0]/(70128.))
        plt.close()
        ax=subplot(1,1,1)
        plot(b[0:-1],lines[0],label="Load", linewidth=2.5, color ='r')
        plot(b[0:-1],lines[1],label=r"Mismatch", linewidth=1.5, color=sugar_colours[2]) 
        plot(b[0:-1],lines[2],label=r"Mismatch after exp/imp", linewidth=1.5, color=sugar_colours[4])
        plot(b[0:-1],lines[3],label=r"Mismatch after exp/imp + storage", linewidth=1.5, color=sugar_colours[5])
        plt.ylabel('P($\Delta$)')
        plt.xlabel(r'Mismatch power [normalised]')
        plt.text(-1.9,max(lines[0]*1.00)*0.95, "a)" , size='x-large')
        #plt.title(r'M : '+str(n.label))
        gcf().set_size_inches([ dcolwidth , dcolwidth*0.40])#[9*1.5,3*1.75])
        gcf().set_dpi(400)
        plt.ylim(0,max(lines[0])*1.075)
        plt.yticks(np.arange(0,max(lines[0])*1.075,0.01))
        plt.xlim(-2,3.001)
        plt.xticks([-2,-1,0,1,2,3])
        handles,labels=ax.get_legend_handles_labels()
        leg=ax.legend(handles,labels,title=ISO2LONG(n.label))
        setp(leg.get_texts(),fontsize="small")
        plt.grid(which="major",axis='x')
        plt.tight_layout()
        #show()
        savefig('./figures/ESS_' + str(n.label) +'.pdf', dpi=400)
        plt.close()

    #### Now the same for aggregated EU
    lines=[]
    nostor=[]
    sistor=[]
    EU_load = sum(n.load for n in N)
    EU_mean = mean(EU_load)
    EU_load/= EU_mean
    EU_delta = sum(n.mismatch for n in N)
    b=arange(x0,x1,(x1-x0)/250.)
    lines.append(hist(-EU_load,bins=b)[0]/70128.)
    EU_delta = sum(n.mismatch for n in N)
    EU_cur= sum(n.curtailment for n in N)
    EU_bal= sum(n.balancing for n in N)
    EU_curS= sum(n.curtailment for n in Ns)
    EU_balS= sum(n.balancing for n in Ns)

    for x in (EU_cur-EU_bal):
        if x >1 or x<-1: nostor.append(x)
    for x in (EU_curS-EU_balS):
        if x >1 or x<-1: sistor.append(x)

    x0=-2
    x1=3

    lines.append(hist(EU_delta/EU_mean,bins=b)[0]/70128.)
    lines.append(hist(nostor/EU_mean,bins=b)[0]/70128.)
    lines.append(hist(sistor/EU_mean,bins=b)[0]/70128.)

    plt.close()
    ax=subplot(1,1,1)
    plot(b[0:-1],lines[0],label="Load", linewidth=1.5, color='r')
    plot(b[0:-1],lines[1],label=r"Mismatch", linewidth=1.5, color=sugar_colours[2]) 
    plot(b[0:-1],lines[2],label=r"Mismatch after exp/imp", linewidth=1.5, color=sugar_colours[4])
    plot(b[0:-1],lines[3],label=r"Mismatch after exp/imp + storage", linewidth=1.5, color=sugar_colours[5])
    plt.ylabel('P($\Delta$)')
    plt.xlabel(r'Mismatch power [normalised]')
    plt.text(-1.9,max(lines[0]*1.00)*0.95, "b)" , size='x-large')
    #plt.title(r'M : '+str(n.label))
    gcf().set_size_inches([ dcolwidth , dcolwidth*0.40])#[9*1.5,3*1.75])
    gcf().set_dpi(400)
    plt.ylim(0,max(lines[0])*1.075)
    plt.xlim(-2,3.001)
    plt.xticks([-2,-1,0,1,2,3])
    plt.yticks(np.arange(0,max(lines[0])*1.075,0.01))
    handles,labels=ax.get_legend_handles_labels()
    leg=ax.legend(handles,labels, title="Europe")
    setp(leg.get_texts(),fontsize="small")
    plt.grid(which="major",axis='x')
    plt.tight_layout()
    #show()
    savefig('./figures/ESS_EU.pdf', dpi=300)
    plt.close()

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
    alphas=[0,.4,.6,.8,1.0]
    plt.close()
    ###### Setting colours
    cnorm=matplotlib.colors.Normalize(vmin=0,vmax=1)
    colormap=matplotlib.cm.ScalarMappable(norm=cnorm,cmap=my_cmap)
    ######
    N.set_alphas(0)
    for n in N:
        lo=-n.load/n.mean
        x0=1.10*min(n.mismatch/n.mean)
        x1=1.10*max(n.mismatch/n.mean)
        b=arange(x0,x1,(x1-x0)/250)
        l=hist(lo,normed=True,histtype='step',bins=b,align='mid',color='k',label="Normalized Load: ")
        for a in alphas:
            n.set_alpha(a)
            hist(n.mismatch/n.mean,normed=True,histtype='step',bins=b,linewidth=1.5,color=colormap.to_rgba(a),label="Mismatch, Alpha = "+str(a))
        ax=subplot(1,1,1)    
        plt.ylabel('P(\Delta)')
        plt.xlabel('Power (a.l.h.)')
        #plt.text(0.8*n.load/n.mean,-1,"Mismatch: "  + ISO2LONG(n.id))
        #plt.title(r'M : '+str(n.label))
        gcf().set_size_inches([9*1.5,3*1.75])
        #plt.ylim(0,max(l[0])*1.075)
        handles,labels=ax.get_legend_handles_labels()
        ax.legend(handles,labels)
        plt.grid(True)
        #show()
        savefig('./figures/Multialpha_' + str(n.label) +'.pdf', dpi=300)
        plt.close()

def AtoKh_nonsparse(N,F,G,pathadmat='./settings/admat.txt'):
    Ad=np.genfromtxt(pathadmat,dtype='d')
    L=0
    #G=nx.DiGraph()
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

    for l in range(L):
        if F[l] >= 0 :
            G.add_edge( listFlows[l][0], listFlows[l][1] , weight= F[l] )
        if F[l] <= 0 :
            G.add_edge( listFlows[l][1] , listFlows[l][0] , weight = -F[l] )    
    return K,h,listFlows,G,F

def ttodate(t):
    date0=dt.datetime(2000,1,1)
    dated=dt.timedelta(t/24)
    date1=date0+dated
    return str(date1.year)+'.'+str(date1.month)+'.'+str(date1.day)+' -- '+str(t - t/24*24)+':00'

def show_hist(link,tit,mean,filename='results/copper_flows.npy',e=1,b=500,limit=10000):
    f=np.load(filename)
    flows=[]
#    ax=subplot(1,1,1)
    for i in f:
        flows.append(i)
#        flows.append(i*(i>e))
#        flows.append(i*(i<-e))
#    a=hist(flows[link],bins=b,range=[0.1,flows[link].max()],normed=0,histtype='stepfilled')
#    b=hist(flows[link+1],bins=b,range=[flows[link+1].min(),-0.1],normed=0,histtype='stepfilled')
    a=hist(flows[link]/mean,bins=b,normed=1,histtype='stepfilled',weights=np.zeros((70128))+ 100./70128)
    plt.ylabel('P($F_l$) (%)',size='large')
    plt.xlabel(r'Transmission Magnitude ($\langle L_{DK} \rangle$)',size='large')
    plt.text(-5.3,.47,"Germany to Denmark",size="x-large")
    #plt.title(r'Histogram of Power Flows : '+str(tit))
    vlines(-19785/mean,0,0.025,color='r',linewidth=1.2)
    plt.text(-19785/mean-0.15,0.030,'0.01% Q')
    vlines(-16213/mean,0,0.085,color='r',linewidth=1.2)
    plt.text(-16213/mean-0.45,0.090,'0.1% Q')
    vlines(-11266/mean,0,0.10,color='r',linewidth=1.2)
    plt.text(-11266/mean-0.05,0.105,'1.0% Q')
    vlines(-6045/mean,0,.15,color='r',linewidth=1.2)
    plt.text(-6045/mean-0.05,.155,'5.0% Q')
    vlines(7328/mean,0,.15,color='r',linewidth=1.2)
    plt.text(7328/mean-0.02,.155,'95% Q')
    vlines(12941/mean,0,0.10,color='r',linewidth=1.2)
    plt.text(12941/mean-0.02,0.105,'99% Q')
    vlines(18205/mean,0,0.085,color='r',linewidth=1.2)
    plt.text(18205/mean-0.04,0.09,'99.9% Q')
    vlines(20838/mean,0,0.025,color='r',linewidth=1.2)
    plt.text(20838/mean-0.07,0.03,'99.99% Q')

#    plt.vlines(0,0,1000,linestyle='dashed')
#    plt.axvline(15357.73,0,800,linestyle='dashed',color='r')
#    plt.axvline(7678.86,0,800,linestyle='dashed',color='r')
#    plt.axvline(3993.00,0,800,linestyle='dashed',color='r')
#    plt.axvline(307.15,0,800,linestyle='dashed',color='r')
#    plt.text(15457,400,'99Q')
#    plt.text(7778,500,'95Q')
#    plt.text(4093,600,'90Q')
#    plt.text(407,700,'75Q')
    gcf().set_size_inches([3*colwidth,3*colwidth*0.5])
    plt.axis([-5.35,5.65,0,0.5])#plt.axis([0.5, 1, 0, 38000])
    plt.grid(True)
    plt.tight_layout
    savefig('./figures/' + str(link) +'.pdf', dpi=400)
    plt.close()
    #show()


def savehist(link,name):
    f=np.load('results/copper_flows.npy')
    N=Nodes()
    a,b,c,d=AtoKh_nonsparse(N)
    a=hist(f[link],bins=500,normed=1)
    b=[a[1][:-1],a[0]]
    save('./figures/'+c[link][0],b)
#    g=csv.writer(open(name,'wb'),delimiter=',')
#    for column in f[link]: 
#        g.writerow(column)

def get_transmissions():
    NA=Nodes(load_filename='Case_A_Beta_1.0.npz')
    a,b,ha,d=AtoKh(NA)
    hc=get_quant(0.9999)
    hq=get_quant(0.99)
    hd=0.4*get_quant(0.99)

    HC=biggestpair(hc)
    HQ=biggestpair(hq)
    HD=biggestpair(hd)
    HA=biggestpair(ha)
    i=0
    while i<(44):
        print d[i][0],'&',str(round(hc[2*i]/1000,2)),'&',str(round(hq[2*i]/1000,2)),'&',str(round(hd[2*i]/1000,2)),'&',str(round(ha[2*i]/1000,2)),'\\\\'
        print '$\\blacktriangleleft$','&',str(round(hc[2*i+1]/1000,2)),'&',str(round(hq[2*i+1]/1000,2)),'&',str(round(hd[2*i+1]/1000,2)),'&',str(round(ha[2*i+1]/1000,2)),'\\\\'   
        i+=1 
    print 'Total & ', sum(HC),'&', sum(HQ),'&', sum(HD),'&', sum(HA),'\\\\'


def importexport():
    NC=Nodes(load_filename='copper_nodes.npz')
    NQ=Nodes(load_filename='Case_C_Beta_1.0.npz')
    ND=Nodes(load_filename='Case_C_Beta_0.4.npz')
    NA=Nodes(load_filename='Case_A_Beta_1.0.npz')
    impwanted=np.ones(28)
    imp=np.ones(28)
    i=[]
    imports=[]
    expwanted=np.ones(28)
    expo=np.ones(28)
    e=[]
    exports=[]
    CASES=[NC,NQ,ND,NA]
    imports=np.zeros((28,4))
    exports=np.zeros((28,4))
    for N in CASES:
        impwanted[27]=0
        expwanted[27]=0
        for n in N:
            impwanted[n.id]=sum(get_positive(-1*n.mismatch))
            expwanted[n.id]=sum(get_positive(n.mismatch))
            imp[n.id]=impwanted[n.id]-sum(n.balancing)
            #imp[27]=sum(imp[:-1])
            expo[n.id]=expwanted[n.id]-sum(n.curtailment)
            #expo[27]=sum(expo[:-1])
            x=imp[n.id]/impwanted[n.id]
            y=expo[n.id]/expwanted[n.id]
            #if n.id == 19:
            #    print impwanted,imp,x
            i.append(x)
            e.append(y)
        impwanted[27] = sum(impwanted[:-1])
        expwanted[27] = sum(expwanted[:-1])
        imp[27]=sum(imp[:-1])
        expo[27]=sum(imp[:-1])
        i.append(imp[27]/impwanted[27])
        e.append(expo[27]/expwanted[27])
        if N==NC: p=0    
        if N==NQ: p=1
        if N==ND: p=2
        if N==NA: p=3
        imports[:,p]=i[p*28:(p+1)*28]
        exports[:,p]=e[p*28:(p+1)*28]       

    return imports,exports


#scenario A
def Case_A(betas=[0.0,0.25,0.5,0.75,1.0,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.0,2.25,2.50,2.75,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.5,15.0,17.5,20.0,25.0,30.0]):
    N=Nodes()    
    for b in betas:
        N,F=sdcpf(N,b=b)
        N.save_nodes('Case_A_Beta_'+str(b))
        save('./results/'+'Flows_Case_A_Beta_'+str(b),F)

def Case_B(links=np.arange(0.0,30000.1,1500.0)):
    hopt=get_quant(.99)
    h0=get_quant(.99)
    N=Nodes()
    for l in links:
        for h in range(len(hopt)):
            h0[h]=l
            if hopt[h]<l: h0[h]=hopt[h]
        N,F=sdcpf(N,h0=h0)
        N.save_nodes('Case_B_Link_'+str(l))
        save('./results/'+'Flows_Case_B_Link_'+str(l),F)

#scenario C
def Case_C(h0=None,betas=[0.0,0.05,0.1,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,1.0,1.05,1.1,1.2,1.3,1.4,1.5,2.0,2.5,3.0,4.0,5.0]):
    if h0 == None: h0=get_quant(.99)
    N=Nodes()
    for b in betas:
        N,F=sdcpf(N,b=b,h0=h0)
        N.save_nodes('Case_C_Beta_'+str(b))
        save('./results/'+'Flows_Case_C_Beta_'+str(b),F)


#0.60,0.70,0.75,0.80,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.9920,.9940,.9960,.9980,.9999
def Case_D(quants=[0.75,0.8,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.9920,.9940,.9960,.9980,.9999]):
    N=Nodes()
    for q in quants:
        h=get_quant(q)
#        if sum(h) >=1:
#            for i in range(len(h)):
#                if h[i]==0: h[i]=50
#            print h
        N,F=sdcpf(N,h0=h)
        N.save_nodes('Case_D_Quant_'+str(q))
        save('./results/'+'Flows_Case_D_Quant_'+str(q),F)

def biggestpair(H):
    H0=np.zeros((len(H))/2)
    for i in range(len(H0)):
        H0[i]=max(H[2*i],H[2*i+1])
    return H0

def Plot_A():
    betas=[0.0,0.25,0.5,0.75,1.0,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.0,2.25,2.50,2.75,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.5]
    PlotA=np.zeros((len(betas),2))
    j=0
    for b in betas:
        PlotA[j,0]=b
        N=Nodes(load_filename='Case_A_Beta_'+str(b)+'.npz')
        a=0
        d=0
        for i in N:
            a+=np.sum(i.balancing)
            d+=i.mean*i.nhours
        c=a/d
        PlotA[j,1]=c
        j+=1
    save('./results/PlotA',PlotA)
    return PlotA


def Plot_B():
    links=np.arange(0.0,30000.1,1500.0)
    N=Nodes()
    K,k2,Hac,lF=AtoKh(N)
    Hact=biggestpair(Hac)
    hopt=get_quant(.99)
    h0=get_quant(.99)
    PlotB=np.zeros((len(links),2))
    j=0
    for l in links:
        for h in range(len(hopt)):
            h0[h]=l
            if hopt[h]<l: h0[h]=hopt[h]
        Hopt=biggestpair(h0)
        PlotB[j,0]=sum(Hopt)/sum(Hact)
        N=Nodes(load_filename='Case_B_Link_'+str(l)+'.npz')
        a=0
        d=0
        for i in N:
            a+=sum(i.balancing)
            d+=i.mean*i.nhours
        c=a/d
        PlotB[j,1]=c
        j+=1
    save('./results/PlotB',PlotB)
    return PlotB

def Plot_C():
    betas=[0.0,0.05,0.1,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,1.0,1.05,1.1,1.2,1.3,1.4,1.5,2.0,2.5,3.0,4.0,5.0]
    N=Nodes()
    K,k2,Hac,lF=AtoKh(N)
    Hact=biggestpair(Hac)
    Hop=get_quant(.99)
    Hopt=biggestpair(Hop)
    PlotC=np.zeros((len(betas),2))
    j=0
    for b in betas:
        PlotC[j,0]=b*sum(Hopt)/sum(Hact)
        N=Nodes(load_filename='Case_C_Beta_'+str(b)+'.npz')
        a=0
        d=0
        for i in N:
            a+=sum(i.balancing)
            d+=i.mean*i.nhours
        c=a/d
        PlotC[j,1]=c
        j+=1
    save('./results/PlotC',PlotC)
    return PlotC

def Plot_D():
    quants=[0.60,0.70,0.75,0.80,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.9920,.9940,.9960,.9980,.9999]
    N=Nodes()
    K,k2,Hac,lF=AtoKh(N)
    Hact=biggestpair(Hac)
    PlotD=np.zeros((len(quants),2))
    j=0
    for q in quants:
        Hop=get_quant(q)
        Hopt=biggestpair(Hop)
        PlotD[j,0]=sum(Hopt)/sum(Hact)
        N=Nodes(load_filename='Case_D_Quant_'+str(q)+'.npz')
        a=0
        d=0
        for i in N:
            a+=sum(i.balancing)
            d+=i.mean*i.nhours
        c=a/d
        PlotD[j,1]=c
        j+=1
    save('./results/PlotD',PlotD)
    return PlotD

def drawgrid(N):
    a,b,c,G=AtoKh_nonsparse(N)
    pos=nx.spring_layout(G)
    pos['FI']=[1.0,1.0]
    pos['SE']=[0.75,1.0]
    pos['NO']=[0.5,1.0]
    pos['DK']=[0.5,0.875]
    pos['PL']=[0.75,0.8]
    pos['GR']=[0.7,0.0]
    pos['BG']=[0.9,0.0]
    pos['SK']=[0.90,0.55]
    pos['CZ']=[0.75,0.60]
    pos['HU']=[1.0,0.45]
    pos['RO']=[1.0,0.15]
    pos['RS']=[0.85,0.15]
    pos['BA']=[0.65,0.15]
    pos['HR']=[0.75,0.3]
    pos['IE']=[0.0,0.95]
    pos['GB']=[0.15,0.85]
    pos['FR']=[0.15,0.60]
    pos['ES']=[0.15,0.35]
    pos['PT']=[0.0,0.15]
    pos['BE']=[0.3,0.8]
    pos['NL']=[0.40,0.85]
    pos['LU']=[0.325,0.575]
    pos['DE']=[0.45,0.7]
    pos['CH']=[0.4,0.45]
    pos['IT']=[0.4,0.2]
    pos['AT']=[0.55,0.45]
    pos['SI']=[0.55,0.3]



    nx.draw_networkx_nodes(G,pos,node_size=1600,node_color='b',facecolor=(1,1,1))
    esmall=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']<=500]
    emid=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>500 and d['weight']<=1500]
    elarge=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1500]

    nx.draw_networkx_edges(G,pos,edgelist=esmall,width=2)#1
    nx.draw_networkx_edges(G,pos,edgelist=emid,width=2)#3
    nx.draw_networkx_edges(G,pos,edgelist=elarge,width=2)#6
    nx.draw_networkx_labels(G,pos,font_size=24,font_color='w',font_family='sans-serif')
    axis('off')
    savefig("weighted_graph.png") # save as png
    show() # display

def plotbars(imports,exports): 
    order = [24,11,10,22,3,9,17,13,4,23,12,0,14,16,26,25,15,6,18,2,8,7,19,21,20,5,1]
    N=Nodes()
    names=['Ave.','','']
    for o in order:
        names.append(str(N[o].label))
    gcf().set_dpi(400)
    double_column_width = (2*3.425+0.236)*1.5 #extrawide!
    #gcf().set_size_inches([15,5])
    impcop=[0.38016,0,0]
    expcop=[0.38016,0,0]
    imp99=[0.369589,0,0]
    exp99=[0.369589,0,0]
    imp2x=[0.277789,0,0]
    exp2x=[0.277789,0,0]
    impact=[0.129063,0,0]
    expact=[0.129063,0,0]

    width=0.65
    for o in order:
        impcop.append(imports[o,0])
        imp99.append(imports[o,1])
        imp2x.append(imports[o,2])
        impact.append(imports[o,3])
        expcop.append(exports[o,0])
        exp99.append(exports[o,1])
        exp2x.append(exports[o,2])
        expact.append(exports[o,3])
    ind=np.arange(30)
    ax = subplot(111)
    rects1=bar(ind*1.25,expcop,width,align='center',color=(1,1,1))
    rects2=bar(ind*1.25+.1,exp99,width,align='center',color=(.8,.8,.8))
    rects3=bar(ind*1.25+.2,exp2x,width,align='center',color=(.6,.6,.6))
    rects4=bar(ind*1.25+.3,expact,width,align='center',color=(.4,.4,.4))
    pp = (rects1,rects2,rects3,rects4)
    pp_txtlabels = (r'Copper Plate',r'99th Quantile',r'Double of Actual',r'Actual Capacity')
    leg = legend(pp,pp_txtlabels,loc='upper right',ncol=len(pp));
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    gcf().set_size_inches([double_column_width*1.5,3*1.75])
    axis(ymin=0,ymax=0.75,xmin=amin(ind)-.875,xmax=amax(ind)+.875)
    xticks(ind*1.25+.65,names,rotation=60,ha='right',va='top')
    ylabel(r'Export Capacity Ratio')
    savefig('./figures/exports.pdf', dpi=300)
    #show()

def plot_allcases():
    plt.close()
    plota=np.load('./results/PlotA.npy')
    plotb=np.load('./results/PlotB.npy')
    plotc=np.load('./results/PlotC.npy')
#    plotd=load('./results/PlotD.npy')
    plota[:,0]*=74.83
    plotb[:,0]*=74.83
    plotc[:,0]*=74.83
#    plotd[:,0]*=74.83
    plota[:,1]*=100
    plotb[:,1]*=100
    plotc[:,1]*=100
    ax=subplot(1,1,1)
    plt.axis([0, 900, 10, 30.1])
    plt.xlabel('Total Installed Transmission Capacity (GW)')
    plt.ylabel('Balancing Energy Required (% of Annual Consumption)')
    #plt.title(r'Reduction in required balancing energy with increased transmission capacity', fontsize=11)
#    plt.text(-10,24.440991,'X',fontsize=20,color='r')
#    plt.text(20,24.640991,'Individual Nodes')
#    plt.text(839,14.9595,'X',fontsize=20,color='r')
#    plt.text(74.83*9.20,16.31,'\'Copper Plate\'')
#    plt.text(150,15.6,'Lower limit for Agreggated Nodes')
    m=np.min(plotc[:,1])
    p0=ax.plot([0,1500],[m,m],linestyle='dashed',color='k')
#    xticks([])
#    savefig('./figures/allcases.pdf', dpi=300)
#    plt.close()

#    ax=subplot(1,1,1)
#    plt.axis([0, 900.1, 10, 30.1])
#    plt.xlabel('Total Installed Transmission Capacity (GW)')
#    plt.ylabel('Balancing Energy Required (% of Annual Consumption)')
    #plt.title(r'Reduction in required balancing energy with increased transmission capacity', fontsize=11)
#    plt.text(-10,24.440991,'X',fontsize=20,color='r')
#    plt.text(20,24.640991,'Individual Nodes')
#    plt.text(74.83*5.204590,15.2519634,'X',fontsize=20,color='r')
#    plt.text(74.83*3.20-40,14.31,'99th Q -- Sufficiently Large Capacity')
#    plt.text(839,014.9595,'X',fontsize=20,color='r')
#    plt.text(74.83*9.20,16.31,'\'Copper Plate\'')
#    p0=ax.plot([0,1500],[15.2595,15.2595],linestyle='dashed',color='k')
#    savefig('./figures/allcases2.pdf', dpi=300)
#    plt.close()

 #   ax=subplot(1,1,1)
 #   plt.axis([0, 900.1, 10, 30.1])
 #   plt.xlabel('Total Installed Transmission Capacity (GW)')
 #   plt.ylabel('Balancing Energy Required (% of Annual Consumption)')
    #plt.title(r'Reduction in required balancing energy with increased transmission capacity', fontsize=11)
 #   plt.text(-10,24.440991,'X',fontsize=20,color='r')
#    plt.text(20,.24640991,'Individual Nodes')
 #   plt.text(74.83*5.204590,15.2519634,'X',fontsize=20,color='r')
#    plt.text(74.83*3.20-40,.1631,'99th Q -- Sufficiently Large Capacity')
 #   plt.text(840,14.9595,'X',fontsize=20,color='r')
 #   p0=ax.plot([0,1500],[15.2595,15.2595],linestyle='dashed',color='k')
    ax.vlines(74.83,0.1,25,linestyle='dashed',color='k')
#    ax.vlines(74.83*2,0.1,22,linestyle='dashed',color='k')
#    ax.vlines(74.83*3,0.1,20,linestyle='dashed',color='k')
#    ax.vlines(74.83*4,0.1,19,linestyle='dashed',color='k')
#    ax.vlines(74.83*5,0.1,18.5,linestyle='dashed',color='k')
    plt.text(74.83 + 10 ,13.1,'Actual Transmission Capacity')
#    plt.text(74.83*2,22.1,'x2')
#    plt.text(74.83*3,20.1,'x3')
#    plt.text(74.83*4,19.1,'x4')
#    plt.text(74.83*5,18.6,'x5')
#    plote=np.append(plotd,[[900,m]],0)
    p1,=ax.plot(plota[:,0],plota[:,1],label='Increase of Actual Caps.',linewidth=3)
    p2,=ax.plot(plotb[:,0],plotb[:,1],label='Reduction of 99% Q Caps.',linewidth=3)
    p3,=ax.plot(plotc[:,0],plotc[:,1],label='Implementation of X% Q Caps.',linewidth=3)
    handles,labels=ax.get_legend_handles_labels()
    ax.legend(handles,labels)
    savefig('./figures/allcases3.pdf', dpi=300)
    plt.close()

    #show()

def find_balancing_reduction_quantiles(reduction=0.90,eps=1.e-3,guess=0.98):
    '''Loop over different quantile line capacities until the quantile
is found that leads to a reduction of balancing by <reduction>
percent, with a possible relative uncertainty of <eps>. Guess specifies
your first guess for the quantile.'''

    baltarget=.16197916
    olddist=0
    balreal=0
    N=Nodes()
    #N.set_alphas(0.7)
    #N.set_gammas(1.0)
    quant=guess # initial guess
    step=0.01
    print step
    while True:
        h=get_quant(quant)
        N,F=sdcpf(N,h0=h)
        a=0.; b=0.
        for i in N:
            a+=sum(i.balancing)
            b+=i.mean*i.nhours
        balreal=a/b
        reldist=abs(1.-balreal/baltarget)
        dist=baltarget-balreal
        if (reldist < eps):
            print '%15s %15s %15s %15s %15s' % ('distance','old distance','relative dist.','quantile','stepsize')
            print '%15.8f %15.8f %15.4f %15.4f %15.6f' % (dist, olddist, reldist,quant,step)
            del N
            break
        if (dist*olddist<0.): # sign change = we passed the perfect point! now reduce step size
            step=step/2.
        if dist<0:
            quant +=step
        if dist>0:
            quant -=step
        if (quant>=1.):
            step=step/2.
            quant=1.-step
        print '%15s %15s %15s %15s %15s' % ('distance','old distance','relative dist.','quantile','stepsize')
        print '%15.8f %15.8f %15.4f %15.4f %15.6f' % (dist, olddist, reldist,quant,step)
        olddist=dist
    del N
    return quant


def writetext(data,mean,title,lapse=8784):
    x=data[:8784]/mean
    try:
    # This will create a new file or **overwrite an existing file**.
        f = open(title+'.txt', "w")
        try:
            for y in x:
                f.write(str(round(y,8)))
                f.write('\n')
        finally:
            f.close()
    except IOError:
        pass

def get_squant(series,quant):
    a=hist(series,cumulative=True,bins=100,normed=True)
    for i in range(len(a[0])):
        if (a[0][i]>=quant):
            return a[1][i]
            break
        if i == len(a[0])-1:
            return max(a[1][:])

def get_lquant(N,quant=0.99):
    hs=np.zeros(len(N))
    for i in N:
        a=hist(i.load,cumulative=True,bins=100,normed=True)
        for j in range(len(a[0])):
            if (a[0][j]>=quant):
                hs[i.id]=a[1][j]
                break
            if j == len(a[0]) - 1:
                hs[i.id]=max(a[1][:])
    return hs

def get_bquant(N,quant=0.99):
    hs=np.zeros(len(N))
    for i in N:
        a=hist(i.balancing,cumulative=True,bins=100,normed=True)
        for j in range(len(a[0])):
            if (a[0][j]>=quant):
                hs[i.id]=a[1][j]
                break
            if j == len(a[0]) - 1:
                hs[i.id]=max(a[1][:])
    return hs


def bPlot_A(quant=1.0):
    betas=[0.0,0.25,0.5,0.75,1.0,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.0,2.25,2.50,2.75,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.5]
    PlotA=np.zeros((len(betas),2))
    j=0
    for b in betas:
        PlotA[j,0]=b
        N=Nodes(load_filename='Case_A_Beta_'+str(b)+'.npz')
        PlotA[j,1]=sum(get_bquant(N,quant))
        j+=1
        del N
    save('./results/b4PlotA',PlotA)
    return PlotA

def bPlot_C(quant=1.0):
    betas=[0.0,0.05,0.1,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,1.0,1.05,1.1,1.2,1.3,1.4,1.5,2.0,2.5,3.0,4.0,5.0]
    N=Nodes()
    K,k2,Hac,lF=AtoKh(N)
    Hact=biggestpair(Hac)
    Hop=get_quant(.99)
    Hopt=biggestpair(Hop)
    PlotC=np.zeros((len(betas),2))
    j=0
    for b in betas:
        PlotC[j,0]=b*sum(Hopt)/sum(Hact)
        N=Nodes(load_filename='Case_C_Beta_'+str(b)+'.npz')
        PlotC[j,1]=sum(get_bquant(N,quant))
        j+=1
        del N
    save('./results/b4PlotC',PlotC)
    return PlotC

def bPlot_D(quant=1.0):
    quants=[0.60,0.70,0.75,0.80,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.9920,.9940,.9960,.9980,.9999]
    N=Nodes()
    K,k2,Hac,lF=AtoKh(N)
    Hact=biggestpair(Hac)
    PlotD=np.zeros((len(quants),2))
    j=0
    for q in quants:
        Hop=get_quant(q)
        Hopt=biggestpair(Hop)
        PlotD[j,0]=sum(Hopt)/sum(Hact)
        N=Nodes(load_filename='Case_D_Quant_'+str(q)+'.npz')
        PlotD[j,1]=sum(get_bquant(N,quant))
        j+=1
        del N
    save('./results/b4PlotD',PlotD)   
    return PlotD




def plot_ballcases():
    plt.close()
    ###### 100% #######
    plota=load('./results/bPlotA.npy')
    plotc=load('./results/bPlotC.npy')
    plotd=load('./results/bPlotD.npy')
    plota[:,0]*=74.83
    plotc[:,0]*=74.83
    plotd[:,0]*=74.83
    plota[:,1]/=1000
    plotc[:,1]/=1000
    plotd[:,1]/=1000
    ####### 99% #######
    plota2=load('./results/b2PlotA.npy')
    plotc2=load('./results/b2PlotC.npy')
    plotd2=load('./results/b2PlotD.npy')
    plota2[:,0]*=74.83
    plotc2[:,0]*=74.83
    plotd2[:,0]*=74.83
    plota2[:,1]/=1000
    plotc2[:,1]/=1000
    plotd2[:,1]/=1000
    ####### 95% #######
    plota3=load('./results/b3PlotA.npy')
    plotc3=load('./results/b3PlotC.npy')
    plotd3=load('./results/b3PlotD.npy')
    plota3[:,0]*=74.83
    plotc3[:,0]*=74.83
    plotd3[:,0]*=74.83
    plota3[:,1]/=1000
    plotc3[:,1]/=1000
    plotd3[:,1]/=1000
    ####### 90% #######
    plota4=load('./results/b4PlotA.npy')
    plotc4=load('./results/b4PlotC.npy')
    plotd4=load('./results/b4PlotD.npy')
    plota4[:,0]*=74.83
    plotc4[:,0]*=74.83
    plotd4[:,0]*=74.83
    plota4[:,1]/=1000
    plotc4[:,1]/=1000
    plotd4[:,1]/=1000
    ax=subplot(1,1,1)
    plt.axis([0, 700, 200, 250])
    plt.xlabel('Total Installed Transmission Capacity (GW)')
    plt.ylabel('Total Balancing Power Required (GW)')
    #ax.vlines(74.83,0.01,30,linestyle='dashed',color='k')
    #plt.text(74.83*1,13.51,'  Actual Installed Capacity')
    #plote=np.append(plotd,[[900,15.2595]],0)
    #pa,=ax.plot(plota[:,0],plota[:,1],label='Increase of Actual Caps.',linewidth=2)
    #pb,=ax.plot(plotc[:,0],plotc[:,1],label='Reduction of 99% Q Caps.',linewidth=2)
    #pc,=ax.plot(plotd[:,0],plotd[:,1],label='Implementation of X% Q Caps.',linewidth=2)
    #plt.text(50,500,'100% Q Balancing')
    #plt.text(100,365,'99% Q B.')
    #plt.text(200,300,'95% Q B.')
    #plt.text(400,220,'90% Q B.')
    plt.text(200,230,'90% Q Balancing')
    #pa2,=ax.plot(plota2[:,0],plota2[:,1],linewidth=2,linestyle='dashed')
    #pb2,=ax.plot(plotc2[:,0],plotc2[:,1],linewidth=2,linestyle='dashed')
    #pc2,=ax.plot(plotd2[:,0],plotd2[:,1],linewidth=2,linestyle='dashed')
    #pa3,=ax.plot(plota3[:,0],plota3[:,1],linewidth=2,linestyle='dashed')
    #pb3,=ax.plot(plotc3[:,0],plotc3[:,1],linewidth=2,linestyle='dashed')
    #pc3,=ax.plot(plotd3[:,0],plotd3[:,1],linewidth=2,linestyle='dashed')
    pa4,=ax.plot(plota4[:,0],plota4[:,1],linewidth=2,linestyle='dashed')
    pb4,=ax.plot(plotc4[:,0],plotc4[:,1],linewidth=2,linestyle='dashed')
    pc4,=ax.plot(plotd4[:,0],plotd4[:,1],linewidth=2,linestyle='dashed')
    #handles,labels=ax.get_legend_handles_labels()
    #ax.legend(handles,labels)
    savefig('./figures/b2allcases.pdf', dpi=300)
    plt.close()



def drawg(N,F,t,lod,wind,solar,gmis,framenumber,D=False):
    close()
    G=nx.DiGraph()
    nodelist=[]
    for n in N:
        G.add_node(str(n.label))
        nodelist.append(str(n.label))
    #print G.nodes()
    fig = figure(figsize=(9,6))
    a,b,ListF,G,F=AtoKh_nonsparse(N,F[:,t],G)
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

    ax1= fig.add_axes([0.875,0.05,0.05,0.9])  #Drawing Color bar on the right
    cmap=mpl.cm.RdYlGn
    norm=mpl.colors.Normalize(vmin=-1,vmax=1)
    cbl=mpl.colorbar.ColorbarBase(ax1,cmap,norm,orientation= 'vertical' )
    cbl.add_lines([gmis],colors='k',linewidths='1.5')


    ax2= fig.add_axes([0.05, 0.175, 0.4, 0.125]) # Writing date above time series
    ax2.text(0.05, 0.2, ttodate(int(t)), fontsize=17)
    ax2.axis('off')

    ax3= fig.add_axes([0.05,0.05,0.8,0.125]) # plotting time series
    ax3.plot(lod,color='r',linewidth=1.5)
    ax3.plot(wind,color='b',linewidth=1.5)
    ax3.plot(solar,color='y',linewidth=1.5)
    #ax3.set_xlabel('Time of day (h)')
    ax3.set_xticks(arange(0,168,24))
    ax3.set_yticks([0.5,1.0])
    ax3.grid(axis='x',linestyle='--')
    ax3.axis([0,168,0,1.5])
    

    ax4= fig.add_axes([-0.075,0.10,.975,0.95]) #For displaying graph
    for l in ListF:
        if -200<F[l[2]]<200: continue      # if low flow, no arrow
        if F[l[2]]>=200:                    # assing positions for positive flow
            x0=pos[l[0]][0]
            y0=pos[l[0]][1]
            x1=pos[l[1]][0]
            y1=pos[l[1]][1]
        if F[l[2]]<=-200:                   # reverse positions for negative flow
            x1=pos[l[0]][0]
            y1=pos[l[0]][1]
            x0=pos[l[1]][0]
            y0=pos[l[1]][1]
        dist=sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1))  # calculate distance, for heeadlength
        ax4.arrow(x0,y0,.5*(x1-x0),.5*(y1-y0),fc='k',ec='k',head_width=min(0.03,0.02*(abs(F[l[2]]))/3000.0+0.01),head_length=0.1*dist)              # actually draw the arrows


    node_c=[]
    #ax4.arrow(0.15,0.35,0.0,-0.2,fc='k',ec='k',head_width=0.05,head_length=0.1)
    for n in N:
        node_c.append(cmap(  0.5 +  0.5*n.mismatch[t]/n.mean    ))
    #print G.nodes()    
    nx.draw_networkx_nodes(G,pos,node_size=1200,nodelist=nodelist,node_color=node_c,facecolor=(1,1,1))
        
    e0=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']<=200]
    e1=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>200 and d['weight']<=400]
    e2=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>400 and d['weight']<=600]
    e3=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>600 and d['weight']<=800]
    e4=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>800 and d['weight']<=1000]
    e5=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1000 and d['weight']<=1250]
    e6=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1250 and d['weight']<=1500]
    e7=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1500 and d['weight']<=1750]
    e8=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>1750 and d['weight']<=2000]
    e9=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>2000 and d['weight']<=2250]
    e10=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>2250 and d['weight']<=2500]
    e11=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>2500 and d['weight']<=2750]
    e12=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>2750 and d['weight']<=3000]
    e13=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>3000 and d['weight']<=3250]
    e14=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>3250 and d['weight']<=3500]
    e15=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>3500 and d['weight']<=3750]
    e16=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight']>3750]
    
    nx.draw_networkx_edges(G,pos,edgelist=e0,width=1.00,edgecolor=(0.2,0.2,0.2),arrows=0,style ='dashed')
    nx.draw_networkx_edges(G,pos,edgelist=e1,width=1.25,edgecolor=(0.2,0.2,0.2),arrows=D)
    nx.draw_networkx_edges(G,pos,edgelist=e2,width=1.50,edgecolor=(0.2,0.2,0.2),arrows=D)
    nx.draw_networkx_edges(G,pos,edgelist=e3,width=1.75,edgecolor=(0.2,0.2,0.2),arrows=D)
    nx.draw_networkx_edges(G,pos,edgelist=e4,width=2.00,edgecolor=(0.2,0.2,0.2),arrows=D)
    nx.draw_networkx_edges(G,pos,edgelist=e5,width=2.25,edgecolor=(0.2,0.2,0.2),arrows=D)
    nx.draw_networkx_edges(G,pos,edgelist=e6,width=2.50,edgecolor=(0.2,0.2,0.2),arrows=D)
    nx.draw_networkx_edges(G,pos,edgelist=e7,width=2.75,edgecolor=(0.2,0.2,0.2),arrows=D)
    nx.draw_networkx_edges(G,pos,edgelist=e8,width=3.00,edgecolor=(0.2,0.2,0.2),arrows=D)
    nx.draw_networkx_edges(G,pos,edgelist=e9,width=3.25,edgecolor=(0.2,0.2,0.2),arrows=D)
    nx.draw_networkx_edges(G,pos,edgelist=e10,width=3.50,edgecolor=(0.2,0.2,0.2),arrows=D)
    nx.draw_networkx_edges(G,pos,edgelist=e11,width=3.75,edgecolor=(0.2,0.2,0.2),arrows=D)
    nx.draw_networkx_edges(G,pos,edgelist=e12,width=4.00,edgecolor=(0.2,0.2,0.2),arrows=D)
    nx.draw_networkx_edges(G,pos,edgelist=e13,width=4.25,edgecolor=(0.15,0.15,0.15),arrows=D)
    nx.draw_networkx_edges(G,pos,edgelist=e14,width=4.50,edgecolor=(0.1,0.1,0.1),arrows=D)
    nx.draw_networkx_edges(G,pos,edgelist=e15,width=4.75,edgecolor=(0.05,0.05,0.05),arrows=D)
    nx.draw_networkx_edges(G,pos,edgelist=e16,width=5.00,arrows=D)
#1
    #nx.draw_networkx_edges(G,pos,edgelist=emid,width=2)#3
    #nx.draw_networkx_edges(G,pos,edgelist=elarge,width=2)#6
    nx.draw_networkx_labels(G,pos,font_size=20,font_color='w',font_family='sans-serif')
    axis('off')
    savefig("./movie/frame%06.i.png" %framenumber) # save as png
    print "Generated frame %06.i" %framenumber
    #show() # display

def makeweek(tinit):   ## identify closest monday 00.00
    N=Nodes(load_filename="Gurobi.npz")    
    F=load('./results/gurobiflows.npy')
    t0=tinit/168*168 + 48
    print 'Picking closest Moday 00.00, t0 = ', t0
    lapse=arange(t0,t0+337,1)
    lod=[]
    wind=[]
    solar=[]
    framenumber=0
    for t in lapse:
        lnew=sum(n.load[t] for n in N)/sum(n.mean for n in N)
        wnew=sum(n.get_wind()[t] for n in N)/sum(n.mean for n in N)
        snew=sum(n.get_solar()[t] for n in N)/sum(n.mean for n in N)
        lod.append(lnew)
        wind.append(wnew)
        solar.append(snew)
        if len(lod)>120 and max(lapse)-t0-framenumber>48:
            for i in range(120):
                lod[i]=lod[i+1]
                wind[i]=wind[i+1]
                solar[i]=solar[i+1]
            del(lod[-1])
            del(wind[-1])
            del(solar[-1])
        gmis=sum(n.mismatch[t] for n in N)/sum(n.mean for n in N)
        if gmis >1:gmis=1
        if gmis <-1: gmis=-1
        drawg(N,F,t,lod,wind,solar,gmis,framenumber)
        framenumber+=1
    
def Case_Q(quants=[0.75,0.8,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.9920,.9940,.9960,.9980,.9999],b=1):
    hours=[3,6,12,24]
    N=Nodes()
    for i in hours:
        N.set_ESS(1.0,i)
        for q in quants:
            h=get_quant(q)
    #        if sum(h) >=1:
    #            for i in range(len(h)):
    #                if h[i]==0: h[i]=50
    #            print h
            M,F=sdcpf(N,b=b,h0=h)
            M.save_nodes(str(i)+'_Hours_'+str(q)+"_Quant")
            save('./results/'+str(i)+'_Hours_'+str(q)+'_Quant_Flows',F)

def Plot_Q(hour):
    quants=[0.5,0.75,0.8,0.85,0.86,0.87,0.88,0.89,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.9920,.9940,.9960,.9980,.9999]
    N=Nodes()
    K0, k1, k2, Hac, lF=AtoKh(N)
    Hact=biggestpair(Hac)
    PlotQ=np.zeros((len(quants),2))
    j=0
    for q in quants:
        Hop=get_quant(q)
        Hopt=biggestpair(Hop)
        PlotQ[j,0]=sum(Hopt)/sum(Hact)
        N=Nodes(load_filename=str(hour)+'_Hours_'+str(q)+'_Quant.npz')
        a=0
        d=0
        for i in N:
            a+=sum(i.balancing)
            d+=i.mean*i.nhours
        c=a/d
        PlotQ[j,1]=c
        j+=1
    save('./results/PlotQ_'+str(hour)+'_Hours',PlotQ)

# A after beta=15 is bad    
def Plot_allQases():
    plota=np.load('./results/PlotQ_3_Hours.npy')
    plotb=np.load('./results/PlotQ_6_Hours.npy')
    plotc=np.load('./results/PlotQ_12_Hours.npy')
    plotd=np.load('./results/PlotQ_24_Hours.npy')
    plote=np.load('./results/PlotD.npy')
    plota[:,0]*=74.83
    plotb[:,0]*=74.83
    plotc[:,0]*=74.83
    plotd[:,0]*=74.83
    plote[:,0]*=74.83
    plota[:,1]*=1
    plotb[:,1]*=1
    plotc[:,1]*=1
    plotd[:,1]*=1
    plote[:,1]*=1
    print min(plota[:,1])
    print min(plotb[:,1])
    print min(plotc[:,1])
    print min(plotd[:,1])
    print min(plote[:,1])
    ax=subplot(1,1,1)
    p1=ax.plot([0,1500],[min(plota[:,1]),min(plota[:,1])],linestyle='dashed',color='k')
    p2=ax.plot([0,1500],[min(plotb[:,1]),min(plotb[:,1])],linestyle='dashed',color='k')
    p3=ax.plot([0,1500],[min(plotc[:,1]),min(plotc[:,1])],linestyle='dashed',color='k')
    p4=ax.plot([0,1500],[min(plotd[:,1]),min(plotd[:,1])],linestyle='dashed',color='k')
    p5=ax.plot([0,1500],[min(plote[:,1]),min(plote[:,1])],linestyle='dashed',color='k')
    plt.axis([0, 700, 0.0, 0.26])
    plt.xlabel('Total Installed Transmission Capacity (GW)')
    plt.ylabel('Total Balancing Power Required (GW)')
    ax.vlines(74.83,0.0,0.27,linestyle='dashed',color='k')
    plt.text(74.83*1,.1351,'  Actual Installed Capacity')
    ax.plot(plota[:,0],plota[:,1])
    ax.plot(plotb[:,0],plotb[:,1])
    ax.plot(plotc[:,0],plotc[:,1])
    ax.plot(plotd[:,0],plotd[:,1])
    ax.plot(plote[:,0],plote[:,1])
    savefig('./figures/qases.pdf', dpi=300)
#Case_A()
#Case_B()
def get_optimal_mix_balancing(L, GW, GS, gamma=1., returnall=False, normalized=True):
    L, GW, GS = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2) #Ensure minimum dimension to 2 to alow the weighed sum to be calculated correctly.
    weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))

    l = weighed_sum(L)
    Gw = weighed_sum(GW)	
    Gs = weighed_sum(GS)

    mismatch = lambda alpha_w: gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - l
    res_load_sum = lambda alpha_w: sum(get_positive(-mismatch(alpha_w)))

    alpha_w_opt = fmin(res_load_sum,0.5,disp=False)
    res_load_sum_1p_interval = lambda alpha_w: res_load_sum(alpha_w)-(res_load_sum(alpha_w_opt)+.01*sum(l))

    alpha_w_opt_1p_interval = array([brentq(res_load_sum_1p_interval, 0, alpha_w_opt),brentq(res_load_sum_1p_interval, alpha_w_opt, 1)])

    if normalized:
        mismatch_opt = mismatch(alpha_w_opt)
    else:
        mismatch_opt = mismatch(alpha_w_opt)*mean(sum(L,axis=0))
    res_load_sum_opt = sum(get_positive(-mismatch_opt))

    if returnall:
    #Returns: alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt
        return alpha_w_opt, alpha_w_opt_1p_interval, res_load_sum_opt, mismatch_opt
    else:
        return alpha_w_opt
#Plot_allQases()
def optbars(): 
    order = [24,11,10,22,3,9,17,13,4,23,12,0,14,16,26,25,15,6,18,2,8,7,19,21,20,5,1]
    alphas=[]
    alphas_inter=[]
    N=Nodes()
    names=['EU','']
    for o in order:
        names.append(str(N[o].label))

    #gcf().set_size_inches([15,5])
    width=0.95
    EU_load = sum(n.load for n in N)
    EU_mean = sum(n.mean for n in N)
    EU_load/= EU_mean
    EU_wind = sum(n.normwind for n in N)/27.0
    EU_solar= sum(n.normsolar for n in N)/27.0
    alpha,alpha_inter,wtf1,wtf2=get_optimal_mix_balancing(EU_load,EU_wind,EU_solar,returnall=True)
    alphas.append(alpha)
    alphas.append(0.000001)
    alphas_inter.append(alpha_inter)
    alphas_inter.append(array([0.0000001,0.00001]))
    for o in order:
        alpha,alpha_inter,wtf1,wtf2=get_optimal_mix_balancing(N[o].load/N[o].mean,N[o].normwind,N[o].normsolar,returnall=True)
        alphas.append(alpha)
        alphas_inter.append(alpha_inter)
    ind=np.arange(29)
    ax = subplot(1,1,1)
    rects=bar(ind*1.25,alphas,width=width,align='center',color=(.3,.3,.3))
    axis(ymin=0,ymax=1.0,xmin=amin(ind)-.875,xmax=amax(ind*1.25)+.875)
    gcf().set_size_inches([ 3*colwidth , 3*colwidth*0.50])#[9*1.5,3*1.75])
    gcf().set_dpi(400)
    p1=ax.plot([-2,36],[mean(alphas[2:]),mean(alphas[2:])],linestyle='dashed',color='k')


    errorbar(0,alphas[0],yerr=[ [float(alphas[0]-alphas_inter[0][0]) , float(alphas_inter[0][1]-alphas[0]) ]],color='k',linewidth=1.5)
    for i in range(2,29):
        errorbar(i*1.25,alphas[i],yerr=float(alphas[i]-alphas_inter[i][0]),color='k',linewidth=1.5)
        #errorbar(i*1.25,alphas[i],yerr=[[ float(alphas[i]-alphas_inter[i][0]) , float(alphas_inter[i][1]-alphas[i]) ]],color='k')

    yticks([0.0,0.25,0.5,0.75,1.0])
    xticks(ind*1.25+.35,names,rotation=60,ha='right',va='top')
    ylabel(r'Wind fraction ($\alpha_W^X$)')
    plt.tight_layout()
    savefig('./figures/optbar.pdf', dpi=300)
    plt.close()
    #show()
def balbars(): 
    order = [24,11,10,22,3,9,17,13,4,23,12,0,14,16,26,25,15,6,18,2,8,7,19,21,20,5,1]
    load_today=[]
    load_99Q=[]
    load_NoC=[]
    N=Nodes(load_filename="TodayCap.npz")
    M=Nodes(load_filename="99QCap.npz")
    O=Nodes(load_filename="NoCap.npz")
    names=['EU','']
    for o in order:
        names.append(str(N[o].label))

    #gcf().set_size_inches([15,5])
    width=0.95
    EU_delta_TOD=sum(n.balancing for n in N)
    EU_delta_99Q=sum(m.balancing for m in M)
    EU_delta_NoC=sum(o.balancing for o in O)
    EU_load=sum(n.nhours*n.mean for n in N)
    #load_today.append(sum(get_positive(EU_delta_TOD))/EU_load)
    #load_99Q.append(sum(get_positive(EU_delta_99Q))/EU_load)
    #load_NoC.append(sum(get_positive(EU_delta_NoC))/EU_load)
    load_today.append(sum(EU_delta_TOD)/EU_load)
    load_99Q.append(sum(EU_delta_99Q)/EU_load)
    load_NoC.append(sum(EU_delta_NoC)/EU_load)
    print load_today,load_99Q,load_NoC
    load_today.append(0)
    load_99Q.append(0)
    load_NoC.append(0)
    for o in order:
        load_today.append(sum(get_positive(N[o].balancing))/(N[o].nhours*N[o].mean))
        load_99Q.append(sum(get_positive(M[o].balancing))/(N[o].nhours*N[o].mean))
        load_NoC.append(sum(get_positive(O[o].balancing))/(N[o].nhours*N[o].mean))
    ind=np.arange(29)
    ax = subplot(1,1,1)
    rects0=bar(ind*1.25,load_NoC,width=width,align='center',color=(.85,.85,.85),label="No Transmission")
    rects1=bar(ind*1.25,load_today,width=width,align='center',color=(.50,.50,.50),label="Actual Capacities")
    rects2=bar(ind*1.25,load_99Q,width=width,align='center',color=(.25,.25,.25),label="$99^{\mathrm{th}}$ Q Capacities")
    axis(ymin=0,ymax=.35,xmin=amin(ind)-.875,xmax=amax(ind*1.25)+.875)
    gcf().set_size_inches([ 3*colwidth , 3*colwidth*0.50])#[9*1.5,3*1.75])
    gcf().set_dpi(400)
    yticks([0.05,0.15,0.25,0.35])
    xticks(ind*1.25+.35,names,rotation=60,ha='right',va='top')
    ylabel(r'Residual Load ($\sum L$)')
    pp = (rects0,rects1,rects2)
    pp_txtlabels = (r'No Transmission',r'Actual Capacities',r'$99^{\mathrm{th}}$ Q Capacities')
    leg = legend(pp,pp_txtlabels,loc='upper right',ncol=len(pp));
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize
    plt.tight_layout()
    savefig('./figures/balbar.pdf', dpi=300)
    plt.close()
#Case_C(betas=[0.93,0.94,0.95,1.0,1.05,1.1,1.2,1.3,1.4,1.5,2.0,2.5,3.0,4.0,5.0])
#Case_D(quants=[0.60,0.70,0.75,0.80,0.85,0.86,0.87,0.88,0.89])
#Plot_A()
#Plot_B()
#Plot_C()
#Plot_D()
def findneighbours(i,matr):
    """ matr should be = np.genfromtxt("./settings/admat.txt")"""
    neighbours=[]
    for j in range(len(matr[i])):
        if matr[i][j]>0 and i<>j:
            neighbours.append(j)
    return neighbours

def findlinks(n,matr):
    """ matr should be = np.genfromtxt("./settings/admat.txt")"""
    links=[]
    linkno=-1
    for j in range(len(matr[n])):
        for i in range(len(matr[n])):
            if matr[i][j]>0 and i>j:
                linkno+=1
                if n==i: links.append(-linkno)
                if n==j: links.append(linkno)
    return links

def scatterneigh(m,N,matr):
    close()
    neigh=findneighbours(m,matr)
    ymism=sum(N[i].mismatch*1 for i in neigh)
    ymism/=max(ymism)
    xmism=N[m].mismatch*1
    xmism/=max(xmism)
    xbin=int((max(xmism)-min(xmism))/0.05)
    ybin=int((max(ymism)-min(ymism))/0.05)
    #scatter(xmism,ymism)
    hexbin(xmism,ymism,gridsize=(xbin,ybin),mincnt=1,cmap=matplotlib.cm.OrRd,linewidths=1)
    xlabel(r'Mismatch: '+str(N[m].label))
    ylabel(r'Mismatch: '+str(N[m].label)+"'s Neighbours")
    savefig("./figures/scatter_type_1_"+str(N[m].label)+".pdf")

def scatterEXPneigh(m,N,matr,F):
    close()
    links=findlinks(m,matr)
    ymism=sum(sign(l)*F[abs(l)] for l in links)
    ymism/=max(ymism)
    xmism=N[m].curtailment-N[m].balancing
    xmism/=max(xmism)
    xbin=int((max(xmism)-min(xmism))/0.05)
    ybin=int((max(ymism)-min(ymism))/0.05)
    scatter(xmism,ymism)
    #hexbin(xmism,ymism,gridsize=(xbin,ybin),mincnt=1,cmap=matplotlib.cm.OrRd,linewidths=1)
    xlabel(r'Residual Load: '+str(N[m].label))
    ylabel(str(N[m].label)+"'s Net exports")
    savefig("./figures/scatter_type_3_"+str(N[m].label)+".png")

def scatterEU(m,N,matr):
    close()
    #neigh=findneighbours(m,matr)
    ymism=sum(i.mismatch*1 for i in N)-N[m].mismatch
    ymism/=max(ymism)
    xmism=N[m].mismatch*1
    xmism/=max(xmism)
    xbin=int((max(xmism)-min(xmism))/0.05)
    ybin=int((max(ymism)-min(ymism))/0.05)
    #scatter(xmism,ymism)
    hexbin(xmism,ymism,gridsize=(xbin,ybin),mincnt=1,cmap=matplotlib.cm.OrRd,linewidths=1)
    xlabel(r'Mismatch: '+str(N[m].label))
    ylabel(r'Mismatch: EU - '+str(N[m].label))
    savefig("./figures/scatter_type_2_"+str(N[m].label)+".pdf")

interest=[4,7,10,18,21,24]

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
                hs[i]=max(a[1][:])
    return hs




