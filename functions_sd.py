import numpy as np
from int_sd import *
from scipy import integrate

#-------------Model DS------------------#
def Wmu(k):
    return 2.8*D2*(np.exp(-np.power(k/1360,2)/(1+np.power(k/260,0.3)+k/340))-np.exp(-np.power(k/32,2)))

def Wy(k):
    return D2*np.exp(-np.power(k/32.0,2))/2.0

def Diy(x,y):
    return I0*(np.power(x,4.0)*np.exp(x)/np.power(np.exp(x)-1.0,2.0))*(x*(1/np.tanh(x/2.0))-4.0)*y

def Dimu(x,mu):
    return I0*(np.power(x,4.0)*np.exp(x)/np.power(np.exp(x)-1.0,2.0))*(1.0/beta - 1.0/x)*mu

#------------Model Monodromy-------------------#

def N0(p):
    return Ns + np.power(phied,2)/(2*p)

def phik(k,p,N0):
    return np.sqrt(2*p*(N0(p)-np.log(k/k0)))     # Mpl

def alpha(f,pf,p,N0):
    return (1+pf)*(phi0/(2*f*N0(p)))*np.power(np.sqrt(2*p*N0(p))/phi0,1+pf)

def dns(f,pf,p,alpha,N0):
    return 3*b*np.power(2*np.pi/alpha(f,pf,p,N0),0.5)

def Pkm(k,f,pf,p,ns,Dphi,alpha,phik,N0):
    return ((2*np.power(np.pi,2))/np.power(k,3))*A*np.power(k/k0,ns-1)*(1+dns(f,pf,p,alpha,N0)*np.cos((phi0/f)*np.power(phik(k,p,N0)/phi0,pf+1)+Dphi))
    

#-----------Integration--------------------#

def int_func(x,f,pf,p,ns,Dphi,alpha,phik,N0):
    return (np.power(x,2)/(2.0*np.power(np.pi,2)))*Wmu(x)*Pkm(x,f,pf,p,ns,Dphi,alpha,phik,N0)
    
def get_mu_sd(f,pf,ns,Dphi,int_func,alpha,phik):
    return integrate.quad(int_func,1,np.infty,args=(f,pf,ns,Dphi,alpha,phik),limit=1000)[0]
     
def get_mu_sd2(f,pf,p,ns,Dphi,int_func,alpha,phik,N0):     
    return integrate.quad(lambda x: int_func(x,f,pf,p,ns,Dphi,alpha,phik,N0),1,np.infty,epsabs=1e-30, epsrel=1e-8,limit=300)[0]

#----------------------------#


