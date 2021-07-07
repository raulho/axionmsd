#--------------------------------------------------------------
# Code for getting CMB spectral μ-distortions on
# Axion Monodromy model 
# Autor: Raúl Henríquez
# date: 7/07/2021
#--------------------------------------------------------------




import numpy as np
from int_inflation import *
from functions_sd import *
import scipy.constants as c



#-------------------------------------------------#


N = 6000


pf1,pf2=-0.75,1.0

f=5e-03
pf=np.linspace(pf1,pf2,N)


Dphi=0

ns_limit=0.05

p1,p2,p3=2.0/3.0,1.0,4.0/3.0


#-------------------------------------------#
doc = open('mu_p23_pf_f'+str(f)+'.dat', 'w')

for i in range(N):
         if (dns(f,pf[i],p1,alpha,N0)) <= ns_limit and (dns(f,pf[i],p1,alpha,N0) >= 0):
             mu=get_mu_sd2(f,pf[i],p1,ns,Dphi,int_func,alpha,phik,N0)
             doc.write("%.9f %.12f %.14f\n" % (pf[i],dns(f,pf[i],p1,alpha,N0),mu))

doc.close()

#-------------------------------------------#

doc = open('mu_p1_pf_f'+str(f)+'.dat', 'w')

for i in range(N):
         if (dns(f,pf[i],p2,alpha,N0)) <= ns_limit and (dns(f,pf[i],p2,alpha,N0) >= 0):
             mu=get_mu_sd2(f,pf[i],p2,ns,Dphi,int_func,alpha,phik,N0)
             doc.write("%.9f %.12f %.14f\n" % (pf[i],dns(f,pf[i],p2,alpha,N0),mu))
             
doc.close()


#-------------------------------------------# 

doc = open('mu_p43_pf_f'+str(f)+'.dat', 'w')

for i in range(N):
         if (dns(f,pf[i],p3,alpha,N0)) <= ns_limit and (dns(f,pf[i],p3,alpha,N0) >= 0):
             mu=get_mu_sd2(f,pf[i],p3,ns,Dphi,int_func,alpha,phik,N0)
             doc.write("%.9f %.12f %.14f\n" % (pf[i],dns(f,pf[i],p3,alpha,N0),mu))
             
doc.close()


#-------------------------------------------#  
