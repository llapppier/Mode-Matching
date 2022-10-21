####################İmporting Needed Libraries############################################
import numpy as np; # numpy library for needed mathematical functionality
import matplotlib.pyplot as plt# pyplot library for praphing operations
import pandas as pd
from mmt import *
############################################################################################
############################################################################################
giga=10**9;
mega=10**6;
mili=10**(-3);
fmin=10*giga;
fmax=20*giga;
mode_trnc=20;
s_size=100;
################################################################################
Wmean=300*mili;
Wvar = 10*mili;
Lmean=50*mili;
Lvar = 10*mili;
#np.random.seed(20)
""" w=np.random.normal(Wmean,Wvar,20);
l=np.random.normal(Lmean,Lvar,20); """
####################################################################################
w=mili*np.array([10.7,8.7]);
l=mili*np.array([1,1]);
###################################################################################
"""w=np.array([15.7,25,7.44])*mili;
l=np.array([50,50,50])*mili;"""
###################################################################################
k=sparam(fmin,fmax,mode_trnc,w,l,s_size);
f=np.linspace(fmin,fmax,s_size);
"""cst=pd.read_csv("aa.txt",sep="\t")
cst=np.array(cst)
plt.plot(f,cst[:,1],"r")"""
plt.plot(f,np.abs(k[:,0,0]))
plt.plot(f,np.abs(k[:,20,0]))
plt.show()
###############################################################################
