#----------------------------Python Code for Mode Matching Technique Applied at Cascaded Discontinuities----------

############################Structure of a Single Discontinuity#########################################
"""
<a1/2>-------------------|
                         ------------------------<-a0/2>
      

      REGION II            REGION I


                          -----------------------<-a0/2>
<-a1/2>-------------------|
"""
####################İmporting Needed Libraries############################################
import numpy as np; # numpy library for needed mathematical functionality
import matplotlib.pyplot as plt# pyplot library for praphing operations
import pandas as pd
###################Some Necessary Constants##################################################
pi=np.pi;
eps_r=1;
eps=eps_r*(8.85418782*(10**(-12)));
mu=1.25663706*(10**(-6));

############################################################################################
def wImp(f,mode,width):#function to calculate wave impedance
    omg=2*pi*f;
    return (omg*mu)/(propConst(f,mode,width));
def wNumber(f):#function to calculate wave number
    omg=2*pi*f;
    return (omg**2)*mu*eps;
def propConst(f,mode,width): #function to calculate propagation constant
    """
    f    :frequency
    mode :propagation mode
    width:width of waveguide along discontinuity axis 
    """
    kfree=wNumber(f);
    kguide=(mode*pi/width)**2;
    j=complex(0,1);
    if kfree > kguide:
        return np.sqrt(np.abs(kfree-kguide))
    else:
        return -j*np.sqrt(np.abs(kguide-kfree))
#############################################################################################################
#############################FUNCTION TO CALCULATE SCATTERING MATRIX FOR A DISCONTINUITY#########################################
#############################################################################################################
def s_matrix(fmin,fmax,mode_trnc,w0,w1,s_size):
    #[fmin:fmax];frequency range
    #mode_trnc  ;after which mode to truncate
    #[w0,w1]    ;widt of the first and second region of the waveguide respectively
    #s_size     ;how many frequency sample to include
    #global S11,S12,S21,S22,LE,scat #declared these matrices global just for debugging
    LE=np.empty((s_size,mode_trnc,mode_trnc),dtype=complex);#placeholder for LE matrix
    I=np.eye(mode_trnc) # identity matrix
    #Allocating empty matrices for S parameters
    S11=np.empty((s_size,mode_trnc,mode_trnc),dtype=complex);
    S12=np.empty((s_size,mode_trnc,mode_trnc),dtype=complex);
    S21=np.empty((s_size,mode_trnc,mode_trnc),dtype=complex);
    S22=np.empty((s_size,mode_trnc,mode_trnc),dtype=complex);
    scat=np.empty((s_size,2*mode_trnc,2*mode_trnc),dtype=complex);
   #---------------------------------------------------------------------------------------------------
    fStep=((fmax-fmin)/s_size);#corresponding frequency step for each indices of array
    k=np.empty((s_size,2,mode_trnc),dtype=complex);#array to be populated with propogation constants
    for i in range(s_size):
        """
        this loop will result an array with shape of k[s_size,2,mode_trnc]
        first indice of this array will denote certain frequency sample(each indices will correspond to certain freq.)
        second indice will denote witch region of the waveguide this k number wil be belong to;
            such as 0 for region I
                    1 for region II 
        
        third indice will specify the propogation mod
        """
        f=fmin+i*fStep;
        for j in range(mode_trnc):
            k[i,0,j]=propConst(f,j+1,w0);
            k[i,1,j]=propConst(f,j+1,w1);
#---------------------------------------------------------------------------------------------
    if w1>w0:
        for q in range(s_size):#loop for frequencies
            for i in range(mode_trnc):#loop for mods at second region("m")
                for j in range(mode_trnc):#loop for mods at first region("n")
                    km=k[q,1,i];
                    kn=k[q,0,j];
                    c0=2*(np.sqrt(km/(w0*w1*kn)));
                    c1=(j+1)*pi/w0;
                    c2=(i+1)*pi/w1;
                    c3=(i+1)*pi*(w0+w1)/(2*w1);
                    c4=(i+1)*pi*(w1-w0)/(2*w1);
                    c5=((c2**2)-(c1**2));
                    LE[q,i,j]=c0*(c1*(((-1)**(j+1))*np.sin(c3)-np.sin(c4))/c5)
            
            l_e=LE[q,:,:];
            l_h=np.transpose(l_e);
            temp=(np.linalg.inv(I+l_h@l_e))@(I-l_h@l_e);
            
            S11[q,:,:]=temp;
            S12[q,:,:]=(2*np.linalg.inv(I+l_h@l_e))@l_h;
            S21[q,:,:]=l_e@(I+S11[q,:,:]);
            S22[q,:,:]=(l_e@S12[q,:,:])-I;
    else:
        for q in range(s_size):#loop for frequencies
            for j in range(mode_trnc):#loop for mods for first region("n")
                for i in range(mode_trnc):#loop for mods for second region("m")
                    km=k[q,1,i];
                    kn=k[q,0,j];
                    c0=2*(np.sqrt(kn/(w0*w1*km)));
                    c1=(i+1)*pi/w1;
                    c2=(j+1)*pi/w0;
                    c3=(j+1)*pi*(w0-w1)/(2*w0);
                    c4=(j+1)*pi*(w0+w1)/(2*w0);
                    LE[q,j,i]=c0*(c1*((np.sin(c3)-((-1)**(i+1))*np.sin(c4))/((c1**2)-(c2**2))));
            l_e=LE[q,:,:];
            l_h=np.transpose(l_e);
            temp=(np.linalg.inv(l_e@l_h+I))@(l_e@l_h-I);
            S11[q,:,:]=temp;
            S12[q,:,:]=2*((np.linalg.inv(l_e@l_h+I))@l_e)
            S21[q,:,:]=l_h@(I-S11[q,:,:])
            S22[q,:,:]=I-l_h@S12[q,:,:];
# Unifying s parameters under  scattering matrix
    scat[:,0:mode_trnc, 0:mode_trnc]=S11[:,:,:];
    scat[:,0:mode_trnc, mode_trnc:2*mode_trnc ]=S12[:,:,:];
    scat[:,  mode_trnc:2*mode_trnc,0:mode_trnc]=S21[:,:,:];
    scat[:,mode_trnc:2*mode_trnc, mode_trnc:2*mode_trnc]=S22[:,:,:];
    return scat;
################################END OF FUNCTION#################################################
###################FUNCTION TO CALCULATE SCATTERING MATRIX FOR A HOMOGENIOUS REGION#################################
def homoS(fmin,fmax,mode_trnc,width,length,s_size):
    #width: width of the waveguide
    #length:length of homogenious region 
    fStep=((fmax-fmin)/s_size);#corresponding frequency step for each indices of array
    scat=np.zeros((s_size,2*mode_trnc,2*mode_trnc),dtype=complex);
    D=np.zeros((s_size,mode_trnc,mode_trnc),dtype=complex);
    k=np.empty((s_size,mode_trnc),dtype=complex);
    for z in range(s_size):
        f=fmin+z*fStep;
        for i in range(mode_trnc):
            k[z,i]=propConst(f,i+1,width);
    for z in range(s_size):
        for i in range(mode_trnc):
            D[z,i,i]=np.exp((-1j)*k[z,i]*length);
        scat[z,0:mode_trnc,mode_trnc:]=D[z,:,:];
        scat[z,mode_trnc:,0:mode_trnc]=D[z,:,:];
    return scat;
############################################################################################################
##################################S MATRIX OF TWO CASCADED SECTION ##################################################
def cascade(sL,sR):
    global sT,sR11,sR12,sR21,sR22,h;
    h=int(np.shape(sL)[1]/2);#for choosing halfpoints of scattering matrices
    I=np.eye(h) # identity matrix
    
    sT=np.empty(np.shape(sL),dtype=complex);#S parameters for combined structure
    sR11=sR[:,:h,:h];#S11 for Right waveguide
    sR12=sR[:,:h,h:];#S12 for Right waveguide
    sR21=sR[:,h:,:h];#S21 for Right waveguide
    sR22=sR[:,h:,h:];#S22 for Right waveguide
    sL11=sL[:,:h,:h];#S11 for left waveguide
    sL12=sL[:,:h,h:];#S12 for left waveguide
    sL21=sL[:,h:,:h];#S21 for left waveguide
    sL22=sL[:,h:,h:];#S22 for left waveguide
    sT[:,:h,:h]=sR11+sR12@(np.linalg.inv(I-sL11@sR22))@sL11@sR21;#sT11
    sT[:,:h,h:]=sR12@(np.linalg.inv(I-sL11@sR22))@sL12;          #sT12
    sT[:,h:,:h]=sL21@(np.linalg.inv(I-sR22@sL11))@sR21;          #sT21
    sT[:,h:,h:]=sL22+sL21@(np.linalg.inv(I-sR22@sL11))@sR22@sL12;#sT22
    return sT;


##################################CALCULATING S MATRIX OF A GIVEN STRUCTURE##################################################
def sparam(fmin,fmax,mode_trnc,w,l,s_size):
    """
    This function will give final s parameters for whole structure so user only needs to use this function

    [fmin:fmax] :frequency range
    mode_trnc   :after which mode to truncate
    s_size      :how many frequency sample to include
    w           :array that holds width value of each vaveguide section
    l           :array that holds length value of each vaveguide section
    """
    Slatest=homoS(fmin,fmax,mode_trnc,w[0],l[0],s_size);
    for i in range(0,np.shape(w)[0]-1):
            Sdis=s_matrix(fmin,fmax,mode_trnc,w[i],w[i+1],s_size);
            Slatest=cascade(Sdis,Slatest);
            Safter=homoS(fmin,fmax,mode_trnc,w[i+1],l[i+1],s_size);
            Slatest=cascade(Safter,Slatest);         
    return Slatest;
