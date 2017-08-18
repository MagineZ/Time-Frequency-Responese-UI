"""
import all needed packages
"""
import numpy as np
import math as math
import scipy as sp

"""
functions
"""

def hermf(N,M,tm):
    """
    hermf.m
    computes a set of orthonormal Hermite functions 
    (for use with tfrrsp_h.m) 
    
    input  : - N : number of points (must be odd) 
             - M : maximum order 
             - tm : half time support (>= 6 recommended) 
             
    output : - h : Hermite functions (MxN) 
             - Dh : H' (MxN) 
             - tt : time vector (1xN) 
    """
    
    dt = 2*tm/(N-1) ; 
    tt = np.linspace(-tm,tm,N) ; 
    g = np.exp(-np.square(tt)/2) ; 
    P = [np.ones(N)] 
    P = np.append(P, [2*tt], axis=0)

    for k in range(2,M+1):
        P = np.append(P,[2*np.multiply(tt,P[k-1]) - 2*(k-1)*P[k-2]],axis=0) ; 
     
    Htemp = [np.multiply(P[0],g)/math.sqrt(math.sqrt(math.pi)*math.gamma(1))*math.sqrt(dt)] ; 
    for k in range(1,M+1): 
        Htemp = np.append(Htemp, [np.multiply(P[k],g)/math.sqrt(math.sqrt(math.pi)
                           *math.pow(2,k)*math.gamma(k+1))*math.sqrt(dt)],axis = 0); 
       

    h = Htemp[0:M,:];          
    Dh = [(np.multiply(tt,Htemp[0]) - math.sqrt(2)*Htemp[1])*dt] ;
    for k in range(1,M): 
        Dh = np.append(Dh,[(np.multiply(tt,Htemp[k]) - math.sqrt(2*(k+1))*Htemp[k+1])*dt] ,axis=0);
      
    
    return h,Dh,tt


def dwindow(h):
    """
    %DWINDOW Derive a window.
    %	DH=DWINDOW(H) derives a window H.
    %
    %	Example : 
    %	 plot(dwindow(tftb_window(210,'hanning')))
    %
    %	See also WINDOW.
    
    %	F. Auger, August 1994, July 1995.
    %	Copyright (c) 1996 by CNRS (France).
    %
    %  This program is free software; you can redistribute it and/or modify
    %  it under the terms of the GNU General Public License as published by
    %  the Free Software Foundation; either version 2 of the License, or
    %  (at your option) any later version.
    %
    %  This program is distributed in the hope that it will be useful,
    %  but WITHOUT ANY WARRANTY; without even the implied warranty of
    %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %  GNU General Public License for more details.
    %
    %  You should have received a copy of the GNU General Public License
    %  along with this program; if not, write to the Free Software
    %  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
    """

    [hrow,hcol]=h.shape 
    if not (hrow==1):
        raise ValueError('h must have only one row')
    
    Lh = (hcol-1)/2
    step_height = (h[0,0]+h[0,-1])/2
    ramp = (h[0,-1]-h[0,0])/(hcol-1)
    h2 = np.append([0],h-step_height-ramp*np.arange(-Lh,Lh+1))
    h2 = np.append(h2,[0])
    Dh = (h2[2:hcol+2]-h2[0:hcol])/2+ramp
    Dh[0]+=step_height 
    Dh[-1]-=step_height
    return np.array(Dh)

def sqSTFTbase(x, lowFreq, highFreq, alpha, tDS, h, Dh, Smooth, Hemi):
    """
    % Synchrosqueezing modifed from tfrrsp.m, by Hau-tieng Wu, 2013 
    %
    %	computes the STFT and its SST
    %
    %   Example:
    %
    %	Hz = 32 ; t=[1/Hz:1/Hz:32]' ;
    %	x=cos(2*pi*(4*t+cos(t/2))) ;
    %	[h, Dh] = hermf(71, 1, 6) ;
    %		%% get the TF representation of the signal x with the frequency range
    %		%% [0,1, 0.4]*Hz with the frequency resolution 0.001*Hz
    %	[tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFT(x, 0.1, 0.4, 0.001, 1, h', Dh');
    %	imageRTF(t, tfrsqtic*Hz, abs(tfrsq)) ;
    %
    %		%% the first version can be recovered by
    %   [tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFT(x, 0, 0.5, 0.5/length(x), 1, h', Dh');
    %
    %
    %======================================
    %	X     : analysed signal.
    %	[lowFreq, highFreq] : frequency range \subset [0 0.5]. For the sake of computational efficiency.
    %	alpha : the resolution in the frequency axis
    %	tDS   : the time axis in the final TF representation is downsampled by tDS (set it to 1 unless you know what you are doing)
    %	H     : frequency smoothing window, H(0) being forced to 1
    %   DH    : differentiation of H	
    %	TFR   : STFT
    %	TFRSQ : synchrosqueezed STFT 
    %
    %	F. Auger, May-July 1994, July 1995.
    %	Copyright (c) 1996 by CNRS (France).
    %
    %	------------------- CONFIDENTIAL PROGRAM -------------------- 
    %	This program can not be used without the authorization of its
    %	author(s). For any comment or bug report, please send e-mail to 
    %	f.auger@ieee.org 
    """
    x=np.array([x])
    xrow,xcol = x.shape
    t = np.arange(1,xcol+1)
    tLen = len(np.arange(1,xcol+1,tDS))

    """ for tfr """
    N = len(np.arange(-0.5+alpha,0.5,alpha))+1
    """ for tfrsq """
    Lidx = round( (N/2)*(lowFreq/0.5) )+1 
    Hidx = round( (N/2)*(highFreq/0.5))
    fLen = Hidx - Lidx + 1

    "=========================================================================="
    "check input signals"
   
    if (xrow!=1):
        raise ValueError('X must have only one row')
    elif (highFreq > 0.5):
        raise ValueError('TopFreq must be a value in [0, 0.5]')
    elif (tDS < 1) or (tDS%1): 
        raise ValueError('tDS must be an integer value >= 1')
     
    hrow,hcol = h.shape
    Lh = int((hcol-1)/2) 
    if ((hrow!=1)or((hcol%2)==0)):
        raise ValueError('H must be a smoothing window with odd length')

    "=========================================================================="
    "run STFT and reassignment rule"
    tfr = np.zeros((tLen, int(N/2)),dtype=complex)
    tfrtic = np.linspace(0, 0.5, math.ceil(N/2)) 
    tfrsq = np.zeros((tLen, fLen),dtype=complex)
    tfrsqtic = np.linspace(lowFreq, highFreq, fLen)


    """Ex = mean(abs(x(min(t):max(t))).^2);"""
    Ex = np.mean(np.square(np.absolute(x)))
    Threshold = math.pow(10,-8)*Ex  
    """% originally it was 1e-6*Ex"""

    Mid = round(len(tfrsqtic)/2)
    Delta = 20*np.square(tfrsqtic[1]-tfrsqtic[0]) 
    weight = np.exp(-np.square(tfrsqtic[Mid-11:Mid+10]-tfrsqtic[Mid-1])/Delta)
    weight = weight/np.sum(weight)
    weightIDX = np.arange(Mid-10,Mid+11,1) - Mid

    for tidx in range(1,tLen+1):
        ti = t[(tidx-1)*tDS] 
        tau = (np.arange(-np.amin([round(N/2)-1,Lh,ti-1]),np.amin([round(N/2)-1,Lh,xcol-ti])+1)).astype(int)
        norm_h=np.linalg.norm(h[:,(Lh+tau).astype(int)])
        """%norm_h = h(Lh+1)""" 
        indices= ((N+tau)%N+1).astype(int)
        tf0 = np.zeros((1, N),dtype=complex)
        tf1 = np.zeros((1, N),dtype=complex) 
        tf0[:,indices-1] = x[:,ti+tau-1]*np.conjugate(h[:,Lh+tau])/norm_h
        tf1[:,indices-1] = x[:,ti+tau-1]*np.conjugate(Dh[:,Lh+tau])/norm_h
        tf0 = np.fft.fft(tf0);tf0 = tf0[:,0:int(N/2)];tf1 = np.fft.fft(tf1);tf1 = tf1[:,0:int(N/2)]
        """% get the first order omega"""
        omega = np.zeros(tf1.shape)
        _,avoid_warn = np.nonzero(tf0)
        omega[:,avoid_warn]= np.round(np.imag(N*np.divide(tf1[:,avoid_warn],tf0[:,avoid_warn])/(2.0*math.pi)))
        sst = np.zeros((1,fLen),dtype=complex) 

        for jcol in range(1,int(N/2)+1):
            if abs(tfr[0,jcol-1]) > Threshold:
                jcolhat = int(jcol - omega[:,jcol-1])
                """%jcolhat = rem(rem(jcolhat-1,N)+N,N)+1;"""
                
                """if(jcolhat < Hidx + 1) & (jcolhat >= Lidx):
                sst(jcolhat-Lidx+1) = sst(jcolhat-Lidx+1) + tf0(jcol) ;
	    	      end"""     
                if (jcolhat <= Hidx) & (jcolhat >= Lidx):   
                    """%IDXa = unique(min(Hidx, max(Lidx, jcolhat-Lidx+1+weightIDX))) ; """
                    if Smooth:
                        IDXb = np.where((jcolhat-Lidx+1+weightIDX <= Hidx) & (jcolhat-Lidx+1+weightIDX >= Lidx))
                        IDXa = (jcolhat-Lidx+1+weightIDX[IDXb]).astype(int)
                        
                        if Hemi:
                            if np.real(tf0[0,jcol-1])>0:
                                sst[0,IDXa-1]=sst[0,IDXa-1]+tf0[:,jcol-1]*weight[IDXb]
                                """%sst(jcolhat-Lidx+1) = sst(jcolhat-Lidx+1) + tf0(jcol) ;
                                %sst(jcolhat,icol) = rtfr(jcolhat,icol) + tfr(jcol,icol) ;"""
                            else:
                                sst[0,IDXa-1]=sst[0,IDXa-1]-tf0[:,jcol-1]*weight[IDXb]
                                """%sst(jcolhat-Lidx+1) = sst(jcolhat-Lidx+1) - tf0(jcol) ;
                                %sst(jcolhat,icol) = rtfr(jcolhat,icol) - tfr(jcol,icol) ;"""
                        else:
                            sst[0,IDXa-1]=sst[0,IDXa-1]+tf0[:,jcol-1]*weight[IDXb]
                	
                    else:
                        if Hemi:
                            if np.real(tf0[0,jcol-1])>0:
                                sst[0,jcolhat-Lidx]=sst[0,jcolhat-Lidx]+tf0[0,jcol-1]
                                """%sst(jcolhat,icol) = rtfr(jcolhat,icol) + tfr(jcol,icol)"""
                            else:
                                sst[0,jcolhat-Lidx]=sst[0,jcolhat-Lidx]-tf0[0,jcol-1]
                                """%sst(jcolhat,icol) = rtfr(jcolhat,icol) - tfr(jcol,icol)"""
                        else:
                            sst[0,jcolhat-Lidx]=sst[0,jcolhat-Lidx]+tf0[0,jcol-1] 
		   			
        tfr[tidx-1, :] = tf0[0,0:int(N/2)]
        tfrsq[tidx-1,:]=sst
    return np.transpose(tfr), tfrtic, np.transpose(tfrsq), tfrsqtic  





def sqSTFTbase2nd(x, lowFreq, highFreq, alpha, tDS, h, Dh, DDh, online):
    """
    % Synchrosqueezing modifed from tfrrsp.m, by Hau-tieng Wu, 2013 
    %
    %	computes the STFT and its SST
    %
    %   Example:
    %
    %	Hz = 32 ; t=[1/Hz:1/Hz:32]' ;
    %	x=cos(2*pi*(4*t+cos(t/2))) ;
    %	[h, Dh] = hermf(71, 1, 6) ;
    %		%% get the TF representation of the signal x with the frequency range
    %		%% [0,1, 0.4]*Hz with the frequency resolution 0.001*Hz
    %	[tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFT(x, 0.1, 0.4, 0.001, 1, h', Dh');
    %	imageRTF(t, tfrsqtic*Hz, abs(tfrsq)) ;
    %
    %		%% the first version can be recovered by
    %   [tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFT(x, 0, 0.5, 0.5/length(x), 1, h', Dh');
    %
    %
    %======================================
    %	X     : analysed signal.
    %	[lowFreq, highFreq] : frequency range \subset [0 0.5]. For the sake of computational efficiency.
    %	alpha : the resolution in the frequency axis
    %	tDS   : the time axis in the final TF representation is downsampled by tDS (set it to 1 unless you know what you are doing)
    %	H     : frequency smoothing window, H(0) being forced to 1
    %   DH    : differentiation of H	
    %	TFR   : STFT
    %	TFRSQ : synchrosqueezed STFT 
    %
    %	F. Auger, May-July 1994, July 1995.
    %	Copyright (c) 1996 by CNRS (France).
    %
    %	------------------- CONFIDENTIAL PROGRAM -------------------- 
    %	This program can not be used without the authorization of its
    %	author(s). For any comment or bug report, please send e-mail to 
    %	f.auger@ieee.org 
    """
    x=np.array([x])
    xrow,xcol = x.shape
    t = np.arange(1,xcol+1)
    tLen = len(np.arange(1,xcol+1,tDS))

    """ for tfr """
    N = len(np.arange(-0.5+alpha,0.5,alpha))+1
    """ for tfrsq """
    Lidx = int(np.ceil( (N/2)*(lowFreq/0.5) )+1) 
    Hidx = int(np.floor( (N/2)*(highFreq/0.5)))
    fLen = Hidx - Lidx + 1
    """===================================================================="""
    """check input signals"""
    if (xrow!=1):
        raise ValueError('X must have only one row')
    elif (highFreq > 0.5):
        raise ValueError('TopFreq must be a value in [0, 0.5]')
    elif (tDS < 1) or (tDS%1): 
        raise ValueError('tDS must be an integer value >= 1')

    hrow,hcol = h.shape
    Lh = int((hcol-1)/2) 
    if ((hrow!=1)or((hcol%2)==0)):
        raise ValueError('H must be a smoothing window with odd length')

    ht = np.arange(-Lh,Lh+1)
    """===================================================================="""
    """run STFT and reassignment rule"""
    if online: 
        tfr = np.zeros((100, int(N/2)),dtype=complex)     
        """ for h"""
        tfrsq = np.zeros((100, fLen),dtype=complex)
        tfrsq2nd = np.zeros((100, fLen),dtype=complex)
    else: 
        tfr = np.zeros((tLen, int(N/2)),dtype=complex) 
        """for h"""
        tfrsq = np.zeros((tLen, fLen),dtype=complex)
        tfrsq2nd = np.zeros((tLen, fLen),dtype=complex) 

    tfrtic = np.linspace(0, 0.5, N/2)
    tfrsqtic = np.linspace(lowFreq, highFreq, fLen)
    
    """Ex = mean(abs(x(min(t):max(t))).^2);"""
    Ex = np.mean(np.square(np.absolute(x[0,np.amin(t):np.amax(t)+1])))
    Threshold = math.pow(10,-8)*Ex  
    """% originally it was 1e-6*Ex"""

    for tidx in range(1,tLen+1):
        """ti is the current time"""
        ti = t[(tidx-1)*tDS]
        """tau is the relevant index associated with ti"""
        tau = (np.arange(-np.amin([round(N/2)-1,Lh,ti-1]),np.amin([round(N/2)-1,Lh,xcol-ti])+1)).astype(int)
        """indices is the absolute index in the "evaluation window" """
        indices= ((N+tau)%N+1).astype(int)
        norm_h=np.linalg.norm(h[:,(Lh+tau).astype(int)])
        tf0 = np.zeros((1, N),dtype=complex)
        tf1 = np.zeros((1, N),dtype=complex)
        tf2 = np.zeros((1, N),dtype=complex)
        tfx0 = np.zeros((1, N),dtype=complex)
        tfx1 = np.zeros((1, N),dtype=complex)
        tf0[:,indices-1] = x[:,ti+tau-1]*np.conjugate(h[:,Lh+tau])/norm_h
        tf1[:,indices-1] = x[:,ti+tau-1]*np.conjugate(Dh[:,Lh+tau])/norm_h
        tf2[:,indices-1] = x[:,ti+tau-1]*np.conjugate(DDh[:,Lh+tau])/norm_h
        tfx0[:,indices-1] = x[:,ti+tau-1]*np.conjugate(h[:,Lh+tau])*ht[Lh+tau]/norm_h
        tfx1[:,indices-1] = x[:,ti+tau-1]*np.conjugate(Dh[:,Lh+tau])*ht[Lh+tau]/norm_h
        tf0 = np.fft.fft(tf0);tf0 = tf0[:,0:int(N/2)]
        tf1 = np.fft.fft(tf1);tf1 = tf1[:,0:int(N/2)]
        tf2 = np.fft.fft(tf2);tf2 = tf2[:,0:int(N/2)]
        tfx0 = np.fft.fft(tfx0);tfx0 = tfx0[:,0:int(N/2)]
        tfx1 = np.fft.fft(tfx1);tfx1 = tfx1[:,0:int(N/2)]	
        """% get the first order omega"""
        omega = np.round(N*np.imag(tf1/tf0)/(2.0*math.pi))
        """% get the 2nd order omega"""
        omega2nd = np.round(N * np.imag(tf1/tf0 - (tf0*tf2-tf1*tf1)/(tfx1*tf0-tfx0*tf1)*tfx0/tf0)/(2.0*math.pi))
        
        
        sst = np.zeros((1,fLen),dtype=complex)
        sst2nd = np.zeros((1,fLen),dtype=complex)
        
        for jcol in range(1,int(N/2)+1):
            if abs(tfr[0,jcol-1]) > Threshold:
                
                jcolhat = int(jcol - omega[:,jcol-1])
                jcolhat2nd = int(jcol - omega2nd[:,jcol-1])
                
                if (jcolhat < Hidx+1) and (jcolhat >= Lidx):
                    sst[0,jcolhat-Lidx] = sst[0,jcolhat-Lidx]+tf0[0,jcol-1]
                if (jcolhat2nd < Hidx+1) and (jcolhat2nd >= Lidx):
                    sst2nd[0,jcolhat2nd-Lidx] = sst2nd[0,jcolhat2nd-Lidx]+tf0[0,jcol-1]
                
                
        if online:
            tfr[0:99,:] = tfr[1:100,:]
            tfrsq[0:99,:] = tfrsq[1:100,:]
            tfrsq2nd[0:99,:] = tfrsq2nd[1:100,:]
            tfr[99,:] = tf0[0,0:int(N/2)]
            tfrsq[99,:] = sst
            tfrsq2nd[99,:] = sst2nd
        
        else:
            tfr[tidx-1,:] = tf0[0,0:int(N/2)]
            tfrsq[tidx-1,:] = sst
            tfrsq2nd[tidx-1,:] = sst2nd

    return np.transpose(tfr), tfrtic, np.transpose(tfrsq), np.transpose(tfrsq2nd), tfrsqtic



def ConceFT_sqSTFT_C(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp, MT, Second, Smooth, Hemi):
    """
    Usage: 
        [tfrsq, ConceFT, tfrsqtic] = sqSTFT(t, x, lowFreq, highFreq, alpha, WinLen, dim, supp, MT)
        MT = 1: ordinary SST; MT > 1: ConceFT
        alpha: resolution in the frequency axis
        WinLen, dim, supp: for hermf.m
    Example:
        [tfrsq, ConceFT, tfrsqtic] = ConceFt_sqSTFT_C([1:length(y)]', y, 0,0.5, 0.0002, 121, 4, 6, 10);
    """ 
    """       
    Multitapering
    generate the window for short time Fourier transform (STFT)
    """
    h,Dh,_ = hermf(WinLen, dim, supp)
                  

    print('Run ordinary STFT-SST (Smooth = ',Smooth,', Hemi = ',Hemi,')\n')
    tfr, tfrtic, tfrsq, tfrsqtic = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, np.conj(np.array([h[0,:]])), np.conj(np.array([Dh[0,:]])), Smooth, Hemi)
    ConceFT = tfrsq 


    if MT > 1:
        print('Complex sphere\n') 
        print('STFT-ConceFT total (Smooth = ',Smooth,', Hemi = ',Hemi,'): ',MT,'; now:     ')
        
        for ii in range(1, MT+1):
            print('\b\b\b\b') 
            print('%4d' % ii) 
            rv = np.random.randn(1, dim) + np.sqrt(-1+0j)*np.random.randn(1, dim)
            rv = rv/np.linalg.norm(rv)
            rh = np.dot(rv , h)  
            rDh = np.dot(rv , Dh)
            
            if not Second:
                _, _, tfrsqX, tfrsqtic = sqSTFTbase(x, lowFreq, highFreq, alpha, hop, np.conj(rh), np.conj(rDh), Smooth, Hemi)
            else:
                _, _, _, tfrsqX, tfrsqtic = sqSTFTbase2nd(x, lowFreq, highFreq, alpha, hop, np.conj(rh), np.conj(rDh), np.array([dwindow(np.conj(rDh))]), 0);
		
            ConceFT = ConceFT + tfrsqX 
                 
    

        ConceFT = ConceFT/(MT+1)
        print('\n') 

    return tfr, tfrtic, tfrsq, ConceFT, tfrsqtic

def imageSQ(fig, axx, clb, t, ytic, M, Qv, scale_selection, cm = 'Greys'):
    axx.clear()
    clb.remove()
    import matplotlib.colors as colors
    q = np.percentile(M,Qv)
    """ truncate the upper bound"""
    M[np.where(M>q)]=q
    """M = M/q"""
    """
    truncate the lower bound
    %m = quantile(Q, 0.002) ;
    %M(find(M<m)) = m ;
    """
    
    M = M[:-1,:-1]
    
    if scale_selection == 0:
        cax = axx.pcolorfast(t, ytic, M, cmap = cm)
        #clbb = fig.colorbar(cax, ax=axx, extend='max')
    
    elif scale_selection == 1:
        """
        if M.min() == 0:
            M+=1
            
        cax= axx.pcolorfast(t, ytic, M, norm=colors.LogNorm(vmin=M.min(), vmax=M.max()), cmap=cm)
        """
        cax= axx.pcolorfast(t, ytic, M, norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,vmin=M.min(), vmax=M.max()), cmap=cm)
        #clb.set_clim(vmin=M.min(),vmax=M.max()) 
        #clb.draw_all()
    clb = fig.colorbar(cax, ax=axx, extend='max')
    
    axx.set_xlim([np.min(t),np.max(t)])
    axx.set_ylim([np.min(ytic),np.max(ytic)])
    return clb
    

  
def wwhat(f,beta,gam,k):
    """% This function calculates morse complex wavelets (in frequency domain)
    % It is very similar to morsefreqs"""
    from scipy.special import gamma
    from math import sqrt, pi
    from numpy import power
    c=-1+((2*beta+1)/gam)
    con=sqrt(gam*power(2,((2.*beta+1)/gam))*gamma(k+1))
    con2=sqrt(pi*gamma(k+((2.*beta+1)/gam)))
    con=sqrt(pi/2)*(con/con2)
    arr=2*(power(2*pi*abs(f),gam))
    m=con*power(2*pi*abs(f),beta)*np.exp(-power(2*pi*abs(f),gam))*laggen(arr,k,c)
    m=sqrt(4*pi)*m
    
    for j in range(0,len(f)):
        if (f[j]<0):m[:,j]=0
    
    return m

def laggen(x,k,c):
    """
    % compute generalized Laguerre poly L_k^c(x)
    %
    %
    """
    from scipy.special import gammaln
    from numpy import power
    l=len(x)
    sn=-1;
    s=np.zeros((1,l))
    for m in range(0,k+1):
        sn = -sn
        ga = gammaln(k+c+1)-gammaln(k-m+1)-gammaln(c+m+1)-gammaln(m+1)
        ga = np.exp(ga)
        s=s+sn*ga*power(x,m)
    return s

    
class init:
    motherwavelet = 'Cinfc'
    CENTER = 1
    FWHM = 0.3
    
def CWT(t, x, opts = init):
    """%
    % Continuous wavelet transform ver 0.1
    % Modified from wavelab 85
    %
    % INPUT:
    % OUTPUT:
    % DEPENDENCY:
    %
    % by Hau-tieng Wu 2011-06-20 (hauwu@math.princeton.edu)
    % by Hau-tieng Wu 2012-12-23 (hauwu@math.princeton.edu)
    %
    """

    """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    """prepare for the input"""
    from math import sqrt, pi
    from numpy import power
    nvoice = 32;
    scale = 2;
    Oct = 1;
    
    """++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    """start to do CWT"""
    
    n = len(x)
    """%% assume the original signal is on [0,L]."""
    """%% assume the signal is on [0,1]. Frequencies are rescaled to xi/L"""
    
    
    if (n%2==0):
        xi = np.concatenate((np.arange(0,n/2+1) , np.arange(-n/2+1,0)))
    else:
        xi = np.concatenate((np.arange(0,n/2) , np.arange(-int(n/2),0)))
    
    xhat = np.fft.fft(x)
    
    noctave = np.floor(np.log2(n)) - Oct
    tfr = np.zeros((n,int(nvoice*noctave)),dtype=complex)
    kscale = 1
    tfrtic = np.zeros(int(nvoice*noctave))
    for jj in range(1 , int(nvoice*noctave)+1):
        tfrtic[jj-1] = scale*power(2,jj/nvoice)



    if (opts.motherwavelet=='morse-b' or opts.motherwavelet=='morse-c'):
        dim = opts.dim
        uFix = np.random.randn(dim,1)
        uFix = -sp.linalg.orth(uFix)


    for jo in range(1,int(noctave)+1): 
        """# of scales"""
        for jv in range(1,nvoice+1):
            qscale = scale * power(2,jv/nvoice)
            omega =  xi/qscale            
             
            if (opts.motherwavelet == 'morse'):
                windowq = wwhat(np.conj(omega), opts.beta, opts.gam, opts.k)
                windowq = np.conj(windowq.T)
                
            elif (opts.motherwavelet == 'morse-b'):
                u = np.random.randn(opts.dim,1)
                u = -sp.linalg.orth(u)
                W = np.zeros((len(omega), opts.dim))
                for ki in range(1,opts.dim+1):
                    W[:,ki-1] = wwhat(np.conj(omega), opts.beta, opts.gam, ki-1)
                
                windowq = np.dot(W ,u)
            
            elif (opts.motherwavelet == 'morse-c'):
                W = np.zeros(len(omega), opts.dim)
                for ki in range(1,opts.dim+1):
                    W[:,ki-1] = wwhat(np.conj(omega), opts.beta, opts.gam, ki-1)
            
                windowq = np.dot(W,uFix)
            
            elif (opts.motherwavelet == 'morse-a'):
                if(opts.k == 0):
                    windowq1 = wwhat(np.conj(omega),opts.beta,opts.gam,0)
                    windowq2 = wwhat(np.conj(omega),opts.beta,opts.gam,1)
                    windowq = ( windowq1 + windowq2 ) / sqrt(2)
                elif (opts.k == 1):
                    windowq1 = wwhat(np.conj(omega),opts.beta,opts.gam,0);
                    windowq2 = wwhat(np.conj(omega),opts.beta,opts.gam,1);
                    windowq = ( windowq1 - windowq2 ) / sqrt(2);
            
                windowq = np.conj(windowq.T)
            
            elif (opts.motherwavelet == 'Cinfc'):
                tmp0 = (omega-opts.CENTER)/opts.FWHM
                tmp1 = power(tmp0,2)-1
                windowq = np.exp(1/tmp1)
                windowq[omega >= (opts.CENTER+opts.FWHM)] = 0
                windowq[omega <= (opts.CENTER-opts.FWHM)] = 0
                windowq = np.array([windowq]).T

            elif (opts.motherwavelet == 'morlet'):
                windowq = 4*sqrt(pi)*np.exp(-4*power(omega-0.69*pi,2))-4.89098e-4*4*sqrt(pi)*np.exp(-4*power(omega,2))
                windowq = np.array([windowq]).T

            elif (opts.motherwavelet == 'gaussian'):
                psihat = lambda f: np.exp( -np.log(2)*power( 2*(f-opts.CENTER)/opts.FWHM,2) )
                windowq = psihat(omega)
                windowq = np.array([windowq]).T


            elif (opts.motherwavelet == 'meyer'):
                windowq = np.zeros(len(omega))
                int1 = np.logical_and(omega>=5/8*0.69*pi , omega<0.69*pi)
                int2 = np.logical_and(omega>=0.69*pi , omega<7/4*0.69*pi)
                meyeraux = lambda f: 35*power(f,4)-84*power(f,5)+70*power(f,6)-20*power(f,7)
                windowq[int1] = np.sin(pi/2*meyeraux((omega[int1]-5/8*0.69*pi)/(3/8*0.69*pi)))
                windowq[int2] = np.cos(pi/2*meyeraux((omega[int2]-0.69*pi)/(3/4*0.69*pi)))
                windowq = np.array([windowq]).T


            elif (opts.motherwavelet == 'BL3'):
                phihat = power(2*pi,-0.5)*power((np.sin(omega/4)/(omega/4)),4)
                phihat[0] = power(2*pi,-0.5)
                aux1 = 151/315 + 397/840*np.cos(omega/2) + 1/21*np.cos(omega) + 1/2520*np.cos(3*omega/2)
                phisharphat = phihat*power(aux1,-0.5)

                aux2 = 151/315 - 397/840*np.cos(omega/2) + 1/21*np.cos(omega) - 1/2520*np.cos(3*omega/2)
                aux3 = 151/315 + 397/840*np.cos(omega) + 1/21*np.cos(2*omega) + 1/2520*np.cos(3*omega)
                msharphat = power(np.sin(omega/4),4)*power(aux2,0.5)*power(aux3,-0.5)
                windowq = phisharphat*msharphat*np.exp(1j*omega/2)*(omega>=0)
                windowq = np.array([windowq]).T
                


            windowq = windowq/sqrt(qscale)
            what = windowq*np.array([xhat]).T
            w = np.fft.ifft(what.T)
            tfr[:,kscale-1] = w
            kscale+=1



        scale*=2




    """%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    """%% calculate the constant for reconstruction
    %% TODO: calculate Rpsi for other mother wavelets"""
    xi = np.arange(0.05,10+1/10000,1/10000)
    
    if (opts.motherwavelet == 'gaussian'):   
        """%% Gaussian (not really wavelet)"""

        psihat = lambda f: np.exp( -np.log(2)*power( 2*(f-opts.CENTER)/opts.FWHM,2))
        windowq = psihat(xi)
        Rpsi = sum(windowq/xi)/10000


    elif (opts.motherwavelet == 'morlet'):
        windowq = 4*sqrt(pi)*np.exp(-4*power(xi-0.69*pi,2))-4.89098e-4*4*sqrt(pi)*np.exp(-4*power(xi,2))
        Rpsi = sum(windowq/xi)/10000

    elif (opts.motherwavelet == 'Cinfc'):
        tmp0 = (xi - opts.CENTER)/opts.FWHM
        tmp1 = power(tmp0,2)-1
        windowq = np.exp(1/tmp1)
        windowq[xi >= (opts.CENTER+opts.FWHM)] = 0
        windowq[xi <= (opts.CENTER-opts.FWHM)] = 0
        Rpsi = sum(windowq/xi)/10000
 
    else:
        """%% Normalization is not implemented for Other mother wavelets"""
        Rpsi = 1 
  
    tfr = tfr/Rpsi
    
    return tfr.T, tfrtic


def sqCWTbase(t, x, freqlow, freqhigh, alpha, opts, Smooth, Hemi):
    """
    %
    % Synchrosqueezing transform ver 0.5 (2015-03-09)
    % You can find more information in 
    %	http://sites.google.com/site/hautiengwu/
    %
    % Example: [~, tfrsq, ~, tfrsqtic] = sqCWT(time, xm, lowfreq, highfreq, alpha, opts);
    %	time: 	time of the signal
    %	xm: 	the signal to be analyzed
    %	[lowfreq, highfreq]: the frequency range in the output time-frequency representation. For the sake of computational efficiency.
    %	alpha:	the frequency resolution in the output time-frequency representation
    %	opts:	parameters for the CWT analysis. See below
    %	tfr/tfrtic:	the CWT and its scale tic
    %	tfrsq/tfrsqtic: the SST-CWT and its frequency tic
    %
    % by Hau-tieng Wu v0.1 2011-06-20 (hauwu@math.princeton.edu)
    %		  v0.2 2011-09-10
    %		  v0.3 2012-03-03
    %		  v0.4 2012-12-12
    %		  v0.5 2015-03-09
    %% you can play with these 4 parameters, but the results might not be
    %% that different if the change is not crazy
    """
    Ex = np.mean(np.power(abs(x),2))
    Gamma = 1.0e-8*Ex  
    """% originally it was 1e-6*Ex"""    

    """%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    dt = t[1] - t[0]

    
    """%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    """%% Continuous wavelet transform""" 
    tfr, tfrtic = CWT(t, x, opts)
    Dtfr = (-1j/2/math.pi/dt)*np.concatenate((tfr[:,1:] - tfr[:,0:-1], np.array([tfr[:,-1]-tfr[:,-2]]).T),axis = 1) 
    Dtfr[(abs(tfr) < Gamma)] = float('nan')
    omega = Dtfr/tfr

    """%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    """%% Synchro-squeezing transform"""

    tfrsq, tfrsqtic = SQ(tfr, omega, freqlow, freqhigh, alpha, Smooth, Hemi)
    return tfr, tfrtic, tfrsq, tfrsqtic    




"""%====================================================================="""
"""%% function for CWT-based SST"""
def SQ(tfd, omega, tfrsqticlow, tfrsqtichigh, alpha, Smooth, Hemi):
    
    nvoice = 32
    scale = 2
    omega = abs(omega)
    nscale,n = tfd.shape
    nalpha = np.floor((tfrsqtichigh - tfrsqticlow)/alpha)
    tfrsq = np.zeros((int(nalpha),n),dtype=complex)
    tfrsqtic = np.arange(1,nalpha+1)*alpha + tfrsqticlow
    ntfrsqtic = len(tfrsqtic)
    
    Mid = round(ntfrsqtic/2)-1
    Delta = 20*np.power(tfrsqtic[1]-tfrsqtic[0],2) 
    weight = np.exp(-np.power(tfrsqtic[Mid-10:Mid+11]-tfrsqtic[Mid],2)/Delta)
    weight = weight / sum(weight) 
    weightIDX = np.arange(Mid-10,Mid+11) - Mid

    for b in range(1,n+1):
        for kscale in range(1,nscale+1):
            qscale = scale*np.power(2,(kscale/nvoice))
            
            if (np.isfinite(omega[kscale-1,b-1]) and (omega[kscale-1,b-1]>0)):
                k = int(np.floor( ( omega[kscale-1,b-1] - tfrsqticlow )/ alpha )+1)

                if (np.isfinite(k) and (k > 0) and (k < ntfrsqtic-1)):
                    ha = tfrsqtic[k]-tfrsqtic[k-1]


                    if Smooth:
                        IDXb = np.where(np.logical_and((k+weightIDX < ntfrsqtic-1) , (k+weightIDX > 0)))
                        IDXa = k+weightIDX[IDXb]
                        IDXa = IDXa.astype(int)
                        
                        if Hemi:
                            if (tfd[kscale-1,b-1].real > 0):
                                tfrsq[IDXa-1,b-1] = tfrsq[IDXa-1,b-1] + weight[IDXb]*np.log(2)*tfd[kscale-1,b-1]*np.sqrt(qscale)/ha/nvoice
                            else:
                                tfrsq[IDXa-1,b-1] = tfrsq[IDXa-1,b-1] - weight[IDXb]*np.log(2)*tfd[kscale-1,b-1]*np.sqrt(qscale)/ha/nvoice
                    	
                        else:
                            tfrsq[IDXa-1,b-1] = tfrsq[IDXa-1,b-1] + np.log(2)*tfd[kscale-1,b-1]*np.sqrt(qscale)/ha/nvoice
                	
                    else:
                        if Hemi:
                            if (tfd[kscale-1,b-1].real>0):
                                tfrsq[k-1,b-1] = tfrsq[k-1,b-1] + np.log(2)*tfd[kscale-1,b-1]*np.sqrt(qscale)/ha/nvoice
                            else:
                                tfrsq[k-1,b-1] = tfrsq[k-1,b-1] - np.log(2)*tfd[kscale-1,b-1]*np.sqrt(qscale)/ha/nvoice
                                     
                        else:
                            tfrsq[k-1,b-1] = tfrsq[k-1,b-1] + np.log(2)*tfd[kscale-1,b-1]*np.sqrt(qscale)/ha/nvoice
                                 
    return tfrsq, tfrsqtic


def ConceFT_CWT(t, x, lowfreq, highfreq, alpha, MT, opts, Smooth, Hemi):
     
    """%% ordinary SST"""
    print('Run ordinary CWT-SST (Smooth = ',Smooth,', Hemi = ',Hemi,')\n')
    tfr, tfrtic, tfrsq, tfrsqtic = sqCWTbase(t, x, lowfreq, highfreq, alpha, opts, Smooth, 0)

    """%===========================
    %% get the ConceFT"""	

    ConceFT = tfrsq


    if MT > 1:
        print('CWT-ConceFT total (Smooth = ',Smooth,', Hemi = ',Hemi,'): ',MT,'; now:     ')
        
        for ii in range(1,MT+1):
            print('\b\b\b\b') 
            print('%4d' % ii)
            _, _, tfrsqX, tfrsqtic = sqCWTbase(t, x, lowfreq, highfreq, alpha, opts, Smooth, Hemi)
            
            ConceFT = ConceFT + tfrsqX 
	
        print('\n')

        ConceFT = ConceFT/(MT+1)
        
    return tfr, tfrtic, tfrsq, ConceFT, tfrsqtic

def stran(h,t,ifsq=0,freqlow=0,freqhigh=0.5, alpha=0.001, Smooth=0, Hemi=0):
    """% Compute S-Transform without for loops

    %%% Coded by Kalyan S. Dash %%%
    %%% IIT Bhubaneswar, India %%%"""
    import numpy.matlib
    import scipy.linalg
     
    """% originally it was 1e-6*Ex"""    

    """%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
    N = len(h) 
    """% h is a 1xN one-dimensional series"""

    nhaf=int(np.fix(N/2))
    odvn=1

    if (nhaf*2==N):
        odvn=0
    
    f = np.concatenate((np.arange(0,nhaf+1),np.arange(-nhaf+1-odvn,0)))/N
                      
    Hft=np.fft.fft(h)

    """%Compute all frequency domain Gaussians as one matrix"""

    invfk=np.conj(np.array([1/f[2:nhaf+2]])).T
    
    W=2*np.pi*numpy.matlib.repmat(f,nhaf,1)*numpy.matlib.repmat(invfk,1,N)
    G=np.exp((-np.power(W,2))/2) 
    
    """%Gaussian in freq domain"""
    
    """% End of frequency domain Gaussian computation"""

    """% Compute Toeplitz matrix with the shifted fft(h)"""

    HW=scipy.linalg.toeplitz(np.conj(np.array([Hft[0:nhaf+1]])).T,Hft)

    """% Exclude the first row, corresponding to zero frequency"""

    HW=HW[1:nhaf+1,:]

    """% Compute Stockwell Transform"""

    tfr=np.fft.ifft(HW*G,axis=1) 
    """%Compute voice"""
    
    """%Add the zero freq row"""
    st0=np.mean(h)*np.ones((1,N))
    
    tfr=np.concatenate((st0,tfr),axis = 0)
    
    tfrtic = np.linspace(0, 0.5, tfr.shape[0])
    
    tfrsq = []
    tfrsqtic = []
    
    if ifsq:
        Ex = np.mean(np.power(abs(h),2))
        Gamma = 1.0e-8*Ex
        dt = t[1] - t[0]
        Dtfr = (-1j/2/math.pi/dt)*np.concatenate((tfr[:,1:] - tfr[:,0:-1], np.array([tfr[:,-1]-tfr[:,-2]]).T),axis = 1) 
        Dtfr[(abs(tfr) < Gamma)] = float('nan')
        omega = Dtfr/tfr

        """%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
        """%% Synchro-squeezing transform"""

        tfrsq, tfrsqtic = SQ(tfr, omega, freqlow, freqhigh, alpha, Smooth, Hemi)
    
    return tfr, tfrtic, tfrsq,tfrsqtic


def rsSTFTbase(x, lowFreq, highFreq, alpha, tDS, h, Dh, Causality):
    """
    %
    % Synchrosqueezing modifed from tfrrsp.m, by Hau-tieng Wu, 2013 
    %
    %	computes the STFT and its SST
    %
    %   Example:
    %
    %	Hz = 32 ; t=[1/Hz:1/Hz:32]' ;
    %	x=cos(2*pi*(4*t+cos(t/2))) ;
    %	[h, Dh] = hermf(71, 1, 6) ;
    %		%% get the TF representation of the signal x with the frequency range
    %		%% [0,1, 0.4]*Hz with the frequency resolution 0.001*Hz
    %	[tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFT(x, 0.1, 0.4, 0.001, 1, h', Dh');
    %	imageRTF(t, tfrsqtic*Hz, abs(tfrsq)) ;
    %
    %		%% the first version can be recovered by
    %   [tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFT(x, 0, 0.5, 0.5/length(x), 1, h', Dh');
    %
    %
    %======================================
    %	X     : analysed signal.
    %	[lowFreq, highFreq] : frequency range \subset [0 0.5]. For the sake of computational efficiency.
    %	alpha : the resolution in the frequency axis
    %	tDS   : the time axis in the final TF representation is downsampled by tDS (set it to 1 unless you know what you are doing)
    %	H     : frequency smoothing window, H(0) being forced to 1
    %   DH    : differentiation of H	
    %	TFR   : STFT
    %	TFRSQ : synchrosqueezed STFT 
    %
    %	F. Auger, May-July 1994, July 1995.
    %	Copyright (c) 1996 by CNRS (France).
    %
    %	------------------- CONFIDENTIAL PROGRAM -------------------- 
    %	This program can not be used without the authorization of its
    %	author(s). For any comment or bug report, please send e-mail to 
    %	f.auger@ieee.org 
    """
    x=np.array([x])
    xrow,xcol = x.shape
    t = np.arange(1,xcol+1) 
    tLen = len(np.arange(1,xcol+1,tDS)) 

    """% for tfr"""
    N = len(np.arange(-0.5+alpha,0.5,alpha))+1

    """% for tfrsq"""
    Lidx = round( (N/2)*(lowFreq/0.5)) + 1 
    Hidx = round( (N/2)*(highFreq/0.5)) 
    fLen = Hidx - Lidx + 1 



    """%===================================================================="""
    """%% check input signals"""
    if (xrow!=1):
        raise ValueError('X must have only one row')
    elif (highFreq > 0.5):
        raise ValueError('TopFreq must be a value in [0, 0.5]')
    elif (tDS < 1) or (tDS%1): 
        raise ValueError('tDS must be an integer value >= 1')
    
    
    hrow,hcol = h.shape
    Lh = int((hcol-1)/2) 
    if ((hrow!=1)or((hcol%2)==0)):
        raise ValueError('H must be a smoothing window with odd length')

    Th = h*np.arange(-Lh,Lh+1)
    Dt = 1
    
    "=========================================================================="
    "run STFT and reassignment rule"
    tfr = np.zeros((tLen, int(N/2)),dtype=complex)
    tfrtic = np.linspace(0, 0.5, int(N/2)) 
    tfrrs = np.zeros((tLen, fLen),dtype=complex)
    tfrrstic = np.linspace(lowFreq, highFreq, fLen)
    
    Ex = np.mean(np.square(np.absolute(x[0,np.amin(t):np.amax(t)+1])))
    Threshold = 10e-8*Ex
    """% originally it was 1e-6*Ex"""


    for tidx in range(1,tLen+1):
        ti = t[(tidx-1)*tDS]
        tau = (np.arange(-np.amin([round(N/2)-1,Lh,ti-1]),np.amin([round(N/2)-1,Lh,xcol-ti])+1)).astype(int)
        norm_h=np.linalg.norm(h[:,(Lh+tau).astype(int)])
        indices= ((N+tau)%N+1).astype(int)
        tf0 = np.zeros((1, N),dtype=complex)
        tf1 = np.zeros((1, N),dtype=complex)
        tf2 = np.zeros((1, N),dtype=complex)
        tf0[:,indices-1] = x[:,ti+tau-1]*np.conjugate(h[:,Lh+tau])/norm_h
        tf1[:,indices-1] = x[:,ti+tau-1]*np.conjugate(Dh[:,Lh+tau])/norm_h
        tf2[:,indices-1] = x[:,ti+tau-1]*np.conjugate(Th[:,Lh+tau])/norm_h
        tf0 = np.fft.fft(tf0);tf0 = tf0[:,0:int(N/2)]
        tf1 = np.fft.fft(tf1);tf1 = tf1[:,0:int(N/2)]
        tf2 = np.fft.fft(tf2);tf2 = tf2[:,0:int(N/2)]

        """% get the first order omega"""
        omega = np.zeros(tf1.size)
        omegaT = np.zeros(tf1.size)
        avoid_warn_x,avoid_warn_y = np.where(tf0!=0)
        omega[avoid_warn_y] = np.round(np.imag(N*tf1[avoid_warn_x,avoid_warn_y]/tf0[avoid_warn_x,avoid_warn_y]/(2.0*np.pi)))
        omegaT[avoid_warn_y] = np.round(np.real(tf2[avoid_warn_x,avoid_warn_y]/tf0[avoid_warn_x,avoid_warn_y]/Dt))

        """%rs = np.zeros((1,fLen) ;"""

        for jcol in range(1,int(N/2)+1):
            if abs(tfr[0,jcol-1]) > Threshold:
                """%% to keep the causality, use max(tidxhat,tidx-Lh) instead of max(tidxhat,1)"""
                tidxhat = int(tidx + omegaT[jcol-1])
                if Causality:
                    tidxhat = min(max(max(tidxhat, tidx-Lh), 1), tLen)
                else:
                    tidxhat = min(max(tidxhat, 1), tLen)

                jcolhat = int(jcol - omega[jcol-1]) 
                """jcolhat = rem(rem(jcolhat-1%N)+N,N)+1;"""
                if(jcolhat <= Hidx) and (jcolhat >= Lidx):
                    tfrrs[tidxhat-1,jcolhat-Lidx] = tfrrs[tidxhat-1,jcolhat-Lidx] + tf0[0,jcol-1]


        tfr[tidx-1,:] = tf0[0:int(N/2)] 


    return np.transpose(tfr), tfrtic, np.transpose(tfrrs), tfrrstic


def ConceFT_rsSTFT_C(x, lowFreq, highFreq, alpha, hop, WinLen, dim, supp, MT):
    
    """
    %
    % Usage: 
    % 	tfrsq, ConceFT, tfrsqtic = sqSTFT(t, x, lowFreq, highFreq, alpha, WinLen, dim, supp, MT)
    %
    % MT = 1: ordinary SST; MT > 1: ConceFT
    % alpha: resolution in the frequency axis
    % WinLen, dim, supp: for hermf.m
    %
    % Example:
    % 	tfrsq, ConceFT, tfrsqtic = sqSTFT([1:length(y)]', y, 0,0.5, 0.0002, 121, 4, 6, 10);
    """
    

    """%%%% Multitapering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"""
    """generate the window for short time Fourier transform (STFT)"""
    
    h, Dh, _ = hermf(WinLen, dim, supp) 
    """%======================================="""
    print('Run ordinary RS\n')
    
    tfr, tfrtic, tfrrs, tfrrstic = rsSTFTbase(x, lowFreq, highFreq, alpha, 1, np.conj(np.array([h[0,:]])),  np.conj(np.array([Dh[0,:]])), 1)

    """%======================================="""
    ConceFT = tfrrs

    if MT > 1:
        print('Complex sphere\n') 
        """%% Conceft"""
        ConceFT = np.zeros(tfrrs.shape) 
        print('RS-ConceFT total: ',MT,'; now:     ')
        
        for ii in range(1, MT+1):
            print('\b\b\b\b') 
            print('%4d' % ii) 
            rv = np.random.randn(1, dim) + np.sqrt(-1+0j)*np.random.randn(1, dim)
            rv = rv/np.linalg.norm(rv)
            rh = np.dot(rv , h)  
            rDh = np.dot(rv , Dh)
            _, _, tfrrs, tfrrstic = rsSTFTbase(x, lowFreq, highFreq, alpha, 1, np.conj(rh), np.conj(rDh), 1)
            
            ConceFT = ConceFT + tfrrs
        
        ConceFT = ConceFT/MT 
        print('\n') 
    
    return tfr, tfrtic, tfrrs, ConceFT, tfrrstic




def signal_filter(t,x,sampling_rate, detrend, DS, btype,f_type='butter',N=1,cutoff=[0.1,0.49],rp=None,rs = None):
    from scipy import signal
    nyq = sampling_rate/2
    time = t
    f_x = x
    cutoff*=2
    if detrend:
        f_x = signal.detrend(x)
    
    
    if btype != 'None':
        analog = False
        b,a = signal.iirfilter(N, cutoff, rp, rs, btype, analog, f_type)
        f_x = signal.filtfilt(b,a,f_x)
    
    if DS:
        f_x,time = signal.resample(f_x,int(len(time)/DS),time)
    
    return time, f_x