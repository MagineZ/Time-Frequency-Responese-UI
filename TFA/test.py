def rnd(n):
    import math
    if n%1==0.5:
        return math.ceil(n)
    else:
        return round(n)
    
    
def buffer(x, n, p=0, opt=None):
    '''Mimic MATLAB routine to generate buffer array

    MATLAB docs here: https://se.mathworks.com/help/signal/ref/buffer.html

    Args
    ----
    x:   signal array
    n:   number of data segments
    p:   number of values to overlap
    opt: initial condition options. default sets the first `p` values
         to zero, while 'nodelay' begins filling the buffer immediately.
    '''
    import numpy

    if p >= n:
        raise ValueError('p ({}) must be less than n ({}).'.format(p,n))

    # Calculate number of columns of buffer array
    cols = int(numpy.ceil(len(x)/(n-p)))

    # Check for opt parameters
    if opt == 'nodelay':
        # Need extra column to handle additional values left
        cols += 1
    elif opt != None:
        raise SystemError('Only `None` (default initial condition) and '
                          '`nodelay` (skip initial condition) have been '
                          'implemented')

    # Create empty buffer array
    b = numpy.zeros((n, cols))

    # Fill buffer by column handling for initial condition and overlap
    j = 0
    for i in range(cols):
        # Set first column to n values from x, move to next iteration
        if i == 0 and opt == 'nodelay':
            b[0:n,i] = x[0:n]
            continue
        # set first values of row to last p values
        elif i != 0 and p != 0:
            b[:p, i] = b[-p:, i-1]
        # If initial condition, set p elements in buffer array to zero
        else:
            b[:p, i] = 0

        # Get stop index positions for x
        k = j + n - p

        # Get stop index position for b, matching number sliced from x
        n_end = p+len(x[j:k])

        # Assign values to buffer array from x
        b[p:n_end,i] = x[j:k]

        # Update start index location for next iteration of x
        j = k

    return b
def tftb_window(N = None, name = None, param = None, param2 = None):
    # tftb_window	Window generation.
    #	H=tftb_window(N,NAME,PARAM,PARAM2)
    #	yields a window of length N with a given shape.
    #
    #	N      : length of the window
    #	NAME   : name of the window shape (default : Hamming)
    #	PARAM  : optional parameter
    #	PARAM2 : second optional parameters
    #
    #	Possible names are :
    #	'Hamming', 'Hanning', 'Nuttall',  'Papoulis', 'Harris',
    #	'Rect',    'Triang',  'Bartlett', 'BartHann', 'Blackman'
    #	'Gauss',   'Parzen',  'Kaiser',   'Dolph',    'Hanna'.
    #	'Nutbess', 'spline',  'Flattop'
    #
    #	For the gaussian window, an optionnal parameter K
    #	sets the value at both extremities. The default value is 0.005
    #
    #	For the Kaiser-Bessel window, an optionnal parameter
    #	sets the scale. The default value is 3*pi.
    #
    #	For the Spline windows, h=tftb_window(N,'spline',nfreq,p)
    #	yields a spline weighting function of order p and frequency
    #	bandwidth proportional to nfreq.
    #
    #       Example: 
    #        h=tftb_window(256,'Gauss',0.005); 
    #        plot(0:255, h); axis([0,255,-0.1,1.1]); grid
    #
    #	See also DWINDOW.

    #	F. Auger, June 1994 - November 1995.
    #	Copyright (c) 1996 by CNRS (France).
    #
    #  This program is free software; you can redistribute it and/or modify
    #  it under the terms of the GNU General Public License as published by
    #  the Free Software Foundation; either version 2 of the License, or
    #  (at your option) any later version.
    #
    #  This program is distributed in the hope that it will be useful,
    #  but WITHOUT ANY WARRANTY; without even the implied warranty of
    #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #  GNU General Public License for more details.
    #
    #  You should have received a copy of the GNU General Public License
    #  along with this program; if not, write to the Free Software
    #  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
    #
    #	References : 
    #	- F.J. Harris, "On the use of windows for harmonic
    #	analysis with the discrete Fourier transform",
    #	Proceedings of the IEEE, Vol 66, No 1, pp 51-83, 1978.
    #	- A.H. Nuttal, "A two-parameter class of Bessel weighting 
    #	functions for spectral analysis or array processing", 
    #	IEEE Trans on ASSP, Vol 31, pp 1309-1311, Oct 1983.
    #	- Y. Ho Ha, J.A. Pearce, "A New window and comparison to
    #	standard windows", Trans IEEE ASSP, Vol 37, No 2, 
    #	pp 298-300, February 1989.
    #	- C.S. Burrus, Multiband Least Squares FIR Filter Design,
    #	Trans IEEE SP, Vol 43, No 2, pp 412-421, February 1995.
    import numpy as np
    from scipy import special
    nargin = np.sum([N != None, name != None, param != None, param2 != None])
    
    if (nargin == 0):
        print('at least 1 parameter is required');return
    if (N<=0):
        print('N should be strictly positive.');return
    if (nargin==1): name = 'Hamming'
    name=name.upper()
    if (name=='RECTANG') or (name=='RECT'): 
        h=np.ones([N,1])
    elif (name == 'HAMMING'):
        h=0.54 - 0.46*np.cos(2.0*np.pi*np.arange(1,N+1)/(N+1))
    elif (name == 'HANNING') or (name=='HANN'):
        h=0.50 - 0.50*np.cos(2.0*np.pi*np.arange(1,N+1)/(N+1))
    elif (name =='KAISER'):       
        if (nargin==3): beta=param
        else: beta=3.0*np.pi
        ind=np.arange(-(N-1)/2,(N-1)/2+1)*2/N
        beta=3.0*np.pi
        h=np.real(special.jv(0,1j*beta*np.sqrt(1.0-ind*ind))/np.real(special.jv(0,1j*beta)))
    elif (name=='NUTTALL'):
        ind=np.arange(-(N-1)/2,(N-1)/2+1)*2.0*np.pi/N
        h=+0.3635819 + 0.4891775*np.cos(ind) + 0.1363995*np.cos(2.0*ind) + 0.0106411*np.cos(3.0*ind)
    elif (name=='BLACKMAN'):
        ind=np.arange(-(N-1)/2,(N-1)/2+1)*2.0*np.pi/N
        h= +0.42 + 0.50*np.cos(ind) + 0.08*np.cos(2.0*ind)
    elif (name == 'BLACKMANHARRIS'):
        ind=(-(N-1)/2,(N-1)/2+1)*2.0*np.pi/N
        h= +0.35875 - 0.48829*np.cos(ind) + 0.14128*np.cos(2.0*ind) - 0.01168*np.cos(3.0*ind)
    elif (name == 'HARRIS'):
        ind=np.arange(1,N+1)*2.0*np.pi/(N+1)
        h=+0.35875 - 0.48829 *np.cos(ind)+0.14128 *np.cos(2.0*ind)-0.01168 *np.cos(3.0*ind)
    elif (name == 'BARTLETT') or (name == 'TRIANG'):
        h=2.0*np.minimum(np.arange(1,N+1),np.arange(N,0,-1))/(N+1)
    elif (name=='BARTHANN'):
        h=  0.38 * (1.0-np.cos(2.0*np.pi*np.arange(1,N+1)/(N+1)))+ 0.48 * np.minimum(np.arange(1,N+1),np.arange(N,0,-1))/(N+1)
    elif (name=='PAPOULIS'):
        ind=np.arange(1,N+1)*np.pi/(N+1)
        h=np.sin(ind)
    elif (name=='GAUSS'):
        if (nargin==3): K=param
        else: K=0.005       
        h= np.exp(np.log(K) * np.power(np.linspace(-1,1,N),2) )
    elif (name == 'PARZEN'):
        ind=np.abs(np.arange(-(N-1)/2,(N-1)/2+1))*2/N
        temp=2*np.power((1.0-ind),3)
        h= np.minimum(temp-np.power((1-2.0*ind),3),temp)
    elif (name == 'HANNA'):
        if (nargin==3): L=param
        else: L=1
        ind=np.arange(0,N)
        h=np.power(np.sin((2*ind+1)*np.pi/(2*N)),(2*L))
    elif (name == 'DOLPH') or (name == 'DOLF'):
        if ((N%2)==0): oddN=1; N=2*N+1;
        else: oddN=0
        if (nargin==3): A=np.power(10,(param/20))
        else: A=1e-3
        
        K=N-1; Z0=np.cosh(np.arccosh(1.0/A)/K); x0=np.arccos(1/Z0)/np.pi; x=np.arange(0,K+1)/N; 
        indices1=np.where(((x<x0)+(x>1-x0)))
        indices2=np.where((x>=x0)*(x<=1-x0))
        h = np.zeros(len(x))+0j
        h[indices1]= np.cosh(K*np.arccosh(0j+Z0*np.cos(np.pi*x[indices1])))
        h[indices2]= np.cos(K*np.arccos(Z0*np.cos(np.pi*x[indices2])))
        h=np.fft.fftshift(np.real(np.fft.ifft(A*np.real(h)))); h=h/h[int(K/2)];
        if oddN: h=h[np.arange(1,K,2)]
    elif (name =='NUTBESS'):
        if (nargin==3): beta=param; nu=0.5 
        elif (nargin==4): beta=param; nu=param2
        else: beta=3*np.pi; nu=0.5
 
        ind=np.arange(-(N-1)/2,(N-1)/2+1) *2/N 
        h=np.sqrt(1-ind**2)**nu*np.real(special.jv(nu,1j*beta*np.sqrt(1.0-ind**2)))/np.real(special.jv(nu,1j*beta))
    elif (name == 'SPLINE'):
        if (nargin < 3):
            print('Three or four parameters required for spline windows');return
        elif (nargin==3):
            nfreq=param; p=np.pi*N*nfreq/10.0
        else: nfreq=param; p=param2
        
        ind=np.arange(-(N-1)/2,(N-1)/2+1) 
        h=np.sinc((0.5*nfreq/p)*ind)**p
    elif (name=='FLATTOP'):
        ind=np.arange(-(N-1)/2,(N-1)/2+1)*2.0*np.pi/(N-1)
        h=+0.2810639 + 0.5208972*np.cos(ind) +0.1980399*np.cos(2.0*ind)
    # Poisson window
    elif (name == 'POISSON'):
        if (nargin==3): K=param
        else: K=4
        ind=abs(np.arange(-(N-1)/2,(N-1)/2+1))*2/N
        h=np.exp(-K*ind)
    # Hanning-Poisson window
    elif (name=='HANNINGPOISSON'):
        if (nargin==3): K=param
        else: K=2
        ind=np.arange(-(N-1)/2,(N-1)/2+1)*2/N
        h=0.50*(1+ np.cos(np.pi*ind))*np.exp(-K*abs(ind))
    # Cauchy window
    elif (name == 'CAUCHY'):
        if (nargin==3): K=param 
        else: K=4
        ind=np.arange(-(N-1)/2,(N-1)/2+1)*2/N
        h=1/(1+abs(K*ind)**2)
    else: print('unknown window name') 
    return h


def STFT_IFD_fast(x, alpha, Hop, h, Dh):
    # Synchrosqueezing by Li Su, 2015
    
        
    import numpy as np
    from numpy import matlib
    import math
    
    
    def half(n):
        half_n = round(n/2)
        if (n%2) == 1: half_n = math.ceil(n/2)
        return half_n
        
    
    if len(h.shape) == 1:
        h = np.array([h])
    if len(Dh.shape) == 1:
        Dh = np.array([Dh])
        
    hrow, hcol = h.shape
    Dhrow, Dhcol = Dh.shape
    
    if (hcol > hrow):
        h = h.T;hrow, hcol = h.shape
    if (Dhcol > Dhrow):
        Dh = Dh.T; Dhrow, Dhcol = Dh.shape
    

	# for tfr
    N = len(np.arange(-0.5+alpha,0.5+alpha,alpha))
    Win_length = np.max(h.shape)
    TH = 7*N/hrow
    tfrtic = np.linspace(0, 0.5, half(N))
    # for tfrsq
    # Lidx = ceil( (N/2)*(lowFreq/0.5) ) ; 
    # Hidx = floor( (N/2)*(highFreq/0.5) ) ; 
    # fLen = Hidx - Lidx + 1 ;

    Overlap = Win_length-Hop
    x_Frame = buffer(x, Win_length, Overlap)
    x_Frame = x_Frame[:,math.ceil(Overlap/Hop/2):]
    x_Frame = x_Frame - matlib.repmat(np.mean(x_Frame,axis=0),x_Frame.shape[0], 1)
    Stime = int(Hop-half(Win_length)+1+math.ceil(Overlap/Hop/2)*Hop)

    tfr = x_Frame*matlib.repmat(h, 1, x_Frame.shape[1]) 	# for h
    tfr = np.fft.fft(tfr, N, 0)
    tfr = tfr[0:half(N),:]

    tf2 = x_Frame*matlib.repmat(Dh, 1, x_Frame.shape[1])	# for h
    tf2 = np.fft.fft(tf2, N, 0)
    tf2 = tf2[0:half(N),:]

	 #get the first order omega
    omega = np.zeros(tf2.shape)
    avoid_warn = np.where(tfr!=0)
    omega[avoid_warn] = np.imag(N*tf2[avoid_warn]/tfr[avoid_warn]/(2.0*np.pi))
    ifd = omega
    """
    % omega = round(omega);
    % % omega(abs(omega)>TH/2)=0;
    % 
    % OrigIndex = repmat((1:round(N/2))', [1 size(tfr,2)]);
    % omega(OrigIndex - omega < 1 | OrigIndex - omega > round(N/2))=0;
    % % ReasIndex = OrigIndex - omega;
    % 
    % Ex = mean(abs(x).^2);
    % Threshold = 1.0e-8*Ex;	% originally it was 1e-6*Ex
    % tfr(abs(tfr) < Threshold) = 0;
    % 
    % totLength = size(tfr,1)*size(tfr,2);
    % rtfr = accumarray((1:totLength)'-omega(:),tfr(:));
    % rtfr = [rtfr; zeros(totLength-length(rtfr),1)];
    % 
    % rtfr = reshape(rtfr,size(tfr,1),size(tfr,2));
    """
    return tfr, ifd, tfrtic, Stime
def  STFT_IFD_lpc_fast(x, alpha, Hop, h, Dh, lpc_num):
    # Synchrosqueezing by Li Su, 2015
    import numpy as np
    from numpy import matlib
    import math
    from scipy import signal
    
    
    def half(n):
        half_n = round(n/2)
        if (n%2) == 1: half_n = math.ceil(n/2)
        return half_n
        
    
    if len(h.shape) == 1:
        h = np.array([h])
    if len(Dh.shape) == 1:
        Dh = np.array([Dh])

    hrow, hcol = h.shape
    Dhrow, Dhcol = Dh.shape
    
    if (hcol > hrow):
        h = h.T;hrow, hcol = h.shape
    if (Dhcol > Dhrow):
        Dh = Dh.T; Dhrow, Dhcol = Dh.shape
    # for tfr
    N = len(np.arange(-0.5+alpha,0.5+alpha,alpha))
    Win_length = np.max(h.shape)
    TH = 7*N/hrow
    tfrtic = np.linspace(0, 0.5, half(N))
    # for tfrsq
    # Lidx = ceil( (N/2)*(lowFreq/0.5) ) 
    # Hidx = floor( (N/2)*(highFreq/0.5) )  
    # fLen = Hidx - Lidx + 1 



    # for tfrsq
    # Lidx = ceil( (N/2)*(lowFreq/0.5) ) 
    # Hidx = floor( (N/2)*(highFreq/0.5) ) 
    # fLen = Hidx - Lidx + 1 
    Overlap = Win_length-Hop
    x_Frame = buffer(x, Win_length, Overlap)
    x_Frame = x_Frame[:,math.ceil(Overlap/Hop/2):]
    x_Frame = x_Frame - matlib.repmat(np.mean(x_Frame,axis=0),x_Frame.shape[0], 1)
    Stime = int(Hop-half(Win_length)+1+math.ceil(Overlap/Hop/2)*Hop)

    for xi in range(0,x_Frame.shape[1]):
        a = lpc_ref(x_Frame[:,xi],lpc_num)
        est_x = signal.lfilter(a,1,x_Frame[:,xi])
        x_Frame[:,xi] = est_x

    
    tfr = x_Frame*matlib.repmat(h, 1, x_Frame.shape[1]) 	# for h
    tfr = np.fft.fft(tfr, N, 0)
    tfr = tfr[0:half(N),:]

    tf2 = x_Frame*matlib.repmat(Dh, 1, x_Frame.shape[1])	# for h
    tf2 = np.fft.fft(tf2, N, 0)
    tf2 = tf2[0:half(N),:]



	 # get the first order omega
    omega = np.zeros(tf2.shape)
    avoid_warn = np.where(tfr!=0)
    omega[avoid_warn] = np.imag(N*tf2[avoid_warn]/tfr[avoid_warn]/(2.0*np.pi))
    ifd = omega
    """
    % omega = round(omega);
    % % omega(abs(omega)>TH/2)=0;
    % 
    % OrigIndex = repmat((1:round(N/2))', [1 size(tfr,2)]);
    % omega(OrigIndex - omega < 1 | OrigIndex - omega > round(N/2))=0;
    % % ReasIndex = OrigIndex - omega;
    % 
    % Ex = mean(abs(x).^2);
    % Threshold = 1.0e-8*Ex;	% originally it was 1e-6*Ex
    % tfr(abs(tfr) < Threshold) = 0;
    % 
    % totLength = size(tfr,1)*size(tfr,2);
    % rtfr = accumarray((1:totLength)'-omega(:),tfr(:));
    % rtfr = [rtfr; zeros(totLength-length(rtfr),1)];
    % 
    % rtfr = reshape(rtfr,size(tfr,1),size(tfr,2));
    """
    return tfr, ifd, tfrtic, Stime

def cepstrum_convert(tfr, tfrtic, g, fs, Tc, num_s, HighFreq, LowFreq):
    import numpy as np
    from scipy import signal, interpolate
    from copy import deepcopy


    if (g!=0):
        ceps = np.real(np.fft.ifft(abs(tfr)**g,2*tfr.shape[0],0))
    else:
        ceps = np.real(np.fft.ifft(np.log(abs(tfr)),2*tfr.shape[0],0))

    for mi in range(0,ceps.shape[1]):
        tra = np.min(np.where(ceps[2:rnd(1/LowFreq),mi]<0))-1+3
        if tra>=0:
            ceps[0:np.max([rnd(1/HighFreq), tra+1]),mi]=0
    

    """
    ceps[0:round(1/HighFreq),:]=0
    % ceps(ceps<Tc)=0; %ceps(ceps>0)=1;
    % % 
    % if num_s>1
    % for kk = 2:num_s
    %     ceps_temp = zeros(size(ceps));
    %     ceps_size = length(kk+1:kk:size(ceps,1));
    %     ceps_temp(2:ceps_size+1,:) = ceps(kk+1:kk:size(ceps,1),:);
    %     ceps = ceps.*ceps_temp;
    % end
    % end
    """
    UpSample = 20

    ceps = ceps[0:rnd(1/LowFreq),:]
    ceps0 = deepcopy(ceps)
    f =  interpolate.interp1d(np.arange(1,ceps.shape[0]+1), ceps,axis=0)
    ceps = f(np.append(np.arange(1,ceps.shape[0],1/UpSample),[ceps.shape[0]]))
    tceps = np.zeros((len(tfrtic), ceps.shape[1]))
    # cepstrum quefrency scale
    freq_scale = UpSample*fs/(np.arange(1,ceps.shape[0]))
    empt = list([])


    for ii in range(1,len(tfrtic)-1):
        # index in quefency
        p_index = np.where((freq_scale > (tfrtic[ii-1]+tfrtic[ii])*fs/2) * (freq_scale <= (tfrtic[ii+1]+tfrtic[ii])*fs/2))
        if len(p_index[0])==0:
            empt.append(ii)
        else:
            #tceps(ii,:)=sum(ceps(p_index,:),1)
            weight = abs(1/p_index[0])
            tceps[ii,:]=np.sum(np.dot(np.diag(weight),ceps[p_index[0],:]),axis=0)


    tceps[tceps<Tc]=0
    tceps_new = deepcopy(tceps)

    # how many peaks in cepstrum (num_s)
    if num_s>1:
        
        for fi in range( rnd(LowFreq*fs*num_s),tceps.shape[0]):
            for kk in range( 2,num_s+1):
                tceps_new[fi,:]=tceps[fi,:]*tceps[max(1,rnd((fi+1)/kk))-1,:]
        
        # for kk = 2:num_s
        #     tceps_temp = zeros(size(tceps))
        #     tceps_size = length(kk+1:kk:size(tceps,1))
        #     tceps_temp(2:ceps_size+1,:) = ceps(kk+1:kk:size(ceps,1),:)
        #     ceps = ceps+ceps_temp
        # end

    tceps = deepcopy(tceps_new)
    """
    % ceps_low = interp1(2:50, ceps(2:50,:), 2:0.001:50)
    % freq_scale_low=1000.*fs./(1:49)
    % for ii=2:length(tfrtic)-1
    %     p_index = find(freq_scale_low > (tfrtic(ii-1)+tfrtic(ii))*fs/2 & freq_scale_low < (tfrtic(ii+1)+tfrtic(ii))*fs/2)
    %     if isempty(p_index)
    % %         empt = [empt; ii]
    %     else
    %         tceps(ii,:)=sum(ceps_low(p_index,:),1);%./(ii)
    %     end
    % end
    % 
    % ceps_mid = interp1(51:200, ceps(51:200,:), 51:0.01:200)
    % freq_scale_mid = 100.*fs./(50:199)
    % for ii=2:length(tfrtic)-1
    %     p_index = find(freq_scale_mid > (tfrtic(ii-1)+tfrtic(ii))*fs/2 & freq_scale_mid < (tfrtic(ii+1)+tfrtic(ii))*fs/2)
    %     if isempty(p_index)
    % %         empt = [empt; ii]
    %     else
    %         tceps(ii,:)=sum(ceps_mid(p_index,:),1);%./(ii)
    %     end
    % end
    % 
    % 
    % ceps_high = interp1(201:size(ceps,1), ceps(201:end,:), 201:0.1:size(ceps,1))
    % freq_scale_high = 10.*fs./(200:size(ceps,1)-1)
    % for ii=2:length(tfrtic)-1
    %     p_index = find(freq_scale_high > (tfrtic(ii-1)+tfrtic(ii))*fs/2 & freq_scale_high < (tfrtic(ii+1)+tfrtic(ii))*fs/2)
    %     if isempty(p_index)
    % %         empt = [empt; ii]
    %     else
    %         tceps(ii,:)=sum(ceps_high(p_index,:),1);%./(ii)
    %     end
    % end
    %tceps(tceps>0)=1
    """
    return ceps0, tceps

def  STFT_IFD_lpc_fast(x, alpha, Hop, h, Dh, lpc_num):
    # Synchrosqueezing by Li Su, 2015
    import numpy as np
    from numpy import matlib
    import math
    from scipy import signal
    
    
    def half(n):
        half_n = round(n/2)
        if (n%2) == 1: half_n = math.ceil(n/2)
        return half_n
        
    
    if len(h.shape) == 1:
        h = np.array([h])
    if len(Dh.shape) == 1:
        Dh = np.array([Dh])

    hrow, hcol = h.shape
    Dhrow, Dhcol = Dh.shape
    
    if (hcol > hrow):
        h = h.T;hrow, hcol = h.shape
    if (Dhcol > Dhrow):
        Dh = Dh.T; Dhrow, Dhcol = Dh.shape
    # for tfr
    N = len(np.arange(-0.5+alpha,0.5+alpha,alpha))
    Win_length = np.max(h.shape)
    TH = 7*N/hrow
    tfrtic = np.linspace(0, 0.5, half(N))
    # for tfrsq
    # Lidx = ceil( (N/2)*(lowFreq/0.5) ) 
    # Hidx = floor( (N/2)*(highFreq/0.5) )  
    # fLen = Hidx - Lidx + 1 



    # for tfrsq
    # Lidx = ceil( (N/2)*(lowFreq/0.5) ) 
    # Hidx = floor( (N/2)*(highFreq/0.5) ) 
    # fLen = Hidx - Lidx + 1 
    Overlap = Win_length-Hop
    x_Frame = buffer(x, Win_length, Overlap)
    x_Frame = x_Frame[:,math.ceil(Overlap/Hop/2):]
    x_Frame = x_Frame - matlib.repmat(np.mean(x_Frame,axis=0),x_Frame.shape[0], 1)
    Stime = int(Hop-half(Win_length)+1+math.ceil(Overlap/Hop/2)*Hop)

    for xi in range(0,x_Frame.shape[1]):
        a = lpc_ref(x_Frame[:,xi],lpc_num)
        est_x = signal.lfilter(a,1,x_Frame[:,xi])
        x_Frame[:,xi] = est_x

    
    tfr = x_Frame*matlib.repmat(h, 1, x_Frame.shape[1]) 	# for h
    tfr = np.fft.fft(tfr, N, 0)
    tfr = tfr[0:half(N),:]

    tf2 = x_Frame*matlib.repmat(Dh, 1, x_Frame.shape[1])	# for h
    tf2 = np.fft.fft(tf2, N, 0)
    tf2 = tf2[0:half(N),:]



	 # get the first order omega
    omega = np.zeros(tf2.shape)
    avoid_warn = np.where(tfr!=0)
    omega[avoid_warn] = np.imag(N*tf2[avoid_warn]/tfr[avoid_warn]/(2.0*np.pi))
    ifd = omega
    """
    % omega = round(omega);
    % % omega(abs(omega)>TH/2)=0;
    % 
    % OrigIndex = repmat((1:round(N/2))', [1 size(tfr,2)]);
    % omega(OrigIndex - omega < 1 | OrigIndex - omega > round(N/2))=0;
    % % ReasIndex = OrigIndex - omega;
    % 
    % Ex = mean(abs(x).^2);
    % Threshold = 1.0e-8*Ex;	% originally it was 1e-6*Ex
    % tfr(abs(tfr) < Threshold) = 0;
    % 
    % totLength = size(tfr,1)*size(tfr,2);
    % rtfr = accumarray((1:totLength)'-omega(:),tfr(:));
    % rtfr = [rtfr; zeros(totLength-length(rtfr),1)];
    % 
    % rtfr = reshape(rtfr,size(tfr,1),size(tfr,2));
    """
    return tfr, ifd, tfrtic, Stime

def synchrosqueeze1win(tfr, ifd, alpha, h, num_s, fr, HighFreq, fs, ths, *args):
    import numpy as np
    import math
    from numpy import matlib
    from scipy import signal
    Smooth = 0  
    Reject = 0
    for var_i in range( 0,len(args)):
        if (args[var_i]== 'Smooth'):
            Smooth = args[var_i + 1]
        if (args[var_i] == 'Reject'):
            Reject = args[var_i + 1]


    M, N = tfr.shape
    K = len(np.arange(-0.5+alpha,0.5+alpha,alpha))
    if Reject>0:
        TH = Reject/fr #*K/length(h); % should be M
        omega = ifd
        if num_s ==1:
            tfr[abs(omega)>TH/2]=0
        else:
            tfr[abs(omega)>TH/2]=0
            for kk in range (2,num_s+1):
                omega_temp = np.empty(omega.shape)
                omega_temp[:] = np.inf
                omega_size = len(np.arange(kk+1,M+1,kk))
                omega_temp[1:omega_size+1,:] = omega[np.arange(kk,M,kk),:]
                tfr[abs(omega_temp)>TH/2]=0
        
    
    else:
        omega = ifd

    tfr = tfr[0:rnd(HighFreq*fs/fr),:]
    omega = np.round(omega[0:rnd(HighFreq*fs/fr),:])

    M, N = tfr.shape
    OrigIndex = matlib.repmat(np.arange(1,M+1),  tfr.shape[1],1).T
    omega[(OrigIndex - omega < 1+2*Smooth) + (OrigIndex - omega > M-2*Smooth)]=0

    Ex = np.mean(sum(abs(tfr)))
    Threshold = ths*Ex	# originally it was 1e-6*Ex
    tfr[abs(tfr) < Threshold] = 0

    totLength = tfr.shape[0]*tfr.shape[1]
    new_idx = np.arange(1,totLength+1)-omega.flatten(order='F')
    new_idx = new_idx.astype(int)
    
    if Smooth == 0:
        rtfr = np.bincount(np.concatenate(([1], new_idx[1:totLength-1], [totLength]))-1,tfr.flatten(order='F'))
    else:
        SmoothWin = signal.triang(1+2*Smooth)/sum(signal.triang(1+2*Smooth))
        rtfr = np.bincount(np.concatenate(([1], new_idx[1:totLength-1], [totLength]))-1, SmoothWin[0]*tfr.flatten(order='F'))
        for ii in range( 1,Smooth):
            rtfr = rtfr + np.bincount(np.concatenate((np.arange(1,ii+1), new_idx[ii:totLength-ii]-ii, np.arange(totLength-ii+1,totLength+1)))-1, SmoothWin[ii]*tfr.flatten(order='F'))
            rtfr = rtfr + np.bincount(np.concatenate((np.arange(1,ii+1), new_idx[ii:totLength-ii]+ii, np.arange(totLength-ii+1,totLength+1)))-1, SmoothWin[ii]*tfr.flatten(order='F'))
        
    rtfr = np.concatenate((rtfr, np.zeros((totLength-len(rtfr)))))

    rtfr = np.reshape(rtfr,(M,N), order='F')
    return tfr, rtfr 