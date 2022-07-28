function  [fy,f]=PSD(y,Fs)

%==========================================================================
%
%  Usage :  
%  Function  [fy]=PSD(y,Fs)
%
%  1)computes the Power spectral density and Amplitude spectrum (P(f),F(f))
%  of 1d signal y(t) with sample rate  Fs (Nyquist rate) which is known
%  apriori. The results are plotted in 3 figures which correspond to simple
%  PSD,logarithmic PSD (dB) and  Amplitude Specturm respectively.
%                             _______________
%  such that  Ampitude(f) = \/    PSD(f)
%  
%  2)The usefulness of this function is the adjustment of the frequency axis
%
%  3)The fast Fourier transform is computed with Matlab built-in function
%  fft, but for signals whose lengths <1000 points, one can use the nested
%  function   y=Fast_Fourier_Transform(X,N) .
%
% (c) KHMOU Youssef , Signal Processing 2013
%==========================================================================

% In case that the input vector is matrix :  Maping with vect{} .
y=y(:).';
L=length(y);

%  (2^N) :Number of points for computing the FFT 
N=ceil(log2(length(y)));

% FFT 
fy=fft(y,2^N)/(L/2);
%------------------------------------------------
% for length<1000 one can replace fft with function :
% fy=Fast_Fourier_Transform(y,2^N)/(L/2); (line 84)
%------------------------------------------------

% Amplitude adjustment by checking for complex input y 
if isreal(y)==0
    fy=fy/2;
end

% PSD
Power=fy.*conj(fy);
%Phase Angle
phy=angle(fy);

%  Frequency axis
f=(Fs/2^N)*(0:2^(N-1)-1);