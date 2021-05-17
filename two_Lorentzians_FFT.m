function [path_difference_half,interferogram_half,path_length_difference_in,interferogram_out,energy_vector,lineshape,full_path,full_interferogram]=two_Lorentzians_FFT(path_length_difference_in, a1,a2,E01,E02,gamma1,gamma2,c)

% create some fake zeta
zeta=[-100:0.001:100]; % in meV. 

[zeta_in,spectral_corr,energy_vector,lineshape]=two_Lorentzians(zeta, a1,a2,E01,E02,gamma1,gamma2,c);

eV2cm=8065.54429;
cm2eV=1/eV2cm;

N=length(zeta);
delta=(max(zeta)-min(zeta))/N;

%get reciprocal space (wavenumbers).
increment=1/delta;
path_difference=linspace(-0.5*increment,0.5*increment,N)*cm2eV*1000; %converted to meV


%% reversed from the PCFS code -> interferogram to spectral correlation.

%%take the FFT of the interferogram to get the spectral correlation.All
%%that shifting is to shift the zero frequency component to the middle
%%of the FFT vector. We take the real part of the FFT because the
%%interferogram is by definition entirely symmetric.
interferogram = real(fftshift(fft(ifftshift(spectral_corr),N)));
NN=((N-1)/2);%taking the positive path length differences only.
interferogram_half=interferogram(NN:end); 
path_difference_half=path_difference(NN:end);
interferogram_out=interp1(path_difference_half,interferogram_half,path_length_difference_in);
interferogram_out=interferogram_out/max(interferogram_out);

full_path=path_difference;
full_interferogram=interferogram;

end