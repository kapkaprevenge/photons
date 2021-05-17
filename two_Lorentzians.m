function [zeta_in,spectral_corr,energy_vector,lineshape]=two_Lorentzians(zeta_in, a1,a2,E01,E02,gamma1,gamma2,c)
n=(length(zeta_in)-1)/2;

delta=zeta_in(2)-zeta_in(1); % energy difference. 
energy_vector=[(-n/2*delta):delta:n/2*delta];

gamma2=0;
Lorentz1=1./((energy_vector-E01).^2+(0.5*gamma1).^2);
Lorentz2=1./((energy_vector-E02).^2+(0.5*gamma1).^2);

% gamma is the FWHM. 
lineshape=a1*Lorentz1+a2*Lorentz2;

spectral_corr=xcorr(lineshape,lineshape);
spectral_corr=spectral_corr/max(spectral_corr)+c;
spectral_corr=spectral_corr/max(spectral_corr);

end