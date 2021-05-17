function cost=two_Lorentzian_cost(zeta,spectral_corr,params)
a1=params(1);
a2=params(2);
E01=params(3);
E02=params(4);
gamma1=params(5);
gamma2=params(6);
c=params(7);

[zeta,spectral_corr_model,energy_vector,lineshape]=two_Lorentzians(zeta, a1,a2,E01,E02,gamma1,gamma2,c);
beep = spectral_corr(2:end);
cost=sum((spectral_corr_model-beep.^2));
end