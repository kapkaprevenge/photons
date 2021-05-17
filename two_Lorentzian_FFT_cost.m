function cost=two_Lorentzian_FFT_cost(path_diff,interferogram,params)

a1=params(1);
a2=params(2);
E01=params(3);
E02=params(4);
gamma1=params(5);
gamma2=params(6);
c=params(7);

[path_difference_half,interferogram_half,path_length_difference_in,interferogram_out,energy_vector,lineshape]=two_Lorentzians_FFT(path_diff, a1,a2,E01,E02,gamma1,gamma2,c);


cost=sum((interferogram_out(2:end)-interferogram(2:end)).^2);
%cost=interpolated_cost(path_difference_half, interferogram_half,path_diff,interferogram);

end