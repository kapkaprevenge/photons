%% rename the files that have small letters at the begining
% f=dir('*.stream')
% for id=1:length(f)
% [~,fid]=fileparts(f(id).name);
% if strcmp(fid(1:4),'dotb')==1;
%     new_file=strcat('DotB',fid(4:end),'.stream')
%     movefile(f(id).name,new_file)
% end
% end

% path = ('/Users/andrewproppe/Bawendi_Lab/PCFS/Dots/2020/dotB');
% data_path = ('/Users/andrewproppe/Bawendi_Lab/PCFS/PCFS_for_practice original/DotPCFS');
% data_path = ('/Users/andrewproppe/Bawendi_Lab/PCFS/PCFS_practice two');
data_path = ('/Users/andrewproppe/Bawendi_Lab/PCFS/Dots/2020/dotB');

% DotPCFS = PCFS1('E:/Dropbox (MIT)/codes/PCFS_for_practice/DotPCFS',2E6,'t2')
% DotPCFS = PCFS1(data_path, 2E6, 't2'); % folder, buffer size, mode (t2 or t3)
load('DotPCFS_analyzed.mat');
DotPCFS = DotPCFS;
% tic
% DotPCFS.get_photons_all()  
% toc
% DotPCFS.get_sum_signal_all()
% tic
% DotPCFS.get_intensity_correlations(2,[1E5,1E12],5);
% toc
DotPCFS.get_blinking_corrected_PCFS_interferogram();

DotPCFS.plot_interferogram([1E8,1E11]);
% 
fringe_index = 2; % you must manually set this parameter

if ~isprop(DotPCFS,'white_fringe')
    dummy=addprop(obj,'white_fringe')
end
%semi-automate getting the photon white fringe
DotPCFS.white_fringe.index = fringe_index;
DotPCFS.white_fringe.pos = DotPCFS.stage_positions(fringe_index);
%% analyze the spectral diffusion. 
% need a function that plots the interferogram decay and the spectral
% correlations as a function of tau.  

DotPCFS.get_blinking_corrected_PCFS_interferogram()
% 
DotPCFS.plot_spectral_diffusion([300E6,1000E6],DotPCFS.stage_positions(2));

DotPCFS.get_low_T_spectral_correlation(DotPCFS.white_fringe.pos,DotPCFS.white_fringe.index);
DotPCFS.plot_low_T_spectral_corr([50E6,300E6,10000E6,100000E6],[-3,3])
%THIS NOW GIVES US Tau Indices so we can find the appropriate spectral
%corrs below

%% Data analysis (framework prototype by Alex)
t_index = DotPCFS.tau_indices(2); %from plot_low_T_spectral_corr
zeta_bounds = sort([-1,1]); %set this here, automates the bounds below. in meV.
zeta_ind = zeros(length(zeta_bounds),1);
for i=1:length(zeta_bounds)
    [dummy,index] = min((DotPCFS.mirrored_intf_spectral_corr.zeta - zeta_bounds).^2);
    zeta_ind(i) = index;
end

mean_spec_corr = mean(DotPCFS.mirrored_intf_spectral_corr.spectral_corr(t_index-5:t_index+5,:));
mean_spec_corr = abs(mean_spec_corr./max(mean_spec_corr));
mean_spec_corr = mean_spec_corr(zeta_ind(1):zeta_ind(2));
%curates an appropriate avg'd corr around the tau desired, and focused on
%the zeta range set above. Test to make sure it looks good below
figure()
plot(DotPCFS.mirrored_intf_spectral_corr.zeta,mean_spec_corr,'Linewidth',3)

%if all looks good create an object to house this data (along with the
%fitting)
if ~isprop(DotPCFS,'spectral_correlation')
    dummy=addprop(obj,'spectral_correlation')
end

DotPCFS.spectral_correlation.corr = mean_spec_corr;
DotPCFS.spectral_correlation.zeta = DotPCFS.mirrored_intf_spectral_corr.zeta(zeta_ind(1):zeta_ind(2))

num_Lor = 2; %how many lorentzians should this be fit by?
FSS = [];
FSS(1) = 400/1000; %you have to manually approx the FSS, in meV, in the future you might have more than one.

DotPCFS.fit_spectral_corr(FSS,num_Lor); %this fit saves a lot of space compared to the old paradigm, check it out in PCFS.m!
%after running this, "spectral correlation" object should have the zeta,
%the correlation, and the fitting parameters.

%I downsized the amount of space this took up but in the process haven't
%completely overhauled everything. this framework is easy to expand and
%copy, however.
%%%%%



