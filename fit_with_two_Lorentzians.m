% figure()
% for i=1:2:30
% zeta=DotD.mirrored_intf_spectral_corr.zeta;
% spectral_corr=mean(DotD.mirrored_intf_spectral_corr.spectral_corr(40+i,:),1);
% spectral_corr=spectral_corr/min(spectral_corr);
% plot(zeta,spectral_corr,'color',[0.5,0.5,0.03*i])
% hold on
% end

%% curated means that we take two runs and stitch them together
% to increase the resolution of the measurement. 



mean_spec_corr=mean(DotD.mirrored_intf_spectral_corr.spectral_corr(45:55,:),1);
% mean_spec_corr=abs(mean_spec_corr/min((mean_spec_corr)));
mean_spec_corr = abs(mean_spec_corr);
mean_spec_corr = mean_spec_corr./max(mean_spec_corr);

figure()
plot(DotD.mirrored_intf_spectral_corr.zeta,mean_spec_corr,'Linewidth',3)
xlim([-2,2])


zeta=DotD.mirrored_intf_spectral_corr.zeta;

%% here we define a function that takes the spectral correlation 
% between -1 and 1 meV and fits it with two offset Lorentzians. 

fun = @(params) two_Lorentzian_cost(zeta(935:1140),mean_spec_corr(935:1140),params);
params0=[1,1,0,1.5,0.6,0,2];
lb=[0,0,-3,-3,0.1,0,0];
ub=[200,200,3,3,1,0,1E7];

gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',params0,...
    'objective',fun,'lb',lb,'ub',ub);
optim_params = run(gs,problem);





%% fitting the PCFS interferogram.
% this is with the non-interpolated data which puts 
% strong bias on the early time data of the mdoel.

a1=optim_params(1);
a2=optim_params(2);
E01=optim_params(3);
E02=optim_params(4);
gamma1=optim_params(5);
gamma2=optim_params(6);
c=optim_params(7);


% get path_lentgh_vector
wf_pos=DotD.stage_positions(2);
path_diff_DotD=2*(DotD.stage_positions(2:end)-wf_pos)/10;
interferogram_DotD=mean(DotD.PCFS_interferogram(44:55,2:end),1)/min(mean(DotD.PCFS_interferogram(44:55,2:end),1));

[path_difference_half,interferogram_half,path_length_difference_in,interferogram_out,energy_vector,lineshape]=two_Lorentzians_FFT(path_diff_DotD, a1,a2,E01,E02,gamma1,gamma2,c)

figure()
plot(path_difference_half,interferogram_half)

%interpolate the interferogram first. 
path_diff_DotD_interp=linspace(min(path_diff_DotD),max(path_diff_DotD),100)
interferogram_DotD_interp=interp1(path_diff_DotD,interferogram_DotD,path_diff_DotD_interp)
path_diff_DotD_interp(isnan(interferogram_DotD_interp))=[];
interferogram_DotD_interp(isnan(interferogram_DotD_interp))=[];


fun = @(params) two_Lorentzian_FFT_cost(path_diff_DotD_interp,interferogram_DotD_interp,params);
params0=[optim_params(1:6),0.01];
lb=[0,0,-1.5,-1.5,0,0,0];
ub=[200,200,1.5,1.5,1,0,1E7];

gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',params0,...
    'objective',fun,'lb',lb,'ub',ub);
optim_params_interferogram = run(gs,problem);

%% calculate and plot the fit results. 
a1=optim_params_interferogram(1);
a2=optim_params_interferogram(2);
E01=optim_params_interferogram(3);
E02=optim_params_interferogram(4);
gamma1=optim_params_interferogram(5);
gamma2=optim_params_interferogram(6);
c=optim_params_interferogram(7);
[path_difference_half,interferogram_half,path_length_difference_in,interferogram_out,energy_vector,lineshape,full_path,full_interferogram]=two_Lorentzians_FFT(path_diff_DotD, a1,a2,E01,E02,gamma1,gamma2,c);

% calculate an exponential decay guide for the eye for the coherence
% decay. 
eV2cm=8065.54429;

exp_decay_const=eV2cm/1000*gamma1*2*pi/10
exp_coherence_decay=0.505*exp(-exp_decay_const*path_diff_DotD*10);

% calcualte the sqrt property. 
sqrt_interferogram_out(interferogram_out>0)=sqrt(interferogram_out(interferogram_out>0));
sqrt_interferogram_out(interferogram_out<0)=-sqrt(-interferogram_out(interferogram_out<0));
sqrt_interferogram_DotD(interferogram_DotD<0)=-sqrt(-interferogram_DotD(interferogram_DotD<0));
sqrt_interferogram_DotD(interferogram_DotD>0)=sqrt(interferogram_DotD(interferogram_DotD>0));


%% making the final plot of the results.
zeta_fit=[-6:0.001:6];
[a,b,c,d]=two_Lorentzians(zeta_fit,optim_params(1),optim_params(2),optim_params(3),optim_params(4),optim_params(5),optim_params(6),optim_params(7));

meV2GHz=241.8; % for secondary axis in GHz. 
speed_of_light=299792458000/1E12; % in mm per ps.  

figure('pos',[10,10,700,600])
subplot(2,1,1)
plot(zeta*meV2GHz,mean_spec_corr,'o-','Linewidth',1)
hold on
plot(a*meV2GHz,b,'Linewidth',2)
plot(a,zeros(length(a),1),'--','color','black')
xlim([-1,1]*meV2GHz)
ylim([-0.2,1.2])
box off
xlabel('\zeta [GHz]')
ylabel('p(\zeta) [a.u.]')
set(gca,'fontsize',16)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'xminorgrid','on')
legend('data','two Lorentzian fit')

ax1=gca;
ax1_pos=ax1.Position;
ax2=axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
ax2.XColor='red';
line(a,b,'Parent',ax2,'Color','red')
xlim([-1,1])
ylim([-0.2,1.2])
xlabel('\zeta [meV]')
ylabel('p(\zeta) [a.u.]')
set(gca,'fontsize',16)

subplot(2,1,2)
plot(path_diff_DotD(1:100)*10/speed_of_light,interferogram_DotD(1:100),'o-')
hold on
plot(path_diff_DotD(1:100)*10/speed_of_light,interferogram_out(1:100),'Linewidth',3)
plot(path_diff_DotD(1:100)*10/speed_of_light,exp_coherence_decay(1:100),'--','Linewidth',2,'color','black')
plot(path_diff_DotD(1:100)*10/speed_of_light,zeros(length(path_diff_DotD(1:100)),1),'--','color','black')
box off
set(gca,'xminorgrid','on')

xlim([-0.5,30]/speed_of_light)
xlabel('\delta [ps]')

ylabel('Normalized g^{2}_{cross}-g^{2}_{auto}')
set(gca,'fontsize',16)
ylim([-0.2,1])
legend('data','two Lorentzian fit','exponential component')

ax1=gca;
ax1_pos=ax1.Position;
ax2=axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
ax2.XColor='red';
line(path_diff_DotD*10,interferogram_out,'Parent',ax2,'Color','red')
xlim([-0.5,30])
xlabel('\delta [mm]')
ylabel('Normalized g^{2}_{cross}-g^{2}_{auto}')
set(gca,'fontsize',16)
ylim([-0.2,1])

% 
% subplot(3,1,3)
% plot(path_diff_DotD*10,sqrt_interferogram_DotD,'o-')
% hold on
% plot(path_diff_DotD*10,sqrt_interferogram_out,'Linewidth',3)
% plot(path_diff_DotD*10,zeros(length(path_diff_DotD),1),'--','color','black')
% box off
% set(gca,'xminorgrid','on')
% 
% xlim([-0.5,80])
% xlabel('\delta [mm]')
% ylabel('Normalized g^{2}_{cross}-g^{2}_{auto}')
% set(gca,'fontsize',16)
% ylim([-0.2,1.2])
% 
% ax1=gca;
% ax1_pos=ax1.Position;
% ax2=axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% ax2.XColor='red';
% line(path_diff_DotD*10/speed_of_light,sqrt(interferogram_out),'Parent',ax2,'Color','red')
% xlim([-0.5,80]/speed_of_light)
% xlabel('\delta [ps]')
% ylabel('Normalized g^{2}_{cross}-g^{2}_{auto}')
% set(gca,'fontsize',16)
% ylim([-0.2,1.2])
%%c=linspace(-2,2,1000)

%% a plot of the underlying spectral lineshape. 
figure()
plot(c-optim_params(4),d/max(d),'Linewidth',2)
xlim([-1,1])
title('DotD reconstructred lineshape from fit')
xlabel('Energy -E0 [meV]')
ylabel('intensity [a.u.]')
set(gca,'fontsize',14)
% 

%% make a plot of the temporal evolution of DotD. 
figure()
plot(DotD.mirrored_intf_spectral_corr.zeta,DotD.mirrored_intf_spectral_corr.spectral_corr(40,:)/min(DotD.mirrored_intf_spectral_corr.spectral_corr(40,:)),'o-','Linewidth',2)
hold on
%plot(DotD.mirrored_intf_spectral_corr.zeta,DotD.mirrored_intf_spectral_corr.spectral_corr(50,:)/min(DotD.mirrored_intf_spectral_corr.spectral_corr(50,:)),'Linewidth',2)
plot(DotD.mirrored_intf_spectral_corr.zeta,DotD.mirrored_intf_spectral_corr.spectral_corr(60,:)/min(DotD.mirrored_intf_spectral_corr.spectral_corr(60,:)),'o-','Linewidth',2)
plot(DotD.mirrored_intf_spectral_corr.zeta,DotD.mirrored_intf_spectral_corr.spectral_corr(81,:)/min(DotD.mirrored_intf_spectral_corr.spectral_corr(81,:)),'o-','Linewidth',2)
plot(DotD.mirrored_intf_spectral_corr.zeta,DotD.mirrored_intf_spectral_corr.spectral_corr(92,:)/min(DotD.mirrored_intf_spectral_corr.spectral_corr(92,:)),'o-','Linewidth',2)
plot(DotD.mirrored_intf_spectral_corr.zeta,DotD.mirrored_intf_spectral_corr.spectral_corr(100,:)/min(DotD.mirrored_intf_spectral_corr.spectral_corr(100,:)),'o-','Linewidth',2)
%plot(DotD.mirrored_intf_spectral_corr.zeta,zeros(length(DotD.mirrored_intf_spectral_corr.zeta),1),'--','color','black')
set(gca,'xminorgrid','on')

xlim([-2,2])
xlabel('\zeta [meV]')
ylabel('p(\zeta) [a.u.]')
legend('\tau = 0.03 ms','\tau = 0.4 ms','\tau = 8 ms','\tau = 36 ms','\tau = 115 ms')
set(gca,'fontsize',16)
ylim([0,1])



% 
% %%% not relevant because it is clear we have to use the curated dot. 
% %% doing the same with the uncurated DotD data. 
% zeta=DotD_4K_run_one.mirrored_intf_spectral_corr.zeta;
% spectral_corr=mean(DotD_4K_run_one.mirrored_intf_spectral_corr.spectral_corr(60:70,:),1);
% spectral_corr=spectral_corr/min(spectral_corr);
% 
% fun = @(params) two_Lorentzian_cost(zeta,spectral_corr,params);
% params0=[1,1,0,0.4,0.1,0,0.01];
% lb=[0,0,0,0,0,0,0];
% ub=[2,2,0.3,20,1,0,0.2];
% 
% gs = GlobalSearch;
% problem = createOptimProblem('fmincon','x0',params0,...
%     'objective',fun,'lb',lb,'ub',ub);
% optim_params = run(gs,problem);

%% plotting the results;
[a,b]=two_Lorentzians(zeta,optim_params(1),optim_params(2),optim_params(3),optim_params(4),optim_params(5),optim_params(6),optim_params(7));

figure()
plot(zeta,mean_spec_corr,'o-','Linewidth',2)
hold on
plot(zeta,b,'Linewidth',2)
xlim([-2,2])
title('Cleaned DotD Spectral Correlation and Lorentzian Fit')
set(gca,'fontsize',14)

Gamma_optim_uncurated=optim_params(5)





