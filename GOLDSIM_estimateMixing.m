clc
clear variables
close all

rng('default')
%%  Load optimized data
load("GOLDSIM_optimal.mat")
times = compare_data.times;

%%  Solve at idealized conditions
ic = [0.0001, 0.0003];
plot_reduct = 0.05;
[sol, mySettings] = GOLDSIM_simulateGoldParticles(prmGlobalOpt,plot_reduct, ic);
[sim_diams, sim_PSDs] = PBElib_getPSDs(sol,mySettings,times);
[n_diams, n_times] = size(sim_PSDs);

ideal_diams = sim_diams;
ideal_PSDs = sim_PSDs;

%%  Vary the initial conditions
n_sims = 1000;
alldiams = zeros(1+n_sims, n_diams, n_times);
alldiams(1,:,:) = sim_diams;

allPSDs = zeros(1+n_sims, n_diams, n_times);
allPSDs(1,:,:) = sim_PSDs;

%   Geometric std
sigma = log(10);
%   Mean values
mu1 = log(0.0001);
mu2 = log(0.0003);

tic
parfor iii=1:n_sims
    %   perturb the initial conditions
    ic1 = lognrnd(mu1,sigma);
    ic2 = lognrnd(mu2,sigma);
    ic = [ic1,ic2];
    %   simulate
    
    [sol, mySettings] = GOLDSIM_simulateGoldParticles(prmGlobalOpt,plot_reduct, ic);
    [sim_diams, sim_PSDs] = PBElib_getPSDs(sol,mySettings,times);
    alldiams(1+iii,:,:) = sim_diams;
    allPSDs(1+iii,:,:) = sim_PSDs;
end
toc

%%
avg_psd = squeeze(mean(allPSDs,1));

%%
%--------------------------------------------------------------------------
%   Plot PSDs over time
%--------------------------------------------------------------------------
figure
ax = gca;
h = animatedline('LineStyle','--',...
    'Marker','none',...
    'LineWidth',2,...
    'DisplayName',"Ideal Simulation",...
    'Color',[0 0.4470 0.7410]);
h2 = animatedline('LineStyle','--',...
    'Marker','none',...
    'LineWidth',2,...
    'DisplayName',"Data",...
    'Color',[0.8500 0.3250 0.0980]);
h3 = animatedline('LineStyle','--',...
    'Marker','none',...
    'LineWidth',2,...
    'DisplayName',"Mixing Estimate",...
    'Color',[0.9290 0.6940 0.1250]);
axis([0,15,0,Inf])
plotpts = times;
anim_seconds = 5;

diams = ideal_diams;
for iii=1:length(plotpts)
    clearpoints(h)
    clearpoints(h2)
    clearpoints(h3)

    psdideal = ideal_PSDs(:,iii);
    areaideal = trapz(diams(:,iii),psdideal);
    addpoints(h,diams(:,iii), psdideal/areaideal);
    
    psdmixed = avg_psd(:,iii);
    areamixed = trapz(diams(:,iii), psdmixed);
    addpoints(h3,diams(:,iii), psdmixed/areamixed);

    
    
    addpoints(h2,compare_data.diam, compare_data.PSDs(:,iii))


    title(ax,sprintf("PSD at %.f seconds",plotpts(iii)))
    legend('Location','northeast')
    drawnow
    pause(anim_seconds/length(plotpts))
end