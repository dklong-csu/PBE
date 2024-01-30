clc
clear variables
% close all

rng('default')

%%  Save log
overwrite_log = true;
if overwrite_log
    fileID = fopen("GOLDSIM_estimateMixing.log",'w');
    fclose(fileID);
end
diary GOLDSIM_estimateMixing.log

%%  Time how long the entire script takes
tstart = datetime;

%%  Load optimized data
load("DATA_GOLDSIM_optimizeParameters.mat")
times = compare_data.times;

%%  Solve at idealized conditions
ic = [0.0001, 0.0003];
plot_reduct = 0.1;
[sol, mySettings] = GOLDSIM_simulateGoldParticles(prmGlobalOpt,plot_reduct, ic, 300.0 );
[sim_diams, sim_PSDs] = PBElib_getPSDs(sol,mySettings,times);
[n_diams, n_times] = size(sim_PSDs);

ideal_diams = sim_diams;
ideal_PSDs = sim_PSDs;

%%  Vary the initial conditions
n_sims = 2;

%   Want to target geometric standard deviation 
target_std = 2;  % = exp(sigma)
%   Geometric std
sigma = log(target_std);
%   Mean values
mu1 = log(0.0001) - sigma^2/2;  %   gives mean of 0.0001 in lognormal
mu2 = log(0.0003) - sigma^2/2;  %   gives mean of 0.0003 in lognormal

%   Time for quenching to take effect varies
sigmaMixT = log(1.5);
% muMixT = log(prmGlobalOpt(5)) - sigmaMixT^2/2;
% mixT = prmGlobalOpt(5);

manySims = zeros(1+n_sims, n_diams, n_times);
ic_draws = zeros(3,n_sims);
tic
parfor iii=1:n_sims
    %   perturb the initial conditions
    ic1 = lognrnd(mu1,sigma);
    ic2 = lognrnd(mu2,sigma);
    % mixT = lognrnd(muMixT,sigmaMixT);
    ic = [ic1;ic2];
    ic_draws(:,iii) = [ic;mixT];
    %   simulate
    
    [sol, mySettings] = GOLDSIM_simulateGoldParticles(prmGlobalOpt,plot_reduct, ic, 300);
    [sim_diams, sim_PSDs] = PBElib_getPSDs(sol,mySettings,times);
    %   Volume weight the PSDs
    sim_PSDs = (sim_diams.^3) .* sim_PSDs;
    %   Normalize to area=1
    area = trapz(sim_diams(:,1),sim_PSDs);
    sim_PSDs = sim_PSDs ./ area;
    manySims(iii,:,:) = sim_PSDs;
end
toc

%%  Scatter plot of ic estimates
figure
scatter(ic_draws(1,:),ic_draws(2,:),20,[0 0.4470 0.7410],"filled",'o','DisplayName','Perturbed Initial Concentration')
hold on
scatter(0.0001, 0.0003,50,[0.8500 0.3250 0.0980],'o','filled','DisplayName','Average Initial Concentration')
xlabel('NaAuCl_4 / mol/m^{-3}')
ylabel('NaBH_4 / mol/m^{-3}')
legend
hold off
%%  Average many sims
mixPSDs_type2 = squeeze(mean(manySims,1));

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
% h3 = animatedline('LineStyle','--',...
%     'Marker','none',...
%     'LineWidth',2,...
%     'DisplayName',"Mixing Estimate",...
%     'Color',[0.9290 0.6940 0.1250]);
h4 = animatedline('LineStyle','--',...
    'Marker','none',...
    'LineWidth',2,...
    'DisplayName',"Mixing Estimate",...
    'Color',[0.4940 0.1840 0.5560]);
axis([0,15,0,1.2])
plotpts = times;
anim_seconds = 15;

diams = ideal_diams;
for iii=1:length(plotpts)
    clearpoints(h)
    clearpoints(h2)
    % clearpoints(h3)
    clearpoints(h4)

    psdideal = ideal_PSDs(:,iii);
    %   Data is volume weighted
    psdideal = diams(:,iii).^3 .* psdideal;
    areaideal = trapz(diams(:,iii), psdideal);
    addpoints(h,diams(:,iii), psdideal/areaideal);
    
    % psdmixed = mixPSDs(:,iii);
    % %   Data is volume weighted
    % psdmixed = diams(:,iii).^3 .* psdmixed;
    % areamixed = trapz(diams(:,iii), psdmixed);
    % addpoints(h3,diams(:,iii), psdmixed/areamixed);

    psdmixed2 = mixPSDs_type2(:,iii);
    % psdmixed2 = diams(:,iii).^3 .* psdmixed2;
    areamixed2 = trapz(diams(:,iii), psdmixed2);
    addpoints(h4, diams(:,iii), psdmixed2/areamixed2);

    
    
    addpoints(h2,compare_data.diam, compare_data.PSDs(:,iii))


    title(ax,sprintf("PSD at %.f seconds",plotpts(iii)))
    legend('Location','northeast')
    drawnow
    pause(anim_seconds/length(plotpts))
end

%%  Save workspace data
save("DATA_GOLDSIM_estimateMixing.mat")

%%  End time for script
tend = datetime;

dur = tend-tstart;
fprintf("Script time (hh:mm:ss): %s\n",dur)
%%  end of diary
diary off