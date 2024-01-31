clc
clear variables
close all

rng('default')

%%
filename_log = "GOLDSIM_estimateMixing.log";
filename_Data = "DATA_GOLDSIM_estimateMixing.mat";
filename_LoadData = "DATA_GOLDSIM_optimizeParameters_5percent.mat";

save_folder = fullfile('.','figures','mixing');
save_file_root = 'mixing';


%%  Save log
overwrite_log = true;
if overwrite_log
    fileID = fopen(filename_log,'w');
    fclose(fileID);
end
diary(filename_log);



%%  Load optimized data
load(filename_LoadData)
times = compare_data.times;

%%
n_sims = 100;
target_std = 1.5;

%%  Plot the PSD that was optimized or the more accurate version?
plot_accurate = true;
if plot_accurate
    sol = solAccurate;
    mySettings = mySettingsAccurate;
    plot_reduct = plotting_reduct;
else
    sol = solOptimize;
    mySettings = mySettingsOptimize;
    plot_reduct = optimization_reduct;
end

%%  Time how long the entire script takes
tstart = datetime;

%%  Solve at idealized conditions
ic = [0.0001, 0.0003];
[sol, mySettings] = GOLDSIM_simulateGoldParticles(prmGlobalOpt,plot_reduct, ic, 300.0 );
[sim_diams, sim_PSDs] = PBElib_getPSDs(sol,mySettings,times);
[n_diams, n_times] = size(sim_PSDs);

ideal_diams = sim_diams;
ideal_PSDs = sim_PSDs;

%%  Vary the initial conditions
%   Geometric std
sigma = log(target_std);
%   Mean values
mu1 = log(0.0001) - sigma^2/2;  %   gives mean of 0.0001 in lognormal
mu2 = log(0.0003) - sigma^2/2;  %   gives mean of 0.0003 in lognormal

%   Time for quenching to take effect varies
% sigmaMixT = log(1.5);
% muMixT = log(prmGlobalOpt(5)) - sigmaMixT^2/2;
% mixT = prmGlobalOpt(5);

manySims = zeros(1+n_sims, n_diams, n_times);
ic_draws = zeros(2,n_sims);
tic
parfor iii=1:n_sims
    %   perturb the initial conditions
    ic1 = lognrnd(mu1,sigma);
    ic2 = lognrnd(mu2,sigma);
    % mixT = lognrnd(muMixT,sigmaMixT);
    ic = [ic1;ic2];
    ic_draws(:,iii) = ic;
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
% xlabel('NaAuCl_4 / mol/m^{-3}')
% ylabel('NaBH_4 / mol/m^{-3}')
% legend
hold off
filenamePic = fullfile(save_folder, 'perturbations');
print('-depsc','-image',filenamePic)
cleanfigure;
filenameTex = strcat(filenamePic,'.tex');
matlab2tikz(filenameTex)
close all

%%  Average many sims
mixPSDs_type2 = squeeze(mean(manySims,1));


%%
%--------------------------------------------------------------------------
%   Plot PSDs over time
%--------------------------------------------------------------------------

for iii=1:length(times)
    figure
    psdideal = ideal_PSDs(:,iii);
    %   Data is volume weighted
    psdideal = sim_diams(:,iii).^3 .* psdideal;
    areaideal = trapz(sim_diams(:,iii), psdideal);
    psdideal = psdideal/areaideal;
    % cdfideal = cumtrapz(sim_diams(:,iii),psdideal);


    psdmixed2 = mixPSDs_type2(:,iii);
    % psdmixed2 = diams(:,iii).^3 .* psdmixed2;
    areamixed2 = trapz(sim_diams(:,iii), psdmixed2);
    psdmixed = psdmixed2/areamixed2;
    % cdfmixed = cumtrapz(sim_diams(:,iii), psdmixed);

    % cdf_data = cumtrapz(data_diam, data_PSDs(:,iii));

    hold on
    plot(sim_diams(:,iii), psdideal, ...
        'Color',[0 0.4470 0.7410],...
        'LineStyle','--',...
        'LineWidth',2,...
        'DisplayName','Simulated PSD')
    plot(data_diam, data_PSDs(:,iii),...
        'Color',[0.8500 0.3250 0.0980],...
        'LineStyle','--',...
        'LineWidth',2,...
        'DisplayName','Measured PSD')

    plot(sim_diams(:,iii), psdmixed,...
        'Color',[0.6350 0.0780 0.1840],...
        'LineStyle','--',...
        'LineWidth',2,...
        'DisplayName','Simulated PSD Perturbed')


    text(8,0.2,sprintf("%.f seconds",times(iii)))
    box on
    set(gca,'YColor','k')

    hold off
    fontsize(gcf,scale=1.75)
    filenamePic = fullfile(save_folder, strcat(save_file_root,num2str(iii)));
    print('-depsc','-image',filenamePic)
    cleanfigure;
    filenameTex = strcat(filenamePic,'.tex');
    matlab2tikz(filenameTex)
    close all
end





%%  Save workspace data
close all
save(filename_Data)

%%  End time for script
tend = datetime;

dur = tend-tstart;
fprintf("Script time (hh:mm:ss): %s\n",dur)
%%  end of diary
diary off