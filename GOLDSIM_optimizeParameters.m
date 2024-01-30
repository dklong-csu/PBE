clc
clear variables
close all

rng('default')
%%  Save log
overwrite_log = true;
if overwrite_log
    fileID = fopen("GOLDSIM_optimzation.log",'w');
    fclose(fileID);
end
diary GOLDSIM_optimzation.log

%%  Time how long the entire script takes
tstart = datetime;

%%

% maxEvals = 2000;
maxEvals = 2000;
% n_starts = 10;
%--------------------------------------------------------------------------
%   Set initial guess at parameters
%--------------------------------------------------------------------------
a = 0.15;
b = 0.89;
c = -4;
K = 2.7;

mixT = 1.0;


% p0 = [K,a,b,c,mixT];
p0 = [K,a,b,c];
%--------------------------------------------------------------------------
%   Some parameter bounds
%--------------------------------------------------------------------------
% lb = [0,0,0,-6,0.0];
% ub = [8, 1, 1.5, 0, 5.0];
lb = [0,0,0,-6];
ub = [8, 1, 1.5, 0];

%--------------------------------------------------------------------------
%   How much model reduction is done to speed up solves?
%   0 -> no model reduction (too slow)
%   0.01 -> bin particles with 1% diameter gaps 
%           (Fine for 1 plot but not optimization)
%   0.05 -> bin partciles with 5% diameter gaps 
%           (Good sweet spot for fast but fairly accurate)
%--------------------------------------------------------------------------
reduct = 0.1;

%--------------------------------------------------------------------------
%   Data to compare to
%--------------------------------------------------------------------------
[data_diam, data_PSDs] = GOLDSIM_cleanRawData("Au_quench_corediameter_hplc.txt",1);
data_times = [2,5,7,10,20,30,60,120,180,240,300];
compare_data.diam = data_diam;
compare_data.PSDs = data_PSDs;
compare_data.times = data_times;

%--------------------------------------------------------------------------
%   Function to minimize
%--------------------------------------------------------------------------

F = @(p) GOLDSIM_objectiveFunction(p, compare_data, reduct);

%--------------------------------------------------------------------------
%   Optimization
%       Use multistart in case there are many local minima
%--------------------------------------------------------------------------

%   options in optimization
% opts = optimoptions(@fmincon, ...
%     'MaxFunctionEvaluations',maxEvals,...
%     'display','iter');
% opt_problem = createOptimProblem('fmincon',...
%     'objective', F,...
%     'x0',p0,...
%     'lb',lb,...
%     'ub',ub,...
%     'options',opts);

%   multistart
% 
% ms = MultiStart("Display","final",...
%     'PlotFcn','gsplotbestf');
% [prmGlobalOpt, fbest] = run(ms,opt_problem,n_starts);
hybridopts = optimoptions('fmincon',...
    'Display','iter',...
    'MaxFunctionEvaluations',maxEvals,...
    'PlotFcn','optimplotfval');
opts = optimoptions(@simulannealbnd,...
    'PlotFcn','saplotbestf',...
    'HybridFcn',{@fmincon,hybridopts},...
    'MaxFunctionEvaluations',maxEvals,...
    'Display','iter');
prmGlobalOpt = simulannealbnd(F,p0,lb,ub,opts);
fprintf("Optimal parameters are:\n")
for iii=1:length(prmGlobalOpt)
    fprintf("\t%.5e\n",prmGlobalOpt(iii))
end

%%
%--------------------------------------------------------------------------
%   Plot the optimal solutions
%--------------------------------------------------------------------------
ic = [0.0001; 0.0003];
plot_reduct = 0.1;
% mixingtime = prmGlobalOpt(5);
tic
[sol, mySettings] = GOLDSIM_simulateGoldParticles(prmGlobalOpt,plot_reduct, ic, data_times(end));
toc

%%  Save data necessary for recreating plots
save("DATA_GOLDSIM_optimizeParameters.mat");

%%
%--------------------------------------------------------------------------
%   Plot the chemical species over time
%--------------------------------------------------------------------------

%   What times to plot
tplot = linspace(0,data_times(end),1000);
species = PBElib_getSpeciesConc(sol, mySettings, tplot);
%   Each row of 'species' corresponds to a different chemical species
figure
plot(tplot,species,'LineWidth',2)
% legend('Precursor','Solvated Precursor','Ligand')
xlabel("Time / s")
ylabel("Conc / mol/L")

%--------------------------------------------------------------------------
%   Plot PSDs over time (animated)
%--------------------------------------------------------------------------
figure
ax = gca;
h = animatedline('LineStyle','--',...
    'Marker','none',...
    'LineWidth',2,...
    'DisplayName',"Simulation",...
    'Color',[0 0.4470 0.7410]);
h2 = animatedline('LineStyle','--',...
    'Marker','none',...
    'LineWidth',2,...
    'DisplayName',"Data",...
    'Color',[0.8500 0.3250 0.0980]);
axis([0,15,0,1.0])
plotpts = data_times;
anim_seconds = 15;
[diams,PSDs] = PBElib_getPSDs(sol,mySettings,plotpts);
%   Volume weight the PSDs
PSDs = PSDs .* (diams.^3);
for iii=1:length(plotpts)
    clearpoints(h)
    clearpoints(h2)
    area = trapz(diams(:,iii),PSDs(:,iii));
    addpoints(h,diams(:,iii), PSDs(:,iii)/area);
    addpoints(h2,data_diam, data_PSDs(:,iii))
    title(ax,sprintf("PSD at %.f seconds",plotpts(iii)))
    legend('Location','northeast')
    drawnow
    pause(anim_seconds/length(plotpts))
end
%%
%--------------------------------------------------------------------------
%   Plot PSDs over time
%       FIXME
%       (a) Time as text box on northeast corner instead of title
%       (b) Shared axes so only display numbers and labels on outside
%       (c) Thicker box around border
%       (d) Also display CDF?
%--------------------------------------------------------------------------
tlfig = figure;
tl = tiledlayout('flow','TileSpacing','compact','Padding','compact');
for iii=1:length(plotpts)
    nexttile
    area = trapz(diams(:,iii),PSDs(:,iii));
    psd = PSDs(:,iii)/area;
    cdf_sim = cumtrapz(diams(:,iii), psd);
    cdf_data = cumtrapz(data_diam, data_PSDs(:,iii));
    hold on
    plot(diams(:,iii), psd, ...
        'Color',[0 0.4470 0.7410],...
        'LineStyle','--',...
        'LineWidth',2,...
        'DisplayName','Simulated PSD')
    plot(data_diam, data_PSDs(:,iii),...
        'Color',[0.8500 0.3250 0.0980],...
        'LineStyle','--',...
        'LineWidth',2,...
        'DisplayName','Measured PSD')
    yyaxis right
    plot(diams(:,iii), cdf_sim,...
        'Color',[0.4940 0.1840 0.5560],...
        'LineStyle','--',...
        'LineWidth',2,...
        'DisplayName','Simulated CDF')
    plot(data_diam, cdf_data,...
        'Color',[0.6350 0.0780 0.1840],...
        'LineStyle','--',...
        'LineWidth',2,...
        'DisplayName','Measured CDF')
    ylim([0 1.1])
    % legend('Location','northeast')
    s = sprintf("%.f seconds",plotpts(iii));
    text(10,0.8,s)
    box on
    % xlabel('Diameter / nm')
    % ylabel('Number Density / nm^{-1}')
    hold off
end
leg = legend();
leg.Layout.Tile = length(plotpts)+1;
xlabel(tl,'Diameter / nm')
ylabel(tl, 'Number Density / nm^{-1}')

% yyaxis(tl,'right')
% ylabel(tl,'Cumulative Number Density / -')

%%  End time for script
tend = datetime;

dur = tend-tstart;
fprintf("Script time (hh:mm:ss): %s\n",dur)
%%  end of diary
diary off