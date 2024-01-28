clc
clear variables
close all

%%  Save log
overwrite_log = true;
if overwrite_log
    fileID = fopen("GOLDSIM_optimzation.log",'w');
    fclose(fileID);
end
diary GOLDSIM_optimzation.log

%%
%--------------------------------------------------------------------------
%   Save all of the data needed to reproduce plots
%       sol --> the solution to the ODE which can be evaluated with deval
%               at any time
%       mySettings --> gives information about vector indices and so forth
%       diams --> gives particle diameters
%--------------------------------------------------------------------------
save_data = false;   %   False if you don't need to save data for some reason
saveFileName = "DATA_gold_2500.mat";

%--------------------------------------------------------------------------
%   Set initial guess at parameters
%--------------------------------------------------------------------------
a = 0.5;
b = 10.1;
K = 1e6;
p0 = [K,a,b];

%--------------------------------------------------------------------------
%   Some parameter bounds
%--------------------------------------------------------------------------
lb = [0,0,0];
ub = [1e10, 1, 100];

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
%   Global optimization
%       This is done first because the general values of the parameters are
%       unknown a priori. This will give a rough idea of what the optimal
%       parameters are but is not extremely accurate
%   Local optimization
%       With the hybrid setting, a local optimization is performed after
%       the global optimization to ``polish'' the solution
%--------------------------------------------------------------------------
maxEvals = 1000;
%   options in optimization
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

%%
%--------------------------------------------------------------------------
%   Plot the optimal solutions
%--------------------------------------------------------------------------
ic = [0.0001, 0.0003];
plot_reduct = 0.01;
tic
[sol, mySettings] = GOLDSIM_simulateGoldParticles(prmGlobalOpt,plot_reduct, ic, data_times(end));
toc

%%  Save data necessary for recreating plots
save("GOLDSIM_optimal.mat","sol","mySettings","prmGlobalOpt","compare_data");
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
xlabel("Time / h")
ylabel("Conc / mol/L")

%--------------------------------------------------------------------------
%   Plot PSDs over time
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
axis([0,15,0,Inf])
plotpts = data_times;
anim_seconds = 5;
[diams,PSDs] = PBElib_getPSDs(sol,mySettings,plotpts);
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

%%  end of diary
diary off