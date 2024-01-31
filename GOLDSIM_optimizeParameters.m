clc
clear variables
close all

rng('default')
%%
filename_log = "GOLDSIM_optimzation_5percent.log";
filename_MatData = "DATA_GOLDSIM_optimizeParameters_5percent.mat";

optimization_reduct = 0.05;
plotting_reduct = 0.01;

%%  Save log
overwrite_log = true;
if overwrite_log
    fileID = fopen(filename_log,'w');
    fclose(fileID);
end
diary(filename_log)

%%  Time how long the entire script takes
tstart = datetime;

%%

maxEvals = 10000;
%--------------------------------------------------------------------------
%   Set initial guess at parameters
%--------------------------------------------------------------------------
a = 0.15;
b = 0.89;
c = -4;
K = 2.7;

p0 = [K,a,b,c];
%--------------------------------------------------------------------------
%   Some parameter bounds
%--------------------------------------------------------------------------

lb = [0,0,0,-6];
ub = [8, 1, 1.5, 0];


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

F = @(p) GOLDSIM_objectiveFunction(p, compare_data, optimization_reduct);

%--------------------------------------------------------------------------
%   Optimization
%--------------------------------------------------------------------------

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

tic
[solOptimize, mySettingsOptimize] = GOLDSIM_simulateGoldParticles(prmGlobalOpt,optimization_reduct, ic, data_times(end));
toc

tic
[solAccurate, mySettingsAccurate] = GOLDSIM_simulateGoldParticles(prmGlobalOpt,plotting_reduct, ic, data_times(end));
toc

%%  Save data necessary for recreating plots
save(filename_MatData);



%%  End time for script
tend = datetime;

dur = tend-tstart;
fprintf("Script time (hh:mm:ss): %s\n",dur)
%%  end of diary
diary off