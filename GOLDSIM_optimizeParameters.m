clc
clear variables
close all

rng('default')
%%
filename_log = "GOLDSIM_optimzation_surface_boric1.log";
filename_MatData = "DATA_GOLDSIM_optimizeParameters_surface_boric1.mat";

optimization_reduct = 0.015;% 0.05
plotting_reduct = 0.01; % 0.01

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

maxEvals = 1000; %10000
%--------------------------------------------------------------------------
%   Set initial guess at parameters
%--------------------------------------------------------------------------

% K = 7;
% E = -1;
% M = 0;
K = log10(8.573367e+05);
E = log10(1.448856e-03);
M = log10(9.838769e-01);

p0 = [K, E, M];
lb = [0, -10, -4];
ub = [6.89, 10, 0];


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
opt_strat = 2;
switch opt_strat
    case 1
        hybridopts = optimoptions('fmincon',...
            'Display','iter',...
            'MaxFunctionEvaluations',maxEvals,...
            'PlotFcn','optimplotfval',...
            'UseParallel',false);
        opts = optimoptions(@simulannealbnd,...
            'PlotFcn','saplotbestf',...
            'HybridFcn',{@fmincon,hybridopts},...
            'MaxFunctionEvaluations',maxEvals,...
            'Display','iter');
        
        prmGlobalOpt = simulannealbnd(F,p0,lb,ub,opts);
    case 2
        hybridopts = optimoptions('fmincon',...
            'Display','iter',...
            'MaxFunctionEvaluations',maxEvals,...
            'PlotFcn','optimplotfval',...
            'UseParallel',true);
        opts = optimoptions(@surrogateopt,...
            'PlotFcn','surrogateoptplot',...
            'MaxFunctionEvaluations',maxEvals,...
            'Display','final',...
            'UseParallel',true,...
            'InitialPoints',p0);
        prmGlobalOpt = surrogateopt(F,lb,ub,opts);
        [prmGlobalOpt,fval] = fmincon(F,prmGlobalOpt,[],[],[],[],lb,ub,[],hybridopts);
end

%%
fprintf("Optimal parameters are:\n")
fprintf("\t%.6e",10.^prmGlobalOpt)
fprintf("\n")
% fprintf("K = %.10e;\n",prmGlobalOpt(1))
% fprintf("E = %.10e;\n",prmGlobalOpt(2))
fprintf("F(opt prm) = %f\n", fval)
%%
%--------------------------------------------------------------------------
%   Plot the optimal solutions
%--------------------------------------------------------------------------
ic = [0.0001; 0.0003; 0.0096];

tic
[solOptimize, mySettingsOptimize] = GOLDSIM_simulateGoldParticles(prmGlobalOpt,optimization_reduct, ic, data_times(end));
toc
%%
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