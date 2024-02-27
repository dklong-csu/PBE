clc
clear variables
close all

rng('default')
%%
filename_log = "GOLDSIM_optimzation_test.log";
filename_MatData = "DATA_GOLDSIM_optimizeParameters_test.mat";

optimization_reduct = 0.01;% 0.05
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

%   x .* y --> ~70
% K = 5.44432e+00;
% gamma = -1.07601e+00;
% a0 = -7.96271e-02;
% b0 = -8.74913e-01;
% K = 4.45061e+00;
%   x + y --> ~64
K = 4.57595e+00;
gamma = -1.05754e+00; 
a0 = -9.40497e-03;
b0 = -1.15462e+00;
p0 = [K,gamma,a0,b0];
%--------------------------------------------------------------------------
%   Some parameter bounds
%--------------------------------------------------------------------------

% lb = [0,0,0,-6];
% ub = [8, 1, 1.5, 0];
%   [reduct K, ligand equilib, a0, b0]
lb = [0,-10,-10,-10];
ub = [10,10,0,10];


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

% hybridopts = optimoptions('fmincon',...
%     'Display','iter',...
%     'MaxFunctionEvaluations',maxEvals,...
%     'PlotFcn','optimplotfval');
% opts = optimoptions(@simulannealbnd,...
%     'PlotFcn','saplotbestf',...
%     'HybridFcn',{@fmincon,hybridopts},...
%     'MaxFunctionEvaluations',maxEvals,...
%     'Display','iter');
% 
% prmGlobalOpt = simulannealbnd(F,p0,lb,ub,opts);


opts = optimoptions(@surrogateopt,...
    'PlotFcn','surrogateoptplot',...
    'MaxFunctionEvaluations',maxEvals,...
    'Display','final',...
    'UseParallel',true,...
    'InitialPoints',p0);
prmGlobalOpt = surrogateopt(F,lb,ub,opts);

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
%%
tic
[solAccurate, mySettingsAccurate] = GOLDSIM_simulateGoldParticles(prmGlobalOpt,plotting_reduct, ic, data_times(end));
toc

%%
b = -10^prmGlobalOpt(end);
decayfcn = @(x,y) exp(b*(x+y));
figure
fsurf(decayfcn,[0 10 0 10])
view(2)
colorbar
%%  Save data necessary for recreating plots
save(filename_MatData);



%%  End time for script
tend = datetime;

dur = tend-tstart;
fprintf("Script time (hh:mm:ss): %s\n",dur)
%%  end of diary
diary off