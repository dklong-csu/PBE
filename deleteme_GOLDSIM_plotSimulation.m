clc
clear variables
close all

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
%   Set parameters
%--------------------------------------------------------------------------
a = 0.5;
b = 10.1;
K = 1e6;
reduct = 0.05;
parameters = [K,a,b];
T = 300;


%--------------------------------------------------------------------------
%   Simulation
%--------------------------------------------------------------------------

tic
[sol, mySettings] = GOLDSIM_simulateGoldParticles(parameters,reduct);
toc

%%
%--------------------------------------------------------------------------
%   Plot the chemical species over time
%--------------------------------------------------------------------------

%   What times to plot
tplot = linspace(0,T,1000);
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
h = animatedline('LineStyle','none',...
    'Marker','o',...
    'LineWidth',2);
axis([0,15,0,Inf])
n_times = 10;
plotpts = linspace(0,T,n_times+1);
plotpts = plotpts(2:end);
anim_seconds = 5;
[diams,PSDs] = PBElib_getPSDs(sol,mySettings,plotpts);
for iii=1:length(plotpts)
    clearpoints(h)
    area = trapz(diams(:,iii),PSDs(:,iii));
    addpoints(h,diams(:,iii), PSDs(:,iii)/area);
    title(ax,sprintf("PSD at %.2f hours",plotpts(iii)))
    drawnow
    pause(anim_seconds/n_times)
end

%%
%--------------------------------------------------------------------------
%   Save all of the data needed to reproduce figures
%--------------------------------------------------------------------------

if save_data
    save(saveFileName,"sol","mySettings","diams","PSDs",'-mat');
end