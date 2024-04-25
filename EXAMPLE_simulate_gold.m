clc
clear variables
close all

%%  Example script solving an Iridium system
use_jac = true;     %   Use analytic Jacobian (much faster)
time_solve = false;  %   Time the ODE solve?
use_emom = false;

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
%   Maximum particle size in number of atoms
%       Parameters below are fit to maxsize=2500
%   Conversion function from # atoms -> diameter
%   Final time to simulate to
%--------------------------------------------------------------------------
maxsize = 60000;
M = 59999;
firstsize = 1;
T = 3e2;

%--------------------------------------------------------------------------
%   Will the number of particles be reduced via binning?
%--------------------------------------------------------------------------
reductionAmount = 0.05;
fcn_reduceVal = @atoms2diam; % If reducing based on a different criteria, change this
fcn_reduceValInv = @diam2atoms; % inverse function to above fcn


%--------------------------------------------------------------------------
%   Growth and agglomeration kernels at bottom of file
%--------------------------------------------------------------------------
a = 0.5;
b = 10.1;
fcn_gKernel = @(size) growthKernel(size);
fcn_aKernel = @(x,y) aggregationKernel(x,y,a,b);

%--------------------------------------------------------------------------
%   Set up chemical reactions
%--------------------------------------------------------------------------
A = zeros(3, 1);
%   A + B -> P1
A(1:3,1) = [1;1;-1];

%   Reaction rates
K = [1e6];

%--------------------------------------------------------------------------
%   Particle Growth due to monomer addition
%       The mechanism states that
%           A + Pi -> Pi+1 + L
%       We need to define the index of A and L so the solver simulates this
%           gidx --> index of A
%           gRxnIdx --> indices of A,L
%           gRxnCoeff --> Stoichiometry for A,L, negative means a product
%           gKernel --> Reaction rate for each A + Pi -> Pi+1 + L
%--------------------------------------------------------------------------
gidx = 1;
gRxnIdx = [];
gRxnCoeff = [];

%--------------------------------------------------------------------------
%   emom
%--------------------------------------------------------------------------
emomInflowRate = 0.0;
emomDelx = 1.0;
emomGrowthRate = 0.0;



%--------------------------------------------------------------------------
%   Consolidate necessary settings
%--------------------------------------------------------------------------

mySettings = PBElib_MakeSettings(A=A,...
    K=K,...
    atoms=firstsize:maxsize,...
    reductionAmount=reductionAmount,...
    fcn_reduceVal=fcn_reduceVal,...
    fcn_reduceValInv=fcn_reduceValInv,...
    useEMoM=use_emom,...
    fcn_gkern=fcn_gKernel,...
    gidx=gidx,...
    gRxnIdx=gRxnIdx,...
    gRxnCoeff=gRxnCoeff,...
    cutoff=M,...
    fcn_akern=fcn_aKernel,...
    emomInflowRate=emomInflowRate,...
    emomDelx=emomDelx,...
    emomGrowthRate=emomGrowthRate,...
    fcn_atoms2size=@atoms2diam);



%--------------------------------------------------------------------------
%   Solve the ODE
%--------------------------------------------------------------------------

%   Initial condition 
y0 = zeros(mySettings.vecSize,1);
y0(1) = 0.0001;
y0(2) = 0.0003;

%   ODE solve

sol = PBElib_solveODEs(mySettings, y0, T, use_jac);

if time_solve
    fcn_time = @() PBElib_solveODEs(mySettings, y0, T, use_jac);
    time_to_solve = timeit(fcn_time);
    fprintf("The ODEs take %.2fs to solve on average\n",time_to_solve)
end

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

%%  Helper functions
function d = atoms2diam(a)
    d = 0.3186854328 * nthroot(a,3);
end

function a = diam2atoms(d)
    a = (d / 0.3186854328).^3;
end

function G = growthKernel(s)
    G = 0*s;
end

function A = aggregationKernel(size1,size2,a,b)
    %   Diameters from # atoms
    d1 = atoms2diam(size1);
    d2 = atoms2diam(size2);

    %   Adhesion probability
    adprob = a * exp(-b * min(d1.^2,d2.^2));

    %   Brownian kernel
    brown = ( d1.^(-1) + d2.^(-1) ) .* (d1 + d2);

    %   Physical constants
    kB = 1.380649e-23;
    Temp = 273.15 + 20;
    viscosity = 1.0016e-3;
    AvoNumber = 6.02214076e23;
    
    %   Final result
    const = kB * Temp * AvoNumber / viscosity;
    A = const * brown * min(1.0, adprob);
end




