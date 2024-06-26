clc
clear variables
% close all

%%  Example script solving an Iridium system
use_jac = true;     %   Use analytic Jacobian (much faster)
time_solve = true;  %   Time the ODE solve?
use_emom = false;

%--------------------------------------------------------------------------
%   Save all of the data needed to reproduce plots
%       sol --> the solution to the ODE which can be evaluated with deval
%               at any time
%       mySettings --> gives information about vector indices and so forth
%       diams --> gives particle diameters
%--------------------------------------------------------------------------
save_data = false;   %   False if you don't need to save data for some reason
saveFileName = "DATA_iridium_2500.mat";

%--------------------------------------------------------------------------
%   Mechanism modeled
%       FIXME more detail
%   A ->[kf*S*S] As + L
%   As + L ->[kb] A
%   A + 2As ->[k1] P3 + L
%   A + Pi ->[G(k2,k3,M)*r(i)] Pi+1 + L
%   (agglom FIXME)
%--------------------------------------------------------------------------

kb = 1.37e5;
kf = 5e-7*kb;
k1 = 7.69e4;
k2 = 1.4e4;
k3 = 7.15e3;
k4 = 1.74e3;
M = 111; 
S = 11.3;

%--------------------------------------------------------------------------
%   Maximum particle size in number of atoms
%       Parameters below are fit to maxsize=2500
%   Conversion function from # atoms -> diameter
%   Final time to simulate to
%--------------------------------------------------------------------------
if use_emom
    maxsize = 2*M;
else 
    maxsize = 2500;
end
firstsize = 3;
atoms2diam = @(a) 0.3000805 * nthroot(a,3);
T = 10.0;

%--------------------------------------------------------------------------
%   Will the number of particles be reduced via binning?
%--------------------------------------------------------------------------
reductionAmount = 0.0;
fcn_reduceVal = atoms2diam; % If reducing based on a different criteria, change this
fcn_reduceValInv = @(d) (d/0.3000805) .^ 3; % inverse function to above fcn



%--------------------------------------------------------------------------
%   Growth and agglomeration kernels at bottom of file
%--------------------------------------------------------------------------

fcn_gKernel = @(size) growthKernel(size,k2,k3,M);
fcn_aKernel = @(x,y) aggregationKernel(x,y,k4,M);

%--------------------------------------------------------------------------
%   Set up chemical reactions
%       Matrix A has a column for each reaction
%       The rows correspond to the rows of the vector
%       A(i,j) = stoichiometric coefficient of species j in reaction i
%           A(i,j) > 0 --> a reactant
%           A(i,j) < 0 --> a product
%           A(i,j) = 0 --> not in the reaction
%   The matrix should be as small as possible for computational performance
%   Most likely, the chemical reactions create the smallest tracked
%   particle but do not interact with any of the larger particles. Hence,
%   matrix A can be limited to number of rows equal to nSpecies+1.
%
%       Vector K is a row vector where each entry is the reaction rate for
%       the corresponding reaction
%--------------------------------------------------------------------------
A = zeros(4, 3);
%   A -> As + L
A(1:3,1) = [1;-1;-1];
%   As + L -> A
A(1:3,2) = [-1;1;1];
%   A + 2As -> P + L
A(1:4,3) = [1;2;-1;-1];

%   Reaction rates
K = [kf*S*S, kb, k1];

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
gRxnIdx = [1,3];
gRxnCoeff = [1,-1];

%--------------------------------------------------------------------------
%   emom information
%--------------------------------------------------------------------------
emomInflowRate=k3;
emomGrowthRate=k3;
emomDelx=0.3000805;

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
    fcn_atoms2size=atoms2diam);



%--------------------------------------------------------------------------
%   Solve the ODE
%--------------------------------------------------------------------------

%   Initial condition --> A(0) = 0.0012, all others 0
y0 = zeros(mySettings.vecSize,1);
y0(1) = 0.0012;

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
legend('Precursor','Solvated Precursor','Ligand')
xlabel("Time / h")
ylabel("Conc / mol/L")

%--------------------------------------------------------------------------
%   Plot PSDs over time
%--------------------------------------------------------------------------
figure
ax = gca;
h = animatedline('LineStyle','-',...
    'Marker','none',...
    'LineWidth',2);
axis([0,max(5,atoms2diam(maxsize)),0,Inf])
n_times = 10;
plotpts = linspace(0,T,n_times+1);
plotpts = plotpts(2:end);
anim_seconds = 5;
[diams,PSDs] = PBElib_getPSDs(sol,mySettings,plotpts);
for iii=1:length(plotpts)
    clearpoints(h)
    addpoints(h,diams(:,iii), PSDs(:,iii));
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
function r = rFcn(size)
    r = 8/3 * (size .^ (2/3) );
end

function G = growthKernel(s,k1,k2,M)
    G = k1 * ones(size(s)) .* (s <= M) + k2 * ones(size(s)) .* (s > M);
    G = rFcn(s) .* G;
end

function A = aggregationKernel(size1,size2,k, M)
    A = rFcn(size1) * rFcn(size2) * k * (size1 <= M) * (size2 <= M);
end