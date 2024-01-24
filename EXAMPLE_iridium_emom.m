clc
clear variables
close all

%%  Example script solving an Iridium system
use_jac = true;     %   Use analytic Jacobian (much faster)
use_mex = true;     %   Compile functions (a bit faster)
time_solve = true;  %   Time the ODE solve?
regen_mex = true;   %   Does the mex function need to be recompiled? (Needed if vector/matrix sizes change)
use_emom = true;    %   Use emom for large particles
%--------------------------------------------------------------------------
%   Save all of the data needed to reproduce plots
%       sol --> the solution to the ODE which can be evaluated with deval
%               at any time
%       mySettings --> gives information about vector indices and so forth
%       diams --> gives particle diameters
%--------------------------------------------------------------------------
save_data = true;   %   False if you don't need to save data for some reason
saveFileName = "DATA_iridium_emom.mat";

%--------------------------------------------------------------------------
%   Maximum particle size in number of atoms
%       Parameters below are fit to maxsize=2500
%   Conversion function from # atoms -> diameter
%   Final time to simulate to
%--------------------------------------------------------------------------
maxsize = 222;
delx = 0.3000805;
atoms2diam = @(a) delx * nthroot(a,3);
T = 10.0;

%--------------------------------------------------------------------------
%   Growth and agglomeration kernels at bottom of file
%       Growth(i atoms) = 8/3 * i^(2/3) * k2  --> i<=M
%       Growth(i atoms) = 8/3 * i^(2/3) * k3  --> i>M
%       Agglom(i,j atoms) = 8/3 * i^(2/3) * 8/3 * j^(2/3) * k4 --> i,j<=M
%       Agglom(i,j atoms) = 0 --> i or j > M
%--------------------------------------------------------------------------

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
%   Set up vectors
%       nRxns --> known from mechanism
%       nSpecies --> known from mechanism
%       firstSize --> known from mechanism
%--------------------------------------------------------------------------
nRxns = 3;
nSpecies = 3;
firstSize = 3;

particles = firstSize:maxsize;
sizes = atoms2diam(particles);
dsizes = sizes - atoms2diam(particles-1);
pstart = nSpecies+1;                    %   start index for particles
pend = pstart + length(particles) - 1;  %   end index for particles
vecSize = pend;                         %   emom adds 4
if use_emom
    vecSize = vecSize + 4;
end

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
A = zeros(pstart, nRxns);
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
%               gKernel(end) = 0 so that mass is not lost
%                   if eMoM is used, do not make this edit!
%--------------------------------------------------------------------------
gidx = 1;
gRxnIdx = [1,3];
gRxnCoeff = [1,-1];
gKernel = growthFcn(particles, k2, k3, M);

%--------------------------------------------------------------------------
%   Particle Agglomeration
%       aKernel --> agglomeration kernel
%           aKernel(i,j) --> agglomeration rate between particle(i) and
%           particle(j)
%--------------------------------------------------------------------------

aKernel = agglomFcn(M,k4);


%--------------------------------------------------------------------------
%   Solve the ODE
%--------------------------------------------------------------------------

%   Consolidate settings
mySettings = PBElib_Settings(vecSize, ...
        A=A,...
        K=K,...
        pstart=pstart,...
        pend=pend,...
        gidx=gidx,...
        gKernel=gKernel,...
        gRxnIdx=gRxnIdx,...
        gRxnCoeff=gRxnCoeff,...
        cutoff=M,...
        aKernel=aKernel,...
        useEMoM=use_emom,...
        emomDelx=delx,...
        emomInflowRate=k3,...
        emomGrowthRate=k3,...
        emomBigParticleDiam=atoms2diam(maxsize),...
        sizes=sizes,...
        dsizes=dsizes);

%   Initial condition --> A(0) = 0.0012, all others 0
y0 = zeros(vecSize,1);
y0(1) = 0.0012;

%   Recompile functions with new vector size
if regen_mex
    fprintf("Compiling right hand side function\n\t")
    codegen -report PBElib_RHS.m -args {0.0, y0, mySettings}
    fprintf("Compiling Jacobian function\n\t")
    codegen -report PBElib_Jacobian.m -args {0.0, y0, mySettings}
end

%   ODE solve
sol = PBElib_solveODEs(mySettings, y0, T, use_mex, use_jac);
if time_solve
    fcn_time = @() PBElib_solveODEs(mySettings, y0, T, use_mex, use_jac);
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
axis([0,5,0,Inf])
n_times = 10;
plotpts = linspace(0,T,n_times+1);
plotpts = plotpts(2:end);
[diams, PSDs] = PBElib_getPSDs(sol,mySettings,plotpts);
anim_seconds = 5;
for iii=1:length(plotpts)
    ttt = plotpts(iii);
    clearpoints(h)
    addpoints(h,diams(:,iii), PSDs(:,iii));
    title(ax,sprintf("PSD at %.2f hours",ttt))
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

function G = growthFcn(s,k1,k2,M)
    G = k1 * ones(size(s)) .* (s <= M) + k2 * ones(size(s)) .* (s > M);
    G = rFcn(s) .* G;
end

function A = agglomFcn(M,k)
    sizes = 1:M;
    A = rFcn(sizes)' * rFcn(sizes) * k;
end