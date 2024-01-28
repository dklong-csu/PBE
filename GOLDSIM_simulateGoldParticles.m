function [sol, mySettings] = GOLDSIM_simulateGoldParticles(parameters, reduction, ic)

%%  Example script solving an Iridium system
use_jac = true;     %   Use analytic Jacobian (much faster)
use_emom = false;


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
reductionAmount = reduction;
fcn_reduceVal = @atoms2diam; % If reducing based on a different criteria, change this
fcn_reduceValInv = @diam2atoms; % inverse function to above fcn


%--------------------------------------------------------------------------
%   Growth and agglomeration kernels at bottom of file
%--------------------------------------------------------------------------
a = parameters(2);
b = parameters(3);
fcn_gKernel = @(size) growthKernel(size);
fcn_aKernel = @(x,y) aggregationKernel(x,y,a,b);

%--------------------------------------------------------------------------
%   Set up chemical reactions
%--------------------------------------------------------------------------
A = zeros(3, 1);
%   A + B -> P1
A(1:3,1) = [1;1;-1];

%   Reaction rates
K = [parameters(1)];

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
y0(1) = ic(1); % 0.0001
y0(2) = ic(2); % 0.0003

%   ODE solve

sol = PBElib_solveODEs(mySettings, y0, T, use_jac);

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