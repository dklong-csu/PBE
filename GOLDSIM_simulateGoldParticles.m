function [sol, mySettings] = GOLDSIM_simulateGoldParticles(pre_parameters, reduction, ic, Tend)

parameters = 10.^pre_parameters;
K = parameters(1);
% A = parameters(2);
% B = parameters(3);
% C = parameters(4);
% parameters(end) = -parameters(end);
% fprintf("%.20f\n",parameters)
% parameters(1) = 10^pre_parameters(1); % reduction rxn constant
% parameters(2) = 10^pre_parameters(2); % adhesion equilibrium

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
% maxsize = 1000;
% M = 999;
firstsize = 1;
T = Tend; % 300

%--------------------------------------------------------------------------
%   Will the number of particles be reduced via binning?
%--------------------------------------------------------------------------
reductionAmount = reduction;
fcn_reduceVal = @atoms2diam; % If reducing based on a different criteria, change this
fcn_reduceValInv = @diam2atoms; % inverse function to above fcn


%--------------------------------------------------------------------------
%   Growth and agglomeration kernels at bottom of file
%--------------------------------------------------------------------------
% gamma = parameters(2);
% b = parameters(3);
% c = parameters(4);
fcn_gKernel = @(size) growthKernel(size);
fcn_aKernel = @(x,y) aggregationKernel(x,y);
% fcn_aProb = @(x,y,vec) adhesionProbability(x,y,vec(3),a,b,c);
% adhesionLangmuir(size1,vec,A,B,C)
% FIXME
% fcn_aProb = @(x,dx,sumint23,sumint1,vec) adhesionProbability(x,dx,sumint23,sumint1,vec,pre_parameters);
fcn_aProb = @(dx,sumint23,sumint1,vec) adhesionModel(dx, vec, parameters, sumint1, sumint23);
% fcn_aProb = @(size, L) adhesionLangmuir(size,L,A,B,C);

%--------------------------------------------------------------------------
%   Set up chemical reactions
%--------------------------------------------------------------------------
A = zeros(4, 1);
%   A + B -> P1
A(1:4,1) = [1;1;-1;-1];

%   Reaction rates
% K = [parameters(1)];

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
    fcn_atoms2size=@atoms2diam,...
    fcn_agglomProbability=fcn_aProb);



%--------------------------------------------------------------------------
%   Solve the ODE
%--------------------------------------------------------------------------

%   Initial condition 
y0 = zeros(mySettings.vecSize,1);
sizeIC = length(ic);
y0(1:sizeIC) = ic;
% y0(1) = ic(1); % 0.0001
% y0(2) = ic(2); % 0.0003

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

% function A = aggregationKernel(size1,size2,a,b)
function A = aggregationKernel(size1,size2)
    %   Diameters from # atoms
    d1 = atoms2diam(size1);
    d2 = atoms2diam(size2);

    %   Brownian kernel
    brown = ( d1.^(-1) + d2.^(-1) ) .* (d1 + d2);

    %   Physical constants
    kB = 1.380649e-23;
    Temp = 273.15 + 20;
    viscosity = 1.0016e-3;
    AvoNumber = 6.02214076e23;
    
    %   Final result
    const = 2/3 * kB * Temp * AvoNumber / viscosity;
    A = const * brown;
    %   fixme
    % if (size1+size2)<10
    %     A = A;
    % else
    %     A = 0;
    % end
end


function P = adhesionModel(binSizes, x, parameters, sumint1, sumint23)
    concBin = x(4:end);
    concBinAvg = concBin ./ diff(binSizes);
    nAuAtoms = sum( concBinAvg .* sumint1 );
    nAuSurfAtoms = sum( concBinAvg .* sumint23 );

    % Boric acid molecules = nAuAtoms
    B = x(3);
    % if B > 0
    %     fprintf("\nB/SA =%e\nAu/SA=%e\n",B/nAuSurfAtoms,nAuAtoms/nAuSurfAtoms)
    % end
    % R = nAuAtoms / nAuSurfAtoms;
    R = B / nAuSurfAtoms;
    equib = parameters(2);
    coverage = max( 0, min(1, equib * R) );

    %   What model to relate surface coverage and adhesion probability?
    P = coverage;
    % P = coverage ^2;
    % P = sqrt(coverage);
    % P = 3*coverage^2 - 2*coverage^3;
    P = P * parameters(3);
    P = 1 - P;
end