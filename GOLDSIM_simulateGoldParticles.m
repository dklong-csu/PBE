function [sol, mySettings] = GOLDSIM_simulateGoldParticles(pre_parameters, reduction, ic, Tend)

parameters = 10.^pre_parameters;
parameters(end) = -parameters(end);
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
gamma = parameters(2);
% b = parameters(3);
% c = parameters(4);
fcn_gKernel = @(size) growthKernel(size);
fcn_aKernel = @(x,y) aggregationKernel(x,y);
% fcn_aProb = @(x,y,vec) adhesionProbability(x,y,vec(3),a,b,c);
fcn_aProb = @(x,dx,sumint23,sumint1,vec) adhesionProbability(x,dx,sumint23,sumint1,vec,gamma,parameters(3),parameters(4));

%--------------------------------------------------------------------------
%   Set up chemical reactions
%--------------------------------------------------------------------------
A = zeros(4, 1);
%   A + B -> P1
A(1:4,1) = [1;1;-1;-1];

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

function P = adhesionProbability(size1,binSizes,sumint23,sumint1,vec,ligandEquilib,a0,b0)
    %   Boric acid concentration = vec(3)
    B = vec(3);
    %   Percentage of atoms on surface
    % A = vec(4:end-1); % conc particle bins mol/L
    %   A*avo = number of atoms in bin I
    %   SA = % on surface * (A*avo)
    %   surf atoms ~ i^(2/3)
    %   perc on surface = surf atoms/ total atoms (per cluster)
    %       ~ i^(-1/3)

    %   V = 4/3 pi (size1/2)^3
    %   SA = 4 pi (size1/2)^2
    %   SA/V = 3/(size1/2)

    %   

    %   size1 = # atoms in (average bin size)
    % expandatoms = (1:binSizes(end))';
    concBin = vec(4:end);
    binwidth = diff(binSizes);
    concmid = concBin ./ binwidth;
    surfatoms = 0;
    masscheck = 0;
    sum0 = 0;
    sum0_m = 0;
    % z23 = zeta(-2/3);
    % z1 = zeta(-1);
    
    masscheck = sum(concmid .* sumint1);
    massfix = B / masscheck;
    concmid = concmid * massfix;

    surfatoms = sum(concmid .* sumint23);
    % for eee=1:length(binSizes)-1
    %     up = binSizes(eee+1)-1;
    %     if eee == (length(binSizes)-1)
    %         up = up+1;
    %     end
    % 
    %     sum1 = z23 - hurwitzZeta(-2/3,up+1);
    %     psum = sum1 - sum0;
    %     surfatoms = surfatoms + concmid(eee)*psum;
    %     sum0 = sum1;
    % 
    %     sum1m = z1 - hurwitzZeta(-1,up+1);
    %     psumm = sum1m - sum0_m;
    %     masscheck = masscheck + concmid(eee)*psumm;
    %     sum0_m = sum1m;
    % end
    % gi = griddedInterpolant(size1,concmid,'nearest','nearest');
    % expandconc = gi(expandatoms);
    % expandconc = 0*expandatoms;
    % for eee=1:length(binSizes)-1
    %     idx1=binSizes(eee);
    %     idx2=binSizes(eee+1)-1;
    %     expandconc(idx1:idx2) = concmid(eee);
    %     if eee == (length(binSizes)-1)
    %         expandconc(end) = concmid(eee);
    %     end
    % end


    % surfatoms = expandatoms .^ (2/3);
    % totalsurf = sum( surfatoms .* expandconc );
    totalsurf = surfatoms;




    % surfatoms = [size1;mean(60000,2*59999)] .^ (2/3);
    % totalsurf = sum(surfatoms .* [A;vec(end)]);
    % SA = sum( min(size1,(size1 .^ (-1/3))) .* A );
    % SA = 
    %   Ratio between boric acid and surface atoms
    if B==0
        R=0;
    else
        R = B / totalsurf;
    end
    % fprintf("R=%e \t \t (B-atoms)/B=%f%%\n",R,100*(B-masscheck)/B)
    % 100*(B-sum(expandconc.*expandatoms))/B
    %   Average coverage
    %   prm1 is an equilibrium giving an average percentage of ligand in
    %   the solution vs attached to particles and also covers the unknown
    %   factor that gives surface atoms
    % coverage = min(1, ligandEquilib*R/(1+ligandEquilib*R));
    % coverage = (ligandEquilib*R)/(1+ligandEquilib*R);
    coverage = min(1, ligandEquilib*R);
    % fprintf("Coverage=%f\n",coverage)

    %   Hopefully correct
    [a,b] = adhesionModel(coverage,a0,b0);

    %   a*exp(-b*(x*y))
    d1 = atoms2diam(size1(1:end-1));
    d1 = d1 * ones(1,length(d1));

    % d2 = atoms2diam(size2);
    D = (d1 + d1').^1;
    P = a*exp(b * D);
    P = min(1,P);
    P = max(0,P);


end