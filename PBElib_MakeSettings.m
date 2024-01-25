function settings = PBElib_MakeSettings(options)
%PBELIB_MAKESETTINGS TODO WRITE DOCS

%--------------------------------------------------------------------------
%   Required inputs
%       options.A
%       options.K
%       options.atoms
%       options.reductionAmount
%       options.fcn_reduceVal
%       options.fcn_reduceValInv
%       options.useEMOM
%       options.fcn_gkern
%       options.gidx
%       options.gRxnIdx
%       options.gRxnCoeff
%       options.cutoff
%       options.fcn_akern
%       options.emomInflowRate
%       options.emomDelx
%       options.emomGrowthRate
%       options.fcn_atoms2size
%--------------------------------------------------------------------------
arguments
    options.A
    options.K
    options.atoms
    options.reductionAmount
    options.fcn_reduceVal
    options.fcn_reduceValInv
    options.useEMoM
    options.fcn_gkern
    options.gidx
    options.gRxnIdx
    options.gRxnCoeff
    options.cutoff
    options.fcn_akern
    options.emomInflowRate
    options.emomDelx
    options.emomGrowthRate
    options.fcn_atoms2size
end

%--------------------------------------------------------------------------
%   Resulting settings:
%       A           --> coefficient matrix for chemical reactions 
%       K           --> reactions rate vector for chemical reactions
%       vecSize     --> number of elements in the vector
%       pstart      --> index of the first particle
%       pend        --> index of the last particle
%       pBins       -->
%       gKernel     --> gKernel(i)=reaction rate for particle in index
%                       pstart+i-1
%       gidx        --> index of the precursor that causes particle growth
%       gRxnIdx     --> indices of vector that are involved in particle 
%                       growth
%       gRxnCoeff   --> coefficients of chemical species in particle growth
%                       reaction corresponding to gRxnIdx
%       cutoff      --> index of the particle that is the last to aggregate
%       aggrInfo    --> {i,j} gives the indices of particles created when
%                       particle i and particle j aggregate as well as the 
%                       amount created.
%       aKernel     --> (i,j) gives the aggregation rate between particle i
%                       and particle j.
%       useEMoM     --> (true/false) whether emom is used
%       emomInflowRate --> 
%       emomDelx
%       emomGrowthRate
%       emomBigParticleDiam
%--------------------------------------------------------------------------

    
%--------------------------------------------------------------------------
%   Chemical species information
%--------------------------------------------------------------------------

%   Assume chemical reactions create first particle
nSpecies = size(options.A,1)-1;
settings.A = options.A;
settings.K = options.K;

%--------------------------------------------------------------------------
%   Particle information
%--------------------------------------------------------------------------

%   Reduce particle sizes
pBins = zeros(length(options.atoms),1);
pBins(1:2) = [options.atoms(1);options.atoms(1)+1];
binIdx = 2;
binMax = pBins(2);
tol = 1 + options.reductionAmount;
while binMax < options.atoms(end) && binIdx < length(options.atoms)
    %   Calculate next
    binMax = floor(options.fcn_reduceValInv( ...
                    tol * options.fcn_reduceVal( pBins(binIdx) )...
                    ) );
    %   Next starts at one more
    binMax = binMax+1;
    binMax = max(binMax, pBins(binIdx)+1);
    binMax = min(max(options.atoms),binMax);
    binIdx = binIdx + 1;
    pBins(binIdx) = binMax;
end
pBins(binIdx+1) = pBins(binIdx)+1;
pBins = pBins(1:binIdx+1);
avgPSize = round( 0.5*pBins(1:end-1) + 0.5*(pBins(2:end)-1) );

%   Vector indexing
nParticles = length(pBins)-1;
vecSize = nSpecies+nParticles + 4*options.useEMoM;

settings.vecSize = vecSize;
settings.pstart = nSpecies+1;
settings.pend = vecSize - 4*options.useEMoM;
settings.pBins = pBins;


%--------------------------------------------------------------------------
%   Particle growth
%--------------------------------------------------------------------------
settings.gKernel = options.fcn_gkern(avgPSize);
settings.gidx = options.gidx;
settings.gRxnIdx = options.gRxnIdx;
settings.gRxnCoeff = reshape(options.gRxnCoeff, length(options.gRxnCoeff),1);

%--------------------------------------------------------------------------
%   Particle aggregation
%--------------------------------------------------------------------------

biggerThanCutoff = find(pBins > options.cutoff);
settings.cutoff = settings.pstart + biggerThanCutoff(1) - 1;


%   This function is really slow if you don't compile it
% fprintf("Compiling binning function\n\t")
% codegen -report PBElib_CalcAggregationBins.m -args {1, 1, pBins}
tic
nAggrBins = biggerThanCutoff(1);
aggregationMap = cell(nAggrBins,nAggrBins);
aggregationKernel = zeros(nAggrBins,nAggrBins);
for iii=1:nAggrBins
    for jjj=iii:nAggrBins
        [createdIdx, createdPercent] = PBElib_CalcAggregationBins(iii,jjj,pBins,settings.pstart);
        aggregationMap{iii,jjj} = [createdIdx, createdPercent];
        aggregationMap{jjj,iii} = [createdIdx, createdPercent];

        aggregationKernel(iii,jjj) = options.fcn_akern(avgPSize(iii),avgPSize(jjj));
        aggregationKernel(jjj,iii) = aggregationKernel(iii,jjj);
    end
end
toc
settings.aggrInfo = aggregationMap;
settings.aKernel = aggregationKernel;

%--------------------------------------------------------------------------
%   eMOM 
%--------------------------------------------------------------------------
settings.useEMoM = options.useEMoM;
settings.emomInflowRate = options.emomInflowRate;
settings.emomDelx = options.emomDelx;
settings.emomGrowthRate = options.emomGrowthRate;
settings.emomBigParticleDiam = options.fcn_atoms2size(pBins(end));
end