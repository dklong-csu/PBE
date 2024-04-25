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
    options.fcn_agglomProbability
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
%       avgPSize
%       dSize
%       aProb
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
    %   Don't jump over aggregation cutoff point
    if pBins(binIdx) <= options.cutoff
        binMax = min(binMax, options.cutoff);
    end
    %   Next starts at one more
    binMax = binMax+1;
    binMax = max(binMax, pBins(binIdx)+1);
    binMax = min(max(options.atoms),binMax);
    binIdx = binIdx + 1;
    pBins(binIdx) = binMax;
end
biggerThanCutoff = find(pBins > options.cutoff);
maxAgglomBin = biggerThanCutoff(1);
pBins(binIdx+1) = max(pBins(binIdx)+1, (pBins(maxAgglomBin)-1)*2+1);
pBins = pBins(1:binIdx+1);
avgPSize = round( 0.5*pBins(1:end-1) + 0.5*(pBins(2:end)-1) );
settings.avgPSize = avgPSize;
settings.sizes = options.fcn_atoms2size(avgPSize);
dsizes = settings.sizes - options.fcn_atoms2size(avgPSize-1);
settings.dsizes = dsizes;

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
settings.cutoff = settings.pstart + biggerThanCutoff(1) - 2;

nAggrBins = biggerThanCutoff(1)-1;
aggregationMap = dictionary;
aggregationKernel = zeros(nAggrBins,nAggrBins);
for iii=1:nAggrBins
    for jjj=iii:nAggrBins
        [createdIdx, createdPercent] = PBElib_CalcAggregationBins(iii,jjj,pBins,settings.pstart);
        for k=1:length(createdIdx)
            key=createdIdx(k);
            if isConfigured(aggregationMap) && isKey(aggregationMap,key)
                val = aggregationMap{key};
            else
                val = sparse(nAggrBins,nAggrBins);
            end
            val(jjj,iii) = createdPercent(k);   %   fixme?
            aggregationMap(key) = {val};
        end
        aggregationKernel(iii,jjj) = options.fcn_akern(avgPSize(iii),avgPSize(jjj));
        aggregationKernel(jjj,iii) = aggregationKernel(iii,jjj);
    end
end

settings.aggrInfo = aggregationMap;
settings.aKernel = aggregationKernel;

% diam = avgPSize(1:nAggrBins) * ones(1,nAggrBins);
diam = avgPSize(1:nAggrBins);
% @(x,dx,sumint23,sumint1,vec)
z1 = zeta(-1);
z23 = zeta(-2/3);
pUB = pBins-1;
tau = -2/3;
sumint23 = z23 - hurwitzZeta(tau,pUB(1:end-1)+1);
sumint23 = z23 - hurwitzZeta(tau,pUB(2:end)+1) - sumint23;
% sumint23 = 
% sumint23 = z23 - hurwitzZeta(pUB(1:end-1)+1,tau);
% sumint23 = z23 - hurwitzZeta(pUB(2:end)+1,tau) - sumint23;

sumint1 = z1 - hurwitzZeta(-1,pUB(1:end-1)+1);
sumint1 = z1 - hurwitzZeta(-1,pUB(2:end)+1) - sumint1;
settings.aProb = @(y) options.fcn_agglomProbability(pBins, sumint23, sumint1, y);
% settings.aProb = @(y) options.fcn_agglomProbability(avgPSize,pBins,sumint23, sumint1,y);
% settings.aProb = @(y) options.fcn_agglomProbability(avgPSize,y(3));

%--------------------------------------------------------------------------
%   eMOM 
%--------------------------------------------------------------------------
settings.useEMoM = options.useEMoM;
settings.emomInflowRate = options.emomInflowRate;
settings.emomDelx = options.emomDelx;
settings.emomGrowthRate = options.emomGrowthRate;
settings.emomBigParticleDiam = options.fcn_atoms2size(pBins(end));

end