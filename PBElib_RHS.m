function dy = PBElib_RHS(~,y,settings) %#codegen
%PBELIB_RHS FIXME WRITE DOCS

dy = zeros(length(y),1);

%--------------------------------------------------------------------------
%   Chemical reactions
%--------------------------------------------------------------------------
[n_rxnrows, ~] = size(settings.A);
dy(1:n_rxnrows) = sum( -settings.A .* ( prod( y(1:n_rxnrows) .^ max(0,settings.A) ) .* settings.K ) ,2);


%--------------------------------------------------------------------------
%   Particle Growth
%--------------------------------------------------------------------------

dGrowth = y(settings.gidx) ...
    * y(settings.pstart:settings.pend) ...
    .* settings.gKernel; 
binSizes = diff(settings.pBins);
w1 = (binSizes - 1) ./ binSizes;
w2 = 1 ./ binSizes;
%   Loss from growth
dy(settings.pstart:settings.pend) = dy(settings.pstart:settings.pend) ...
                                    - dGrowth;
%   Particle growth
%   If binned, w1% stay in same bin, w2% leave bin
dy(settings.pstart+1:settings.pend) = dy(settings.pstart+1:settings.pend) ...
                                        + w2(1:end-1) .* dGrowth(1:end-1);
dy(settings.pstart:settings.pend) = dy(settings.pstart:settings.pend) ...
                                    + w1(1:end) .* dGrowth(1:end);
%   Other chemical species involved in growth
dy(settings.gRxnIdx) = dy(settings.gRxnIdx) - settings.gRxnCoeff .* sum(dGrowth);


%--------------------------------------------------------------------------
%   Particle aggregation
%--------------------------------------------------------------------------
aggregIdx = settings.pstart:settings.cutoff;
aggrRates = y(aggregIdx) * y(aggregIdx)';
aggrRates = aggrRates .* settings.aKernel;
aggrRates = tril(aggrRates);
%---------------------------------------------------
%   aggrRates = [aij*yi*yj] (on lower triangle only)
%   --> gives rxn rate for every Pi+Pj agglomeration
%   --> sum accross columns gives total for each Pi
%---------------------------------------------------
dy(aggregIdx) = dy(aggregIdx) - sum(aggrRates,2);
dy(aggregIdx) = dy(aggregIdx) - (1./binSizes(1:length(aggregIdx))) .* diag(aggrRates);
%---------------------------------------------------
%   For created particles, need to retrieve information
%   of how many of each bin are created for each Pi+Pj
%---------------------------------------------------
for key=keys(settings.aggrInfo)
    weights = settings.aggrInfo{key};
    dy(key) = dy(key) + sum(weights.*aggrRates,"all");
end


%--------------------------------------------------------------------------
%   emom
%--------------------------------------------------------------------------
if settings.useEMoM == true
    idx_mu0 = length(y)-3;
    idx_mu1 = length(y)-2;
    idx_mu2 = length(y)-1;
    idx_mu3 = length(y);
    idxs_emom = idx_mu0:idx_mu3;

    k_inflow = settings.emomInflowRate;
    delx = settings.emomDelx; 
    k_emom = settings.emomGrowthRate; 
    diamM = settings.emomBigParticleDiam; 
    
    flux = nthroot(9*pi / 1, 3) * k_inflow * ( (settings.pBins(end)-1) ^ (2/3) ) * w2(end)*y(idx_mu0-1) / delx / k_emom;
    
    %   Effect on chemical species due to growth reaction
    %   FIXME 8/9 and so on hardcoded
    growthRate = 8/9 * delx * y(settings.gidx) * k_emom ...
        * diamM^3 * flux ...
        + k_emom * 8 / (delx^2 * 3) * y(settings.gidx) * y(idx_mu2);
    dy(settings.gRxnIdx) = dy(settings.gRxnIdx) ...
        - growthRate * settings.gRxnCoeff;


    %   Effect on moments
    diamMpow = (0:3)';
    dy(idxs_emom) = dy(idxs_emom) ...
        + 8 * delx / 9 * y(settings.gidx) * k_emom * flux * diamM.^diamMpow ...
        + 8 * k_emom * delx / 9 * y(settings.gidx) * diamMpow .* y(idxs_emom-1);
end

end









