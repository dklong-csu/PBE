function dfdy = PBElib_Jacobian(~,y,settings) %#codegen
%PBELIB_JACOBIAN FIXME WRITE DOCS

dfdy = zeros(settings.vecSize);

binSizes = diff(settings.pBins);
w1 = (binSizes-1) ./ binSizes;
w2 = 1 ./ binSizes;

[n_rxnrows, ~] = size(settings.A);
for iii=1:settings.vecSize
    noniii = setdiff(1:n_rxnrows,iii);
    %----------------------------------------------------------------------
    %   Chemical reactions
    %----------------------------------------------------------------------
    if iii<=n_rxnrows
    dfdy(1:n_rxnrows,iii) = dfdy(1:n_rxnrows,iii) +...
        sum(...
        -settings.A .* ...
        ( ...
            y(iii) .^ max(settings.A(iii,:)-1,0) ...
            .* prod( y(noniii) .^ max(0, settings.A(noniii,:)) ) ...
            .* settings.K ...
            .* max(0, settings.A(iii,:)) ...
        ), 2);
    end


    %----------------------------------------------------------------------
    %   Particle Growth
    %----------------------------------------------------------------------
    if iii==settings.gidx
        %   A * sum(Pi*G) -> sum(Pi*G)
        dGrowthdy = y(settings.pstart:settings.pend) .* settings.gKernel(:);
        %   Loss in particles
        dfdy(settings.pstart:settings.pend, iii) = dfdy(settings.pstart:settings.pend, iii) - dGrowthdy;
        %   Gain in particles
        dfdy(settings.pstart+1:settings.pend, iii) = dfdy(settings.pstart+1:settings.pend, iii) +   w2(1:end-1).*dGrowthdy(1:end-1);
        dfdy(settings.pstart:settings.pend-1, iii) = dfdy(settings.pstart:settings.pend-1, iii) +   w1(1:end-1).*dGrowthdy(1:end-1);
        %   Gain/loss in chemical species
        dfdy(settings.gRxnIdx, iii) = dfdy(settings.gRxnIdx, iii) - settings.gRxnCoeff .* sum(dGrowthdy);
    elseif (iii >= settings.pstart) && (iii <= settings.pend)
        %   A * Piii*Giii -> A * Giii
        %   A (neg)
        %   Piii (neg)
        %   Piii+1 (pos)
        %   L (pos)
        dGrowthdy = y(settings.gidx) * settings.gKernel(iii-settings.pstart+1);
        dfdy(iii,iii) = dfdy(iii,iii) - dGrowthdy;
        if iii < settings.pend
            dfdy(iii+1,iii) = dfdy(iii+1,iii) + dGrowthdy;
        end
        dfdy(settings.gRxnIdx, iii) = dfdy(settings.gRxnIdx, iii) - settings.gRxnCoeff * dGrowthdy;
    end

    %----------------------------------------------------------------------
    %   Particle aggregation
    %----------------------------------------------------------------------
    
    if (iii >= settings.pstart) && (iii <= settings.cutoff)
        %   Pjjj
        aggregIdx = settings.pstart:settings.cutoff;
        dPiii = zeros(1,length(aggregIdx));
        idxiii = iii-settings.pstart+1;
        dPiii(idxiii) = 1;
        aggrRates = y(aggregIdx) * dPiii;
        aggrRates = aggrRates .* settings.aKernel;
        aggrRates = tril(aggrRates);
        aggrRates(idxiii,idxiii) = 2* aggrRates(idxiii,idxiii);

        dfdy(aggregIdx,iii) = dfdy(aggregIdx,iii) - sum(aggrRates,2);
        dfdy(iii,iii) = dfdy(iii,iii) - aggrRates(idxiii,idxiii)/binSizes(idxiii);

        akeys = keys(settings.aggrInfo);
        for k=1:length(akeys)
            key = akeys(k);
            weights = settings.aggrInfo{key};
            dfdy(key,iii) = dfdy(key,iii) + sum(weights.*aggrRates,"all");
        end
    end
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
    
    % flux = 3 * k_inflow * ( settings.cutoff ^ (2/3) ) * y(idx_mu0-1) / delx / k_emom;
    flux = nthroot(9*pi / 1, 3) * k_inflow * ( (settings.pBins(end)-1) ^ (2/3) ) * w2(end)*y(idx_mu0-1) / delx / k_emom;
    % dfluxdP = 3 * k_inflow * ( settings.cutoff ^ (2/3) ) / delx / k_emom;
    dfluxdP = nthroot(9*pi / 1, 3) * k_inflow * ( (settings.pBins(end)-1) ^ (2/3) ) * w2(end) / delx / k_emom;

    %   Chemical species
    %   FIXME
    dCdA = 8 * delx / 9 * k_emom * diamM^3 * flux ...
        + k_emom * 8 / (delx^2 * 3) * y(idx_mu2);
    dCdP = 8 * delx / 9 * k_emom * diamM^3 * y(settings.gidx) * dfluxdP;
    dCdMu = k_emom / delx^2 * 8 / 3 * y(settings.gidx);
    
    dfdy(settings.gRxnIdx, settings.gidx) = dfdy(settings.gRxnIdx, settings.gidx) ...
        - dCdA * settings.gRxnCoeff;
    dfdy(settings.gRxnIdx, idx_mu0-1) = dfdy(settings.gRxnIdx, idx_mu0-1) ...
        - dCdP * settings.gRxnCoeff;
    dfdy(settings.gRxnIdx, idx_mu2) = dfdy(settings.gRxnIdx, idx_mu2) ...
        - dCdMu * settings.gRxnCoeff;

    %   Moments
    diamMpow = (0:3)';
    dCdA = 8 * delx / 9 * k_emom * flux * diamM.^diamMpow ...
        + 8 * k_emom * delx / 9 * diamMpow .* y(idxs_emom-1);
    dfdy(idxs_emom, settings.gidx) = dfdy(idxs_emom, settings.gidx) + dCdA;

    dCdP = 8 * delx / 9 * y(settings.gidx) * k_emom * dfluxdP * diamM.^diamMpow;
    dfdy(idxs_emom, idx_mu0-1) = dfdy(idxs_emom, idx_mu0-1) + dCdP;

    dCdMu = 8 * k_emom * delx / 9 * y(settings.gidx) * diamMpow;
    dfdy(idx_mu3, idx_mu2) = dfdy(idx_mu3, idx_mu2) + dCdMu(end);
    dfdy(idx_mu2, idx_mu1) = dfdy(idx_mu2, idx_mu1) + dCdMu(end-1);
    dfdy(idx_mu1, idx_mu0) = dfdy(idx_mu1, idx_mu0) + dCdMu(end-2);
end




dfdy = sparse(dfdy);
end