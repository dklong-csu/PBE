function dfdy = PBElib_Jacobian(~,y,settings) %#codegen
%PBELIB_JACOBIAN FIXME WRITE DOCS
    dfdy = zeros(settings.vecSize);

    binSizes = reshape(settings.pBins,numel(settings.pBins),1); %  make sure it's a column vector 
    w1 = (binSizes-1) ./ binSizes;
    w2 = 1 ./ binSizes;

    %   Chemical reactions
    [n_rxnrows, ~] = size(settings.A);
    for iii=1:settings.vecSize
        noniii = setdiff(1:n_rxnrows,iii);
        %   Chemical reactions
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


        %   Particle Growth
        %   A + Pi -> Pi+1 + L
        %       For A and Piii
        
        if iii==settings.gidx
            %   A * sum(Pi*G) -> sum(Pi*G)
            dGrowthdy = y(settings.pstart:settings.pend) .* settings.gKernel';
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

        %   Particle Agglomeration
        %   Piii + P* -> Piii+*
        %   Deriv wrt yiii
        %       Aggl(iii,*) * Piii * P* -> Aggl(iii,*) * P*
        %   dfiiidiii = -sum(...)
        %   dfxdiii = - Aggl(iii,x) * 
        if (iii >= settings.pstart) && (iii <= settings.cutoff)
            jjj = [settings.pstart:iii-1, iii+1:settings.cutoff];   % don't include iii because equation is special
            sizeiii = iii-1; %fixme
            sizejjj = jjj-1; %fixme
            created = sizejjj + sizeiii + 1; %fixme
            dAggldy = y(jjj) .* settings.aKernel(sizejjj, sizeiii);
            
            dfdy(jjj,iii) = dfdy(jjj,iii) - dAggldy;
            dfdy(created,iii) = dfdy(created,iii) + dAggldy;
            dfdy(iii,iii) = dfdy(iii,iii) - sum(dAggldy);
    
            %   Piii + Piii -> P2iii
            %   Special because it's a square and thus has differnt deriv
            %   fiii -= 2 * Piii^2 * Agg(iii)
            %   f2iii += Piii^2 * Agg(iii)
            %   dfiiidiii -= 4 * Piii * Agg(iii)
            %   df2iiidiii += 2 * Piii * Agg(iii)
            dAggl2dy = 2 * y(iii) * settings.aKernel(sizeiii,sizeiii);
            dfdy(iii,iii) = dfdy(iii,iii) - 2 * dAggl2dy;
            iii2idx = sizeiii*2 + 1; %fixme
            dfdy(iii2idx, iii) = dfdy(iii2idx, iii) + dAggl2dy;
        end
    end

    %   emom
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
        flux = nthroot(9*pi / 1, 3) * k_inflow * ( (settings.pend-1) ^ (2/3) ) * w2(end)*y(idx_mu0-1) / delx / k_emom;
        % dfluxdP = 3 * k_inflow * ( settings.cutoff ^ (2/3) ) / delx / k_emom;
        dfluxdP = nthroot(9*pi / 1, 3) * k_inflow * ( (settings.pend-1) ^ (2/3) ) * w2(end) / delx / k_emom;

        %   Chemical species
        %   FIXME
        dCdA = 8 * delx / 9 * k_emom * diamM^3 * flux ...
            + k_emom * 8 / (delx^2 * 3) * y(idx_mu2);
        dCdP = 8 * delx / 9 * k_emom * diamM^3 * y(settings.gidx) * dfluxdP;
        dCdMu = k_emom / delx^2 * 8 / 3 * y(settings.gidx);
        
        dfdy(settings.gRxnIdx, settings.gidx) = dfdy(settings.gRxnIdx, settings.gidx) ...
            - dCdA * settings.gRxnCoeff';
        dfdy(settings.gRxnIdx, idx_mu0-1) = dfdy(settings.gRxnIdx, idx_mu0-1) ...
            - dCdP * settings.gRxnCoeff';
        dfdy(settings.gRxnIdx, idx_mu2) = dfdy(settings.gRxnIdx, idx_mu2) ...
            - dCdMu * settings.gRxnCoeff';

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