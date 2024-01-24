function dy = PBElib_RHS(~,y,settings) %#codegen
    dy = zeros(length(y),1);
    %   Chemical reactions
    [n_rxnrows, ~] = size(settings.A);
    dy(1:n_rxnrows) = sum( -settings.A .* ( prod( y(1:n_rxnrows) .^ max(0,settings.A) ) .* settings.K ) ,2);



    %   Particle Growth
    pstart = settings.pstart;
    pend = settings.pend;
    dGrowth = y(settings.gidx) * y(pstart:pend) .* settings.gKernel';
    dy(pstart:pend) = dy(pstart:pend) - dGrowth;
    dy(pstart+1:pend) = dy(pstart+1:pend) + dGrowth(1:end-1);
    dy(settings.gRxnIdx) = dy(settings.gRxnIdx) - settings.gRxnCoeff' .* sum(dGrowth);

    %   Particle Agglomeration
    for iii=pstart:settings.cutoff
        sizeiii = iii-1; %fixme
        jjj = iii:settings.cutoff;
        sizejjj = jjj-1; %fixme
        created = sizejjj + sizeiii + 1; %fixme
        dAggl = y(iii) * y(jjj) .* settings.aKernel(sizejjj,sizeiii);

        dy(iii) = dy(iii) - sum(dAggl);
        dy(jjj) = dy(jjj) - dAggl;
        dy(created) = dy(created) + dAggl;
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
        
        % flux = 3 * k_inflow * ( (settings.pend-1) ^ (2/3) ) * y(idx_mu0-1) / delx / k_emom;
        flux = nthroot(9*pi / 1, 3) * k_inflow * ( (settings.pend-1) ^ (2/3) ) * y(idx_mu0-1) / delx / k_emom;
        
        %   Effect on chemical species due to growth reaction
        %   FIXME 8/9 and so on hardcoded
        growthRate = 8/9 * delx * y(settings.gidx) * k_emom ...
            * diamM^3 * flux ...
            + k_emom * 8 / (delx^2 * 3) * y(settings.gidx) * y(idx_mu2);
        dy(settings.gRxnIdx) = dy(settings.gRxnIdx) ...
            - growthRate * settings.gRxnCoeff';


        %   Effect on moments
        diamMpow = (0:3)';
        dy(idxs_emom) = dy(idxs_emom) ...
            + 8 * delx / 9 * y(settings.gidx) * k_emom * flux * diamM.^diamMpow ...
            + 8 * k_emom * delx / 9 * y(settings.gidx) * diamMpow .* y(idxs_emom-1);
    end
end









