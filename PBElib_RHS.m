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


% pstart = settings.pstart;
% pend = settings.pend;
% dGrowth = y(settings.gidx) * y(pstart:pend) .* settings.gKernel';
% binSizes = reshape(settings.pBins,numel(settings.pBins),1); %  make sure it's a column vector 
% w1 = (binSizes-1) ./ binSizes;
% w2 = 1 ./ binSizes;
% %   Loss due to reaction
% dy(pstart:pend) = dy(pstart:pend) - dGrowth;
% %   Gain due to reaction
% dy(pstart+1:pend) = dy(pstart+1:pend) + w2(1:end-1) .* dGrowth(1:end-1);
% dy(pstart:pend) = dy(pstart:pend) + w1(1:end) .* dGrowth(1:end);
% dy(settings.gRxnIdx) = dy(settings.gRxnIdx) - settings.gRxnCoeff' .* sum(dGrowth);


%--------------------------------------------------------------------------
%   Particle aggregation
%--------------------------------------------------------------------------
aggregIdx = settings.pstart:settings.cutoff;
for iii=1:length(aggregIdx)
    for jjj=iii:length(aggregIdx)
        createdIdx = settings.aggrInfo{iii,jjj}(:,1);
        createdWeight = settings.aggrInfo{iii,jjj}(:,2);
        %   Loss due to Piii+Pjjj
        if iii==jjj
            rate = y( aggregIdx(iii) )^2 * settings.aKernel(iii,jjj);
            dy( aggregIdx(iii) ) = dy( aggregIdx(iii) ) ...
                - ( 1 + 1/binSizes(iii)) * rate;
        else
            rate = y( aggregIdx(iii) ) * y( aggregIdx(jjj) ) * settings.aKernel(iii,jjj);
            dy( aggregIdx(iii) ) = dy( aggregIdx(iii) ) - rate;
            dy( aggregIdx(jjj) ) = dy( aggregIdx(jjj) ) - rate;
        end
        %   Created particles
        dy( createdIdx ) = dy( createdIdx ) + rate * createdWeight;
    end
end
% for iii=pstart:settings.cutoff
%     idxiii = iii-pstart+1;
%     for jjj=iii:settings.cutoff
%         idxjjj = jjj-pstart+1;
%         %   Pi + Pj ->
%         sizeiii = iii-1; %fixme
%         sizejjj = jjj-1; %fixme
%         agglomInfo = settings.aMapping{idxiii,idxjjj};
%         idxCreated = agglomInfo(:,1);
%         weightCreated = agglomInfo(:,2);
%         created = sizejjj + sizeiii + 1; %fixme
%         dAggl = y(iii) * y(jjj) .* settings.aKernel(idxiii,idxjjj);
%         %   Loss Pi
%         dy(iii) = dy(iii) - dAggl;
%         %   Loss Pj
%         if jjj==iii
%             dy(jjj) = dy(jjj) - binSizes(idxiii)*dAggl;
%         else
%             dy(jjj) = dy(jjj) - dAggl;
%         end
%         %   Gain wkPk
%         dy(idxCreated) = dy(idxCreated) + weightCreated.*dAggl;
%         % dy(iii) = dy(iii) - sum(dAggl);
%         % dy(jjj) = dy(jjj) - dAggl;
%         % dy(created) = dy(created) + dAggl;
%     end
% end


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
    
    % flux = 3 * k_inflow * ( (settings.pend-1) ^ (2/3) ) * y(idx_mu0-1) / delx / k_emom;
    flux = nthroot(9*pi / 1, 3) * k_inflow * ( (settings.pBins(end)) ^ (2/3) ) * w2(end)*y(idx_mu0-1) / delx / k_emom;
    
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









