function [diams, PSDs] = PBElib_geteMoMPSDs(sol, mySettings, times)
%GETEMOMPSDS blah
%   FIXME write docs

if ~mySettings.useEMoM
    warning("EMOM PSDs requested but EMOM was not used in the simulation.\nReturning empty arrays.")
    diams = [];
    PSDs = [];
else
    n_xi = 1000; %FIXME? number of tau as arg?
    diams = zeros( length(times), n_xi );
    PSDs = zeros( length(times), n_xi );
    for iii=1:length(times)
        ttt=times(iii);

        %   Get diameters
        %   xi = x0 + int_ti^T G( A(s) ) ds
        taui = linspace(ttt,0.0,n_xi+1 ); %
        taui = taui(2:end);
        prefactor = 8/9 * mySettings.emomDelx * mySettings.emomGrowthRate;
        G = @(t) prefactor * deval(sol,t,mySettings.gidx);
        for jjj=1:length(taui)
            ti = taui(jjj);
            diams(iii,jjj) = mySettings.emomBigParticleDiam ...
                + integral(G,ti,ttt);
        end

        %   Get density
        %   FIXME size of max particle in atoms
        PSDs(iii,:) = nthroot(9*pi, 3) / mySettings.emomDelx ...
            * mySettings.emomInflowRate / mySettings.emomGrowthRate ...
            * (mySettings.pend-1)^(2/3) * deval(sol,taui,mySettings.pend);
    end
end
end