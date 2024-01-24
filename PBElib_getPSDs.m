function [diams, PSDs] = PBElib_getPSDs(sol, mySettings, times)
    [discreteDiams, discretePSDs] = PBElib_getDiscretePSDs(sol,mySettings,times);
    if mySettings.useEMoM
        [emomDiams, emomPSDs] = PBElib_geteMoMPSDs(sol, mySettings, times);
        emomDiams = emomDiams';
        emomPSDs = emomPSDs';
        
        diams = [discreteDiams;emomDiams];
        PSDs = [discretePSDs;emomPSDs];
    else
        diams = discreteDiams;
        PSDs = discretePSDs;
    end
    
end