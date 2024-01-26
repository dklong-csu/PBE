function [diams, PSDs] = PBElib_getDiscretePSDs(sol, mySettings, times)
%GETPSDS blah
%   FIXME write docs
PSDs = deval(sol,times,mySettings.pstart:mySettings.pend);
sizes = mySettings.sizes;
dsizes = mySettings.dsizes;
PSDs = PSDs ./ reshape(mySettings.dsizes,numel(dsizes),1);
binSize = diff(mySettings.pBins);
PSDs = PSDs ./ binSize;
ntimes = numel(times);
diams = reshape(sizes,numel(sizes),1) * ones(1,ntimes);
end