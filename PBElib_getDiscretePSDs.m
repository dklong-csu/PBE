function [diams, PSDs] = PBElib_getDiscretePSDs(sol, mySettings, times)
%GETPSDS blah
%   FIXME write docs
PSDs = deval(sol,times,mySettings.pstart:mySettings.pend);
PSDs = PSDs ./ reshape(mySettings.dsizes,numel(mySettings.dsizes),1);
ntimes = numel(times);
diams = reshape(mySettings.sizes,numel(mySettings.sizes),1) * ones(1,ntimes);
end