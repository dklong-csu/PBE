function species = PBElib_getSpeciesConc(sol, mySettings, times)
%GETSPECIESCOC blah
%   TODO write docs

species = deval(times,sol,1:mySettings.pstart-1);

end