function [diams, PSDs] = cleanRawData(filename,nHeaderLines)
%cleanRawData   Cleans raw PSD data so it is better compared to simulations
%
%   FIXME DOCUMENTATION

raw_data = readmatrix(filename, ...
    "NumHeaderLines",nHeaderLines);

%   Column 1 contains diameters
diams = raw_data(:,1);
positive_rows = find(diams > 0);

%   Cleaned diameters
diams = diams(positive_rows);

%   Column 2 is time 0 which does not mean much
PSDs = raw_data(:,3:end);
PSDs = PSDs(positive_rows,:);
PSDs = max(0, PSDs);

%   Normalize PSDs so area=1
areas = trapz(diams, PSDs);
PSDs = PSDs./areas;
end