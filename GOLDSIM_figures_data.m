clc
clear variables
close all
%--------------------------------------------------------------------------
%   Create a plot with PSDs for each time a measurement is conducted
%--------------------------------------------------------------------------

%%
m2tpath = fullfile('\\home.rrze.uni-erlangen.de','bi70vahu','Documents','Programming','matlab2tikz','src');
addpath(m2tpath)

%%
[data_diam, data_PSDs] = GOLDSIM_cleanRawData("Au_quench_corediameter_hplc.txt",1);
data_times = [2,5,7,10,20,30,60,120,180,240,300];
%%
figure
hold on
for iii=1:length(data_times)
    plot(data_diam, data_PSDs(:,iii),...,
        'DisplayName',sprintf("%ds",data_times(iii)),...
        'LineStyle','--',...
        'LineWidth',2)
end
legend
xlabel('Diameter')
ylabel('Number Density')
box on
set(gca,'linewidth',2)
hold off

%%
cleanfigure
matlab2tikz(fullfile('.','figures','psd_data.tex'))