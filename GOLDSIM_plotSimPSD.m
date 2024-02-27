clear variables
close all
clc


%%
% load("DATA_GOLDSIM_optimizeParameters_5percent.mat");
load("DATA_GOLDSIM_optimizeParameters_test.mat");
save_folder = fullfile('.','figures','psdfit_test');
save_file_root = 'psdfit';

%%  Plot the PSD that was optimized or the more accurate version?
plot_accurate = false;
if plot_accurate
    sol = solAccurate;
    mySettings = mySettingsAccurate;
else
    sol = solOptimize;
    mySettings = mySettingsOptimize;
end

%%
%--------------------------------------------------------------------------
%   Plot the chemical species over time
%--------------------------------------------------------------------------

%   What times to plot
tplot = linspace(0,data_times(end),1000);
species = PBElib_getSpeciesConc(sol, mySettings, tplot);
%   Each row of 'species' corresponds to a different chemical species
% figure
% plot(tplot,species,'LineWidth',2)
% legend('Precursor','Solvated Precursor','Ligand')
% xlabel("Time / s")
% ylabel("Conc / mol/L")


plotpts = data_times;
[diams,PSDs] = PBElib_getPSDs(sol,mySettings,plotpts);
%   Volume weight the PSDs
PSDs = PSDs .* (diams.^3);

%%
%--------------------------------------------------------------------------
%   Plot PSDs over time
%--------------------------------------------------------------------------

for iii=1:length(plotpts)
    figure
    area = trapz(diams(:,iii),PSDs(:,iii));
    psd = PSDs(:,iii)/area;
    cdf_sim = cumtrapz(diams(:,iii), psd);
    cdf_data = cumtrapz(data_diam, data_PSDs(:,iii));
    hold on
    plot(diams(:,iii), psd, ...
        'Color',[0 0.4470 0.7410],...
        'LineStyle','--',...
        'LineWidth',2,...
        'DisplayName','Simulated PSD')
    plot(data_diam, data_PSDs(:,iii),...
        'Color',[0.8500 0.3250 0.0980],...
        'LineStyle','--',...
        'LineWidth',2,...
        'DisplayName','Measured PSD')
    ylim([0 1])


    yyaxis right
    plot(diams(:,iii), cdf_sim,...
        'Color',[0.4940 0.1840 0.5560],...
        'LineStyle','--',...
        'LineWidth',2,...
        'DisplayName','Simulated CDF')
    plot(data_diam, cdf_data,...
        'Color',[0.6350 0.0780 0.1840],...
        'LineStyle','--',...
        'LineWidth',2,...
        'DisplayName','Measured CDF')
    ylim([0 1.1])
    text(8,0.2,sprintf("%.f seconds",plotpts(iii)))
    box on
    set(gca,'YColor','k')
    legend('Location','northoutside','NumColumns',2)


    fontsize(gcf,scale=1.75)
    hold off
    filenamePic = fullfile(save_folder, strcat(save_file_root,num2str(iii)));
    print('-depsc','-image',filenamePic)
    saveFigureData(gcf, fullfile(save_folder,"psdfit_"), sprintf("_%s.dat",num2str(iii)));
    % cleanfigure;
    % filenameTex = strcat(filenamePic,'.tex');
    % figData = findobj(fig,'-property','YData');
    % matlab2tikz(filenameTex)
    % close all
end


%%  write curves in matlab figure to a data file
function saveFigureData(fig, fileNamePrefix, fileNameSuffix)
    figData = findobj(fig, '-property','YData'); %  find everything with y-values
    %   Loop through all relevant data 
    for iii=1:numel(figData)
        %   Get display name of curve to give meaningful name to data file
        ID = figData(iii).DisplayName;
        ID = erase(ID," ");
        %   Get the plotted x-y data
        xvals = figData(iii).XData;
        yvals = figData(iii).YData;
        %   Construct a meaningful filename
        filename = strcat(fileNamePrefix,ID,fileNameSuffix);
        %   Print to the file
        fileID = fopen(filename,'w');   %   write mode deletes previous contents
        fprintf(fileID,"%.10e\t%.10e\n",[xvals(:)';yvals(:)']);
        fclose(fileID);
    end

end