clear variables
close all
clc

load("DATA_GOLDSIM_optimizeParameters.mat");

%%
%--------------------------------------------------------------------------
%   Plot the chemical species over time
%--------------------------------------------------------------------------

%   What times to plot
tplot = linspace(0,data_times(end),1000);
species = PBElib_getSpeciesConc(sol, mySettings, tplot);
%   Each row of 'species' corresponds to a different chemical species
figure
plot(tplot,species,'LineWidth',2)
% legend('Precursor','Solvated Precursor','Ligand')
xlabel("Time / s")
ylabel("Conc / mol/L")

%--------------------------------------------------------------------------
%   Plot PSDs over time (animated)
%--------------------------------------------------------------------------
figure
ax = gca;
h = animatedline('LineStyle','--',...
    'Marker','none',...
    'LineWidth',2,...
    'DisplayName',"Simulation",...
    'Color',[0 0.4470 0.7410]);
h2 = animatedline('LineStyle','--',...
    'Marker','none',...
    'LineWidth',2,...
    'DisplayName',"Data",...
    'Color',[0.8500 0.3250 0.0980]);
axis([0,15,0,1.0])
plotpts = data_times;
anim_seconds = 1;
[diams,PSDs] = PBElib_getPSDs(sol,mySettings,plotpts);
%   Volume weight the PSDs
PSDs = PSDs .* (diams.^3);
for iii=1:length(plotpts)
    clearpoints(h)
    clearpoints(h2)
    area = trapz(diams(:,iii),PSDs(:,iii));
    addpoints(h,diams(:,iii), PSDs(:,iii)/area);
    addpoints(h2,data_diam, data_PSDs(:,iii))
    title(ax,sprintf("PSD at %.f seconds",plotpts(iii)))
    legend('Location','northeast')
    drawnow
    pause(anim_seconds/length(plotpts))
end
%%
%--------------------------------------------------------------------------
%   Plot PSDs over time
%       FIXME
%       (a) Time as text box on northeast corner instead of title
%       (b) Shared axes so only display numbers and labels on outside
%       (c) Thicker box around border
%       (d) Also display CDF?
%   11 plots so:
%       Plot 1,5,9 are on left edge
%       Plot 1-4 on top
%       Plot 4,8,11 on right edge
%--------------------------------------------------------------------------
tlfig = figure;
set(tlfig,'Position',[-1647 72 1561 891]);
tl = tiledlayout('flow','TileSpacing','tight','Padding','loose');
for iii=1:length(plotpts)
    nexttile
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

    if ~(iii==1 || iii==5 || iii==9)
        yticklabels('')
    end

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
    % legend('Location','northeast')
    % title(sprintf("%.f seconds",plotpts(iii)))
    text(8,0.2,sprintf("%.f seconds",plotpts(iii)))
    box on
    % xlabel('Diameter / nm')
    % ylabel('Number Density / nm^{-1}')
    set(gca,'YColor','k')

    if ~(iii==4 || iii==8 || iii==11)
        yticklabels('')
    end

    if iii<8
        xticklabels('')
    end


    hold off
end
leg = legend();
leg.Layout.Tile = length(plotpts)+1;
xlabel(tl,'Diameter / nm')
ylabel(tl, 'Number Density / nm^{-1}')

ax = axes(tlfig);
han = gca;
han.Visible = 'off';
yyaxis(ax,'right')
ylabel('Cumulative Number Density / -')
han.YLabel.Visible = 'on';
han.YLabel.Color = 'k';


fontsize(gcf,scale=1.75)

%%
print(tlfig,'-depsc','-image',"psd_fit");