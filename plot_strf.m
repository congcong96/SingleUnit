function [tticks, fticks] = plot_strf(rf, n0, pval, mdb, dur, taxis, faxis, varargin)

p = inputParser;
addParameter(p, 'sigstrf', 0)
addParameter(p, 'timelabels', 0:25:100)
addParameter(p, 'freqlabels', 2.^(-2:4))
%addParameter(p, 'cboundary', n0/dur)
addParameter(p, 'contur', 'on')
addParameter(p, 'sig', 0)
addParameter(p, 'properties', 0)
addParameter(p,'save_plot', 0);
addParameter(p,'figure_folder', 0);
addParameter(p,'figure_name_base', 0);
parse(p,varargin{:});
sigstrf = p.Results.sigstrf;
timelabels =  p.Results.timelabels;
freqlabels =  p.Results.freqlabels;
save_plot = p.Results.save_plot;
figure_folder = p.Results.figure_folder;
figure_name_base = p.Results.figure_name_base;
%cboundary = p.Results.cboundary;
contur =  p.Results.contur;
sig =  p.Results.sig;
properties =  p.Results.properties;

fticks = [];
for ii = 1:numel(freqlabels)
    freqlabel = freqlabels(ii); % frequency in kHz
    fticks = [fticks find(faxis/1000 >= freqlabel, 1)];
end
tticks = [];
for ii = 1:numel(timelabels)
    timelabel = timelabels(ii); % time in ms
    tticks = [tticks find(taxis*1000 >= timelabel, 1)];
end

if sigstrf
    rf = significant_strf(rf, pval, n0, mdb, dur);
end
%rf = rf(fticks(1):fticks(end), tticks(1):tticks(end));
%faxis = faxis(fticks(1):fticks(end));
%rf = rf(:, tticks(1):tticks(end));
%taxis = taxis(tticks(1):tticks(end));
%cboundary = max(cboundary/20, max(abs(rf(:))));
cboundary = max(abs(rf(:)));
imagesc(rf);
cmap = cschemes('spectral', 21);
colormap(gca, cmap)

set(gca,'ydir', 'normal');
set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca, 'clim', [-1.05*cboundary-eps 1.05*cboundary+eps]);

xlabel('Time preceding spike (ms)')
xticks(tticks)
xticklabels(timelabels)
xlim(tticks([1 end]))

ylabel('Frequency (kHz)')
yticks(fticks)
yticklabels(freqlabels)
ylim(fticks([1 end]))

if strcmp(contur, 'on')
    rfsig = significant_strf(rf, pval, n0, mdb, dur);
    rfsigabs = abs(rfsig);
    rfsigabs(rfsigabs > 0) = 1;
    hold on
    contour(rfsigabs, 1, 'color', 0.2*[1 1 1 ])
end

if sig
    subfield = find(properties.sig == 1);
    for ii = 1:numel(subfield)
        ax1 = gca;
        fig = figure();
        copyobj(ax1, fig);
        isub = subfield(ii);
        plot_strf_properties(properties, isub, taxis, faxis, tticks, fticks)
        if save_plot
            saveas(gcf, fullfile(figure_folder, sprintf('%s-unit%d-subfield%d.jpg', figure_name_base, properties.unit, isub)))
        end
        close
    end
elseif save_plot
     
        saveas(gcf, fullfile(figure_folder, sprintf('%s-unit%d.jpg', figure_name_base, properties.unit)))
    close
end