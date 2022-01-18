function plot_strf(rf, n0, pval, mdb, dur, taxis, faxis, varargin)

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
rf = rf(:, tticks(1):tticks(end));
taxis = taxis(tticks(1):tticks(end));
%cboundary = max(cboundary/20, max(abs(rf(:))));
cboundary = max(abs(rf(:)));
imagesc(rf);
cmap = bluewhitered(256);
colormap(gca, cmap)

set(gca,'ydir', 'normal');
set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca, 'clim', [-1.05*cboundary-eps 1.05*cboundary+eps]);

xlabel('Time preceding spike (ms)')
xticks(tticks-tticks(1)+1)
xticklabels(timelabels)

ylabel('Frequency (kHz)')
yticks(fticks-fticks(1)+1)
yticklabels(freqlabels)

if strcmp(contur, 'on')
    rfsig = significant_strf(rf, pval, n0, mdb, dur);
    rfsigabs = abs(rfsig);
    rfsigabs(rfsigabs > 0) = 1;
    hold on
    contour(rfsigabs, 1, 'color', 'y')
end

if sig
    subfield = find(properties.sig == 1);
    for ii = 1:numel(subfield)
        ax1 = gca;
        fig = figure();
        copyobj(ax1, fig);
        isub = subfield(ii);
        % draw arrow for best frequency and Latency
        BF = properties.BF(isub);
        Latency = properties.Latency(isub);
        idx_f = find(faxis>=BF*1000, 1);
        idx_t = find(taxis>=Latency/1000, 1);
        arrow([idx_t idx_f], [0 idx_f], 'width', 2, 'color', [.5 .5 .5])
        text(idx_t/2, idx_f+10, sprintf('BF = %.1fkHz', BF))
        arrow([idx_t idx_f], [idx_t 0], 'width', 2, 'color', [.5 .5 .5])
        text(idx_t+5, idx_f/2, sprintf('Latency = %.1fms', Latency))
        % draw contour of response
        c = properties.contour{isub};
        c(1,:) = c(1,:)-tticks(1)+1;
        %c(2,:) = c(2,:)-fticks(1)+1;
        plot(c(1,:),c(2,:), 'k', 'linewidth', 2)
        % draw bandwidth
        BW_max = properties.BW_max(isub);
        BW_min = properties.BW_min(isub);
        [~, top_idx] = max(c(2,:));
        top_t = c(1, top_idx);
        top_f = find(faxis>=BW_max*1000, 1);
        [~, bot_idx] = min(c(2,:));
        bot_t = c(1, bot_idx);
        bot_f = find(faxis>=BW_min*1000, 1);
        right_t = max(c(1,:));
        plot([top_t, right_t + 20], [top_f, top_f], 'k--', 'linewidth', 2)
        plot([bot_t, right_t + 20], [bot_f, bot_f], 'k--', 'linewidth', 2)
        arrow([right_t + 20 top_f+10], [right_t + 20 top_f], 'width', 2, 'color', 'k')
        arrow([right_t + 20 bot_f-10], [right_t + 20 bot_f], 'width', 2, 'color', 'k')
        text(right_t+20, (top_f + bot_f)/2, sprintf('BW = %.1foct', log2(BW_max/BW_min)))
        if isub == properties.ex
            scatter(idx_t, idx_f, '*g')
        end
        if save_plot
            saveas(gcf, fullfile(figure_folder, sprintf('%s-unit%d-subfield%d.jpg', figure_name_base, properties.unit, isub)))
        end
        close
    end
else
    if save_plot
        saveas(gcf, fullfile(figure_folder, sprintf('%s-unit%d.jpg', figure_name_base, properties.unit)))
    end
    close
end