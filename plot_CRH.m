function plot_CRH(crh, tmf, smf, colorflag)

assert(size(crh, 1)*size(crh,2) == (length(tmf)) * (length(smf)))
if ~exist('colorflag', 'var')
    colorflag = 'r';
end
if size(crh, 1) == 1
    crh = reshape(crh, [length(smf), length(tmf)]);
end

imagesc(tmf, smf, crh)
if strcmp(colorflag, 'r')
    cmap = cbrewer('seq', 'Reds', 9);
elseif strcmp(colorflag, 'b')
    cmap = cbrewer('seq', 'Blues', 9);
elseif strcmp(colorflag, 'p')
    cmap = cbrewer('seq', 'Purples', 9);
end
colormap(gca, cmap)
clim = max(abs(crh(:)));
caxis([0 clim]);
axis xy

yticks(1:4)

xticks(-60:30:60)

xlabel('TMF (Hz)')
ylabel('SMF (oct/cyc)')