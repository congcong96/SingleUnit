function crh = batch_crh_RI(crh, sprtmf, sprsmf, spktrain, varargin)

p = inputParser;
addParameter(p,'plot_flag', 0);
addParameter(p,'fig_basename', 0);
addParameter(p,'niter', 1000);
addParameter(p,'nblocks', 10);
parse(p,varargin{:});
plot_flag = p.Results.plot_flag;
fig_basename = p.Results.fig_basename;
niter = p.Results.niter;
nblocks = p.Results.nblocks;

for ii = 1:length(crh)
    [RI, z, p, crhcorr, ncrhcorr] = calc_CRH_RI(sprtmf, sprsmf, spktrain(ii,:),...
        'plot_flag', plot_flag, 'niter', niter, 'nblocks', nblocks);
    crh(ii).RI = RI;
    crh(ii).z = z;
    crh(ii).p = p;
    crh(ii).crhcorr = crhcorr;
    crh(ii).ncrhcorr = ncrhcorr;
    if plot_flag
        saveas(gcf, sprintf('%s-unit%d.jpg', fig_basename, crh(ii).unit))
        close
    end
end
