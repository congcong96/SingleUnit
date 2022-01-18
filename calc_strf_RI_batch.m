function strf = calc_strf_RI_batch(strf, spktrain, stimulus, nlags, varargin)

p = inputParser;
addParameter(p,'niter', 1000);
addParameter(p,'plot_flag', 0);
addParameter(p,'figure_folder', 0);
addParameter(p,'figure_name_base', 0);
parse(p,varargin{:});
niter = p.Results.niter;
plot_flag = p.Results.plot_flag;
figure_folder = p.Results.figure_folder;
figure_name_base = p.Results.figure_name_base;

for ii = 1:length(strf)
    figure_name = sprintf('%s-unit%d.jpg', figure_name_base, strf(ii).unit);
    [RI_sig_z, RI_sig_p, RI, z, p, strfcorr, nstrfcorr] = calc_strf_RI(stimulus, spktrain(ii,:), nlags,...
        'plot_flag', plot_flag, 'figure_folder', figure_folder, 'figure_name', figure_name, 'niter', niter);
    strf(ii).RI_sig_z = RI_sig_z;
    strf(ii).RI_sig_p = RI_sig_p;
    strf(ii).RI = RI;
    strf(ii).z = z;
    strf(ii).p = p;
    strf(ii).strfcorr = strfcorr;
    strf(ii).nstrfcorr = nstrfcorr;
    

end



