function strf = strf_properties(strf, dur,varargin)

p = inputParser;
addParameter(p,'plot_flag', 0);
addParameter(p,'save_plot', 0);
addParameter(p,'figure_folder', 0);
addParameter(p,'figure_name_base', 0);
parse(p,varargin{:});
plot_flag = p.Results.plot_flag;
save_plot = p.Results.save_plot;
figure_folder = p.Results.figure_folder;
figure_name_base = p.Results.figure_name_base;

TM = strf.tm;
SM = strf.sm;
taxis = strf.taxis;
faxis = strf.faxis;
mdb = strf.mdb;
pval = 0.002;

for ii = 1:length(strf)
    rf = strf(ii).rfcontra;
    n0 = strf(ii).n0contra;
    rfsig = significant_strf(rf, pval, n0, mdb, dur);
    [c, BF, BW_max, BW_min, Latency,Q, tw, sig, ex] = strf_bf_bw_latency(rfsig,taxis,faxis, TM, SM);
    
    strf(ii).contour = c;
    strf(ii).BF = BF;
    strf(ii).BW_max = BW_max;
    strf(ii).BW_min = BW_min;
    strf(ii).Latency = Latency;
    strf(ii).Q = Q;
    strf(ii).tw = tw;
    strf(ii).sig = sig;
    strf(ii).ex = ex;
    
    if plot_flag
        if any(sig)
            plot_strf(rf, n0, pval, mdb, dur, taxis, faxis, ...
                'timelabels', 0:25:100, 'freqlabels', 2.^(-2:4), ...
                'sig', 1, 'properties', strf(ii), ...
                'save_plot', save_plot, 'figure_folder', figure_folder, 'figure_name_base', figure_name_base)  
        else
            plot_strf(rf, n0, pval, mdb, dur, taxis, faxis, ...
                'properties', strf(ii),...
                'timelabels', 0:25:100, 'freqlabels', 2.^(-2:4), ...
                'save_plot', save_plot, 'figure_folder', figure_folder, 'figure_name_base', figure_name_base)
        end
        
        close all
    end
end