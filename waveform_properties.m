function waveform = waveform_properties(waveform, varargin)

ip = inputParser;
addParameter(ip, 'flag_plot', 0)
addParameter(ip, 'fig_basename', [])
addParameter(ip, 'index0', 1)
parse(ip, varargin{:})
flag_plot = ip.Results.flag_plot;
fig_basename = ip.Results.fig_basename;
index0 = ip.Results.index0;

ws = waveform(1).taxis(end) - waveform(1).taxis(1);%window size in ms
for ii = 1:length(waveform)
    % waveform shape
    chan = waveform(ii).chan;
    if index0
        average = waveform(ii).avgwaveform(chan+1,:);
        stdw = waveform(ii).stdwaveform(chan+1,:);
    else
        average = waveform(ii).avgwaveform(chan,:);
        stdw = waveform(ii).stdwaveform(chan,:);
    end
    average_interp = interp1(0:0.05:ws,average, 0:0.01:ws, 'spline');
    waveform(ii).average = average_interp;
    stdw_interp = interp1(0:0.05:ws,stdw, 0:0.01:ws, 'spline');
    waveform(ii).std = stdw_interp;
    waveform(ii).taxisu = -ws/2:0.01:ws/2;
    
    % get trough peak difference
    taxis = waveform(ii).taxisu;
    [~,trough_idx] = min(waveform(ii).average);
    t_trough = taxis(trough_idx);
    waveform(ii).t_trough = t_trough;
    [~,peak_idx] = max(waveform(ii).average(trough_idx:end));
    t_peak = taxis(peak_idx+trough_idx-1);
    waveform(ii).t_peak = t_peak;
    tpd = waveform(ii).t_peak - waveform(ii).t_trough;
    waveform(ii).tpd = tpd;
    
    if flag_plot
        figure('visible', 'off')
        plot_waveform(waveform(ii))
        saveas(gcf, sprintf('%s-%d-unit%d.jpg', fig_basename, ii, waveform(ii).unit))
        close
    end
end