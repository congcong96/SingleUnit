function waveform = waveform_properties(waveform, flag_plot, fig_basename)

if nargin == 1
    flag_plot = 0;
end

for ii = 1:length(waveform)
    % waveform shape
    chan = waveform(ii).chan;
    average = waveform(ii).avgwaveform(chan+1,:);
    average_interp = interp1(0:0.05:2,average, 0:0.01:2, 'spline');
    waveform(ii).average = average_interp;
    stdw = waveform(ii).stdwaveform(chan+1,:);
    stdw_interp = interp1(0:0.05:2,stdw, 0:0.01:2, 'spline');
    waveform(ii).std = stdw_interp;
    waveform(ii).taxisu = -1:0.01:1;
    
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
        plot_waveform(waveform(ii))
        saveas(gcf, sprintf('%s-unit%d.jpg', fig_basename, waveform(ii).unit))
        close
    end
end