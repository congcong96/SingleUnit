function plot_waveform(waveform)

taxis = waveform.taxisu;
t_peak = waveform.t_peak;
t_trough = waveform.t_trough;
tpd = waveform.tpd;

curve1 = waveform.average + waveform.std;
curve2 = waveform.average - waveform.std;
if waveform.tpd <= 0.35
    c1 = [0.4, 0.8, 1];
    c2 = 'b';
else
    c1 = [1, 0.8, 0.8];
    c2 = 'r';
end
a = fill([taxis, fliplr(taxis)], [curve1, fliplr(curve2)], 'r');
a.FaceColor = c1;
a.EdgeColor = [1 1 1];
hold on
plot(taxis, waveform.average, c2, 'linewidth', 3)
xlim([taxis(1) taxis(end)])
y = [min(waveform.average), max(waveform.average)];

% plot(t_trough*[1 1], y, 'k--', 'LineWidth', 2)
% plot(t_peak*[1 1], y, 'k--', 'LineWidth', 2)
% arrow([(t_peak+t_trough)/2, y(1)], [t_trough, y(1)])
% arrow([(t_peak+t_trough)/2, y(1)], [t_peak, y(1)])
% text((t_peak+t_trough)/2- tpd*.25, y(1)-10, sprintf('%.2fms', tpd), 'fontsize', 13)

box off