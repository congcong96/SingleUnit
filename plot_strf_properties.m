function plot_strf_properties(properties, isub, taxis, faxis, tticks, fticks)
% draw arrow for best frequency and Latency
BF = properties.BF(isub);
Latency = properties.Latency(isub);
idx_f = find(faxis>=BF*1000, 1);
idx_t = find(taxis>=Latency/1000, 1);
arrow([idx_t idx_f], [tticks(1) idx_f], 'width', 2, 'color', [.5 .5 .5])
text((idx_t-tticks(1))/2, idx_f+30, sprintf('%.1fkHz', BF))
arrow([idx_t idx_f], [idx_t fticks(1)], 'width', 2, 'color', [.5 .5 .5])
text(idx_t+5, idx_f/2, sprintf('%.1fms', Latency))
% draw contour of response
c = properties.contour{isub};
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
plot([top_t, right_t + 20], [top_f, top_f], 'k--', 'linewidth', 1)
plot([bot_t, right_t + 20], [bot_f, bot_f], 'k--', 'linewidth', 1)
arrow([right_t + 20 top_f+10], [right_t + 20 top_f], 'width', 1, 'color', 'k')
arrow([right_t + 20 bot_f-10], [right_t + 20 bot_f], 'width', 1, 'color', 'k')
text(right_t+20, (top_f + bot_f)/2, sprintf('BW = %.1foct', log2(BW_max/BW_min)))
if isub == properties.ex
    scatter(idx_t, idx_f, '*g')
end