function plot_su_strf_rtf_crh(probe, waveform, strf, rtfparams, crh, trigger, basename)

%% get probe layout
x = 0:0.02:0.8;
if contains(probe, 'Tetrode')
    % determine the location of waveforms in the plot
    if strcmp(probe, 'TetrodeB1x64')
        prbfile = '/data/congcong/SqMoPhys_Josh/prb/prbfiles/TetrodeB1x64.csv';
    elseif strcmp(probe, 'Tetrode1x64')
        prbfile = '/data/congcong/SqMoPhys_Josh/prb/prbfiles/Tetrode1x64.csv';
    end
    prblayout = readmatrix(prbfile);
    plotx = prblayout(:,1)/100;
    idx = mod(plotx, 1) > 0.3;
    plotx(idx) = floor(plotx(idx))+ 1;
    idx = mod(plotx, 1) > 0.1;
    plotx(idx) = floor(plotx(idx))+ 0.5;
    ploty = (prblayout(:,2)-50)/75;
    idx = mod(ploty, 1) > 0.4;
    ploty(idx) = floor(ploty(idx))+ 1;
    idx = mod(ploty, 1) > 0.2;
    ploty(idx) = floor(ploty(idx))+ 0.5;
elseif strcmp(probe, 'poly21x48')
    prbfile = '/data/congcong/SqMoPhys_Josh/prb/prbfiles/poly21x48.csv';
    prblayout = readmatrix(prbfile);
    prby =  prblayout(:,2);
    prby = sort(prby);
    ymid = prby(24);
    prblayout(prblayout(:,2) <= ymid, 1) = prblayout(prblayout(:,2) <= ymid, 1) + 2*86;
    prblayout(prblayout(:,2) > ymid, 2) = prblayout(prblayout(:,2) > ymid, 2) - ymid;
    plotx = prblayout(:,1)/86;
    ploty = prblayout(:,2)/50/4;
    prblayout = readmatrix(prbfile);
end
%%
mdb = strf(1).mdb;
taxis = strf(1).taxis;
faxis = strf(1).faxis;
fs = strf(1).fs;
dur = (trigger(end) - trigger(1))/fs;
pval = 0.002;
for ii = 1:length(waveform)
    figure
    figuresetup2savepdf(25, 25)
    %% plot waveform of signal from each channel
    subplot(4,3,1)
    avgwaveform = waveform(ii).avgwaveform;
    avgwaveform = avgwaveform/max(max(abs(avgwaveform)))/2;
    hold on
    for kk = 1:length(prblayout)
        %plot(x + plotx(kk), avgwaveform(kk,:) + ploty(kk))
        plot(x + plotx(kk), avgwaveform(kk,:) + ploty(kk))
    end
    position = abs(waveform(ii).position);
    positionidx = find(prblayout(:,1) == position(1) & prblayout(:,2) == position(2));
    positionx = plotx(positionidx);
    positiony = ploty(positionidx);
    text(positionx, positiony, '*', 'color', 'r')
    title(sprintf('%d spikes', strf(ii).n0contra))
    axis off
    hold off
    % waveform shape
    subplot(4,3,2)
    plot_waveform(waveform(ii))
    axis off
    
    % ACG
    subplot(4,3,3)
    bar(-50:0.5:50, waveform(ii).ACG, 1, 'k')
    xticks([-50 50])
    ymax = max(waveform(ii).ACG*1.1);
    if ymax < 1
        ymax = 1;
    end
    ylim([0 ymax])
    hold on
    plot([1.5 1.5], ylim, 'r--', 'linewidth', 1)
    plot([-1.5 -1.5], ylim, 'r--', 'linewidth', 1)
    title(sprintf('ISI violation: %.3f%%', waveform(ii).ISI_vio))
    
    %% plot STRF
    subplot(4,3,4)
    rf = strf(ii).rfcontra;
    n0 = strf(ii).n0contra;
    [tticks, fticks] = plot_strf(rf, n0, pval, mdb, dur, taxis, faxis, ...
        'timelabels', 0:25:100, 'freqlabels', 2.^(-2:4));
    isub = find(strf(ii).sig == 1, 1);
    % plot properties of dignificant feature
    if ~isempty(isub)
        plot_strf_properties(strf(ii), isub, taxis, faxis, tticks, fticks)
    end
    
    if isfield(strf, 'RI')
        title(sprintf('RI z:%d-p:%d', strf(ii).z, strf(ii).p))
    else
        title('STRF')
    end
    
    % plot temporal correlation index
    subplot(4,3,5)
    tci = rtfparams(ii).tci;
    plot(tci(:,1), tci(:,2))
    xlim([tci(1,1), tci(end,1)])
    xlabel('time shift (s)')
    title('temporal correlation')
    idx = find(tci(:,2)<0.368, 1);
    hold on
    if ~isempty(idx)
    plot(tci(idx,1)*[1 1], [0 tci(idx,2)], '--k')
    plot(tci(idx,1)*[0 1], tci(idx,2)*[1 1], '--k')
    text(0.015, 0.5, sprintf('time at 1/e : %.0fms', tci(idx,1)*1000))
    end
    % plot spectral correlation index
    subplot(4,3,6)
    sci = rtfparams(ii).sci;
    plot(sci(:,1), sci(:,2))
    xlim([sci(end,1), sci(1,1)])
    y = ylim;
    xlabel('frequency shift (oct)')
    title('spectral correlation')
    idx = [find(sci(:,2)>0.368, 1, 'first'), find(sci(:,2)>0.368, 1, 'last')];
    hold on
    if ~isempty(idx)
    plot(sci(idx(1),1)*[1 1], sci(idx(1),2)*[y(1), 1], '--k')
    plot(sci(idx(2),1)*[1 1], sci(idx(2),2)*[y(1), 1], '--k')
    plot(sci(idx,1), sci(idx(1),2)*[1 1], '--k')
    text(.5, 0.5, sprintf('%.2foct', sci(idx(1),1)-sci(idx(2),1)))
    end
    %% plot RTF
    RTF = squeeze(rtfparams(ii).rtf(:,:, end));
    Fm = rtfparams(ii).tmf;
    RD = rtfparams(ii).xmf;
    RTF4plot = RTF./sum(sum(RTF));
    cmap = cschemes('spectral', 21);
    Max = max(max(RTF4plot));
    subplot(4,3,7)
    imagesc(Fm,RD,RTF4plot,[0 Max]),shading flat,colormap(cmap) %axis square
    hold on
    tBMF = rtfparams(ii).RTFparam.tBMF;
    sBMF = rtfparams(ii).RTFparam.sBMF;
    plot(tBMF,sBMF,'k+','linewidth',2)
    axis xy
    xlabel('TMF (Hz)')
    ylabel('SMF (cyc/oct)')
    title('RTF')
    % plot tMTFu
    subplot(4,3,8)
    hold on;
    Fmu = rtfparams(ii).RTFparam.Fmu;
    tMTFu = rtfparams(ii).RTFparam.tMTFu;
    plot(Fmu,tMTFu,'k')
    line([tBMF tBMF],[0 1.2],'Color','k','LineWidth',1,'LineStyle',':');
    plot(tBMF,1,'r+','linewidth',2)
    xlim([-Fm(end) Fm(end)])
    ylim([0 1])
    title(sprintf('TBMF %.2f(Hz)', tBMF))
    xlabel('TMF (Hz)')
     % plot sMTFu
    subplot(4,3,9)
    hold on;
    RDu = rtfparams(ii).RTFparam.RDu;
    sMTFu = rtfparams(ii).RTFparam.sMTFu;
    plot(RDu,sMTFu,'k')
    line([sBMF sBMF],[0 1.2],'Color','k','LineWidth',1,'LineStyle',':');
    plot(sBMF,1,'r+','linewidth',2)
    xlim([0 RD(end)])
    ylim([0 1])
    title(sprintf('SBMF %.2f(cyc/oct)', sBMF))
    xlabel('SMF (cyc/oct)')
    
    %% plot CRH
    mtfhist = crh(ii).mtfhist;
    tmfaxis = crh(ii).tmfaxis;
    smfaxis = crh(ii).smfaxis;
    if size(mtfhist, 1) == 1
        mtfhist = reshape(mtfhist, [length(smfaxis), length(smfaxis)]);
    end
    Max = max(max(mtfhist));
    
    subplot(4,3,10)
    imagesc(tmfaxis, smfaxis, mtfhist,[0 Max]),shading flat,colormap(cmap) %axis square
    hold on
    tBMF = crh(ii).tBMF;
    plot(tBMF,sBMF,'k+','linewidth',2)
    axis xy
    xlabel('TMF (Hz)')
    ylabel('SMF (cyc/oct)')
    if isfield(crh, 'RI')
        title(sprintf('RI  z:%.2f  p:%.2f', crh(ii).z, crh(ii).p))
    else
        title('CRH')
    end
    % plot tMTFu
    subplot(4,3,11)
    hold on;
    tMTF = crh(ii).tMTF;
    sBMF = crh(ii).sBMF;
    tmfaxis = crh(ii).tmfaxis;
    tMTFu = crh(ii).tMTFu;
    tmfaxisu = crh(ii).tmfaxisu;
    
    plot(tmfaxis,tMTF,'ko', 'MarkerFaceColor', 'k')
    plot(tmfaxisu,tMTFu,'k')
    line([tBMF tBMF],[0 1.2],'Color','k','LineWidth',1,'LineStyle',':');
    plot(tBMF,1,'r+','linewidth',2)
    xlim([tmfaxis(1) tmfaxis(end)])
    ylim([0 1])
    title(sprintf('TBMF %.2f(Hz)', tBMF))
    xlabel('TMF (Hz)')
     % plot sMTFu
    subplot(4,3,12)
    hold on;
    sMTF = crh(ii).sMTF;
    smfaxis = crh(ii).smfaxis;
    sMTFu = crh(ii).sMTFu;
    smfaxisu = crh(ii).smfaxisu;
    plot(smfaxisu,sMTFu,'k')
    plot(smfaxis,sMTF,'ko', 'MarkerFaceColor', 'k')
    line([sBMF sBMF],[0 1.2],'Color','k','LineWidth',1,'LineStyle',':');
    plot(sBMF,1,'r+','linewidth',2)
    xlim([0 smfaxis(end)])
    ylim([0 1])
    title(sprintf('SBMF %.2f(cyc/oct)', sBMF))
    xlabel('SMF (cyc/oct)')
    
  
    saveas(gcf, sprintf('%s-unit%d.jpg', basename, waveform(ii).unit))
    close
    
end