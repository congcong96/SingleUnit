function  [tBMF, sBMF, tMTF, sMTF, tmfaxisu, smfaxisu, tMTFu, sMTFu] = crh_parameters(mtfhist, tmfaxis, smfaxis, flag_plot)

if size(mtfhist, 1) == 1
    mtfhist = reshape(mtfhist, [length(smfaxis), length(smfaxis)]);
end

% Temporal and spectral MTFs
tMTF = nrm1(sum(mtfhist));
sMTF = nrm1(sum(mtfhist,2));
% Upsample tMTF, sMTF
tmfaxisu = linspace(min(tmfaxis),max(tmfaxis),1000);
tMTFu = interp1(tmfaxis,tMTF,tmfaxisu,'spline');
smfaxisu = linspace(min(smfaxis),max(smfaxis),1000);
sMTFu = interp1(smfaxis,sMTF,smfaxisu,'spline');
% get tBMF and sBMF
Max_tMTFu = max(tMTFu);
tBMF = tmfaxisu(tMTFu==Max_tMTFu);
Max_sMTFu = max(sMTFu);
sBMF = smfaxisu(sMTFu==Max_sMTFu);
% normalize MTF
tMTF = tMTF/Max_tMTFu;
tMTFu = tMTFu/Max_tMTFu;
sMTF = sMTF/Max_sMTFu;
sMTFu = sMTFu/Max_sMTFu;
if flag_plot
    % plot RTF
    figure('visible', 'off');
    figuresetup2savepdf(30, 10)
    cmap = cschemes('spectral', 21);
    Max = max(max(mtfhist));
    subplot(1,3,1)
    imagesc(tmfaxis, smfaxis, mtfhist,[0 Max]),shading flat,colormap(cmap) %axis square
    hold on
    plot(tBMF,sBMF,'g*','linewidth',2)
    axis xy
    xlabel('TMF (Hz)')
    ylabel('SMF (cyc/oct)')
    % plot tMTFu
    subplot(1,3,2)
    hold on;
    plot(tmfaxis,tMTF,'ko', 'MarkerFaceColor', 'k')
    plot(tmfaxisu,tMTFu,'k')
    line([tBMF tBMF],[0 1.2],'Color','k','LineWidth',1,'LineStyle',':');
    plot(tBMF,1,'r+','linewidth',2)
    xlim([tmfaxis(1) tmfaxis(end)])
    ylim([0 1])
    title(sprintf('TBMF %.2f(Hz)', tBMF))
    xlabel('TMF (Hz)')
     % plot sMTFu
    subplot(1,3,3)
    hold on;
    plot(smfaxisu,sMTFu,'k')
    plot(smfaxis,sMTF,'ko', 'MarkerFaceColor', 'k')
    line([sBMF sBMF],[0 1.2],'Color','k','LineWidth',1,'LineStyle',':');
    plot(sBMF,1,'r+','linewidth',2)
    xlim([0 smfaxis(end)])
    ylim([0 1])
    title(sprintf('SBMF %.2f(cyc/oct)', sBMF))
    xlabel('SMF (cyc/oct)')
end