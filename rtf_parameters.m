%% [RTFparam]=rtf_parameters04(RTF)
% this calculates parameters for RTF
% the original is a part of function rtf from James
% Natsumi 2Oct17
% updated [RTFparam]=rtf_parameters(RTF)
% to analyze more closely
% classify based on the shape of mtf.
% calculate best SM and TM from the center 
% and BW from the bottom of the plotting
% Natsumi 25Feb18 
% updated rtf_parameters02 for removing 2nd calculatin and use 1st min and
% 2nd max and calculate bw from all the points above thershold 0.7
% also use the largest peak for the best s/tMTF 
% Natsumi 26May2018
% updated rtf_parameters03 to calcualte bTMF and bSMF acording to ﻿Atencio
% and Schreiner 2016 Neurosciene
% Jan042019 Natsumi

function  [RTFparam]=rtf_parameters(RTF,RD,Fm, flag_plot)

if nargin == 3
    flag_plot = 0;
end

% Fold across temporal axis
% RTFn = flip(RTF(:,1:find(Fm==0)),2);
% RTFp=RTF(:,find(Fm==0):end);
% RTF = (RTFn + RTFp)/2;
% Fm = Fm(find(Fm==0):end);
% RD = RD(find(RD==0):end);

% Temporal and spectral MTFs
tMTF = nrm1(sum(RTF));
sMTF = nrm1(sum(RTF,2));

% Upsample tMTF, sMTF
Fmu = linspace(min(Fm),max(Fm),1000);
tMTFu = interp1(Fm,tMTF,Fmu,'spline');
RDu = linspace(min(RD),max(RD),1000);
sMTFu = interp1(RD,sMTF,RDu,'spline');

% Temporal BP, BW, BMF ----------------------------------------------------

% tMTF_BP
idx07 = find(tMTFu>.7);
diffidx07 = diff(idx07);
f = find(diffidx07>1);
Max_tMTFu = max(tMTFu);

if isempty(f)
    Max_tMTFu1 = max(tMTFu);
    Max_tMTFu2 = NaN;% there's no second best temporal modulation frequency
    Max_tMTFu12 = 1;
    
    if min(tMTFu(1:find(tMTFu==Max_tMTFu1)))<=.7 && min(tMTFu(find(tMTFu==Max_tMTFu1):end))<=.7
        %on both sides of the best temporal modulation, there're data points with value <0.7 
        tMTF_LP = false;
        tMTF_BP = true;
        tMTF_HP = false;
    elseif min(tMTFu(1:find(tMTFu==Max_tMTFu1))) > .7 && min(tMTFu(find(tMTFu==Max_tMTFu1):end))<=.7
        tMTF_LP = true;
        tMTF_BP = false;
        tMTF_HP = false;
    elseif min(tMTFu(1:find(tMTFu==Max_tMTFu1))) <= .7 && min(tMTFu(find(tMTFu==Max_tMTFu1):end)) > .7
        tMTF_LP = false;
        tMTF_BP = false;
        tMTF_HP = true;
    end
    
elseif length(f) == 1
    Max_tMTFu1 = max(tMTFu(idx07(1):idx07(f))); 
    Max_tMTFu2 = max(tMTFu(idx07(f+1):idx07(end)));
    
    if (sum(tMTFu(idx07(1):idx07(f)))-.7*(idx07(f)-idx07(1)+1)) > (sum(tMTFu(idx07(f+1):idx07(end)))-.7*(idx07(end)-idx07(f+1)+1))
        Max_tMTFu12 = 1;
    else
        Max_tMTFu12 = 2;
    end
    
    if min(tMTFu(1:find(tMTFu==Max_tMTFu1))) > .7 && min(tMTFu(find(tMTFu==Max_tMTFu1):find(tMTFu==Max_tMTFu2)))<=.7 && min(tMTFu(find(tMTFu==Max_tMTFu2):end))<=.7
        tMTF_LP = true;
        tMTF_BP = true;
        tMTF_HP = false;
    elseif min(tMTFu(1:find(tMTFu==Max_tMTFu1))) <= .7 && min(tMTFu(find(tMTFu==Max_tMTFu1):find(tMTFu==Max_tMTFu2)))<=.7 && min(tMTFu(find(tMTFu==Max_tMTFu2):end)) > .7
        tMTF_LP = false;
        tMTF_BP = true;
        tMTF_HP = true;
    elseif min(tMTFu(1:find(tMTFu==Max_tMTFu1))) <= .7 && min(tMTFu(find(tMTFu==Max_tMTFu1):find(tMTFu==Max_tMTFu2)))<=.7 && min(tMTFu(find(tMTFu==Max_tMTFu2):end)) <= .7
        tMTF_LP = false;
        tMTF_BP = true;
        tMTF_HP = false;
    elseif min(tMTFu(1:find(tMTFu==Max_tMTFu1))) > .7 && min(tMTFu(find(tMTFu==Max_tMTFu1):find(tMTFu==Max_tMTFu2)))<=.7 && min(tMTFu(find(tMTFu==Max_tMTFu2):end)) > .7
        tMTF_LP = true;
        tMTF_BP = false;
        tMTF_HP = true;
    end
    
else
    Max_tMTFu1 = max(tMTFu);
    Max_tMTFu2 = NaN;
    Max_tMTFu12 = NaN;
    
    tMTF_LP = false;
    tMTF_BP = false;
    tMTF_HP = false;
    
end

% tMTF_BW
if ~isnan(Max_tMTFu1)
    tx = find(tMTFu==Max_tMTFu1);
    while tMTFu(tx) > .7
        if tx == length(tMTFu)
            break
        end
        tx=tx+1;
    end
    t2=tx; %get the point when tMTFu >= 0.7 right of BTMF
    % if tMTF_BP % commented out lines 79 and 88 to 90 Natsumi 20July2017
    tx = find(tMTFu==Max_tMTFu1);
    while tMTFu(tx) > .7
        if tx == 1
            break
        end
        tx=tx-1;
    end
    t1=tx;
    % else
    %     t1=1;
    % end
    tMTF_BW1=Fmu(t2)-Fmu(t1);
    tMTF_Max1=Fmu(t2); % added 20Jul2017 Natsumi
    tMTF_Min1=Fmu(t1); % added 20Jul2017 Natsumi
    
    % tBMF
%     tBMF1 = Fmu(find(tMTFu==Max_tMTFu1)); % use this for both case % Natsumi 20Jul2017 commented out lines 95 to 99
    if tMTF_LP
        tBMF1 = Fmu(round(t2/2));
    elseif tMTF_HP
        tBMF1 = Fmu(round(t1+(t2-t1)/2));
    else
        tBMF1 = Fmu(find(tMTFu==Max_tMTFu1));
    end
    
else
    tMTF_BW1=NaN;
    tMTF_Max1=NaN;
    tMTF_Min1=NaN;
    tBMF1 = NaN;
end

% tMTF_BW for the 2nd peak
if ~isnan(Max_tMTFu2)
    tx = find(tMTFu==Max_tMTFu2);
    while tMTFu(tx) > .7
        if tx == length(tMTFu)
            break
        end
        tx=tx+1;
    end
    t2=tx;
    % if tMTF_BP % commented out lines 79 and 88 to 90 Natsumi 20July2017
    tx = find(tMTFu==Max_tMTFu2);
    while tMTFu(tx) > .7
        if tx == 1
            break
        end
        tx=tx-1;
    end
    t1=tx;
    % else
    %     t1=1;
    % end
    tMTF_BW2=Fmu(t2)-Fmu(t1);
    tMTF_Max2=Fmu(t2);
    tMTF_Min2=Fmu(t1);
    
    % tBMF
    %tBMF2 = Fmu(find(tMTFu==Max_tMTFu2)); % use this for both case % Natsumi 20Jul2017 commented out lines 95 to 99
    if  tMTF_HP
        assert(tMTF_LP == true || tMTF_BP == true)
        tBMF2 = Fmu(round(t1+(t2-t1)/2));
    else
        assert(tMTF_BP == true)
        tBMF2 = Fmu(find(tMTFu==Max_tMTFu2));
    end
else
    tMTF_BW2=NaN;
    tMTF_Max2=NaN;
    tMTF_Min2=NaN;
    tBMF2 = NaN;
end

if ~isempty(f)
    % calculate from all above 0.7
    dummy_tMTF_BW = zeros(length(f)+1,1);
    for idxf = 1:length(f)+1
        if idxf == 1
            %         idx07(1)-1
            %         idx07(f(idxf))+1
            try dummy_tMTF_BW(idxf) = Fmu( idx07(f(idxf))+1 ) - Fmu(idx07(1)-1);
            catch dummy_tMTF_BW(idxf) = Fmu( idx07(f(idxf))+1 ) - Fmu(idx07(1));
            end
        elseif idxf > 1 & idxf < length(f)+1
            %         idx07(f(idxf-1)+1)-1
            %         idx07(f(idxf))+1
            dummy_tMTF_BW(idxf) = Fmu( idx07(f(idxf))+1 ) - Fmu( idx07(f(idxf-1)+1)-1 );
        elseif  idxf == length(f)+1
            %         idx07(f(idxf-1)+1)-1
            %         idx07(end)+1
            try dummy_tMTF_BW(idxf) = Fmu( idx07(end)+1 ) - Fmu(  idx07(f(idxf-1)+1)-1 );
            catch dummy_tMTF_BW(idxf) = Fmu( idx07(end) ) - Fmu(  idx07(f(idxf-1)+1)-1 );
            end
        end
    end
    tMTF_BW= sum(dummy_tMTF_BW);
    try tMTF_Max = Fmu( idx07(end)+1 );
    catch tMTF_Max = Fmu( idx07(end) );
    end
    try tMTF_Min= Fmu( idx07(1)-1 );
    catch tMTF_Min = Fmu( idx07(1) );
    end
    tBMF = Fmu(find(tMTFu==Max_tMTFu));
    
else
    
    try tMTF_Max = Fmu( idx07(end)+1 );
    catch tMTF_Max = Fmu( idx07(end) );
    end
    try tMTF_Min= Fmu( idx07(1)-1 );
    catch tMTF_Min = Fmu( idx07(1) );
    end
    tMTF_BW= tMTF_Max-tMTF_Min;
    tBMF = Fmu(find(tMTFu==Max_tMTFu));
    
end

% Spectral BP, BW, BMF ----------------------------------------------------

% tMTF_BP
idx07 = find(sMTFu>.7);
diffidx07 = diff(idx07);
f = find(diffidx07>1);
Max_sMTFu = max(sMTFu);

if isempty(f)
    Max_sMTFu1 = max(sMTFu);
    Max_sMTFu2 = NaN;
    Max_sMTFu12 = 1;
    
    if min(sMTFu(1:find(sMTFu==Max_sMTFu1)))<=.7 && min(sMTFu(find(sMTFu==Max_sMTFu1):end))<=.7
        sMTF_LP = false;
        sMTF_BP = true;
        sMTF_HP = false;
    elseif min(sMTFu(1:find(sMTFu==Max_sMTFu1))) > .7 && min(sMTFu(find(sMTFu==Max_sMTFu1):end))<=.7
        sMTF_LP = true;
        sMTF_BP = false;
        sMTF_HP = false;
    elseif min(sMTFu(1:find(sMTFu==Max_sMTFu1))) <= .7 && min(sMTFu(find(sMTFu==Max_sMTFu1):end)) > .7
        sMTF_LP = false;
        sMTF_BP = false;
        sMTF_HP = true;
    end
    
elseif length(f) == 1
    Max_sMTFu1 = max(sMTFu(idx07(1):idx07(f))); 
    Max_sMTFu2 = max(sMTFu(idx07(f+1):idx07(end)));
    
    if (sum(sMTFu(idx07(1):idx07(f)))-.7*(idx07(f)-idx07(1)+1)) > (sum(sMTFu(idx07(f+1):idx07(end)))-.7*(idx07(end)-idx07(f+1)+1))
        Max_sMTFu12 = 1;
    else
        Max_sMTFu12 = 2;
    end
    
    if min(sMTFu(1:find(sMTFu==Max_sMTFu1))) > .7 && min(sMTFu(find(sMTFu==Max_sMTFu1):find(sMTFu==Max_sMTFu2)))<=.7 && min(sMTFu(find(sMTFu==Max_sMTFu2):end))<=.7
        sMTF_LP = true;
        sMTF_BP = true;
        sMTF_HP = false;
    elseif min(sMTFu(1:find(sMTFu==Max_sMTFu1))) <= .7 && min(sMTFu(find(sMTFu==Max_sMTFu1):find(sMTFu==Max_sMTFu2)))<=.7 && min(sMTFu(find(sMTFu==Max_sMTFu2):end)) > .7
        sMTF_LP = false;
        sMTF_BP = true;
        sMTF_HP = true;
    elseif min(sMTFu(1:find(sMTFu==Max_sMTFu1))) <= .7 && min(sMTFu(find(sMTFu==Max_sMTFu1):find(sMTFu==Max_sMTFu2)))<=.7 && min(sMTFu(find(sMTFu==Max_sMTFu2):end)) <= .7
        sMTF_LP = false;
        sMTF_BP = true;
        sMTF_HP = false;
    elseif min(sMTFu(1:find(sMTFu==Max_sMTFu1))) > .7 && min(sMTFu(find(sMTFu==Max_sMTFu1):find(sMTFu==Max_sMTFu2)))<=.7 && min(sMTFu(find(sMTFu==Max_sMTFu2):end)) > .7
        sMTF_LP = true;
        sMTF_BP = false;
        sMTF_HP = true;
    end
    
else
    Max_sMTFu1 = max(sMTFu);
    Max_sMTFu2 = NaN;
    Max_sMTFu12 =  NaN;
    
    sMTF_LP = false;
    sMTF_BP = false;
    sMTF_HP = false;
end

% sMTF_BW
if ~isnan(Max_sMTFu1)
    sx = find(sMTFu==Max_sMTFu1);
    while sMTFu(sx) > .7
        if sx == length(sMTFu)
            break
        end
        sx=sx+1;
    end
    s2=sx;
    % if sMTF_BP % commented out lines 126 and 135 to 137 Natsumi 20July2017
    sx = find(sMTFu==Max_sMTFu1);
    while sMTFu(sx) > .7
        if sx == 1
            break
        end
        sx=sx-1;
    end
    s1=sx;
    % else
    %     s1=1;
    % end
    sMTF_BW1=RDu(s2)-RDu(s1);
    sMTF_Max1 = RDu(s2); % added 20Jul2017 Natsumi
    sMTF_Min1 = RDu(s1); % added 20Jul2017 Natsumi
    
    % sBMF
    %sBMF1 = RDu(find(sMTFu==Max_sMTFu1)); % use this for both case % Natsumi 20Jul2017 commented out lines 143 to 147
    if sMTF_LP
        sBMF1 = RDu(round(s2/2));
    elseif tMTF_HP
        sBMF1 = RDu(round(s1+(s2-s1)/2));
    else
        sBMF1 = RDu(find(sMTFu==Max_sMTFu1));
    end
    
else
    sMTF_BW1= NaN;
    sMTF_Max1 = NaN;
    sMTF_Min1 = NaN;
    sBMF1 = NaN;
end

% sMTF_BW
if ~isnan(Max_sMTFu2)
    sx = find(sMTFu==Max_sMTFu2);
    while sMTFu(sx) > .7
        if sx == length(sMTFu)
            break
        end
        sx=sx+1;
    end
    s2=sx;
    % if sMTF_BP % commented out lines 126 and 135 to 137 Natsumi 20July2017
    sx = find(sMTFu==Max_sMTFu2);
    while sMTFu(sx) > .7
        if sx == 1
            break
        end
        sx=sx-1;
    end
    s1=sx;
    % else
    %     s1=1;
    % end
    sMTF_BW2=RDu(s2)-RDu(s1);
    sMTF_Max2 = RDu(s2); % added 20Jul2017 Natsumi
    sMTF_Min2 = RDu(s1); % added 20Jul2017 Natsumi
    
    % sBMF
    %sBMF2 = RDu(find(sMTFu==Max_sMTFu2)); % use this for both case % Natsumi 20Jul2017 commented out lines 143 to 147
    if sMTF_HP
        assert(sMTF_LP == true || sMTF_BP == true)
        sBMF2 = RDu(round(s1+(s2-s1)/2));
    else
        assert(sMTF_BP == true)
        sBMF2 = RDu(find(sMTFu==Max_sMTFu2));
    end
    
else
    sMTF_BW2 = NaN;
    sMTF_Max2 = NaN;
    sMTF_Min2 = NaN;
    sBMF2 = NaN;
end

% calculate from all above 0.7
if ~isempty(f)
    dummy_sMTF_BW = zeros(length(f)+1,1);
    for idxf = 1:length(f)+1
        if idxf == 1
            %         idx07(1)-1
            %         idx07(f(idxf))+1
            try dummy_sMTF_BW(idxf) = RDu( idx07(f(idxf))+1 ) - RDu(idx07(1)-1);
            catch dummy_sMTF_BW(idxf) = RDu( idx07(f(idxf))+1 ) - RDu(idx07(1));
            end
        elseif idxf > 1 & idxf < length(f)+1
            %         idx07(f(idxf-1)+1)-1
            %         idx07(f(idxf))+1
            dummy_sMTF_BW(idxf) = RDu( idx07(f(idxf))+1 ) - RDu( idx07(f(idxf-1)+1)-1 );
        elseif  idxf == length(f)+1
            %         idx07(f(idxf-1)+1)-1
            %         idx07(end)+1
            try dummy_sMTF_BW(idxf) = RDu( idx07(end)+1 ) - RDu(  idx07(f(idxf-1)+1)-1 );
            catch dummy_sMTF_BW(idxf) = RDu( idx07(end) ) - RDu(  idx07(f(idxf-1)+1)-1 );
            end
        end
    end
    
    sMTF_BW= sum(dummy_sMTF_BW);
    try sMTF_Max= RDu( idx07(end)+1 );
    catch sMTF_Max= RDu( idx07(end) );
    end
    try sMTF_Min= RDu( idx07(1)-1 );
    catch sMTF_Min= RDu( idx07(1) );
    end
    sBMF = RDu(find(sMTFu==Max_sMTFu));
    
else
    
    try sMTF_Max= RDu( idx07(end)+1 );
    catch sMTF_Max= RDu( idx07(end) );
    end
    try sMTF_Min= RDu( idx07(1)-1 );
    catch sMTF_Min= RDu( idx07(1) );
    end
    sMTF_BW= sMTF_Max - sMTF_Min;
    sBMF = RDu(find(sMTFu==Max_sMTFu));
    
end

% Params ------------------------------------------------------------------
RTFparam.tMTF = tMTF;
RTFparam.sMTF = sMTF;
RTFparam.Fm = Fm;
RTFparam.RD = RD;
            
RTFparam.tBMF = tBMF;
RTFparam.sBMF = sBMF;
RTFparam.tMTF_BW = tMTF_BW;
RTFparam.sMTF_BW = sMTF_BW;
RTFparam.tMTF_Min = tMTF_Min; % added 2Oct2017 Natsumi
RTFparam.sMTF_Min = sMTF_Min; % added 2Oct2017 Natsumi
RTFparam.tMTF_Max = tMTF_Max; % added 20Jul2017 Natsumi
RTFparam.sMTF_Max = sMTF_Max; % added 20Jul2017 Natsumi
RTFparam.tMTF_BP = tMTF_BP;
RTFparam.sMTF_BP = sMTF_BP;

% new in version02
RTFparam.tMTFu = tMTFu;
RTFparam.sMTFu = sMTFu;
RTFparam.Fmu = Fmu;
RTFparam.RDu = RDu;

RTFparam.tMTF_LP = tMTF_LP;
RTFparam.sMTF_LP = sMTF_LP;
RTFparam.tMTF_HP = tMTF_HP;
RTFparam.sMTF_HP = sMTF_HP;

RTFparam.tBMF1 = tBMF1;
RTFparam.sBMF1 = sBMF1;
RTFparam.tMTF_BW1 = tMTF_BW1;
RTFparam.sMTF_BW1 = sMTF_BW1;
RTFparam.tMTF_Min1 = tMTF_Min1;
RTFparam.sMTF_Min1 = sMTF_Min1;
RTFparam.tMTF_Max1 = tMTF_Max1;
RTFparam.sMTF_Max1 = sMTF_Max1;

RTFparam.tBMF2 = tBMF2;
RTFparam.sBMF2 = sBMF2;
RTFparam.tMTF_BW2 = tMTF_BW2;
RTFparam.sMTF_BW2 = sMTF_BW2;
RTFparam.tMTF_Min2 = tMTF_Min2;
RTFparam.sMTF_Min2 = sMTF_Min2;
RTFparam.tMTF_Max2 = tMTF_Max2;
RTFparam.sMTF_Max2 = sMTF_Max2;

% new in version04
RTFparam.Max_tMTFu12 = Max_tMTFu12;
RTFparam.Max_sMTFu12 = Max_sMTFu12;

if flag_plot
    % plot RTF
    figure;
    figuresetup2savepdf(30, 10)
    RTF4plot = RTF./sum(sum(RTF));
    cmap = cschemes('spectral', 21);
    Max = max(max(RTF4plot));
    subplot(1,3,1)
    imagesc(Fm,RD,RTF4plot,[0 Max]),shading flat,colormap(cmap) %axis square
    hold on
    plot(tBMF,sBMF,'g*','linewidth',2)
    axis xy
    xlabel('TMF (Hz)')
    ylabel('SMF (cyc/oct)')
    % plot tMTFu
    subplot(1,3,2)
    hold on;
    plot(Fmu,tMTFu,'k')
    line([tBMF tBMF],[0 1.2],'Color','k','LineWidth',1,'LineStyle',':');
    plot(tBMF,1,'r+','linewidth',2)
    xlim([-Fm(end) Fm(end)])
    ylim([0 1])
    title(sprintf('TBMF %.2f(Hz)', tBMF))
    xlabel('TMF (Hz)')
     % plot sMTFu
    subplot(1,3,3)
    hold on;
    plot(RDu,sMTFu,'k')
    line([sBMF sBMF],[0 1.2],'Color','k','LineWidth',1,'LineStyle',':');
    plot(sBMF,1,'r+','linewidth',2)
    xlim([0 RD(end)])
    ylim([0 1])
    title(sprintf('SBMF %.2f(cyc/oct)', sBMF))
    xlabel('SMF (cyc/oct)')
end



end