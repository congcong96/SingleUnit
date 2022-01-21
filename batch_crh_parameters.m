function crh = batch_crh_parameters(crh, flag_plot, basename)
if nargin == 2
    flag_plot = 0;
end
for ii = 1:length(crh)
    mtfhist = crh(ii).mtfhist;
    tmfaxis = crh(ii).tmfaxis;
    smfaxis = crh(ii).smfaxis;
    
    [tBMF, sBMF, tMTF, sMTF, tmfaxisu, smfaxisu, tMTFu, sMTFu] = crh_parameters(mtfhist, tmfaxis, smfaxis, flag_plot);
    crh(ii).tBMF = tBMF;
    crh(ii).sBMF = sBMF;
    crh(ii).tMTF = tMTF;
    crh(ii).sMTF = sMTF;
    crh(ii).tmfaxisu = tmfaxisu;
    crh(ii).smfaxisu = smfaxisu;
    crh(ii).tMTFu = tMTFu;
    crh(ii).sMTFu = sMTFu;

    if flag_plot
       saveas(gcf, sprintf('%s-unit%d.jpg', basename, crh(ii).unit))
       close 
   end
end
