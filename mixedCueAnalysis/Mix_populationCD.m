
%%
addpath ../Rajeev_Code
Datapath =  '../ForContextSwitchProject/DoubleCueDatabase/SOMCre/DualModalityCue/';

colorMD = [0.49,0.18,0.56];
colorPFC = [0.5,0.5,0.5];

%%
ExptDates = {'2017-12-15','2017-12-16','2017-12-17','2017-12-18','2017-12-20','2017-12-22','2017-12-23','2017-12-24','2017-12-26'};

for d = 1:length(ExptDates)
    
    clear Spfc Smd Spfc_mix Smd_mix Z_C1 D
    load( [Datapath ExptDates{d} '/Somcre_mixed.mat'] )
    
    %%
    [Spfc_mix, Smd_mix] = packageData(Z_C1, Smd, Spfc);
    
    
    %%
    timevec = [0.2, 0.8];
    
    clear MD_nosound_proj_x1 MD_sound_proj_x1 MD_nosound_proj_x2 MD_sound_proj_x2
    clear PFC_nosound_proj_x1 PFC_sound_proj_x1 PFC_nosound_proj_x2 PFC_sound_proj_x2
    
    [MDCD_sound,MD_sound_Sel, MD_sound_proj_x1, MD_sound_proj_x2,...
        MDCD_nosound,MD_nosound_Sel, MD_nosound_proj_x1, MD_nosound_proj_x2,...
        PFCCD_sound,PFC_sound_Sel, PFC_sound_proj_x1, PFC_sound_proj_x2,...
        PFCCD_nosound,PFC_nosound_Sel, PFC_nosound_proj_x1, PFC_nosound_proj_x2, time] = computeCD(Smd_mix, Spfc_mix, timevec);
    
    sMD_NS(:,d) = smooth( MD_nosound_proj_x1{1}(1:end-1)) -smooth( MD_sound_proj_x1{1}(1:end-1)) ;
    sMD_S(:,d)  = smooth( MD_nosound_proj_x2{1}(1:end-1)) -smooth( MD_sound_proj_x2{1}(1:end-1))  ;
    
    sPFC_NS(:,d) = smooth( PFC_nosound_proj_x1{1}(1:end-1)) -smooth( PFC_sound_proj_x1{1}(1:end-1)) ;
    sPFC_S(:,d)  = smooth( PFC_nosound_proj_x2{1}(1:end-1)) -smooth( PFC_sound_proj_x2{1}(1:end-1));
    
    clear foo
    foo =  smooth( MD_nosound_proj_x1{1}(1:end-1)) -smooth( MD_nosound_proj_x2{1}(1:end-1));
    AudMD_NS(d) = nanmean( foo );
    
    clear foo
    foo =  smooth( PFC_nosound_proj_x1{1}(1:end-1)) -smooth( PFC_nosound_proj_x2{1}(1:end-1));
    AudPFC_NS(d) = nanmean( foo );
    
    clear foo
    foo =  smooth( MD_sound_proj_x1{1}(1:end-1)) -smooth( MD_sound_proj_x2{1}(1:end-1));
    AudMD_S(d) = nanmean( foo );
    
    clear foo
    foo =  smooth( PFC_sound_proj_x1{1}(1:end-1)) -smooth( PFC_sound_proj_x2{1}(1:end-1));
    AudPFC_S(d) = nanmean( foo );
    
    m = min(length(PFC_sound_Sel{1}), length(PFC_nosound_Sel{1}));
%      PFC_theta{d} = atand( PFC_sound_Sel{1}(1:m)./PFC_nosound_Sel{1}(1:m) );
%      PFC_R{d}     = sqrt( PFC_sound_Sel{1}(1:m).^2 + PFC_nosound_Sel{1}(1:m).^2 );
% %     
%      MD_theta{d} = atand( MD_sound_Sel{1}./MD_nosound_Sel{1} );
%      MD_R{d}     = sqrt( MD_sound_Sel{1}.^2 + MD_nosound_Sel{1}.^2 );
% %     
%      figure(2);
% %     subplot(1,2,1);
%      plot( MD_sound_Sel{1},  MD_nosound_Sel{1},'o','markerfacecolor',colorMD,'markeredgecolor',colorMD); hold on;
% %     
% %     subplot(1,2,2);
%     plot( PFC_sound_Sel{1}(1:m), PFC_nosound_Sel{1}(1:m),'o','markerfacecolor',colorPFC,'markeredgecolor',colorPFC); hold on;
%     
%     MDSS{d} =MD_sound_Sel{1};
%     MDNSS{d} =MD_nosound_Sel{1};
%     
%     PFCSS{d} =  PFC_sound_Sel{1};
%     PFCNSS{d} = PFC_nosound_Sel{1};
%     
end;

%%
figure(1);
subplot(1,2,1)
stdshadenew( time(1:end-1), nanmean(sMD_S,2)', (nanstd(sMD_S,[],2)')./sqrt(d),0.5,colorMD);
stdshadenew( time(1:end-1), nanmean(sPFC_S,2)', (nanstd(sPFC_S,[],2)')./sqrt(d),0.5,colorPFC);
xlim([-0.2, 1.2])

figure(1);
subplot(1,2,2)
stdshadenew( time(1:end-1), nanmean(sMD_NS,2)', (nanstd(sMD_NS,[],2)')./sqrt(d),0.5,colorMD);
stdshadenew( time(1:end-1), nanmean(sPFC_NS,2)', (nanstd(sPFC_NS,[],2)')./sqrt(d),0.5,colorPFC);
xlim([-0.2, 1.2])


%%
% figure(5);
% [a,b] = hist( cell2mat(PFC_R),20);
% plot( b, cumsum(a./sum(a)) ); hold on;
% [a,b] = hist( cell2mat(MD_R),20);
% plot( b, cumsum(a./sum(a)) ); hold on;
