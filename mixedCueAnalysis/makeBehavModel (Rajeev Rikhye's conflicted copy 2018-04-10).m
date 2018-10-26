% clear all;
% VC1 = []; VC2= [];
Vmega = [];
Datapath =  '../ForContextSwitchProject/DoubleCueDatabase/BPM1_3/ContextSwitchHalo/CuePeriodMD/';
% ExptDates = {'2018-02-11','2018-02-12','2018-02-13', '2018-02-14', '2018-02-16'};
%  ExptDates = {'2018-02-17','2018-02-18','2018-02-20','2018-02-22','2018-02-23','2018-02-24'};
%  ExptDates = {'2018-03-13','2018-03-14', '2018-03-15', '2018-03-16'};
%  ExptDates = {'2018-03-17','2018-03-18'};

% ExptDates = {'2018-03-18','2018-03-19','2018-03-20','2018-03-21','2018-03-22','2018-03-24'};
% ExptDates = {'2018-03-13','2018-03-14','2018-03-15','2018-03-16','2018-03-22'};

ExptDates = {'2018-04-07', '2018-04-08'};

for d = 1:length(ExptDates)
    
    clear Spfc Smd Spfc_mix Smd_mix Z_C1 D
    load( [Datapath ExptDates{d} '/BPM25_uniMDHalo_Session.mat'] )
    %     load( [Datapath ExptDates{d} '/Somcre_mixed_ssfo_C1.mat'] )
    Vmega = cat(1, Vmega, Z_C1 );
    SE(d).Behav = Z_C1;
    
end;
%%
% x   = Vmega(:,9);
% xx  = diff(x);
%
% indSwitch = find( xx ~= 0 );
%
% Vm = [];
%
% for i = 1:length(indSwitch)
%     range = [ indSwitch(i) - 15 : indSwitch(i) + 15 ];
%     Vm    = cat(1, Vm, Vmega( range, :) );
% end;

clearvars -except Vmega SE Vm MD B stats
%%
indC1 = find( Vmega(:,9) == 0 );
Vmega = Vmega(indC1,:);

% for  history = 1:1:15

history = 10;

numCVSplits  = 5;
Indices = crossvalind( 'Kfold', size(Vmega,1) , numCVSplits) ;


% clear B dev stats Y X
[currRule, currRew, currChoice, currContext, currLaser, ....
    R, RW, N, NW, L, C ] = getRegressors(Vmega, history);

%%
Y = categorical(currChoice);
X = [currRule', currRew', (currLaser.*currRule)', R, N, RW, NW, L.*R, L.*N, L.*RW, L.*NW] ;

% X = [currRule', currRew', R, N, RW, NW] ;

% clean up the design matrix by removing constant columns.
v  = var(X, [], 1);
ind = find(v==0);
Xnew = X( :, setdiff(1:size(X,2), ind));

% perform multinomial regression..........................................
[B, dev, stats] = mnrfit( Xnew, Y,'Interactions','on' );

[pihat,dlow,dhi] = mnrval(B,Xnew,stats);

if ~isempty(ind)
    
    if length(ind)==2
        Borig(1:ind(1)-1) = B(1:ind(1)-1);
        Borig(ind(1))     = NaN;
        Borig( ind(1)+1 :  ind(2)-1 ) = B(ind(1):ind(2)-2);
        Borig(ind(2))     = NaN;
        Borig( ind(2)+1 : size(B,1)+length(ind) ) = B( ind(2)-1:end);
    else
        Borig(1:ind(1)-1) = B(1:ind(1)-1);
        Borig(ind(1))     = NaN;
        Borig( ind(1)+1 : size(B,1)+length(ind) ) = B( ind(1):end);
    end
    Borig = abs(Borig)./max(abs(Borig));
    
else
    Borig = B;
    Borig = abs(Borig)./max(abs(Borig));
end
%%
figure(3);

yyaxis right
plot( 1:length(currChoice), smooth( pihat(:,2), 2 ) ); hold on;
yyaxis left
plot( 1:length(currChoice), currChoice ); ylim([-2,2]); hold on

%%
Beta = reshape( Borig(5:end), [history, length(Borig(5:end))/history] );
p = stats.p;

figure(4);
plot(1:history, Beta(:,1),'-bo'); hold on
plot(1:history, Beta(:,5),'-ro'); hold on

plot(1:history, Beta(:,2),'-bo'); hold on
plot(1:history, Beta(:,6),'-ro'); hold on
axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')


figure(5);
plot(1:history, Beta(:,3),'-bo'); hold on
plot(1:history, Beta(:,7),'-ro'); hold on

plot(1:history, Beta(:,4),'-bo'); hold on
plot(1:history, Beta(:,8),'-ro'); hold on
axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')


%%

% MFI = ( sum(Beta(:,1))+sum(Beta(:,2)) ) - ( sum(Beta(:,3))+sum(Beta(:,4)) );
% PI  = ( sum(Beta(:,1))-sum(Beta(:,2)) ) + ( sum(Beta(:,3))-sum(Beta(:,4)) );
%
% MFILas = ( sum(Beta(:,5))+sum(Beta(:,6)) ) - ( sum(Beta(:,7))+sum(Beta(:,8)) );
% PILas  = ( sum(Beta(:,5))-sum(Beta(:,6)) ) + ( sum(Beta(:,7))-sum(Beta(:,8)) );


MFI = ( nansum(Beta(:,1)) -  nansum(Beta(:,3))  ) ;

MFILas = ( nansum(Beta(:,5)) - nansum(Beta(:,6)) ); %  - ( sum(Beta(:,7)) + sum(Beta(:,8)) );

figure(8);
bar(3, MFI ); hold on
bar(4, MFILas ) ;

MD(2) = MFI;


%%

VC1 = [];VC2 =[];

for i = 1:size(SE,2)
    thisSe = SE(i).Behav;
    indnolas = find(thisSe(:,10)==0 & thisSe(:,9)==0);
    indlas = find(thisSe(:,10)==1 & thisSe(:,9)==0);
    Perfnolas = thisSe(indnolas,1);
    Perflas   = thisSe(indlas,1);
    
    EFnolas(i) = length(find(Perfnolas)==0)./length(Perfnolas);
    EFlas(i)   = length(find(Perflas)==0)./length(Perflas);

end;
V = [EFnolas; EFlas ];
figure(200); subplot(1,2,1);
ylim([0.2, 1]);
hold on
plotSpread(V');
axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')
ranksum(EFnolas, EFlas)

VC1 = cat(2, VC1, V);

for i = 1:size(SE,2)
    thisSe = SE(i).Behav;
    indnolas = find(thisSe(:,10)==0 & thisSe(:,9)==3);
    indlas = find(thisSe(:,10)==1 & thisSe(:,9)==3);
    Perfnolas = thisSe(indnolas,1);
    Perflas   = thisSe(indlas,1);
    
    EFnolas(i) = length(find(Perfnolas)==0)./length(Perfnolas);
    EFlas(i)   = length(find(Perflas)==0)./length(Perflas);

end;
V = [EFnolas; EFlas ];
figure(200); subplot(1,2,2)
ylim([0.2, 1]);
plotSpread(V');
axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')

VC2 = cat(2, VC2, V);

%%
figure(200); subplot(1,2,1);
ylim([0.2, 1]);
hold on
plotSpread(VC2');
axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')
ranksum(EFnolas, EFlas)
