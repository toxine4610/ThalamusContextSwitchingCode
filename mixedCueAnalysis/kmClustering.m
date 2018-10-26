rng(1); % For reproducibility


numClust = 4;
cmap = lines(numClust);

x = cell2mat( PFC_CMI_C1_we);
y = cell2mat( PFC_CMI_C2_we );
X = [-x; -y];

% MDC1 = cell2mat( MDFRTS_C1' );
% MDC2 = cell2mat( MDFRTS_C2' );

opts = statset('Display','final');
[idx,C] = kmeans(X',numClust,'Distance','sqeuclidean', 'Replicates',20,'Options',opts);

figure(1);
axis square; box off; set(gca,'tickdir','out','fontsize',16);
hold on;
plot( linspace(-1,1,20), linspace(-1,1,20),'color','k' );
line( [0,0], [-1,1],'color','k');
line( [-1,1],[0,0],'color','k');

for i = 1:numClust
    plot( -x( idx == i), -y(idx == i) ,'o','markersize',12,'markerfacecolor',cmap(i,:),'markeredgecolor','w' ); hold on
end;


for i = 1:numClust
    plot( C(i,1), C(i,2), '+','markersize',22,'markeredgecolor',cmap(i,:),'linewidth',4 );
end


% MDC1 = cell2mat( MDFRTS_C1' );
% MDC2 = cell2mat( MDFRTS_C2' );
% 
% for i = 1:numClust
%     figure(2);
%     subplot(numClust,1, i);
%     plot( time, mean( MDC1(idx==i,:)),'color',cmap(i,:) ); hold on;
%     plot( time, mean( MDC2(idx==i,:)),'color',cmap(i,:),'linestyle','--' ); hold on;
%     xlim([-0.2, 1]);
%     axis square; box off; set(gca,'tickdir','out','fontsize',16);
% end;
