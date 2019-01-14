
clear all;

cmiMD = [];
cmiPFC = [];


Dates = {'2017-12-15', '2017-12-16', '2017-12-17', '2017-12-18' , '2017-12-20', '2017-12-22', '2017-12-23' , ...
    '2017-12-24' , '2017-12-26'};

for d = 1:length(Dates)
    
    dataName = [ Dates{d} filesep 'Somcre_mixed.mat'];
    fprintf('Packaging Data....');
    [glmtrial, unitOfTime, binSize, nTrials, binfun, numMD, numPFC, RejMD, RejPFC, keptPFC, keptMD, CMIPFC, CMIMD, Z] = packageData_for_glm(dataName);
    fprintf('...Complete!\n');
     
    cmiMD = [cmiMD, CMIMD(keptMD)];
    cmiPFC = [cmiPFC, CMIPFC(keptPFC)];
    
    md{d} = CMIMD(keptMD);
    pfc{d} = CMIPFC(keptPFC);
end;



%%

MDComp1 = [];
MDComp2 = [];
PFCComp1 = [];
PFCComp2 = [];
Exp1  = [];
PFC_EV1 = [];
PFC_EV2 = [];

MD_EV1 = [];
MD_EV2 = [];
LL  = [];


for i = 1:9
    clear A
    A = load( ['MX_Sess_' num2str(i) '_C1.mat'] );
    
    [V,MDInput1, MDInput2, PFCInput1, PFCInput2, SNRPFC, SNRPFC, MD_expl1, MD_expl2, PFC_expl1, PFC_expl2] = classifyInputs(A.Lc, A.CouplingFilter, A.CouplingFilterPFC);
    
    LL       = [LL, A.Lc];
        
    MDComp1 = cat(1, MDComp1, MDInput1);
    MDComp2 = cat(1, MDComp2, MDInput2);
    
    PFCComp1 = cat(1, PFCComp1, PFCInput1);
    PFCComp2 = cat(1, PFCComp2, PFCInput2);
    
    PFC_EV1 = [PFC_EV1, PFC_expl1];
    PFC_EV2 = [PFC_EV2, PFC_expl2];
    
    MD_EV1 = [MD_EV1, MD_expl1];
    MD_EV2 = [MD_EV2, MD_expl2];
    
    
end

%%

for i = 1
    clear dist A X
    
    this_md = md{i};
    this_pfc = pfc{i};
    
    for j = 1:length(this_pfc)
        for q = 1:length(this_md)
            if sign( this_pfc(j) ) == sign(this_md(q))
                % same context;
                dist(j,q) = 1.*sqrt( this_pfc(j).^2  + this_md(q) ).^2;
            elseif sign( this_pfc(j) ) ~= sign(this_md(q))
                dist(j,q) = -1.*sqrt( this_pfc(j).^2 + this_md(q) ).^2;
            end;
        end;
    end;
    
    A = load( ['MX_Sess_' num2str(i) '_C1.mat'] );
    X = cell2mat( A.AUC' );
    
    C1.meanMDinput{i} = nanmedian( X,2 );
    C1.logL{i} = A.Lc./1e4;
    
    A = load( ['MX_Sess_' num2str(i) '_C2.mat'] );
    X = cell2mat( A.AUC' );
    C2.meanMDinput{i} = nanmedian( X,2 );
    C2.logL{i} = A.Lc./1e4;
    
    muDist{i} = nanmean( dist, 2);
    
   fprintf('... Done\n');
   
end;
    

%%
ll  = cell2mat( C1.logL );
muMD = cell2mat( C1.meanMDinput' );
ind = find( ll > 0.2 );

mudist = cell2mat( muDist );

ll2  = cell2mat( C2.logL );
muMD2 = cell2mat( C2.meanMDinput' );
ind2 = find( ll2 > 0.2 );

indkeep = intersect( ind, ind2 );

muMD = muMD(indkeep );
muMD2 = muMD2(indkeep );

cmipfc = cmiPFC( indkeep );

figure(300); set(gcf,'color','w');
xax = cmipfc; yax = muMD-muMD2;
yy2 = smooth(xax,yax,.5,'rloess');
[xx,ind] = sort(xax);
plot( xax, yax, 'o','markerfacecolor','r', 'markeredgecolor','y'); hold on;
plot( xx, yy2(ind),'k','linewidth',3 );
axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')

figure(400); set(gcf,'color','w');
xax = mudist(indkeep); yax = muMD-muMD2;
yy2 = smooth(xax,yax,.5,'rloess');
[xx,ind] = sort(xax);
plot( xax, yax, 'o','markerfacecolor','r', 'markeredgecolor','y'); hold on;
plot( xx, yy2(ind),'k','linewidth',3 );
axis square; box off; set(gca,'tickdir','out','fontsize',16);set(gca,'xscale','lin')



%%
ind = find( LL > prctile(LL, 55) );
PFCComp1_filt = PFCComp1( ind, :);
PFCComp2_filt = PFCComp2( ind, :);

clear PTR PTT X

ct = 0;
for i = 1:length(ind)
    ct = ct+1;
    foo = PFCComp1(ind(i),2:20);
    PTR(ct) = max(foo)-min(foo);
    [~,tmax] = max(foo);
    [~,tmin] = min(foo);
    PTT(ct)  = tmin-tmax;
    
    E(ct) = trapz( foo(foo>1) );
    I(ct) = trapz( foo(foo<1) );
    CoupStre(ct) = trapz( foo );
    
    EIR(ct) = E(ct)./I(ct);
end;

X = [PTR; PTT];

rng(1);
opts = statset('Display','final');
[idx,C] = kmeans(X',2,'Distance','correlation','Replicates',25,'Options',opts);

u = unique(idx);

GP1 = idx == u(1);
GP2 = idx == u(2);


figure(2);
plot( PTT(GP1), PTR(GP1), 'o','markerfacecolor','r', 'markeredgecolor','y'); hold on;
plot( PTT(GP2), PTR(GP2), 'o','markerfacecolor','b', 'markeredgecolor','y'); hold on;


figure(3);
xax = cmiPFC(GP1); yax = PTR(GP1);
yy2 = smooth(xax,yax,.8,'rloess');
[xx,ind] = sort(xax);
plot( cmiPFC(GP1), PTR(GP1), 'o','markerfacecolor','r', 'markeredgecolor','y'); hold on;
plot( xx, yy2(ind),'k' );


%%
figure(1);
subplot(1,2,1);
plot( PFCComp1_filt( GP1,2:20)' ); hold on;
plot( median( PFCComp1_filt( GP1,2:20), 1), 'k','linewidth', 4);

subplot(1,2,2);
plot( PFCComp1_filt( GP2,2:20)' ); hold on;
plot( median( PFCComp1_filt( GP2,2:20), 1), 'k','linewidth', 4);
