

for i = 1 : 1000
    x1 = normrnd(10,2);
    x2 = normrnd(-3,6);
    if rand < 0.6
        D(i) = x1;
    else
        D(i) = x2;
    end
end

figure(1); set(gcf,'color','w');
histogram( D, 20,'normalization','probability');
hold on;

%%
GMModel = fitgmdist(D',2,'RegularizationValue',0.1)
x = linspace(-20,20,200);
plot(x, pdf(GMModel,x'))


% x = linspace(-20,20,20000);
% 
% pd1 = makedist('Normal',GMModel.mu(1),GMModel.Sigma(:,:,1));
% comp1 =  pdf(pd1, x');
% pd2 = makedist('Normal',GMModel.mu(2),GMModel.Sigma(:,:,2));
% comp2 =  pdf(pd2, x');
% 
% 
% plot(x, comp1); hold on
% plot(x, comp2);
