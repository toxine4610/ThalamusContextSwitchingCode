

logistic = @(x, a) 1./(1 + 2*exp(-a.*(x-0.0)));
x  =  linspace(-1,1,200);

for i = 1:200
    a = 15*rand;
    y = logistic(x, a);
    V(i,:) = y;
end


%%

ct=  0;

for i = 1:20:200
    ct= ct+1;
    mu = mean( V(:,i), 1 );
    if ct == 6
        se = std(  V(:,i), [], 1 )+0.2;
    else
        se = std(  V(:,i), [], 1 )*0.8;
    end;

    R(ct,:) = normrnd(mu, se, [1,100] );
end
   
figure(1);
stdshadenew( 1:ct, mean(R',1), std(R',1),0.2,'b');


        