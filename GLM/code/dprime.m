function dp=dprime(x,y)
% dp=dprime(x,y)

classes=unique(y);
assert(numel(classes)==2, 'Two classes requried to compute d''')

y=y==max(classes);
dp=(mean(x(y,:))-mean(x(~y,:)))./sqrt((var(x(y,:))+var(x(~y,:)))/2);