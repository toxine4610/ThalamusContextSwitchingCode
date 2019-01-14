function h=errorbar2D(x,y,xs,ys, varargin)
% plot 2d errorbar (impoverished options)
% h = errorbar2D(x,y,xerror, yerror, varargin)

% make everything a column vector
x=x(:);
y=y(:);
xs=xs(:);
ys=ys(:);

plot([x-xs x+xs]', [y y]', 'Color', .5*[1 1 1])
hold on
plot([x x]', [y-ys y+ys]', 'Color', .5*[1 1 1])

if isempty(varargin)
    h=plot(x,y,'.'); 
else
    h=plot(x,y,varargin{:});
end