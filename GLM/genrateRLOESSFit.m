
function  [xx,yy2] = genrateRLOESSFit(xax, yax, smthFactor)

yy2 = smooth(xax,yax,smthFactor,'rloess');
[xx,ind] = sort(xax);
yy2 = yy2(ind);
