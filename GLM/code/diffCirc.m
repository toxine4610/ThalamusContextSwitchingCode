function deltaStim = diffCirc(stimDir, prefDir)
% get the distance between two circular variables (in degrees)
% deltaStim = diffCirc(stimDir, prefDir)

dist=@(x,y) angle(exp(1i*x)./exp(1i*y));
deltaStim = (dist(stimDir/180*pi, prefDir/180*pi))*180/pi;