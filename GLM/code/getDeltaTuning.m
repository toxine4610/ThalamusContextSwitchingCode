function deltaStim = getDeltaTuning(stimDir, prefDir)
% get tuning distance for circular variable (in degrees)
% deltaStim = getDeltaTuning(stimDir, prefDir)
% Inputs
%   stimDir - angle (in degrees) of the stimulus
%   prefDir - angle (in degrees) of the preferred direction
% Outputs
%   deltaStim - distance (in degrees) between the two preferred angles

deltaStim = abs(diffCirc(stimDir, prefDir));