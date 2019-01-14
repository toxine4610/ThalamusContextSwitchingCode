function [contrast, phi, dphi] = getSingleGabor(maskphase, basephase)
% convert from two gabor implementation to one gabor implementation
% [contrast, phi] = getSingleGabor(maskphase, basephase)
% 
% Inputs:
%   <maskphase>  [nTrials x nFrames x nGabors]
%   <basephase>  [nTrials x nFrames x nGabors]
%
% Outputs:
%   <contrast>  [nTrials x nFrames x nGabors] contrast
%   <phase>     [nTrials x nFrames x nGabors] phase
%   <dphase>    [nTrials x nFrames x nGabors] time derivative of phase

contrast=2*cosd( (maskphase-basephase)/2);
phi=(maskphase+basephase)/2;
if nargout>2
    dphi=zeros(size(phi));
    if numel(size(maskphase))==2
        dphi(1:end-1,:)=diff(phi, [], 1);
    elseif numel(size(maskphase))==3
        dphi(:,1:end-1,:)=diff(phi, [], 2);
    end
    dphi(abs(dphi)>45)=0;
end

contrast(maskphase==0 & basephase==0)=0;