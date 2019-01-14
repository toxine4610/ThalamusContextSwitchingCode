function xfilt=upsampleFilters(kFilt, upFactor)
% xfilt=upsampleFilters(kFilt, upFactor)
% kFilt [nBins x nFilters]
% upFactor [1 x 1] upscale factor

s=size(kFilt);

xfilt=nan(s(1)*10,s(2)); 
for k=1:s(2)
    xfilt(1:10:end,k)=kFilt(:,k); 
    xfilt(:,k)=repnan(xfilt(:,k),'spline')/upFactor;
end