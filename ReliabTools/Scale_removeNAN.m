function scaled = Scale_removeNAN(data);

scaled = (data - nanmean(data)) / nanstd(data);
scaled(isnan(data)) = 0;