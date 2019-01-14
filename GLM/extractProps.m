


function properties = extractProps(filter)

filter_crop = filter(1:15);

properties.PTR = max(exp(filter_crop))-min(exp(filter_crop));
properties.InhibSubfield = trapz( filter_crop( filter_crop < 0) );
properties.ExcSubfield = trapz( filter_crop( filter_crop > 0) );
properties.Ratio = ( properties.ExcSubfield-abs(properties.InhibSubfield))./ ( properties.ExcSubfield + abs(properties.InhibSubfield));
properties.Diff  = properties.ExcSubfield - abs(properties.InhibSubfield);

[~,tmax] = max(filter_crop);
[~,tmin] = min(filter_crop);
properties.PTT      = tmin-tmax;

