function parsave(filename, model, vstr)
    % cannot save within a parfor

    try
	if nargin < 3
	    save(filename, 'model');
	else
	    save(filename, 'model', vstr);
	end
    catch
	if nargin < 3
	    save(filename, 'model', '-v7.3');
	else
	    save(filename, 'model', '-v7.3', vstr);
	end
    end
end
