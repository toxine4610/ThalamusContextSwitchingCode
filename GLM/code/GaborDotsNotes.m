function [info, exname] = GaborDotsNotes(ex)
% [info, exname] = GaborDotsNotes(ex)

info = [];

if nargin < 1
    ex = 'dummy';
end

exnames = {...
    'p20121104', ...
    'p20121107', ...
    'p20121114', ...
    'p20121204', ...
    'p20121205', ...
    'p20121206', ...
    'p20121212', ...
    'p20130328', ...
    'p20130329', ...
    'p20130502', ...
    'p20130515', ...
    'p20130516', ...
    'p20130517', ...
    'p20130611', ...
    'p20140213', ...
    'p20140303', ...
    'p20140304', ...
    'p20140305', ...
    'p20140306', ...
    'p20140307', ...
    'p20140310', ...
    'n20150304a', ...
    'n20150304b', ...
    'n20150305a', ...
    'n20150305b', ...
    'n20150306b', ...
    'n20150306c', ...
    'n20150309', ...
    'n20150312b', ...
    'n20150316c', ...
    'n20150324a', ...
    'n20150325', ...
    'n20150326a', ...
    'n20150331', ...
    'n20150401', ...
    'n20150407a', ...
    'n20150407b', ...
    'n20150408a', ...
    'n20150408b', ...
    'n20150416', ...
    'n20150501', ...
    'n20150505', ...
    'n20150518', ...
    'n20150519', ...
    'n20150609', ...
    };

if isnumeric(ex)
    if ex <= numel(exnames)
        exname = exnames{ex};
    else
        exname = 'dummy';
    end
else
    exname = ex;
end


%----------------------------------------------------------------------%
% notes about each experiment
switch exname
    case {'p20121031', '20121031'}
        info.exname = 'p20121031';
        info.monkey = 'Pat';
        info.tags = {'LIP'};
        
    case {'p20121104', '20121104'}
        info.exname = 'p20121104';
        info.monkey = 'Pat';
        info.tags = {'LIP', 'good', 'revco', 'behavior', 'simultaneous'};
        
    case {'p20121107', '20121107'}
        info.exname = 'p20121107';
        info.monkey = 'Pat';
        info.tags = {'LIP', 'good', 'revco', 'LFP', 'behavior', 'simultaneous'};
        
    case {'p20121110', '20121110'}
        info.exname = 'p20121110';
        info.monkey = 'Pat';
        info.tags = {'LIP', 'LFP'};
        
    case {'p20121113', '20121113'}
        info.exname = 'p20121113';
        info.monkey = 'Pat';
        info.tags = {'LIP', 'good', 'revco', 'LFP', 'behavior'};
        
    case {'p20121114', '20121114'}
        info.exname = 'p20121114';
        info.monkey = 'Pat';
        info.tags = {'LIP', 'good', 'synchrony peak', 'LFP', 'behavior', 'simultaneous'};
        
    case {'p20121129', '20121129'}
        info.exname = 'p20121129';
        info.monkey = 'Pat';
        info.tags = {'LIP', 'good', 'LFP'};
        
    case {'p20121130', '20121130'}
        info.exname = 'p20121130';
        info.monkey = 'Pat';
        info.tags = {'LIP'};
        
    case {'p20121204', '20121204'}
        info.exname = 'p20121204';
        info.monkey = 'Pat';
        info.tags = {'LIP', 'good', 'revco', 'LFP', 'behavior', 'simultaneous'};
        
    case {'p20121205', '20121205'}
        info.exname = 'p20121205';
        info.monkey = 'Pat';
        info.tags = {'LIP', 'good', 'revco', 'LFP', 'synchrony peak','behavior', 'simultaneous'};
        
    case {'p20121206', '20121206'}
        info.exname = 'p20121206';
        info.monkey = 'Pat';
        info.tags = {'LIP', 'good', 'revco', 'LFP', 'behavior', 'simultaneous'};
        
    case {'p20121212', '20121212'}
        info.exname = 'p20121212';
        info.monkey = 'Pat';
        info.tags = {'LIP', 'good', 'revco', 'LFP', 'behavior', 'simultaneous'};
        
    case {'p20130326', '20130326'}
        info.exname = 'p20130326';
        info.monkey = 'Pat';
        info.tags = {'LIP', 'LFP', 'synchrony peak'};
        
    case {'p20130328', '20130328'}
        info.exname = 'p20130328';
        info.monkey = 'Pat';
        info.tags = {'LIP', 'good', 'revco', 'LFP', 'behavior', 'simultaneous'};
        
    case {'p20130329', '20130329'}
        info.exname = 'p20130329';
        info.monkey = 'Pat';
        info.tags = {'LIP', 'good', 'revco', 'LFP', 'behavior', 'simultaneous'};
        
    case {'p20130502', '20130502'}
        info.exname = 'p20130502';
        info.monkey = 'Pat';
        info.tags = {'MT', 'LIP', 'good', 'revco', 'LFP', 'simultaneous'};
        
    case {'p20130514', '20130514'}
        info.exname = 'p20130514';
        info.monkey = 'Pat';
        info.tags = {'LIP', 'MT', 'LFP'};
        
    case {'p20130515', '20130515'}
        info.exname = 'p20130515';
        info.monkey = 'Pat';
        info.tags = {'LIP', 'good', 'LFP', 'behavior', 'simultaneous'};
        
    case {'p20130516', '20130516'}
        info.exname = 'p20130516';
        info.monkey = 'Pat';
        info.tags = {'LIP', 'good', 'LFP', 'behavior', 'simultaneous'};
        
    case {'p20130517', '20130517'}
        info.exname = 'p20130517';
        info.monkey = 'Pat';
        info.tags = {'MT', 'LIP', 'good', 'revco', 'LFP', 'behavior', 'simultaneous'};
        
    case {'p20130611', '20130611'}
        info.exname = 'p20130611';
        info.monkey = 'Pat';
        info.tags = {'MT', 'LIP', 'good', 'LFP', 'behavior', 'simultaneous'};
        
    case {'p20140213', '20140213'}
        info.exname = 'p20140213';
        info.monkey = 'Pat';
        info.tags = {'MT', 'good', 'revco', 'LFP', 'behavior'};
        
    case {'p20140218', '20140218'}
        info.exname = 'p20140218';
        info.monkey = 'Pat';
        info.tags = {'MT', 'good', 'revco', 'LFP'};
        
    case {'p20140226', '20140226'}
        info.exname = 'p20140226';
        info.monkey = 'Pat';
        info.tags = {'MT', 'LIP', 'LFP', 'good'};
        
    case {'p20140303', '20140303'}
        info.exname = 'p20140303';
        info.monkey = 'Pat';
        info.tags = {'MT', 'LIP', 'good', 'revco', 'LFP', 'behavior', 'simultaneous'};
        
    case {'p20140304', '20140304'}
        info.exname = 'p20140304';
        info.monkey = 'Pat';
        info.tags = {'MT', 'good', 'revco', 'LFP', 'behavior', 'simultaneous'};
        
    case {'p20140305', '20140305'}
        info.exname = 'p20140305';
        info.monkey = 'Pat';
        info.tags = {'MT', 'good', 'revco', 'LFP', 'behavior', 'simultaneous'};
        
    case 'p20140306'
        info.exname = exname;
        info.monkey = 'Pat';
        info.tags = {'MT', 'LIP', 'good', 'revco', 'LFP', 'simultaneous'};
        
    case {'p20140307', '20140307', 'p20140307b'}
        info.exname = exname;
        info.monkey = 'Pat';
        info.tags = {'MT', 'LIP', 'good', 'revco', 'LFP', 'behavior', 'simultaneous'};
        
    case {'p20140310', '20140310'}
        info.exname = exname;
        info.monkey = 'Pat';
        info.tags = {'MT', 'LIP', 'good', 'revco', 'LFP', 'behavior', 'simultaneous'};
        
        
    case {'n20150304a', '20150304a'}
        exname = 'n20150304a';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT', 'LIP', 'good', 'behavior', 'simultaneous'};
        
        
    case {'n20150304b', '20150304b'}
        exname = 'n20150304b';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT', 'LIP', 'good', 'revco', 'behavior', 'simultaneous'};
        
    case {'n20150305a', '20150305a'}
        exname = 'n20150305a';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT', 'LIP', 'good', 'behavior', 'simultaneous'};
        
    case {'n20150305b', '20150305b'}
        exname = 'n20150305b';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT', 'LIP', 'good', 'behavior', 'simultaneous'};
        
    case {'n20150305c', '20150305c'}
        exname = 'n20150305c';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT', 'LIP', 'FreeChoice'};
        
    case {'n20150306b', '20150306b'}
        exname = 'n20150306b';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT', 'good', 'behavior', 'simultaneous', 'LFP'}; %  'behavior'
        
    case {'n20150306c', '20150306c'}
        exname = 'n20150306c';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT', 'LIP', 'good', 'revco', 'simultaneous'}; % 'behavior'
        
    case {'n20150309', '20150309'}
        exname = 'n20150309';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT', 'LIP', 'good', 'behavior', 'simultaneous'};
        
    case {'n20150310', '20150310'}
        exname = 'n20150310';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT','behavior'};
        
    case {'n20150312a', '20150312a'}
        exname = 'n20150312a';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'LFP', 'SaccadeMappingOnly'};
        
        
    case {'n20150312b', '20150312b'}
        exname = 'n20150312b';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'LIP', 'good', 'behavior', 'simultaneous'};
        
    case {'n20150313', '20150313'}
        exname = 'n20150313';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'LIP'};
        
    case {'n20150316c', '20150316c'}
        exname = 'n20150316c';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT','good', 'revco', 'behavior', 'simultaneous'};
        
    case {'n20150324a', '20150324a'}
        exname = 'n20150324a';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT', 'good', 'revco', 'behavior', 'simultaneous'};
        
    case {'n20150325', '20150325'}
        exname = 'n20150325';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'LIP', 'good', 'behavior', 'simultaneous'};
        
    case {'n20150326a', '20150326a'}
        exname = 'n20150326a';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'LIP', 'good', 'revco', 'behavior', 'simultaneous'};
        
    case {'n20150326b', '20150326b'}
        exname = 'n20150326b';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'LIP'};
        
    case {'n20150331', '20150331'}
        exname = 'n20150331';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT', 'LIP', 'good', 'revco', 'behavior', 'simultaneous'};
        
    case {'n20150401', '20150401'}
        exname = 'n20150401';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'LIP', 'good', 'revco', 'behavior', 'simultaneous'};
        
    case {'n20150402', '20150402'}
        exname = 'n20150402';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT', 'LIP', 'good'};
        
    case {'n20150407a', '20150407a'}
        exname = 'n20150407a';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'LIP', 'good', 'revco', 'behavior', 'simultaneous'};
        
    case {'n20150407b', '20150407b'}
        exname = 'n20150407b';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'LIP', 'good', 'revco', 'behavior'};
        
    case {'n20150408a', '20150408a'}
        exname = 'n20150408a';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT', 'LIP', 'good', 'revco', 'behavior', 'simultaneous'};
        
    case {'n20150408b', '20150408b'}
        exname = 'n20150408b';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'LIP', 'good', 'revco', 'behavior', 'simultaneous'};
        
    case {'n20150416', '20150416'}
        exname = 'n20150416';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'LIP', 'revco', 'good', 'behavior', 'simultaneous'};
        
    case {'n20150429', '20150429'}
        exname = 'n20150429';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'LIP'};
        
    case {'n20150501', '20150501'}
        exname = 'n20150501';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'LIP', 'revco', 'good', 'behavior', 'simultaneous'};
        
    case {'n20150505', '20150505'}
        exname = 'n20150505';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'LIP', 'revco', 'good', 'behavior', 'simultaneous'};
        
    case {'n20150508', '20150508'}
        exname = 'n20150508';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'LIP', 'revco', 'behavior'};
        
        
    case {'n20150518', '20150518'}
        exname = 'n20150518';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT', 'LIP', 'revco', 'good', 'behavior', 'simultaneous'};
        
    case {'n20150519', '20150519'}
        exname = 'n20150519';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT', 'LIP', 'revco', 'good', 'behavior', 'simultaneous'};
        
    case {'n20150609', '20150609'}
        exname = 'n20150609';
        info.exname = exname;
        info.monkey = 'Nancy';
        info.tags = {'MT', 'LIP', 'revco', 'good', 'behavior', 'simultaneous'};
        
        
    otherwise
        if nargout == 0
            fprintf('%d experiments to choose from:\r', numel(exnames))
            fprintf('%s\r', exnames{:})
        end
        exname = exnames;
end
