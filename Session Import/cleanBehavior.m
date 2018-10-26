
function [led,ts] = cleanBehavior(behavior);

i = 1;
ct = 0; ct2 = 0;
while i <= size(behavior,1)
    foo = behavior{i};
    foo2 = str2num(foo);
    if ~isempty(foo2)
        ct = ct+1;
        ts{ct} = foo;
    else
        ct2 = ct2+1;
        led{ct2} = foo;
    end;
    i = i+1;
end;
        