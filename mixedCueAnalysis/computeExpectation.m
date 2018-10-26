
function [ExpH, rup, rlo ] = computeExpectation(Z_context1)

Z = 1.96;
useWindow = 1;
if useWindow
    windowlen = 20;
    
    for i = 1:size(Z_context1,1)
        if i == 1 && Z_context1(i,1) == 1
            p = 1; n = 1;
            numH = 1; numT = 0;
        elseif i == 1 && Z_context1(i,1) == 0
            p = 0; n = 0;
            numH = 0; numT = 1;
        elseif i > 1 && i <= windowlen
            foo = Z_context1( 1:i, 1);
            numH = length(find(foo==1));
            numT = length(find(foo==0));
            p = numH./(numH+numT);
            n = length(foo);
        elseif i > windowlen
            foo = Z_context1( i-(windowlen-1) : i, 1);
            numH = length(find(foo==1));
            numT = length(find(foo==0));
            p = numH./(numH+numT);
            n = length(foo);
        end;
        
        rup (i) = p + Z./sqrt(4*n);
        rlo(i)  = p - Z./sqrt(4*n);
        PP(i)   = p;
        
        ExpH(i) = beta( numH+2, numT+1)./beta( numH+1, numT+1 );
        
    end;
    
else
   
    for i = 1:size(Z_context1,1)
        if i == 1 && Z_context1(i,1) == 1
            p = 1; n = 1;
            numH = 1; numT = 0;
        elseif i == 1 && Z_context1(i,1) == 0
            p = 0; n = 0;
            numH = 0; numT = 1;
        elseif i > 1 
            foo = Z_context1( 1:i, 1);
            numH = length(find(foo==1));
            numT = length(find(foo==0));
            p = numH./(numH+numT);
            n = length(foo);
        end;
        
        rup (i) = p + Z./sqrt(4*n);
        rlo(i)  = p - Z./sqrt(4*n);
        PP(i)   = p;
        
        ExpH(i) = beta( numH+2, numT+1)./beta( numH+1, numT+1 );
        
    end;
    
end;