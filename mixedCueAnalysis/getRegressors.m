
function [currSuccess, currRule, currRew, currChoice, currContext, currLaser, ....
            R, RW, N, NW, L, C ] = getRegressors(Vmega, history)

% Vmega = Z_C1;

% ind = find(Vmega(:,9)== 0);
% Vmega = Vmega(ind,:);


% history = 10;

ct = 0;
% N = zeros(1,length(Vmega)-1);
% R = zeros(1,length(Vmega)-1);

for i = (history+1):size(Vmega,1)
    
    ct = ct+1;
    currSuccess(ct) = 2.*Vmega(i,1) - 1;
    currChoice(ct)  = 2.*Vmega(i, 5) - 1; % 1 right, 0 left
    currRule(ct)    = Vmega(i,2);
    currContext(ct) = (Vmega(i,9)-3)+1;
    currLaser(ct)  = Vmega(i,10);
    currRew(ct) = Vmega(i,1);
    
    ct2 = 0;
    for j = i-1:-1:i-history
        ct2 = ct2+1;
        prevChoice = 2.*Vmega(j,5) - 1;
        prevRule   = abs( Vmega(j,2) - currRule(ct) );
        prevRew    = Vmega(j, 1);
        prevLaser  = Vmega(j, 10);
        prevContext  = Vmega(j, 9);
        
        if prevChoice == 1 && prevRule == 0 && prevRew == 1
            R(ct, ct2) = 1;
        elseif  prevChoice == -1 && prevRule == 0 && prevRew == 1
            R(ct, ct2) = -1;
        elseif prevRule == 0 && prevRew == 0
            R(ct, ct2) = 0;
        else
             R(ct,ct2) = 0;
        end;
        
        if prevChoice == 1 && prevRule == 0 && prevRew == 0
            RW(ct, ct2) = 1;
        elseif  prevChoice == -1 && prevRule == 0 && prevRew == 0
            RW(ct, ct2) = -1;
        elseif prevRule == 0 && prevRew == 1
            RW(ct, ct2) = 0;
        else
             RW(ct,ct2) = 0;
        end;
        
        
        if prevChoice == 1 && prevRule == 1 && prevRew == 1
            N(ct,ct2) = 1;
        elseif  prevChoice == -1 && prevRule == 1 && prevRew == 1
            N(ct,ct2) = -1;
        elseif prevRule == 1 && prevRew == 0
            N(ct,ct2) = 0;
        else
             N(ct,ct2) = 0;
        end;
        
        if prevChoice == 1 && prevRule == 1 && prevRew == 0
            NW(ct,ct2) = 1;
        elseif  prevChoice == -1 && prevRule == 1 && prevRew == 0
            NW(ct,ct2) = -1;
        elseif prevRule == 1 && prevRew == 1
            NW(ct,ct2) = 0;
        else
             NW(ct,ct2) = 0;
        end;
        
        if prevLaser == 1
            L(ct,ct2) = 1;
        else
            L(ct,ct2) = 0;
        end;
        
        if prevContext == 3
            C(ct,ct2) = 1;
        else
            C(ct,ct2) = 0;
        end;
        
    end
end
