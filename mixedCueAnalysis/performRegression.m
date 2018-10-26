

function [B, Borig,dev, stats] = performRegression( Xnew, Y, ind, normFlag )

if nargin < 3; normFlag = 1; end; % normalize result by default

%% ===== perform multinomial regression....................................

[B, dev, stats] = mnrfit( Xnew, Y,'Interactions','on' );

%% ==== add nans to missing data.
if ~isempty(ind)
    if length(ind)==2
        Borig(1:ind(1)-1) = B(1:ind(1)-1);
        Borig(ind(1))     = NaN;
        Borig( ind(1)+1 :  ind(2)-1 ) = B(ind(1):ind(2)-2);
        Borig(ind(2))     = NaN;
        Borig( ind(2)+1 : size(B,1)+length(ind) ) = B( ind(2)-1:end);
    else
        Borig(1:ind(1)-1) = B(1:ind(1)-1);
        Borig(ind(1))     = NaN;
        Borig( ind(1)+1 : size(B,1)+length(ind) ) = B( ind(1):end);
    end
else
    Borig = B;
end;

if normFlag == 1

    
    
    Borig = abs(Borig)./max(abs(Borig));
end;