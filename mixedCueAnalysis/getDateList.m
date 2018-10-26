

function Dates = getDateList(folder)

X = dir(folder);

for i = 3:size(X,1)
    Dates{i-2} = X(i).name;
end;