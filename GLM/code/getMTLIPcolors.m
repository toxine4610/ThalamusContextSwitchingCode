function cmap=getMTLIPcolors(tag)
% cmap=getMTLIPcolors()

if nargin<1
    tag='Nature';
end

switch tag
    case 'jnk'
        clrMT=[70 144 205]/255;
        clrLIP=[85, 184, 74]/255;
        cmap=lines;
        cmap(1,:)=clrMT;
        cmap(2,:)=clrLIP;
    case 'Nature'
        clrMT=[36 165 78]/255;
        clrLIP=[81 128 193]/255;
        cmap=lines;
        cmap(1,:)=clrMT;
        cmap(2,:)=clrLIP;
    otherwise
        clrMT=[70 144 205]/255;
        clrLIP=[85, 184, 74]/255;
        cmap=lines;
        cmap(1,:)=clrMT;
        cmap(2,:)=clrLIP;
end