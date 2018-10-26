%%
all_tt_nums = [];


cd(Se.folder{session_num})
mkdir('AnalyzedFiles')
names = dir('*.cut');

for n = 1:length(names)
    nameSU = names(n).name;
    nameSU(strfind(nameSU,'.cut'):length(nameSU)) = [];
    nameSU(1:2) = [];
    all_tt_nums = [all_tt_nums str2num(nameSU)];
end

all_tt_nums = sort(all_tt_nums); %Automatic selection of electrode sets to plot

%%  Auto Import -- manual clust

fprintf('Extracting Data....');

i = 0;

for n = 1:(length(all_tt_nums))
    curr_tt_num = all_tt_nums(n);
    Sc = spikes(Se,session_num,curr_tt_num);
    is_cluster = 1;
    m = 0;
    while is_cluster == 1
        m = m+1;
        cl_holder = cluster(Sc,m);
        is_full = max(size(cl_holder.timestamp));
        if is_full > 1
            i = i+1;
            eval(sprintf('cl%d = cluster(Sc,m)', i));
            num_seq(i,1:2) = [curr_tt_num m];
        else
            m = m-1;
            is_cluster = 0;
        end
    end
    Sc_unit_count(n,1) = curr_tt_num;
    Sc_unit_count(n,2) = m;
end

clear cl_holder i curr_tt_num is_cluster is_full m n
fprintf('....Done!\n');