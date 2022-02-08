function plot_twoday_data(plane1, plane2, plane3)
% Aashish 
% preselected pairs
%% Correct red channel based on green-red relationship
%alter dd based upon where the selected.mat files are found for each mouse
    dd = dir('000-003\*selected*.mat');
% dd = dir('002-005\*selected*.mat');
split_file_name = strsplit(folder1, '\');
mouse_name = split_file_name{3};
exp_date = split_file_name{4};
for ii=1:3 
        % ii is number of the plane (1 to 3)

    se{ii} = load(plane num2str{});
        %loads file containing information about pre and post cells that are
        %matched across all 3 sessions
end

%% concatenate all data % matched pairs for all 3 planes - all 3 sessions 

%just f-fneu and stat
neur_coeff = 0.7;
pre.F =[];
pre.stat = []; 
post.F = [];
post.stat = [];

for ii=1:3
    pre.F = cat(1, pre.F, se{ii}.F_pre-neur_coeff*se{ii}.Fneu_pre);
    post.F = cat(1, post.F, se{ii}.F_post-neur_coeff*se{ii}.Fneu_post);
    pre.stat = cat(1,pre.stat,se{ii}.stat_pre);
    post.stat = cat(1,post.stat,se{ii}.stat_post);
end

npair = size(pre.F(1));
nframe = size(pre.F,2); 