load('class_Af_Am.mat')
[S,comp] = graphconncomp(sparse(handles.SIM), 'Directed', false);
% [S,comp] = graphconncomp(sparse(Spectral.W), 'Directed', false);
comp = comp';
% sum(comp(1:55)==3)


label = cellstr(regionName);close all;

group = {};
female = {};
male = {};
mean_f = {};
mean_m = {};
mean_both = {};
idx = 1;
%get female for each group 
for i = 1:handles.current_NumCluster
    group{i} = mat4SVM(handles.cidx == i,:);
    mean_region = nanmean(group{i},1);
    
    figure(idx)
    barh(mean_region(1:25));
    set(gca, 'YTick', 1:25);set(gca, 'YTickLabel', label(1:25))
    title({strcat('Group ', num2str(i)),...
                    strcat('Total: ', num2str(sum(handles.cidx == i)), '  Male: ', num2str(sum(handles.cidx(56:end) == i)), ...
                    '  Female: ', num2str(sum(handles.cidx(1:55) == i)))}) 
    idx = idx+1;
    
    female{i} = mat4SVM(handles.cidx(1:55) == i,:);
    mean_f = nanmean(female{i},1);
    figure(idx)
    barh(mean_f(1:25));
    set(gca, 'YTick', 1:25);set(gca, 'YTickLabel', label(1:25))
    title(strcat('Group ', num2str(i), ' Female Average')) 
    idx = idx+1;
    male{i} = mat4SVM(handles.cidx(56:end) == i,:);
    mean_m = nanmean(male{i},1);
    figure(idx)
    barh(mean_m(1:25));
    set(gca, 'YTick', 1:25);set(gca, 'YTickLabel', label(1:25))
    title(strcat('Group ', num2str(i), ' Male Average'))
    idx = idx+1;
end 