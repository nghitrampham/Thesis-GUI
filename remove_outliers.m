function mat4SVM_outliers = remove_outliers(mat4SVM)
mat4SVM_outliers = mat4SVM;
% mat4SVM_outliers = mat4SVM;
for i = 1:size(mat4SVM, 2)
%     
%     def replace(group):
%     mean, std = group.mean(), group.std()
%     outliers = (group - mean).abs() > 3*std
%     group[outliers] = mean        % or "group[~outliers].mean()"
%     return group
      avg = mean(mat4SVM(:,i));
      stand_dev = std(mat4SVM(:,i));
      outliers = abs(mat4SVM(:,i) - avg) > 2.5*stand_dev;
      mat4SVM_outliers(outliers,i) = avg;
end 
% csvwrite('mat4SVM_outliers.csv', mat4SVM_outliers);
% save('mat4SVM_outliers.mat', 'mat4SVM_outliers');