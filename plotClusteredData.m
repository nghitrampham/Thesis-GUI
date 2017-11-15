% function plotClusteredData(M, ClusteredData, label, numCluster)
function plotClusteredData(M, ClusteredData, numCluster)
%only plot the first two dimensions 
% @param ClusteredData, matrix containing the index for cluster groups 
% @M, original matrix, n-by-d matrix 
plotDimensions = 1:size(M,2);
colors = 'brgcmykw';
% colors = 'gmbrcykw';



hold on 
for ii = 1:numCluster
    currentColor = [colors(mod(ii, size(colors, 2))) '*'];
    
    scatter(M(ClusteredData(:, ii) == 1, plotDimensions(1)), ...
            M( ClusteredData(:, ii) == 1, plotDimensions(2)), ...
            7, currentColor);
        hold on 
%     text(M(ClusteredData(:, ii) == 1, plotDimensions(1)), ...
%             M( ClusteredData(:, ii) == 1, plotDimensions(2)),...
%             num2str(label(ClusteredData(:, ii) == 1, 1)))
%         hold on 
end

hold off;