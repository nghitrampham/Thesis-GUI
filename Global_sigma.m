% only need to change sigma, data, and numClusters
close all;
% hold on 

data = normalizeData(mat4SVM');

% numCluster = 3;
cluster_SD= [];
distortion= [];
% cluster_LS = [];
final = {};
clusters = {};
numOfCluster = [2,3,4,5,6];
best_sigma = 0;
best_distortion = 1/eps;
final_distortion = [];
% numOfCluster = [];
% sigma = [0.01, 0.1];
sigma = [0.01, 0.03, 0.06,0.08, 0.1, 0.4, 0.6,0.7,0.8, 0.9,1.2, 1.4, 1.6, 2];
% sigma = [0.9];
for k = 1: length(numOfCluster)
best_sigma = 0;
numCluster = numOfCluster(k);
disp('running new k');
best_distortion = 1/eps;
for j = 1:length(sigma)
disp(['running new sigma', num2str(sigma(j))]);
cluster_SD= [];
distortion =[];
W = SimGraph_NearestNeighbors(data, 10, 2, sigma(j));
[ L, U, eigenvalues, cidx, bestD] = Spectral_test(W, numCluster, 3);
distortion(1) = bestD;

if bestD < best_distortion
    best_distortion = bestD;
    best_sigma = sigma(j);
end


for i = 1:length(cidx)
    cluster_SD(cidx{i},i) = 1;
end
disp('Value stored')
clusters = {};
clusters{1} = cluster_SD;
final{k}{j} = {distortion, clusters}; 
end
final_sigma(k) = best_sigma;
final_distortion(k) = best_distortion;
end 
 



% for ii = 1:length(numOfCluster)
%     for jj = 1:length(sigma)
%     KmeansDist_SD(ii,jj) = first_run{ii}{jj}{1}{1};
%     KmeansDist_LS(ii,jj) = first_run{ii}{jj}{1}{2};
%     end 
% end 