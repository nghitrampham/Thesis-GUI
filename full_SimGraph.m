function W = full_SimGraph(M, sigma)
% @return W, adjacency matrix for full similarity graph 
% @param M, n-by-d matrix containing n d-dimensional data points 
% @param sigma, parameter used to control the size of neiborhoods in
% Gaussian similarity function 

% @return W, adjacency matrix of full similarity graph 

% calculate the distance between points using Euclidean metric 
% W = squareform(pdist(M, 'Euclidean'));
W = distance(M,M);
% Applying Gaussian similarity function
W = exp(-W.^2 ./ (2*sigma^2));