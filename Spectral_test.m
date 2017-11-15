function [ L, U, eigenvalues, cidx, bestD] = Spectral_test(W, k, Type)

% @param W, adjancency matrix of the similarity graph 
% @param k, number of clusters
% @param 'type', defines which type of spectral clustering will be used
%                1 - Unormalized 
%                2 - Normalized according to Shi and Malik algorithm (2000)
%                3 - Normalized according to Jordan and Weiss algorithm (2002)

% @return L, Laplacian matrix 
% @return U, matrix containing the eigenvectors 
% @return C, matrix containing k cluster vectors

% calculate degree matrix
degs = sum(W, 2);
D    = sparse(1:size(W, 1), 1:size(W, 2), degs);

L = D - W;

switch Type
    case 1 
        L = D - W;
    case 2
        degs(degs == 0) = eps;
        D = spdiags(1./degs, 0, size(D, 1), size(D, 2));
        
        L = D * L;
    case 3
        degs(degs == 0) = eps;
        D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
        
        L = D * L * D;
end

% compute the eigenvectors corresponding to the k smallest
% eigenvalues
diff   = eps;
[U, eigenvalues] = eigs(L, k, diff);

% in case of the Jordan-Weiss algorithm, we need to normalize
% the eigenvectors row-wise
if Type == 3
    U = bsxfun(@rdivide, U, sqrt(sum(U.^2, 2)));
end

% now use the k-means algorithm to cluster U row-wise
% C will be a n-by-1 matrix containing the cluster number for
% each data point
% C = kmeans(U, k, 'start', 'cluster', ...
%                  'EmptyAction', 'singleton');
[clusts,bestD]= kmeans(U,k);
cidx = clusts;
         
% now convert C to a n-by-k matrix containing the k indicator
% vectors as columns
% C = sparse(1:size(D, 1), C, 1);
% C = full(C);

end