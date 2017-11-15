function [C,L,U, eigvalue, cidx] = SpectralClustering(W,  k, spectral_type)

% @param W, adjancency matrix of the similarity graph 
% @param k, number of clusters
% @param 'type', defines which type of spectral clustering will be used
%                1 - Unormalized 
%                2 - Normalized according to Shi and Malik algorithm (2000)
%                3 - Normalized according to Jordan and Weiss algorithm (2002)

% @return L, Laplacian matrix 
% @return U, matrix containing the eigenvectors 
% @return C, matrix containing k cluster vectors

% degree matrix 
 degree = sum(W, 2); 
 


D = diag(degree, 0);
L = D-W;
D(D == 0) = eps;
switch spectral_type 
    case 1 
        L = D-W;
    case 2

        D = 1./D;
        
        L = D*L; 
        
    case 3
        D = 1./ (D.^(0.5));
        L = D*L*D;
    otherwise 
        error('Invalid spectral type for spectral clustering');
end 

% compute eigenvectors 
[U, eigvalue] = eigs(L, k, eps);

if spectral_type == 3
    for i = 1:size(U, 1)
        n = sqrt(sum(U(i,:).^2));
        U(i,:) = U(i,:) ./ n;
    end 
end 

% apply K-means 
C = kmeans(U,k, 'start', 'cluster', 'EmptyAction', 'singleton');
cidx = C; 
% retrieve cluster data points 
C = sparse(1:size(D, 1), C, 1);
% C = full(C);

% if size(C, 2) > 1
%     indMatrix = zeros(size(C, 1), 1);
%     for ii = 1:size(C, 2)
%         indMatrix(C(:, ii) == 1) = ii;
%     end
% else
%     indMatrix = sparse(1:size(C, 1), C, 1);
% end