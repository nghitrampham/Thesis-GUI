function [U, eval, cidx] =jordan_spectral(data, sigma, k)
% sigma = 2;
data = normalizeDATA(data');
A = SimGraph_NearestNeighbors(data', 14, 2, 3.2);
% dist = distance(data, data);
% A = Gaussian(dist, sigma);
A = update_diagonal(A, 0);
degree = sum(A, 2);

% degree    = sparse(1:size(A, 1), 1:size(A, 2), degree);
degree(degree==0) = eps;
% D = spdiags(1./(degree.^0.5), 0, size(A, 1), size(A, 2));
D = 1./(degree.^(0.5));
D = diag(D);
L = D*A*D;
[U, eval] = eigs(L, k, eps);
U = bsxfun(@rdivide, U, sqrt(sum(U.^2, 2)));

C = kmeans(U, k, 'start', 'cluster', ...
                 'EmptyAction', 'singleton');
cidx = C;