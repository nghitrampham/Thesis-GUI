function D = distance(M,N)
% return a SQUARED Euclidean distance between data points 
% M, N are n-by-d matrix containing n d dimensional data points  
M_squared = sum(M.*M,2);
N_squared = sum(N'.*N',1);
D = M_squared(:,ones(1,size(N,1))) + N_squared(ones(1,size(M,1)),:) - 2*M*N';

 D = update_diagonal(D, 0);
 D = sqrt(D);