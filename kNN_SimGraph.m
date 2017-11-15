function W = kNN_SimGraph(M, k, type, matrix_type, sigma)


% @param M 
% @param k 
% @param type, 'normal'  or 'mutual' 

% both cases are undirected graph (one weighted and another one is not)
% k = 7;

n = size(M,1); m = size(M,2);

if n ~= m
    error('Invalid input');
end 

switch matrix_type 
    case 1 % similarity matrix 
        M = update_diagonal(M, 0);
        for i =1:n
        [~, index] = sort(M(i,:), 'descend');
        M(i, index(k+1:end)) = 0;
        end 
    case 2 % distance matrix - for unweighted graph 
        M = update_diagonal(M, Inf);
        for j = 1:n
        [~, index] = sort(M(j,:), 'ascend');
        M(j, index(k+1:end)) = 0;
         M(j, index(1:k)) = 1;
        end

    otherwise 
        error('Invalid input matrix type for kNN');
       
end 
W = M;

if strcmp(type, 'normal') == 1
    W = max(W, W');
elseif strcmp(type, 'mutual') == 1
    W = min(W, W');
else 
    error('Invalid type of similarity graph for kNN')
end 

if sigma ~= 0 && matrix_type ==2
   W = spfun(@(W) (Gaussian(W, sigma)), W);
   W = full(W);
end 