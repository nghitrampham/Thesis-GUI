function W = epsilon_SimGraph(M, epsilon, matrix_type)
% @return W, adjacency matrix for epsilon similarity graph 
% if M is distance matrix, it will return an unweighted graph 
% matrix type, 1-similarity matrix 2 - distance matrix
n = size(M,2); m = size(M,1);

if n ~= m 
    error('Matrix needs to be square');
end

switch matrix_type 
    case 1 % similarity matrix 
        M = update_diagonal(M, 0);
        
        index = M < epsilon;
        M(index) = 0; %adjacency matrix, weighted by similarity 
        
    case 2 % distance matrix - for unweighted graph 
        M = update_diagonal(M, Inf);
        M = (M <= epsilon); %unweighted graph
        double(M);
        
    otherwise 
        error('Invalid input matrix');
        return; 
end 
W = M;