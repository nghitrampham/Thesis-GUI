function A = update_diagonal(A, value)

if size(A,1) ~= size(A,2)
    error('Not square matrix')
end 

A(logical(eye(size(A)))) = value;