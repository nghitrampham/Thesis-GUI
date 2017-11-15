function final_eps = select_eps(M,  matrix_type, k)

% choose based on the average of k nearest neighborhood 
% @param matrix_type, 1- distance matrix, 2- similarity matrix 
% M distance matrix 
default_eps = 0.5;
% min_eps = 0; 
% max_eps = 100;

switch matrix_type 
    case 1 %distance matrix 
        sorted = sort(M, 2, 'ascend');
        holder = sorted(:, k+1);
        calculated_eps = mean(holder);
    case 2 % similarity matrix 
        sorted = sort(M, 2, 'descend');
        holder = sorted(:, k+1);
        calculated_eps = mean(holder);
    otherwise 
        calculated_eps = default_eps;
end 

final_eps = calculated_eps;
% if calculated_eps >= min_eps && calculated_eps <= max_eps
%     final_eps = calculated_eps;
% else 
%     final_eps = default_eps;
% end 