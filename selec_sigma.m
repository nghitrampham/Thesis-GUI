function final_sigma = selec_sigma(M, k)

default_sigma = 0.5;
% min_eps = 0; 
% max_eps = 100;
% M is distance matrix 
% k = number of nearest neighbors 


        sorted = sort(M, 2, 'ascend');
        holder = sorted(:, k+1);
        calculated_sigma = mean(holder);
        
if calculated_sigma == 0 || calculated_sigma == Inf
    final_sigma = default_sigma;
else 
    final_sigma = calculated_sigma;
end
% if calculated_eps >= min_eps && calculated_eps <= max_eps
%     final_sigma = calculated_eps;
% else 
%     final_sigma = default_eps;
% end 