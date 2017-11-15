function M = Gaussian(M,sigma)
% M, squared distance matrix 
M = exp(-M.^2 ./ (2*(sigma^2)));
