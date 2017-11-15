% data = X;
% sigma = [0.01, 0.05, 0.08,0.1, 0.4, 0.5,0.7, 0.9,1.2, 1.4, 1.6];
% sigma = [0.01, 0.05,0.08];
sigma = 0.7;
distortion = 1000;
best_distortion = [];
cluster = [];
for i = 1:length(sigma)
    W = SimGraph_NearestNeighbors(X', 8, 2, sigma(i));
    for j = 1:5
        [L, U, eigenvalues, cidx, bestD] = Spectral_test(W, 2, 3);
        if bestD < distortion 
            best_distortion(i) = bestD;
            distortion = bestD;
            cluster{i}= cidx;
            
        end
    end
end