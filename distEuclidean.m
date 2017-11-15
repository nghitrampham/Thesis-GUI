function [ dist ] = distEuclidean( M, N )

dist = sqrt(sum((M - N) .^ 2, 1));

end

