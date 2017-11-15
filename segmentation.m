function segmentation(FileName, k, Neighbors)



roundColors = 0;        % Round color values for less strict uniqueness
roundDigits = 2;        % Precision for Uniqueness
% saveData    = 0;        % Save Dataset
markEdges   = 0;        % Outline edges



FileName = fullfile(FileName);

Img = imread(FileName);
[m, n, d] = size(Img);

% convert into list of data points
Data = reshape(Img, 1, m * n, []);

if d >= 2
    Data = (squeeze(Data))';
end

% convert to double and normalize to [0,1]
Data = double(Data);
Data = normalizing(Data);


% Find unique colors
if isequal(roundColors, 1)
    fac = 10^roundDigits;
    rData = round(Data * fac) / fac;
else
    rData = Data;
end

[~, ind, order] = unique(rData', 'rows', 'R2012a');

% crop data
Data = Data(:, ind);

% now for the clustering
fprintf('Busy Creating Similarity Graph...\n');
% Data = distance(Data, Data);
SimGraph = SimGraph_NearestNeighbors(Data, Neighbors, 1);
% SimGraph = kNN_2(Data, Neighbors, 1);

% try
%     comps = graphconncomp(SimGraph, 'Directed', false);
%     fprintf('- %d connected components found\n', comps);
% end

fprintf('Busy Clustering Data...\n');
C = Spectral_test(SimGraph, k, 2);
fprintf('Done Clustering Data...\n');
% convert and restore full size
D = convertClusterVector(C);
D = D(order);

% reshape indicator vector into m-by-n
S = reshape(D, m, n);

% choose colormap
if k == 2
    map = [0 0 0; 1 1 1];
else
    map = zeros(3, k);
    
    for ii = 1:k
        ind = find(D == ii, 1);
        map(:, ii) = rData(:, ind);
    end
    
    map = map';
end

% plot image
figure(2)
set(gca, 'Position', [0 0 1 1], 'Units', 'Normalized');

if isequal(markEdges, 1)
    imshow(Img, 'Border', 'tight');
    
    lS = label2rgb(S);
    BW = im2bw(lS, graythresh(lS));
    [B, L] = bwboundaries(BW, 'holes');

    hold on;
    
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:, 2), boundary(:, 1), 'r', 'LineWidth', 2)
    end
    
    hold off;
else
    imshow(S, map, 'Border', 'tight');
end

hold on;

axis off;
truesize;
hold off;

clear all;