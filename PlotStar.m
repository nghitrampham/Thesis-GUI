function PlotStar
n = size(mat4SVM,2);
mapColor  = 'Hot';
EdgeCol   = 'b';
VertixCol = 'r';
plotColored = 0;
dim1Val = zeros(2*n, 1);
dim2Val = zeros(2*n, 1);

if ~isequal(plotColored, 0)
    sMin = min(min(nonzeros(W)));
    sMax = max(max(W));
    
    if isequal(sMin, sMax)
        sMin = 0;
    end
    
    map = colormap(mapColor);
end

for ii = 1:n
    nnzs = find(W(ii, :) > 0);
    nnzsLength = length(nnzs);
    
    if nnzsLength > 0
        dim1Val(1:2:2*nnzsLength) = mat4SVM(1, ii);
        dim1Val(2:2:2*nnzsLength) = mat4SVM(1, nnzs);
        dim2Val(1:2:2*nnzsLength) = mat4SVM(2, ii);
        dim2Val(2:2:2*nnzsLength) = mat4SVM(2, nnzs);
        
        if isequal(plotColored, 0)
            plot(...
                dim1Val(1:2*nnzsLength), ...
                dim2Val(1:2*nnzsLength), ...
                ['-' EdgeCol]);
        else
            for jj = 1:2:2*nnzsLength
                tempColor = 1 - ...
                    (W(ii, nnzs((jj+1)/2)) - sMin) / (sMax - sMin);
                currentColor = map(1 + floor((size(map, 1) - 1) * tempColor), :);
                
                plot(...
                    dim1Val(jj:jj+1), dim2Val(jj:jj+1), '-', 'Color', 'b');
            end
        end
    end
end

scatter(...
    mat4SVM(1, :), ...
    mat4SVM(2, :), ...
    getPlotMarkerSize(), VertixCol, 'filled');

hold off;