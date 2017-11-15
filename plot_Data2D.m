function plot_Data2D(M)

plotDimensions = 1:size(M, 1 );

scatter(M(:, plotDimensions(1)), ...
        M(:, plotDimensions(2)), ...
        'o',...
        'MarkerFaceColor', 'k');
%     plot(M);
hold off;
