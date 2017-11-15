function plotData3D(DATA)
DATA = DATA';
% plotDimensions = get(handles.lstPlotDimensions, 'Value');
plotDimensions = 1:size(DATA,2); 


view(3);

scatter3(DATA(plotDimensions(1), :), ...
         DATA(plotDimensions(2), :), ...
         DATA(plotDimensions(3), :), ...
         10, ...
         ['k' '*']);
    
hold off;

% guidata(hObject, handles);