function handles = plotClusterStarCoordinates(hObject, handles)

% plotDimensions = get(handles.lstPlotDimensions, 'Value');
% figure()

% if additional VARIABLES is looaded, display average DQ and gender
if isempty(handles.ArbitraryFeature) == 0
    handles.PerFeature = handles.ArbitraryFeature;
end

% as soon as addiontional variables are loaded in, plot DQ, gender, and
% average DQ for each group

if istable(handles.VARIABLES) == 1
    handles.store = handles.VARIABLES(:,15);
    DQ = table2array(handles.store)';
    
    handles.store1 = handles.VARIABLES(:,1);
    gender = string(table2array(handles.store1))';
    
    gender_count = [];
    if ~isequal(handles.current_NumCluster,0)
        for cluster = 1:handles.current_NumCluster
            idx = handles.cidx == cluster;
            handles.avgCluster_DQ = nanmean(DQ(idx));
            handles.avgperCluster_DQ(idx) = handles.avgCluster_DQ;
            %%%%%keep track of how many male and female%%%%%%
            gender_count(cluster, 1) = sum(gender(idx) == "Female");
            gender_count(cluster, 2) = sum(gender(idx) == "Male");
        end
        PerCluster_DQ = handles.avgperCluster_DQ;
    else
        PerCluster_DQ = 0;
    end
end

%%%%% from Cluster_Data_Button_Callback %%%%%
% total number of points in a cluster
PointInCluster = handles.numberPerCluster;
% average DQ for each cluster

% average for arbitrary feature chosen
%%%%%%%%%%%%%%%%%%%%%%%


plotDimensions = 1:size(handles.DATA,2);

workData = normalizeDATA(handles.DATA');
workData = workData';
workData = workData(plotDimensions, :);

d = size(workData, 1);
unitRoot = exp(2 * pi * 1i / d);

starAxes = zeros(1, d);
for ii = 1:d
    starAxes(ii) = unitRoot^ii;
end

starAxes = repmat(starAxes', 1, size(workData, 2));
plotPoints = sum(workData .* starAxes, 1);

gscatter(real(plotPoints), imag(plotPoints), handles.cidx, ...
    'brgcmykw', [], 10, 'off');
hold on

%%%%%%%%%
% text(real(plotPoints)+ 0.2, imag(plotPoints) + 0.2i, num2str(handles.classification))
%        hold on
%     text(M(ClusteredData(:, ii) == 1, plotDimensions(1)), ...
%             M( ClusteredData(:, ii) == 1, plotDimensions(2)),...
%             num2str(label(ClusteredData(:, ii) == 1, 1)))
%         hold on
%%%%%%%%%%

if handles.contents ~= 0
    handles.dcm_obj = datacursormode(gcf);
    datacursormode on
    group = handles.cidx;
    set(handles.dcm_obj,'UpdateFcn',{@myupdatefcn_variables,  DQ, gender, ...
        PointInCluster, PerCluster_DQ, handles.PerFeature, group, ...
        gender_count})
    
elseif handles.contents == 0 && istable(handles.VARIABLES) == 1
    handles.dcm_obj = datacursormode(gcf);
    datacursormode on
    group = handles.cidx;
    set(handles.dcm_obj,'UpdateFcn',{@myupdatefcn,  DQ, gender, ...
        PointInCluster, PerCluster_DQ, group, gender_count})
elseif istable(handles.VARIABLES) == 0
    handles.dcm_obj = datacursormode(gcf);
    datacursormode on
    
    set(handles.dcm_obj,'UpdateFcn',{@myupdatefcn_normal, PointInCluster ...
        })
else
        handles.dcm_obj = datacursormode(gcf);
    datacursormode on
    
    set(handles.dcm_obj,'UpdateFcn',{@myupdatefcn_normal, PointInCluster ...
        })
    
end


    function txt = myupdatefcn(~,event_obj, DQ, gender, PointInCluster,...
            PerCluster_DQ,group, gender_count)
        % Customizes text of data tips
        pos = get(event_obj,'Position');
        I = get(event_obj, 'DataIndex');
        txt = {['X: ',num2str(pos(1))],...
            ['Y: ',num2str(pos(2))],...
            ['Idx: ',num2str(I)], ...
            ['DQ:',num2str(DQ(I))], ...
            ['Gender:', gender(I)], ...
            ['TotalPoint:', num2str(PointInCluster(I))], ...
            ['Avg DQ:', num2str(  PerCluster_DQ(I))],...
            ['ClusterGroup:' , num2str(group(I))], ...
            ['Num_Female:' , num2str(gender_count(group(I), 1))], ...
            ['Num_Male:' , num2str(gender_count(group(I), 2))]};
    end


    function txt = myupdatefcn_variables(~,event_obj, DQ, gender, PointInCluster, ...
            PerCluster_DQ, PerFeature,...
            group, gender_count)
        % Customizes text of data tips
        pos = get(event_obj,'Position');
        index = find(imag(plotPoints) == pos(2) & real(plotPoints) == pos(1));
        I = index;
        % I = get(event_obj, 'DataIndex');
        
        txt = {['X: ',num2str(pos(1))],...
            ['Y: ',num2str(pos(2))],...
            ['Idx: ',num2str(I)], ...
            ['DQ:',num2str(DQ(I))], ...
            ['Gender:', gender(I)], ...
            ['TotalPoint:', num2str(PointInCluster(I))], ...
            ['Avg DQ:', num2str(PerCluster_DQ(I))], ...
            ['Feature Avg:' , num2str(PerFeature(I))], ...
            ['ClusterGroup:' , num2str(group(I))], ...
            ['Num_Female:' , num2str(gender_count(group(I), 1))], ...
            ['Num_Male:' , num2str(gender_count(group(I), 2))]};
    end

    function txt = myupdatefcn_normal(~, event_obj, PointInCluster)
        pos = get(event_obj,'Position');
        I = get(event_obj, 'DataIndex');
        txt = {['X: ',num2str(pos(1))],...
            ['Y: ',num2str(pos(2))],...
            ['Idx: ',num2str(I)], ...
            ['TotalPoint:', num2str(PointInCluster(I))]};
    end

hold off;

guidata(hObject, handles);
end

%
% ['Gender:', num2str(gender(I))], ...
