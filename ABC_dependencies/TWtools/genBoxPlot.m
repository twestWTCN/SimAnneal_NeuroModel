function a = genBoxPlot(G,cmap,orient)
boxplot(G,'colors',cmap,'boxstyle','filled', 'jitter',0.2,'orientation',orient)
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects

% Make Boxes Wide
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',25); % Set width

% Make Median lines Wide
idx=strcmpi(t,'Median');  % Find Box objects
medlines=a(idx);          % Get the children you need
set(medlines,'linewidth',2); % Set width
for i = 1:size(medlines,1)
    medlines(i).Color = cmap((end+1)-i,:).*0.8;
end

% Make Whisker lines Wide
idx=strcmpi(t,'Whisker');  % Find Box objects
whisklines=a(idx);          % Get the children you need
set(whisklines,'linewidth',2); % Set width

% Outliers Filled Circles
idx=strcmpi(t,'Outliers');  % Find Box objects
ols=a(idx);          % Get the children you need
for i = 1:size(ols,1)
ols(i).Marker = 'none';
ols(i).MarkerEdgeColor = 'none';
    ols(i).MarkerFaceColor = cmap((end+1)-i,:);
end
