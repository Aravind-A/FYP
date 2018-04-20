function out = kmeans_clustering(dataset, truth_labels, centroids)

    dataset_size = size(dataset);
    dimensions = dataset_size(2);

    % PLOT DATASET
    fh=figure(1);
    hold on;
    if dimensions == 3
        plot3(dataset(:,1),dataset(:,2),dataset(:,3),'k*');
        view(3);
    elseif dimensions > 3                % CHANGED
        plot(dataset(:,1),dataset(:,2),'k*');
    end

    % PLOT STUFF .. SETTING UP AXIS IN THE FIGURE
    axis equal;
    axis(reshape([min(dataset(:,[1 2]))-2; max(dataset(:,[1 2]))+2],1,[]));    % CHANGED
    hold off;

    [idx,KMEANS_CENTROIDS] = kmeans(dataset,3, 'dist','sqEuclidean', 'display','iter','start','uniform','onlinephase','off');
    ari = clustereval(idx, truth_labels, 'ari');
    fprintf('ARI :  ');
    disp(ari);
                
    % PLOT THE ASSOCIATIONS WITH RESPECT TO THE CLUSTER 
    hold on; 
    cluster_colors = ['m','g','y','b','r','c','g'];
    for centroid=1:centroids  
        if dimensions == 3
            plot3(dataset(idx==centroid,1),dataset(idx==centroid,2),dataset(idx==centroid,3),'o','color',cluster_colors(centroid));
        elseif dimensions > 3              % CHANGED
            plot(dataset(idx==centroid,1),dataset(idx==centroid,2),'o','color',cluster_colors(centroid));
        end
    
    end
    hold off;
    out = 1;
end