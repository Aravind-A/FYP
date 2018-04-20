function out = project_driver(dataset, algo)
    close all;
    switch dataset
        case 'iris'  
            load fisheriris.mat
            dataset_size = size(meas);
            O = dataset_size(1);
            truth_labels = zeros(O,1);
            for i=1:50
                truth_labels(i) = 1;
            end
            for i=51:100
                truth_labels(i) = 2;
            end
            for i=101:150
                truth_labels(i) = 3;
            end
            if strcmp(algo,'pso') == 1
                out = pso_clustering(meas, truth_labels, 3, 20, 100, true);
            elseif strcmp(algo,'dpso') == 1
                out = density_based_pso(meas, truth_labels, 3, 100);
            elseif strcmp(algo,'kmeans') == 1
                out = kmeans_clustering(meas, truth_labels, 3);
            end
        case 'glass'    
            file = csvread('glass.data.txt');
            s = size(file);
            truth_labels = file(:,s(2));
            file = file(:,2:s(2)-1);
            if strcmp(algo,'pso') == 1
                out = pso_clustering(file, truth_labels, 6, 20, 100, true);
            elseif strcmp(algo,'dpso') == 1
                out = density_based_pso(file, truth_labels, 6, 100);
            elseif strcmp(algo,'kmeans') == 1
                out = kmeans_clustering(file, truth_labels, 3);
            end
        case 'seeds'
            file = csvread('seeds.data.txt');
            s = size(file);
            truth_labels = file(:,s(2));
            file = file(:,1:s(2)-1);
            if strcmp(algo,'pso') == 1
                out = pso_clustering(file, truth_labels, 3, 20, 100, true);
            elseif strcmp(algo,'dpso') == 1
                out = density_based_pso(file, truth_labels, 3, 100);
            elseif strcmp(algo,'kmeans') == 1
                out = kmeans_clustering(file, truth_labels, 3);
            end
        case 'wine'
            file = csvread('wine.data.txt');
            s = size(file);
            truth_labels = file(:,1);
            file = file(:,2:s(2));
            if strcmp(algo,'pso') == 1
                out = pso_clustering(file, truth_labels, 3, 20, 100, true);
            elseif strcmp(algo,'dpso') == 1
                out = density_based_pso(file, truth_labels, 3, 100);
            elseif strcmp(algo,'kmeans') == 1
                out = kmeans_clustering(file, truth_labels, 3);
            end
    end
end