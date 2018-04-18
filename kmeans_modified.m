clear;
close all;

k = 3;
dimensions = 4;
iterations = 100;

load fisheriris.mat
dataset_size = size(meas);
n = dataset_size(1);
c_indices = ceil(rand(k,1)*n);
centroids = meas(c_indices,:);
disp('initial centroid indices : ');
disp(c_indices);
distances = zeros(n,k);
mapping = zeros(n,k);
threshold = zeros(k,1);

for i=1:n
    for j=1:k
        distances(i,j) = norm(meas(i,:) - centroids(j,:));
    end
    [dist, index] = min(distances(i,:));
    mapping(i,index) = 1;
    %distances(i,k+1) = index;
    %distances(i,k+2) = dist;
end

for it=2:iterations
    prev_centroids = centroids;
    for i=1:k
        threshold(i,1) = mean(distances(mapping(:,i) == 1, i));
        %threshold(i,1) = mean(distances(distances(:,k+1) == i, k+2));
        centroids(i,:) = mean(meas(mapping(:,i) == 1, :));
    end
    if sum(sum(prev_centroids == centroids)) == 0 && it > 2
        X = ['Iter :', num2str(it)];
        disp(X);
        break;
    end
    for i=1:n
        flag = 0;
        for j=1:k
            distances(i,j) = norm(meas(i,:) - centroids(j,:));
            if distances(i,j) <= threshold(j,1)
                mapping(i,j) = 1;
                flag = 1;
            else
                mapping(i,j) = 0;
            end
        end
        if flag == 0
            for j=1:k
                mapping(i,j) = 0;
            end
            [dist, index] = min(distances(i,:));
            mapping(i,index) = 1;
        end
    end
end