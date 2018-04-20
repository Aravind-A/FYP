function out = pso_clustering(dataset, truth_labels, centroids, particles, iterations, hybrid_pso)
         
simtime=0.01;           % simulation delay btw each iteration
write_video = false;    % enable to grab the output picture and save a video

% VIDEO GRUB STUFF...
if write_video
    writerObj = VideoWriter('PSO.avi');
    writerObj.Quality=100;
%     writerObj.FrameRate=30;
    open(writerObj);
end

dataset_size = size (dataset);
dimensions = dataset_size(2);

if hybrid_pso
    fprintf('Running Matlab K-Means Version\n');
    [~,KMEANS_CENTROIDS] = kmeans(dataset,centroids, 'dist','sqEuclidean', 'display','iter','start','uniform','onlinephase','off');
    fprintf('\n');
end

% GLOBAL PARAMETERS
w  = 0.72; %INERTIA
c1 = 1.49; %COGNITIVE
c2 = 1.49; %SOCIAL

% PLOT STUFF... HANDLERS AND COLORS
pc = []; txt = [];
cluster_colors_vector = rand(particles, 3);

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

% SETTING UP PSO DATA STRUCTURES
swarm_vel = rand(centroids,dimensions,particles)*0.1;
swarm_pos = rand(centroids,dimensions,particles);
swarm_best = zeros(centroids,dimensions);
c = zeros(dataset_size(1),particles);
ranges = max(dataset)-min(dataset); %%scale
swarm_pos = swarm_pos .* repmat(ranges,centroids,1,particles) + repmat(min(dataset),centroids,1,particles);
swarm_fitness(1:particles)=Inf;

% KMEANS_INIT
if hybrid_pso
    swarm_pos(:,:,1) = KMEANS_CENTROIDS;
end

for iteration=1:iterations
      
    %CALCULATE EUCLIDEAN DISTANCES TO ALL CENTROIDS
    distances=zeros(dataset_size(1),centroids,particles);
    for particle=1:particles
        for centroid=1:centroids
            distance=zeros(dataset_size(1),1);
            for data_vector=1:dataset_size(1)
                distance(data_vector,1)=norm(swarm_pos(centroid,:,particle)-dataset(data_vector,:));
            end
            distances(:,centroid,particle)=distance;
        end
    end
    
    %ASSIGN MEASURES with CLUSTERS    
    for particle=1:particles
        [~, index] = min(distances(:,:,particle),[],2);
        c(:,particle) = index;
    end

    % PLOT STUFF... CLEAR HANDLERS
    delete(pc); delete(txt);
    pc = []; txt = [];
    
    % PLOT STUFF...
    hold on;
    for particle=1:particles
        for centroid=1:centroids
            if any(c(:,particle) == centroid)
                if dimensions == 3
                    pc = [pc plot3(swarm_pos(centroid,1,particle),swarm_pos(centroid,2,particle),swarm_pos(centroid,3,particle),'*','color',cluster_colors_vector(particle,:))];
                elseif dimensions > 3          % CHANGED
                    pc = [pc plot(swarm_pos(centroid,1,particle),swarm_pos(centroid,2,particle),'*','color',cluster_colors_vector(particle,:))];
                end
            end
        end
    end
    set(pc,{'MarkerSize'},{12})
    hold off;
 
    %CALCULATE GLOBAL FITNESS and LOCAL FITNESS:=swarm_fitness
    average_fitness = zeros(particles,1);
    for particle=1:particles
        for centroid = 1 : centroids
            if any(c(:,particle) == centroid)
                local_fitness=mean(distances(c(:,particle)==centroid,centroid,particle));
                average_fitness(particle,1) = average_fitness(particle,1) + local_fitness;
            end
        end
        average_fitness(particle,1) = average_fitness(particle,1) / centroids;
        if (average_fitness(particle,1) < swarm_fitness(particle))
            swarm_fitness(particle) = average_fitness(particle,1);
            swarm_best(:,:,particle) = swarm_pos(:,:,particle);     %LOCAL BEST FITNESS
        end
    end    
    [~, index] = min(swarm_fitness);       %GLOBAL BEST FITNESS
    swarm_overall_pose = swarm_pos(:,:,index);          %GLOBAL BEST POSITION
    
    % SOME INFO ON THE COMMAND WINDOW
    %fprintf('%3d. global fitness is %5.4f\n',iteration,global_fitness);            
    pause(simtime);
        
    % VIDEO GRUB STUFF...
    if write_video
        frame = getframe(fh);
        writeVideo(writerObj,frame);
    end
       
    % SAMPLE r1 AND r2 FROM UNIFORM DISTRIBUTION [0..1]
    r1 = rand;
    r2 = rand;
    
    % UPDATE CLUSTER CENTROIDS
    for particle=1:particles        
        inertia = w * swarm_vel(:,:,particle);
        cognitive = c1 * r1 * (swarm_best(:,:,particle)-swarm_pos(:,:,particle));
        social = c2 * r2 * (swarm_overall_pose-swarm_pos(:,:,particle));
        vel = inertia+cognitive+social;
                
        swarm_pos(:,:,particle) = swarm_pos(:,:,particle) + vel ;   % UPDATE PARTICLE POSE
        swarm_vel(:,:,particle) = vel;                              % UPDATE PARTICLE VEL
    end 
end

ari = clustereval(c(:,index), truth_labels, 'ari');
fprintf('ARI :  ');
disp(ari);

% PLOT THE ASSOCIATIONS WITH RESPECT TO THE CLUSTER 
hold on;
particle=index; %select the best particle (with best fitness) 
cluster_colors = ['m','g','y','b','r','c','g'];
for centroid=1:centroids
    if any(c(:,particle) == centroid)
        if dimensions == 3
            plot3(dataset(c(:,particle)==centroid,1),dataset(c(:,particle)==centroid,2),dataset(c(:,particle)==centroid,3),'o','color',cluster_colors(centroid));
        elseif dimensions > 3              % CHANGED
            plot(dataset(c(:,particle)==centroid,1),dataset(c(:,particle)==centroid,2),'o','color',cluster_colors(centroid));
        end
    end
end
hold off;

% VIDEO GRUB STUFF...
if write_video
    frame = getframe(fh);
    writeVideo(writerObj,frame);
    close(writerObj);
end

%fprintf('\nEnd, global fitness is %5.4f\n',global_fitness);
out=1;
end