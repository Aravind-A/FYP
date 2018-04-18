clear;
close all;

%rng('default') % For reproducibility

% INIT PARTICLE SWARM
centroids = 3;          % == clusters here (aka centroids)
dimensions = 4;         % how many dimensions in each centroid  CHANGED
particles = 70;         % how many particles in the swarm, aka how many solutions
iterations = 50;        % iterations of the optimization alg.
simtime=0.01;           % simulation delay btw each iteration
dataset_subset = 2;     % for the IRIS dataset, change this value from 0 to 2
write_video = false;    % enable to grab the output picture and save a video
hybrid_pso = true;     % enable/disable hybrid_pso
manual_init = false;    % enable/disable manual initialization (only for dimensions={2,3})

% VIDEO GRUB STUFF...
if write_video
    writerObj = VideoWriter('PSO.avi');
    writerObj.Quality=100;
%     writerObj.FrameRate=30;
    open(writerObj);
end

% LOAD DEFAULT CLUSTER (IRIS DATASET); USE WITH CARE!
load fisheriris.mat
%meas = meas(:,1+dataset_subset:dimensions+dataset_subset); % CHANGED
dataset_size = size (meas);

% EXECUTE K-MEANS
if hybrid_pso
    fprintf('Running Matlab K-Means Version\n');
    [idx,KMEANS_CENTROIDS] = kmeans(meas,centroids, 'dist','sqEuclidean', 'display','iter','start','uniform','onlinephase','off');
    fprintf('\n');
end

% GLOBAL PARAMETERS (the paper reports this values 0.72;1.49;1.49)
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
    plot3(meas(:,1),meas(:,2),meas(:,3),'k*');
    view(3);
elseif dimensions == 4                % CHANGED
    plot(meas(:,1),meas(:,2),'k*');
end

% PLOT STUFF .. SETTING UP AXIS IN THE FIGURE
axis equal;
axis(reshape([min(meas(:,[1 2]))-2; max(meas(:,[1 2]))+2],1,[]));    % CHANGED
hold off;

% SETTING UP PSO DATA STRUCTURES
swarm_vel = rand(centroids,dimensions,particles)*0.1;
swarm_pos = rand(centroids,dimensions,particles);
swarm_best = zeros(centroids,dimensions);
c = zeros(dataset_size(1),particles);
ranges = max(meas)-min(meas); %%scale
swarm_pos = swarm_pos .* repmat(ranges,centroids,1,particles) + repmat(min(meas),centroids,1,particles);
swarm_fitness(1:particles)=Inf;

% KMEANS_INIT
if hybrid_pso
    swarm_pos(:,:,1) = KMEANS_CENTROIDS;
else
    disp('random init');
    disp(swarm_pos(:,:,1));
end

% MANUAL INITIALIZATION (only for dimension 2 and 3)
if manual_init
    if dimensions == 3
        % MANUAL INIT ONLY FOR THE FIRST PARTICLE
             swarm_pos(:,:,1) = [6 3 4; 5 3 1]; 
    elseif dimensions == 2
        % KEYBOARD INIT ONLY FOR THE FIRST PARTICLE
             swarm_pos(:,:,1) = ginput(2);
    end
end

for iteration=1:iterations
      
    %CALCULATE EUCLIDEAN DISTANCES TO ALL CENTROIDS
    distances=zeros(dataset_size(1),centroids,particles);
    for particle=1:particles
        for centroid=1:centroids
            distance=zeros(dataset_size(1),1);
            for data_vector=1:dataset_size(1)
                %meas(data_vector,:)
                distance(data_vector,1)=norm(swarm_pos(centroid,:,particle)-meas(data_vector,:));
            end
            distances(:,centroid,particle)=distance;
        end
    end
    
    %ASSIGN MEASURES with CLUSTERS    
    for particle=1:particles
        [value, index] = min(distances(:,:,particle),[],2);
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
                elseif dimensions == 4          % CHANGED
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
    [global_fitness, index] = min(swarm_fitness);       %GLOBAL BEST FITNESS
    swarm_overall_pose = swarm_pos(:,:,index);          %GLOBAL BEST POSITION
    
    % SOME INFO ON THE COMMAND WINDOW
    fprintf('%3d. global fitness is %5.4f\n',iteration,global_fitness);            
    %uicontrol('Style','text','Position',[40 20 180 20],'String',sprintf('Actual fitness is: %5.4f', global_fitness),'BackgroundColor',get(gcf,'Color'));        
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


% PLOT THE ASSOCIATIONS WITH RESPECT TO THE CLUSTER 
hold on;
particle=index; %select the best particle (with best fitness) 
cluster_colors = ['m','g','y','b','r','c','g'];
for centroid=1:centroids
    if any(c(:,particle) == centroid)
        if dimensions == 3
            plot3(meas(c(:,particle)==centroid,1),meas(c(:,particle)==centroid,2),meas(c(:,particle)==centroid,3),'o','color',cluster_colors(centroid));
        elseif dimensions == 4              % CHANGED
            plot(meas(c(:,particle)==centroid,1),meas(c(:,particle)==centroid,2),'o','color',cluster_colors(centroid));
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

% SAY GOODBYE
fprintf('\nEnd, global fitness is %5.4f\n',global_fitness);