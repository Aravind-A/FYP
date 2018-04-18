clear;
close all;
%in pso each particle is an possible solution to the give set of objective
%function; pso finds the best particle after algo completes depending on
%the global and local best and also the individual particle's postion and
%velocity with also influence other particle's parameters
%rng('default') % For reproducibility

variable=csvread('f1.dat.txt'); % dataset with the dimension values
sp=csvread('f2.dat.txt');% file that contains the actual value of centroids each data point belong to
% INIT PARTICLE SWARM
dataset_size = size (variable);
centroids = max(sp) ;          % == clusters here (aka centroids) change at line 77 if any od centroid or dimensions changes values in initial inis
dimensions = dataset_size(2) ;         % how many dimensions in each centroid
particles = 20;         % how many particles in the swarm, aka how many solutions
iterations = 20;        % iterations of the optimization alg.
simtime=0.01;           % simulation delay btw each iteration
dataset_subset = 2;     % for the IRIS dataset, change this value from 0 to 2
write_video = false;    % enable to grab the output picture and save a video
hybrid_pso = true;     % enable/disable hybrid_pso
manual_init = false;    % enable/disable manual initialization (only for dimensions={2,3})


% LOAD DEFAULT CLUSTER (IRIS DATASET); USE WITH CARE!
%variable = variable(:,1+dataset_subset:dimensions+dataset_subset); %RESIZE THE DATASET WITH CURRENT DIMENSIONS; USE WITH CARE!

% EXECUTE K-MEANS
if hybrid_pso
    fprintf('Running Matlab K-Means Version\n');
    [idx,KMEANS_CENTROIDS] = kmeans(variable,centroids, 'dist','sqEuclidean', 'display','iter','start','uniform','onlinephase','off');
    %idx - 150*1 value of the centroid(1,2,3) each datapoint belongs to
    %and k means centroids have coordinates of each centroid
    fprintf('\n');
end

% GLOBAL PARAMETERS (the paper reports this values 0.72;1.49;1.49)
w  = 0.72; %INERTIA
c1 = 1.49; %COGNITIVE
c2 = 1.49; %SOCIAL

% PLOT STUFF... HANDLERS AND COLORS
pc = []; txt = [];
cluster_colors_vector = rand(particles, 3);

% PLOT DATASET all the data points for k  means clustering
fh=figure(1);
hold on;
if dimensions >=3
    
    plot(variable(:,1),variable(:,2),'k*');

    
end

% PLOT STUFF .. SETTING UP AXIS IN THE FIGURE
axis equal;
%axis(reshape([min(variable)-2; max(variable)+2],1,[]));
hold off;

% SETTING UP PSO DATA STRUCTURES
swarm_vel = rand(centroids,dimensions,particles)*0.1;% for every particle, its postion and velocity for each dimesion
swarm_pos = rand(centroids,dimensions,particles);% and for each centroid is maintained
swarm_best = zeros(centroids,dimensions);% refers to the values for the best particle for each centroid and each dimesion
c = zeros(dataset_size(1),particles);% contains for each data point and for each particle to what centroid it points to
% from this after the best particle is found out that 150*pth best particle
% is taken and each data point is allocated to that particular centroid
ranges = max(variable)-min(variable); %%scale
%swarm_pos = swarm_pos .* repmat(ranges,centroids,1,particles) + repmat(min(variable),centroids,1,particles); % adding the postion to minimum postition for scaling purposes
swarm_fitness(1:particles)=Inf;%individual fitness of each particle

% KMEANS_INIT
if hybrid_pso
    swarm_pos(:,:,1) = KMEANS_CENTROIDS;
end

% MANUAL INITIALIZATION (only for dimension 2 and 3)
if manual_init
    if dimensions == 3
        % MANUAL INIT ONLY FOR THE FIRST PARTICLE
             %%swarm_pos(:,:,1) = [6 3 4; 5 3 1];
    elseif dimensions == 4
        % MANUAL INIT ONLY FOR THE FIRST PARTICLE
             %swarm_pos(:,:,1) = [2 3 4 6; 8 6 3 4; 4 5 3 1]; % for 3 clusters and 4 dimensions
    elseif dimensions == 2
        % KEYBOARD INIT ONLY FOR THE FIRST PARTICLE
             swarm_pos(:,:,1) = ginput(2);
    end
end

for iteration=1:iterations
      
    %CALCULATE EUCLIDEAN DISTANCES TO ALL CENTROIDS for 1 particle to each
    %of data points and centroids
    distances=zeros(dataset_size(1),centroids,particles);% datapoints*centroids*particles
    for particle=1:particles % for each particle 
        for centroid=1:centroids% for each centroid
            distance=zeros(dataset_size(1),1);
            for data_vector=1:dataset_size(1)% for each data point
                %variable(data_vector,:)
                distance(data_vector,1)=norm(swarm_pos(centroid,:,particle)-variable(data_vector,:));% calc distance between data point and particle
            end
            distances(:,centroid,particle)=distance;
        end
    end
    
    %ASSIGN variableURES with CLUSTERS    
    for particle=1:particles
        [value, index] = min(distances(:,:,particle),[],2);
        c(:,particle) = index;
        % value is 150 * 1 matrix and contains the minimun of the  distance between each data point and all the three centroids, 
        %effectively finding the nearest centroid for each data point and particle and store the centroid in c[datapoint,particle];     
    end

    % PLOT STUFF... CLEAR HANDLERS
    delete(pc); delete(txt);
    pc = []; txt = [];
    
    % PLOT STUFF...
    hold on;
    %plot the position of each particle on every iteration with
    %indicationg the cluster they(only particles) belong with the respective color
    for particle=1:particles
        for centroid=1:centroids
            if any(c(:,particle) == centroid)
                %if dimensions >= 3 
                 %   pc = [pc plot3(swarm_pos(centroid,1,particle),swarm_pos(centroid,2,particle),swarm_pos(centroid,3,particle),'*','color',cluster_colors_vector(particle,:))];
                if dimensions >= 2
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
                local_fitness=mean(distances(c(:,particle)==centroid,centroid,particle));% means of 
                average_fitness(particle,1) = average_fitness(particle,1) + local_fitness; % avg fitness for every particle
           %on every loop centroid for each datapoint,particle can change
            end
        end
        average_fitness(particle,1) = average_fitness(particle,1) / centroids;
        if (average_fitness(particle,1) < swarm_fitness(particle))
            swarm_fitness(particle) = average_fitness(particle,1);
            swarm_best(:,:,particle) = swarm_pos(:,:,particle);     %LOCAL BEST FITNESS of every particle and every data point
        end
    end    % calc local fitness and swarm fitness for all particles ends
    [global_fitness, index] = min(swarm_fitness); %get the index of the minimum value and store it in swarm overall pos      %GLOBAL BEST FITNESS
    swarm_overall_pose = swarm_pos(:,:,index);          %GLOBAL BEST POSITION
    
    % SOME INFO ON THE COMMAND WINDOW
    fprintf('%3d. global fitness is %5.4f\n',iteration,global_fitness);            
    %uicontrol('Style','text','Position',[40 20 180 20],'String',sprintf('Actual fitness is: %5.4f', global_fitness),'BackgroundColor',get(gcf,'Color'));        
    pause(simtime);
        
   % now comes pso algo ; change the velocity and position of each particle
   % as defined by the algo
       
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
particle=index; %select the best particle (with best fitness) selected at the end of last loop
cluster_colors = ['m','g','y','b','r','c','g'];
for centroid=1:centroids
    if any(c(:,particle) == centroid)
        if dimensions >= 2
            plot(variable(c(:,particle)==centroid,1),variable(c(:,particle)==centroid,2),'O','color',cluster_colors(centroid));
        end
    end
end
    
    

    
 
  
        
        
    cnt=hist(c(:,particle),centroids);
    fprintf('%d\n',cnt);
    finalcount=zeros(centroids,centroids);
    for i=1:dataset_size(1)
        finalcount(c(i,particle),sp(i))=1+finalcount(c(i,particle),sp(i));   
                
   end
    

hold off;



% SAY GOODBYE
fprintf('\nEnd, global fitness is %5.4f\n',global_fitness);