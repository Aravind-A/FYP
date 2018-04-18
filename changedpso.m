clear;
close all;
%in pso each particle is an possible solution to the give set of objective
%function; pso finds the best particle after algo completes depending on
%the global and local best and also the individual particle's postion and
%velocity with also influence other particle's parameters
%rng('default') % For reproducibility

% INIT PARTICLE SWARM
centroids = 3 ;          % == clusters here (aka centroids) 
dimensions = 4 ;         % how many dimensions in each centroid
particles = centroids;         % how many particles in the swarm, aka how many solutions
iterations = 30;        % iterations of the optimization alg.
simtime=0.01;           % simulation delay btw each iteration
dataset_subset = 2;     % for the IRIS dataset, change this value from 0 to 2
write_video = false;    % enable to grab the output picture and save a video
hybrid_pso = false;     % enable/disable hybrid_pso
manual_init = false;    % enable/disable manual initialization (only for dimensions={2,3})
alpha=0.001;
V_max=rand(dimensions);
V_min=rand(dimensions);

% LOAD DEFAULT CLUSTER (IRIS DATASET); USE WITH CARE!
load fisheriris.mat
%meas = meas(:,1+dataset_subset:dimensions+dataset_subset); %RESIZE THE DATASET WITH CURRENT DIMENSIONS; USE WITH CARE!
dataset_size = size (meas);
O=dataset_size(1);
% EXECUTE K-MEANS
if hybrid_pso
    fprintf('Running Matlab K-Means Version\n');
    [idx,KMEANS_CENTROIDS] = kmeans(meas,centroids, 'dist','sqEuclidean', 'display','iter','start','uniform','onlinephase','off');
    %idx - 150*1 value of the centroid(1,2,3) each datapoint belongs to
    %and k means centroids have coordinates of each centroid
    cnt2=hist(idx);
    fprintf('%d\n',cnt2);
end

% GLOBAL PARAMETERS (the paper reports this values 0.72;1.49;1.49)
w  = 0.72; %inertia

% PLOT STUFF... HANDLERS AND COLORS
pc = []; txt = [];
cluster_colors_vector = rand(particles, 3);

% PLOT DATASET all the data points for k  means clustering
fh=figure(1);
hold on;
if dimensions >=3
    
    plot(meas(:,1),meas(:,2),'k*');

    
end

% PLOT STUFF .. SETTING UP AXIS IN THE FIGURE
axis equal;
axis(reshape([min(meas)-2; max(meas)+2],1,[]));
hold off;

% SETTING UP PSO DATA STRUCTURES
% social and cognitive learning factor, the gravity force
force=zeros(dimesions,particles);
swarm_vel = rand(dimensions,particles)*0.1;% for every particle, its postion and velocity for each dimesion

swarm_pos = rand(dimensions,particles);% each particle takes care of one single cluster
swarm_bestpos = zeros(dimensions,particles)+Inf;
% refers to the values for the best particle for each centroid and each dimesion
swarm_bestden = zeros(dimensions,particles);
c = zeros(dataset_size(1),particles);% contains for each data point and for each particle to what centroid it points to
% from this after the best particle is found out that 150*pth best particle
% is taken and each data point is allocated to that particular centroid
ranges = max(meas)-min(meas); %%scale
%changing swarm pos
pone=ceil(rand*dataset_size(1));

swarm_pos(:,1) = meas(pone,:);  % 
for particle=2:particles
swarm_pos(:,particle) = floor(mod(swarm_pos(:,particle-1) +O/centroids -1,O))+1;
end

    %CALCULATE EUCLIDEAN DISTANCES TO ALL CENTROIDS for 1 particle to each
    %of data points and centroids
    c=zeros(dataset_size(1));
    distances=zeros(dataset_size(1),particles);% datapoints*centroids*particles
    for data_vector=1:dataset_size(1)
    for particle=1:particles
    distances(data_vector,particle)=norm(swarm_pos(:,particle)-meas(data_vector,:));
    end
    [value,index]=min(distances(data_vector));
    c(data_vector)=index;
    end
    
for iteration=1:iterations
    %plot postion and color  
    force=zeros(dimesions,particles);
    %CALCULATE EUCLIDEAN DISTANCES TO ALL particles 
    distances=zeros(dataset_size(1),particles);% 
    for particle=1:particles
    for data_vector=1:dataset_size(1)
    for particle1=1:particles
    distances(data_vector,particle1)=norm(swarm_pos(:,particle1)-meas(data_vector,:));
    end
    [~,index]=min(distances(data_vector));
    c(data_vector)=index;
    end
        
        sigma=1/nthroot(dataset_size(1),dimensions);
        % personal best position calculation yet to do
        % stored in swarm_bestpos(dimensions,particles)
       [value,index]=min(distances(:,particle));
       if(value<swarms)
       swarm_bestpos(:,particle)=meas(index);
       end
        
        %personal best density calc yet to do
        % stored in swarm_bestden(dimesions,particles)
        M_i=1;M_j=1;
        for data_vector=1:dataset_size(1)
        if c(data_vector)==particle
        force(:,particle)=force(:,particle)+ (rand* M_i * M_j * (meas(data_vector,:)-swarm_pos(:,particle)))/(distances(data_vector,particle)+alpha);
        end
        end
        r1 = rand;
        r2 = rand;
        swarm_vel(:,particle) = w*swarm_vel(:,particle);
        swarm_vel(:,particle) = swarm_vel(:,particle)+r1*force(:,particle)*(swarm_bestpos(:,particle)-swarm_pos(:,particle));
        swarm_vel(:,particle) = swarm_vel(:,particle)+r2*force(:,particle)*(swarm_bestden(:,particle)-swarm_pos(:,particle));
        if swarm_vel(:,particle) > V_max
            swarm_vel(:,particle) = V_max;
        elseif swarm_vel(:,particle) < V_min
            swarm_vel(:,particle) = V_min;
        end
     swarm_pos(:,particle)=swarm_pos(:,particle)+swarm_vel(:,particle);      
    end
    end
end
cnt=hist(c);
printf('%d',c);

