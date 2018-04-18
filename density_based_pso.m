clear;
close all;

centroids = 10;
dimensions = 3706;
particles = centroids;
iterations = 100;
simtime = 0.01;
alpha = 0.01;
vmax = 5;
hybrid_pso = false;

%load fisheriris.mat
meas = csvread('OP1.txt');
dataset_size = size(meas);
O = dataset_size(1);
C = centroids;
sigma = 1/nthroot(O, dimensions);

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
elseif dimensions > 3                % CHANGED
    plot(meas(:,1),meas(:,2),'k*');
end

% PLOT STUFF .. SETTING UP AXIS IN THE FIGURE
axis equal;
axis(reshape([min(meas(:,[1 2]))-2; max(meas(:,[1 2]))+2],1,[]));    % CHANGED
hold off;

swarm_pos = rand(particles,dimensions);
swarm_vel = rand(particles,dimensions)*0.1;
swarm_best = zeros(centroids,dimensions);
pBest = zeros(particles, dimensions);
pDense = zeros(particles, dimensions);
grav_coeff = zeros(particles, dimensions);
p_dense_values = zeros(particles,1);
dp_particle_mapping = zeros(O,1);
swarm_fitness(1:particles)=Inf;

p_one_index = ceil(rand*O);
swarm_pos(1,:) = meas(p_one_index,:);
prev_index = p_one_index;

for particle = 2:particles
    cur_index = mod(floor(prev_index + (O/C) - 1), O) + 1;
    swarm_pos(particle,:) = meas(cur_index,:);
    prev_index = cur_index;
end


% KMEANS_INIT
if hybrid_pso
    swarm_pos(1,:) = KMEANS_CENTROIDS;
end

for iteration = 1:iterations
    
    distances=zeros(O,particles);

    for particle=1:particles
        distance=zeros(O,1);
        for data_vector=1:O
            distance(data_vector,1)=norm(swarm_pos(particle,:)-meas(data_vector,:));
        end
        [val, index] = min(distance);
        pBest(particle,:) = meas(index,:);
        distances(:,particle)=distance;
    end
    
    % PLOT STUFF... CLEAR HANDLERS
    delete(pc); delete(txt);
    pc = []; txt = [];
    
    % PLOT STUFF...
    hold on;
    for particle=1:particles
        pc = [pc plot(swarm_pos(particle,1),swarm_pos(particle,2),'*','color',cluster_colors_vector(particle,:))];  
    end
    set(pc,{'MarkerSize'},{12})
    hold off;
    
    for i=1:O
        [~, index] = min(distances(i,:));
        dp_particle_mapping(i,1) = index;
    end
    
    for particle=1:particles
        if any(dp_particle_mapping == particle)
            dps = meas(dp_particle_mapping == particle,:);
            density_values = zeros(size(dps,1),1);
            coeff_particle = 0;
            for i=1:size(dps,1)
                temp_sum = 0;
                for j=1:size(dps,1)
                    dist = norm(dps(i,:) - dps(j,:));
                    exponent = -1 * ((dist^2) / (2 * (sigma^2)));
                    val = exp(exponent) / ((2 * pi)^0.5);
                    temp_sum = temp_sum + val;
                end
                density_values(i,1) = temp_sum / (O * sigma);
                temp_coeff = (rand * (dps(i,:) - swarm_pos(particle,:))) / (norm(dps(i,:) - swarm_pos(particle,:)) + alpha);
                coeff_particle = coeff_particle + temp_coeff;
            end
            grav_coeff(particle,:) = coeff_particle;
            if max(density_values) > p_dense_values(particle,1)
                [value, index] = max(density_values);
                pDense(particle,:) = dps(index,:);
                p_dense_values(particle,1) = value;
            end
        end
    end
    
    pause(simtime);
    %temp_vel = (w*swarm_vel) + (rand * grav_coeff.*(pBest - swarm_pos)) + (rand * grav_coeff.*(pDense - swarm_pos));
    temp_vel = (w*swarm_vel) + (rand * (pBest - swarm_pos)) + (rand * (pDense - swarm_pos));
    if any(temp_vel > vmax)
        temp_vel(temp_vel > vmax) = vmax;
    end
    swarm_vel = temp_vel;
    swarm_pos = swarm_pos + swarm_vel;
end

% PLOT THE ASSOCIATIONS WITH RESPECT TO THE CLUSTER 
hold on; 
cluster_colors = ['m','g','y','b','r','c','g','m','r','y'];
for particle=1:particles
    plot(meas(dp_particle_mapping == particle,1),meas(dp_particle_mapping == particle,2),'o','color',cluster_colors(particle));
end
hold off;