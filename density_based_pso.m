function out = density_based_pso(meas, truth_labels, centroids, iterations)

particles = centroids;
simtime = 0.01;
alpha = 0.01;
vmax = 1;
dataset_size = size(meas);
O = dataset_size(1);
dimensions = dataset_size(2);
C = centroids;
w  = 0.72; %INERTIA

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
pBest = zeros(particles, dimensions);
pDense = zeros(particles, dimensions);
grav_coeff = zeros(particles, dimensions);
min_dist = zeros(particles,1);
p_dense_values = zeros(particles,1);
dp_particle_mapping = zeros(O,1);

%p_one_index = ceil(rand*O/3);
p_one_index = 32;
swarm_pos(1,:) = meas(p_one_index,:);
prev_index = p_one_index;

for particle = 2:particles
    cur_index = mod(floor(prev_index + (O/C) - 1), O) + 1;
    swarm_pos(particle,:) = meas(cur_index,:);
    prev_index = cur_index;
end

for iteration = 1:iterations
    
    distances=zeros(O,particles);

    for particle=1:particles
        distance=zeros(O,1);
        for data_vector=1:O
            distance(data_vector,1)=norm(swarm_pos(particle,:)-meas(data_vector,:));
        end
        [val, index] = min(distance);
        if min_dist(particle,1) == 0 || min_dist(particle,1) > val
            min_dist(particle,1) = val;
            pBest(particle,:) = meas(index,:); 
        end
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
            dps_size = size(dps);
            dps_len = dps_size(1);
            dps_dist = Inf(dps_len,dps_len);
            for i=1:dps_len
                for j=i+1:dps_len
                    dps_dist(i,j) = norm(dps(i,:)-dps(j,:));
                end
            end
            sigma = min(min(dps_dist));
            density_values = zeros(size(dps,1),1);
            coeff_particle = zeros(1,dimensions);
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
    temp_vel = (w*swarm_vel) + ((pBest - swarm_pos)) + ((pDense - swarm_pos));
    if any(temp_vel > vmax)
        temp_vel(temp_vel > vmax) = vmax;
    end
    swarm_vel = temp_vel;
    swarm_pos = swarm_pos + swarm_vel;
end

ari = clustereval(dp_particle_mapping, truth_labels, 'ari');
fprintf('ARI :  ');
disp(ari);

% PLOT THE ASSOCIATIONS WITH RESPECT TO THE CLUSTER 
hold on; 
cluster_colors = ['m','g','y','b','r','c','g','m','r','y'];
for particle=1:particles
    plot(meas(dp_particle_mapping == particle,1),meas(dp_particle_mapping == particle,2),'o','color',cluster_colors(particle));
end
hold off;
out=1;
end