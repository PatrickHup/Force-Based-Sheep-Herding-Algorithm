clear; clc;

radius = 60;
use_advanced_vision = false;

[fiedler, distance_to_goal, sheep_trj, dog_trj, goal_pos, K] = simulate(radius, use_advanced_vision);


% Save trajectoy and paramaters for visualization
writematrix(sheep_trj, 'sheep_trj.csv');
writematrix(dog_trj, 'dog_trj.csv');
writematrix(goal_pos, 'goal_pos.csv');
writematrix(K, 'K.csv');

plot_util.animate(goal_pos, 10, sheep_trj, dog_trj, K);

function [fiedler, distance_to_goal, sheep_trj, dog_trj, goal_pos, K] = simulate(radius, use_advanced_vision)
%   SIMULATE Runs the herding simulation with multiple sheepdogs and a sheep flock.
%
%   [fiedler, distance_to_goal, sheep_trj, dog_trj, goal_pos, K] = SIMULATE(radius)
%   simulates the motion of sheep and dogs, calculating graph connectivity and
%   movement trajectories over time.
%
%   Input:
%       radius - Interaction radius for vision-based behaviors.
%
%   Outputs:
%       fiedler          - Time series of the second smallest eigenvalue (connectivity).
%       distance_to_goal - Mean distance of sheep from the goal at each timestep.
%       sheep_trj        - 3D trajectory matrix of sheep positions (2xNxK).
%       dog_trj          - 3D trajectory matrix of dog positions (2xDxK).
%       goal_pos         - 2D goal position vector.
%       K                - Number of simulation steps.

    % --- Simulation Parameters ---
    N = 100;         % Number of sheep
    D = 4;           % Number of dogs
    K = 1000;        % Number of time steps
    dt = 0.05;       % Time step size

    % Sheep behavior gains
    K_s1 = 10;       % Sheep-sheep repulsion
    K_s2 = 0.5;      % Velocity alignment
    K_s3 = 2.0;      % Cohesion
    K_s4 = 5000;     % Sheep-dog repulsion

    % Dog behavior gains
    K_f1 = 10;       % Attraction to target sheep
    K_f2 = 200;      % Repulsion from target sheep
    K_f3 = 8;        % Repulsion from goal
    K_f4 = 3000;     % Repulsion from other dogs

    goal_pos = [0; 20];  % Target goal position

    % --- Initialize State ---
    sheep_pos = spawn_sheep(N, [0; 0], 5);
    sheep_vel = zeros(2, N);

    dog_pos = spawn_dogs(D, 15, 3);
    dog_vel = zeros(2, D);

    % Preallocate storage for results
    sheep_trj = zeros(2, N, K);
    dog_trj = zeros(2, D, K);
    fiedler_value_all = zeros(1, K);
    distance_to_goal = zeros(1, K);
    connectivity_matrix = zeros(N, N);

    % --- Simulation Loop ---
    for k = 1:K
        % Store current positions
        sheep_trj(:, :, k) = sheep_pos;
        dog_trj(:, :, k) = dog_pos;

        % Track average distance of sheep from goal
        distance_to_goal(k) = mean(vecnorm(sheep_pos - goal_pos, 2, 1));

        % --- Sheep Update ---
        sheep_vel_new = zeros(2, N);
        for n = 1:N
            [a_i, b_i, c_i, d_i, conn] = calc_u_i_components(n, sheep_pos, sheep_vel, dog_pos, radius, N, use_advanced_vision);
            connectivity_matrix(n, :) = conn;
            sheep_vel_new(:, n) = K_s1 * a_i + K_s2 * b_i + K_s3 * c_i + K_s4 * d_i;
        end

        % Graph Laplacian and Fiedler value
        degree = sum(connectivity_matrix, 2);
        D_inv_sqrt = diag(1 ./ sqrt(degree + 1e-6));  % Avoid division by zero
        L = eye(N) - D_inv_sqrt * connectivity_matrix * D_inv_sqrt;
        eigen_values = sort(eig(L));
        fiedler_value_all(k) = eigen_values(2);  % Fiedler value

        % Limit sheep speed and update position
        speeds = vecnorm(sheep_vel_new);
        exceeds = speeds > 5;
        sheep_vel_new(:, exceeds) = sheep_vel_new(:, exceeds) ./ speeds(exceeds) * 5;
        sheep_pos = sheep_pos + dt * sheep_vel_new;
        sheep_vel = sheep_vel_new;

        % --- Dog Update ---
        dog_vel_new = zeros(2, D);
        for d = 1:D
            [A_i, B_i, C_i, D_i] = calc_v_i_components(d, dog_pos, sheep_pos, sheep_vel, goal_pos, radius);
            dog_vel_new(:, d) = K_f1 * A_i + K_f2 * B_i + K_f3 * C_i + K_f4 * D_i;
        end

        % Limit dog speed and update position
        speeds = vecnorm(dog_vel_new);
        dog_vel_new(:, speeds > 10) = dog_vel_new(:, speeds > 10) ./ speeds(speeds > 10) * 10;
        dog_pos = dog_pos + dt * dog_vel_new;
        dog_vel = dog_vel_new;

        disp("iteration: " + k);
    end

    % Output Fiedler value over time
    fiedler = fiedler_value_all;
end

function multi_trial(num_experiments, vision_tag, radii)
    % multi_trial runs a set of simulations for various radii and vision modes.
    % 
    % Parameters:
    %   num_experiments - Number of experiments to run for each radius.
    %   vision_tag      - A string indicating the vision mode ("normal_Vision" or "advanced_Vision").
    %   radii           - Array of radius values to test.

    for R = radii    
        % Create a folder name based on the current radius and vision mode
        folder_name = sprintf('results_R%d_%s', R, vision_tag);
        
        % Create the folder if it does not already exist
        if ~exist(folder_name, 'dir')
            mkdir(folder_name);
        end
    
        % Run simulations for the specified number of experiments
        for exp_num = 1:num_experiments
            % Run the simulation with current radius and vision mode
            % NOTE: 'advanced_vision_mode' should be a variable or flag defined elsewhere
            [fiedler, distance_to_goal, ~, ~, ~, ~] = simulate(R, advanced_vision_mode);
    
            % Smooth the Fiedler values using a moving average window of 20 samples
            fiedler_smooth = smoothdata(fiedler, 'movmean', 20);
            
            % Normalize the smoothed Fiedler values by their maximum
            fiedler_normalized = fiedler_smooth / max(fiedler_smooth);
    
            % Prepare data matrix: iteration index, raw Fiedler values, and distance to goal
            K = length(fiedler_normalized);
            data = [(1:K)', fiedler(:), distance_to_goal(:)];
            
            % Convert data to a table with appropriate column names
            T = array2table(data, ...
                'VariableNames', {'Iteration', 'Fiedler', 'DistanceToGoal'});
    
            % Build the output file name including radius, experiment number, and vision tag
            filename = fullfile(folder_name, ...
                sprintf('fiedler_R%d_E%d_%s.csv', R, exp_num, vision_tag));
    
            % Save the table as a CSV file
            writetable(T, filename);
        end
    end
    
    % Indicate completion
    disp("finished");
end

function [A_i, B_i, C_i, D_i] = calc_v_i_components(i, dog_pos, sheep_pos, sheep_vel, goal_pos, radius)
%   CALC_V_I_COMPONENTS Computes velocity components for dog i based on nearby agents.
%
%   [A_i, B_i, C_i, D_i] = CALC_V_I_COMPONENTS(...)
%   returns four 2D vectors guiding dog i's behavior toward sheep and away from others.
%
%   Inputs:
%       i         - Index of the current dog.
%       dog_pos   - 2xD matrix of all dog positions.
%       sheep_pos - 2xN matrix of all sheep positions.
%       sheep_vel - 2xN matrix of all sheep velocities (unused here).
%       goal_pos  - 2x1 vector of goal position.
%       radius    - Scalar interaction radius.
%
%   Outputs:
%       A_i - Attraction to target sheep (farthest from goal).
%       B_i - Repulsion from that target sheep (inverse-cube).
%       C_i - Repulsion from the goal (to avoid crowding it).
%       D_i - Repulsion from other dogs (inverse-cube).

    dog_pos_i = dog_pos(:, i);  % Position of current dog
    dog_pos_j = dog_pos(:, [1:i-1, i+1:end]);  % Positions of other dogs

    % Find nearby sheep in interaction radius
    [sheep_pos_in_range, ~, N_s] = get_sheep_in_range(dog_pos_i, sheep_pos, sheep_vel, radius);

    if N_s == 0
        % No sheep in range: minimal movement strategy
        A_i = zeros(2, 1);  % No target to approach
        B_i = zeros(2, 1);  % No repulsion
        diff_goal = goal_pos - dog_pos_i;
        dist_goal = norm(diff_goal);
        C_i = 0.1 * diff_goal / (dist_goal + 1e-6);  % Weak attraction to stay near goal
        D_i = zeros(2, 1);  % No dog-dog repulsion
        return;
    end

    % --- Target the farthest sheep within range ---
    target_sheep = get_target_sheep(goal_pos, sheep_pos_in_range);
    diff = dog_pos_i - target_sheep;
    dist = norm(diff);

    % Component A: Move toward target sheep (unit vector)
    A_i = -diff / (dist + 1e-6);

    % Component B: Repel from target sheep (inverse-cube law)
    B_i = diff / (dist + 1e-6)^3;

    % Component C: Repel from goal to prevent clustering
    diff_goal = dog_pos_i - goal_pos;
    dist_goal = norm(diff_goal);
    C_i = diff_goal / (dist_goal + 1e-6);

    % --- Dog-dog repulsion ---
    [dog_pos_in_range, N_d] = get_dog_in_range(dog_pos_i, dog_pos_j, radius);
    if N_d == 0
        D_i = zeros(2, 1);  % No nearby dogs to repel from
    else
        diff = dog_pos_i - dog_pos_in_range;
        dist = vecnorm(diff);
        D_i = 1 / N_d * sum(diff ./ dist.^3, 2);  % Inverse-cube repulsion from each dog
    end
end

function farthest_sheep_pos = get_target_sheep(goal_pos, sheep_pos)
%   GET_TARGET_SHEEP Finds the sheep farthest from the goal.
%
%   farthest_sheep_pos = GET_TARGET_SHEEP(goal_pos, sheep_pos)
%   returns the position of the sheep that is farthest away from
%   the goal position.
%
%   Inputs:
%       goal_pos  - 2x1 vector specifying the goal location.
%       sheep_pos - 2xN matrix of sheep positions (each column is a sheep).
%
%   Output:
%       farthest_sheep_pos - 2x1 position vector of the farthest sheep.

    % Compute position differences between each sheep and the goal
    diff = sheep_pos - goal_pos;

    % Compute Euclidean distances for each sheep
    dist = vecnorm(diff);

    % Find the index of the maximum distance
    [~, idx] = max(dist);

    % Return the position of the farthest sheep
    farthest_sheep_pos = sheep_pos(:, idx);
end

function [a_i, b_i, c_i, d_i, connectivity_] = calc_u_i_components(i, sheep_pos, sheep_vel, dog_pos, radius, N, use_advanced_vision)
%   CALC_U_I_COMPONENTS Calculates the components of sheep i's movement.
%
%   [a_i, b_i, c_i, d_i, connectivity_] = CALC_U_I_COMPONENTS(...)
%   computes the repulsion, alignment, cohesion, and dog-avoidance vectors
%   for a single sheep agent, based on visible neighbors and nearby dogs.
%
%   Inputs:
%       i          - Index of the current sheep.
%       sheep_pos  - 2xN matrix of all sheep positions.
%       sheep_vel  - 2xN matrix of all sheep velocities.
%       dog_pos    - 2xD matrix of all dog positions.
%       radius     - Vision radius for neighbor detection.
%       N          - Total number of sheep.
%
%   Outputs:
%       a_i         - Repulsion vector from nearby sheep.
%       b_i         - Velocity alignment vector.
%       c_i         - Cohesion (attraction) vector to flock center.
%       d_i         - Repulsion vector from nearby dogs.
%       connectivity_ - 1xN vector of connectivity weights to other sheep.

    sheep_pos_i = sheep_pos(:, i);  % Position of sheep i

    % Get visible neighboring sheep within radius, accounting for occlusion
    if use_advanced_vision
        [pos_in_range, vel_in_range, N_s, index] = get_sheep_in_range_blocked(sheep_pos_i, sheep_pos, sheep_vel, radius);
    else
        [pos_in_range, vel_in_range, N_s, index] = get_sheep_in_range(sheep_pos_i, sheep_pos, sheep_vel, radius);
    end

    % Compute distance-based weights for graph connectivity
    dist = vecnorm(sheep_pos_i - pos_in_range);
    weights = 1 ./ (dist + 1e-6);  % Add epsilon to avoid division by zero

    % Store weights in connectivity vector (sparse, only for visible neighbors)
    connectivity_ = zeros(1, N);
    connectivity_(index) = weights;

    % --- Calculate sheep-sheep interaction terms ---
    if N_s == 0
        % No neighbors: zero contribution from flocking behavior
        a_i = zeros(2,1);  % Repulsion
        b_i = zeros(2,1);  % Velocity alignment
        c_i = zeros(2,1);  % Cohesion
    else
        % Relative vectors to neighbors
        diff = sheep_pos_i - pos_in_range;
        dist = vecnorm(diff);

        % Repulsion (stronger from closer neighbors)
        a_i = 1/N_s * sum(diff ./ dist.^2, 2);

        % Cohesion (pulls toward center of mass)
        c_i = -1/N_s * sum(diff ./ dist, 2);

        % Alignment (average normalized velocity of moving neighbors)
        speeds = vecnorm(vel_in_range);
        valid = speeds > 0;
        if any(valid)
            b_i = 1/sum(valid) * sum(vel_in_range(:, valid) ./ speeds(valid), 2);
        else
            b_i = zeros(2,1);
        end
    end

    % --- Calculate dog avoidance term ---
    [dog_in_range, N_d] = get_dog_in_range(sheep_pos_i, dog_pos, radius);
    if N_d == 0
        d_i = zeros(2,1);  % No dogs nearby
    else
        diff = sheep_pos_i - dog_in_range;
        dist = vecnorm(diff);
        d_i = 1/N_d * sum(diff ./ dist.^3, 2);  % Stronger repulsion from closer dogs
    end
end

function [pos_in_range, vel_in_range, N_s, index] = get_sheep_in_range(pos_i, pos_j, vel_j, radius)
%   GET_SHEEP_IN_RANGE Finds all sheep within a certain radius (excluding self).
%
%   [pos_in_range, vel_in_range, N_s, index] = GET_SHEEP_IN_RANGE(pos_i, pos_j, vel_j, radius)
%   identifies nearby sheep to a given position `pos_i`, filtering by Euclidean
%   distance and excluding the self position.
%
%   Inputs:
%       pos_i  - 2x1 vector representing the position of the querying sheep.
%       pos_j  - 2xN matrix of all sheep positions.
%       vel_j  - 2xN matrix of all sheep velocities.
%       radius - Scalar defining the neighborhood interaction radius.
%
%   Outputs:
%       pos_in_range - 2xM matrix of nearby sheep positions.
%       vel_in_range - 2xM matrix of corresponding sheep velocities.
%       N_s          - Number of sheep found within range.
%       index        - 1xN logical vector indicating which sheep are in range.

    % Compute difference vectors and corresponding Euclidean distances
    diff = pos_i - pos_j;         % 2xN vector differences from pos_i
    dist = vecnorm(diff);         % 1xN vector of Euclidean distances

    % Identify neighbors that are within the radius and not at zero distance (exclude self)
    in_range = (dist > 0) & (dist < radius);

    % Extract positions and velocities of nearby sheep
    pos_in_range = pos_j(:, in_range);
    vel_in_range = vel_j(:, in_range);

    % Count how many neighbors are in range
    N_s = size(pos_in_range, 2);

    % Return logical index for reference
    index = in_range;
end

function [pos_in_range, vel_in_range, N_s, index] = get_sheep_in_range_blocked(pos_i, pos_j, vel_j, radius)
%   GET_SHEEP_IN_RANGE_BLOCKED Returns visible neighbors within radius, excluding blocked ones.
%
%   [pos_in_range, vel_in_range, N_s, index] = GET_SHEEP_IN_RANGE_BLOCKED(pos_i, pos_j, vel_j, radius)
%   determines which sheep are both within a certain radius and not occluded
%   (vision-blocked) by other sheep from the perspective of the reference agent.
%
%   Inputs:
%       pos_i  - 2x1 vector, position of the querying sheep.
%       pos_j  - 2xN matrix, positions of all sheep (columns).
%       vel_j  - 2xN matrix, velocities of all sheep (columns).
%       radius - Scalar interaction radius.
%
%   Outputs:
%       pos_in_range - 2xM matrix of positions of visible sheep within range.
%       vel_in_range - 2xM matrix of corresponding velocities.
%       N_s          - Number of visible sheep in range.
%       index        - 1xN logical vector indicating which sheep are visible.

    N = size(pos_j, 2);
    index = false(1, N);  % Logical mask for visible neighbors

    for j = 1:N
        if all(pos_j(:, j) == pos_i)
            continue;  % Skip self
        end

        % Compute distance to candidate neighbor
        diff = pos_j(:, j) - pos_i;
        dist = norm(diff);
        if dist >= radius
            continue;  % Outside interaction range
        end

        % --- Check for vision occlusion ---
        is_blocked = false;
        for b = 1:N
            if b == j || all(pos_j(:, b) == pos_i)
                continue;  % Skip current neighbor and self
            end

            blocker = pos_j(:, b);

            % Project blocker onto line between pos_i and neighbor
            t = dot(blocker - pos_i, diff) / (norm(diff)^2);
            t = min(max(t, 0), 1);  % Clamp projection to line segment
            proj = pos_i + t * diff;

            % Compute perpendicular distance from blocker to projection
            d_block = norm(proj - blocker);

            % Vision block if too close to line and in front of neighbor
            block_threshold = 1.0;
            if d_block < block_threshold && norm(blocker - pos_i) < dist
                is_blocked = true;
                break;
            end
        end

        if ~is_blocked
            index(j) = true;
        end
    end

    % Filter positions and velocities by visibility mask
    pos_in_range = pos_j(:, index);
    vel_in_range = vel_j(:, index);
    N_s = sum(index);
end

function [dog_in_range, N_d] = get_dog_in_range(pos_i, pos_j, radius)
%   GET_DOG_IN_RANGE Identifies neighboring dogs within a given radius.
%
%   [dog_in_range, N_d] = GET_DOG_IN_RANGE(pos_i, pos_j, radius)
%   returns the positions of dogs that are within the specified radius
%   of the reference dog located at pos_i.
%
%   Inputs:
%       pos_i  - 2x1 position vector of the reference dog.
%       pos_j  - 2xM matrix of positions of all other dogs (each column is a dog).
%       radius - Scalar radius defining the neighborhood.
%
%   Outputs:
%       dog_in_range - 2xK matrix of dog positions within range.
%       N_d          - Number of neighboring dogs found.

    % Compute distance from pos_i to each dog in pos_j
    diff = pos_i - pos_j;               % 2xM difference vectors
    dist = vecnorm(diff);              % 1xM Euclidean distances

    % Logical index of dogs within range but not at zero distance (exclude self)
    in_range = (dist > 0) & (dist < radius);

    % Extract positions of neighboring dogs
    dog_in_range = pos_j(:, in_range);

    % Count number of neighbors
    N_d = size(dog_in_range, 2);
end

function dog_pos = spawn_dogs(D, radius, spacing)
%   SPAWN_DOGS Generates initial positions for D dog agents aligned on a tangent line.
%
%   dog_pos = SPAWN_DOGS(D, radius, spacing) returns a 2xD matrix where
%   each column represents the (x, y) position of a dog. The dogs are placed
%   along a line tangent to a circle centered at the origin, at a distance 
%   'radius' from the origin, and spaced evenly along that tangent.
%
%   Inputs:
%       D       - Number of dog agents.
%       radius  - Distance from the origin to the center of the tangent line.
%       spacing - Distance between adjacent dogs along the tangent.
%
%   Output:
%       dog_pos - 2xD matrix of dog positions.

    % Random angle within a lower semicircle (between 225° and 315°)
    theta_min = 5 * pi / 4;
    theta_max = 7 * pi / 4;
    theta = theta_min + rand() * (theta_max - theta_min);

    % Compute center of tangent line on circle at angle theta
    center_pos = radius * [cos(theta); sin(theta)];

    % Compute unit tangent vector at angle theta
    tangent = [-sin(theta); cos(theta)];

    % Compute symmetric offsets along tangent line centered at center_pos
    offsets = ((0:D-1) - (D-1)/2) * spacing;
    
    % Generate dog positions along the tangent line
    dog_pos = center_pos + tangent * offsets;
end

function sheep_pos = spawn_sheep(N, center, radius)
%   SPAWN_SHEEP Generates initial positions for N sheep within a circular region.
%
%   sheep_pos = SPAWN_SHEEP(N, center, radius) returns a 2xN matrix where
%   each column is the (x, y) position of a sheep randomly placed within a
%   circle of given radius centered at the specified location.
%
%   Inputs:
%       N      - Number of sheep to generate.
%       center - 2x1 vector specifying the center of the circular region [x; y].
%       radius - Radius of the circle within which to distribute the sheep.
%
%   Output:
%       sheep_pos - 2xN matrix of sheep positions.

    % Generate random angles and distances (polar coordinates)
    theta = 2 * pi * rand(1, N);         % Uniform angle from 0 to 2π
    rho = radius * sqrt(rand(1, N));     % Radial distance adjusted for uniform area

    % Convert polar coordinates to Cartesian and shift by the center
    sheep_pos = center + [rho .* cos(theta); rho .* sin(theta)];
end

