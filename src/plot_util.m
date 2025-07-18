classdef plot_util
    methods(Static)
        function animate(goal_pos, goal_rad, sheep_trj, dog_trj, K)
        %   ANIMATE Visualizes the herding simulation over K time steps.
        %
        %   animate(goal_pos, goal_rad, sheep_trj, dog_trj, K) creates a 2D
        %   animation of sheep and dog trajectories guiding toward a circular goal.
        %
        %   Inputs:
        %       goal_pos  - 2x1 vector specifying the goal center position.
        %       goal_rad  - Scalar radius of the goal area.
        %       sheep_trj - 2xN xK matrix of sheep trajectories.
        %       dog_trj   - 2xD xK matrix of dog trajectories.
        %       K         - Total number of simulation time steps.
        
            % --- Set up figure and plot window ---
            close(gcf);        % Close any existing figure window
            figure; hold on;   % Open a new figure and hold all plots
        
            % Compute overall bounds for axes from all positions
            all_pos = cat(2, reshape(sheep_trj, 2, []), reshape(dog_trj, 2, []));
            pos_min = min(all_pos, [], 2);
            pos_max = max(all_pos, [], 2);
            padding = 5;
        
            % --- Initial plot setup ---
            h_sheep = plot(sheep_trj(1,:,1), sheep_trj(2,:,1), 'ko', 'MarkerSize', 3);  % Sheep markers
            h_dog   = plot(dog_trj(1,:,1),   dog_trj(2,:,1),   'rs', 'MarkerSize', 5);  % Dog markers
        
            % Plot goal as a circle
            plot_util.plot_circle(goal_pos, goal_rad, 'k-');  % Custom utility to draw goal area
        
            axis equal;
            xlim([pos_min(1) - padding, pos_max(1) + padding]);
            ylim([pos_min(2) - padding, pos_max(2) + padding]);
        
            % --- Animation loop ---
            for k = 1:K
                % Update sheep and dog positions at step k
                set(h_sheep, 'XData', sheep_trj(1,:,k), 'YData', sheep_trj(2,:,k));
                set(h_dog,   'XData', dog_trj(1,:,k),   'YData', dog_trj(2,:,k));
        
                title(['Step ', num2str(k)]);
                drawnow;        % Force MATLAB to update the plot
                pause(0.01);    % Pause for smooth animation
            end
        end
        
        function h = plot_circle(center, radius, varargin)
            % center: [x; y] or [x, y]
            % radius: scalar
            % varargin: additional plot options, e.g., 'r--' for style
        
            theta = linspace(0, 2*pi, 100);
            x = center(1) + radius * cos(theta);
            y = center(2) + radius * sin(theta);
        
            h = plot(x, y, varargin{:});
        end
    end
end