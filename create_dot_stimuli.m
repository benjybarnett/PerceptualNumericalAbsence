function create_dot_stimuli()
    %% This script creates dot pattern images ranging from 0 - 5 dots.
    % There are three different background colours, and the dots will be white and black,
    % alternating on each successive stimulus.  The total area covered by dot, the density
    % of dots and the weber contrast are held constant across numerosity.
    % You can edit the number of exemplars per numerosity, but this should be a multiple of 6
    % (6 = (3 * backgorund colours) * (2 * dot colours). This will give you an even number of 
    % dot colour/background colour pairings per each numerosity.
    % Each image has a random size, between a min and max bound specified somewhere in the script.
    % At the end, the function spits out a correlation matrix that shows you the relationship
    % between numerosity and total dot area, density, luminance, mean dot area, and convex hull area
    
    % Benjy Barnett, 2024
    
    mkdir('numerosity_0')
    mkdir('numerosity_1')
    mkdir('numerosity_2')
    mkdir('numerosity_3')
    mkdir('numerosity_4')
    mkdir('numerosity_5')

    
    % Number of stimuli (numerosities) from 0 to 5
    numerosities = 0:5;
    
    % Define 3 different background greys
    background_greys = [0.3, 0.5, 0.7];  % 3 different luminance levels (greyscale, between 0 and 1)
    
    % Total dot area (constant across numerosities) - Reduced to 800
    total_dot_area = 800;  % in pixels, adjust depending on display size
    
    % Dot spacing factor (to avoid overlap)
    min_dot_spacing = 1;  % Minimum distance between dots to prevent overlap
    
    % Set the largest frame size so we can reuse the same window
    max_frame_size = 250;
    
    % Number of exemplars to create per numerosity
    num_exemplars_per_numerosity = 120;

    % Initialize Psychtoolbox once
    PsychDebugWindowConfiguration([], 0.5)
    Screen('Preference', 'SkipSyncTests', 1);
    [win, rect] = Screen('OpenWindow', 0, [128 128 128], [0 0 max_frame_size max_frame_size], [], [], [], 6);
    ShowCursor;
    
    % Initialize data structure to store image metadata
    image_metadata = struct('numDots', [], 'totalDotArea', [], 'dotPositions', [], ...
                            'dotDensity', [], 'luminance', [], 'meanDotArea', [], ...
                            'convexHullArea', []);

    % Loop through each numerosity
    for numDots = numerosities
        
        exemplar_count = 0;
        
        % Loop to create num_exemplars_per_numerosity for each numerosity
        while exemplar_count < num_exemplars_per_numerosity
            
            % Randomly vary the screen size (frame size) between 150 and 250
            screen_size = randi([150, 250]);

            % Loop through each background grey and each color index
            for bgIdx = 1:length(background_greys)
                background_luminance = background_greys(bgIdx);  % Select current background grey
                
                for colorIdx = 1:2
                    % Increment the exemplar count and break the loop if we have reached max
                    exemplar_count = exemplar_count + 1;
                    if exemplar_count > num_exemplars_per_numerosity
                        break;
                    end
                    
                    % Set dot color: first three patterns white, next three black
                    if colorIdx == 1
                        dot_color = [255 255 255];  % White dots
                        color_code = 'w';
                    else
                        dot_color = [0 0 0];  % Black dots
                        color_code = 'b';
                    end

                    % Set background color for this grey
                    Screen('FillRect', win, background_luminance * 255, [0 0 screen_size screen_size]);

                    if numDots > 0
                        % Calculate the dot radius for the total dot area to be constant
                        dot_area = total_dot_area/numDots;
                        dot_radius = sqrt(dot_area / pi);  % For circular dots
                        
                        % Add some small random variation to dot size
                        a = 0.98;
                        b = 1.02;
                        factor = (b-a).*rand(1,1) + a; % Slight random variation in dot size
                        dot_radius = dot_radius * factor; % Control area covered by dots
                    
                        % Generate non-overlapping positions for the dots with edge padding
                        [dot_positions] = generate_non_overlapping_positions(numDots, dot_radius, min_dot_spacing, screen_size);

                        % Draw dots with anti-aliasing enabled (smooth edges)
                        Screen('DrawDots', win, dot_positions', dot_radius*2, dot_color, [], 2); % Anti-aliasing enabled with '2'
                    
                        % Store metadata information
                        dot_density = calculate_density(dot_positions); % Density based on median pairwise distance
                        dot_luminance = calculate_weber_contrast(background_luminance * 255, dot_color); % Luminance
                        
                        final_dot_area = (pi * (dot_radius^2)) * numDots;

                        % Mean area of an individual dot (all dots share the same radius per exemplar)
                        mean_dot_area = final_dot_area / numDots;

                        % Convex hull area spanned by the dot positions (needs >=3 dots to form a polygon)
                        hull_area = calculate_convex_hull_area(dot_positions);
                        
                        % Append the information to the image metadata structure
                        image_metadata(numel(image_metadata)+1) = struct('numDots', numDots, ...
                            'totalDotArea', final_dot_area, ...
                            'dotPositions', dot_positions, ...
                            'dotDensity', dot_density, ...
                            'luminance', dot_luminance, ...
                            'meanDotArea', mean_dot_area, ...
                            'convexHullArea', hull_area);
                    end

                    % Flip the screen to display the stimulus
                    Screen('Flip', win);
                    
                    % Save the image as a .png file, including numerosity, background index, dot color, and exemplar index
                    imageArray = Screen('GetImage', win, [0 0 screen_size screen_size]);
                    imwrite(imageArray, fullfile(strcat('numerosity_',num2str(numDots)),[ 'bg_' num2str(bgIdx) '_' color_code '_exemplar_' num2str(exemplar_count) '.png']));
                    
                    % Clear the screen for the next stimulus
                    Screen('FillRect', win, background_luminance * 255, [0 0 screen_size screen_size]);
                    
                end
            end
        end
    end
    
    % Close Psychtoolbox window only after all stimuli are created
    Screen('CloseAll');
    
    % Call the function to calculate the correlation matrix
    calculate_correlation_matrix(image_metadata);
    savefig('corr_matrix.fig')
end

function [dot_positions] = generate_non_overlapping_positions(numDots, dot_radius, min_spacing, screen_size)
    % Generate positions for dots such that they don't overlap and don't touch the edge
    dot_positions = zeros(numDots, 2);
    max_attempts = 100;  % Limit the number of attempts to find non-overlapping positions
    
    for i = 1:numDots
        placed = false;
        attempts = 0;
        
        while ~placed && attempts < max_attempts
            % Generate random position within the bounds (considering dot radius and padding from edges)
            x = randi([round(dot_radius*2), round(screen_size - (dot_radius*2))]);
            y = randi([round(dot_radius*2), round(screen_size - (dot_radius*2))]);
            new_position = [x, y];
            
            % Check if this new dot overlaps with any previous dot
            if i == 1 || all(sqrt(sum((dot_positions(1:i-1, :) - new_position).^2, 2)) >= 2 * dot_radius + min_spacing)
                dot_positions(i, :) = new_position;
                placed = true;
            end
            
            attempts = attempts + 1;
        end
        
        if attempts == max_attempts
            error('Failed to find non-overlapping positions for dots');
        end
    end
end

function weber_contrast = calculate_weber_contrast( bg_color, dot_color)
    % Calculate Weber contrast between dot color and background
    dot_luminance = dot_color(1);  % (assuming all dots have the same color)
    bg_luminance = bg_color(1); 
    
    weber_contrast = (dot_luminance - bg_luminance) / bg_luminance;
    
end

function dot_density = calculate_density(dot_positions)
    % Compute the negative median Euclidean distance between dots
    numDots = size(dot_positions, 1);
    
    if numDots > 1
        % Calculate pairwise distances between dots
        pairwise_distances = pdist(dot_positions); % Euclidean distances
        median_distance = median(pairwise_distances); % Median of distances
        dot_density = -median_distance; % Negative median distance
    else
        % If only one dot, set density to NaN (since no distance between dots)
        dot_density = NaN;
    end
end

function hull_area = calculate_convex_hull_area(dot_positions)
    % Compute the area of the convex hull spanning the dot centre positions.
    % A convex hull polygon requires at least 3 non-collinear points, so
    % numerosities of 0, 1, or 2 (or degenerate/collinear layouts) return NaN.
    numDots = size(dot_positions, 1);
    
    if numDots >= 3
        try
            hull_indices = convhull(dot_positions(:,1), dot_positions(:,2));
            hull_area = polyarea(dot_positions(hull_indices,1), dot_positions(hull_indices,2));
        catch
            % convhull errors out on collinear/degenerate point sets
            hull_area = NaN;
        end
    else
        hull_area = NaN;
    end
end

function calculate_correlation_matrix(image_metadata)
    % Extract data from the image_metadata structure
    numerosities = [image_metadata.numDots];
    total_dot_area = [image_metadata.totalDotArea];
    dot_density = [image_metadata.dotDensity];
    luminance = [image_metadata.luminance];
    mean_dot_area = [image_metadata.meanDotArea];
    convex_hull_area = [image_metadata.convexHullArea];

    % Calculate correlation matrix using pairwise method (ignores NaN, e.g.
    % from numerosities < 3 where the convex hull is undefined)
    data_matrix = [numerosities; total_dot_area; dot_density; luminance; ...
                   mean_dot_area; convex_hull_area]';
    save('data_matrix','data_matrix')
    correlation_matrix = corr(data_matrix, 'Rows', 'pairwise'); % 'pairwise' ignores NaN
    
    % Display the correlation matrix
    disp('Correlation Matrix between Numerosity and Image Properties:');
    disp(correlation_matrix);
    
    % Plot the correlation matrix as an image
    figure;
    imagesc(correlation_matrix);  % Create a heatmap of the correlation matrix
    colorbar;
    
    % Set axis labels
    labels = {'Numerosity', 'Total Dot Area', 'Density', 'Luminance', ...
              'Mean Dot Area', 'Convex Hull Area'};
    xticks(1:6);
    yticks(1:6);
    xticklabels(labels);
    yticklabels(labels);
    
    % Display correlation values on the heatmap cells
    for i = 1:size(correlation_matrix, 1)
        for j = 1:size(correlation_matrix, 2)
            text(j, i, sprintf('%.4f', correlation_matrix(i, j)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k',...
                'FontSize',14);
        end
    end
    colormap summer
    set(gca,'FontSize',16);
    xtickangle(45);

    
end