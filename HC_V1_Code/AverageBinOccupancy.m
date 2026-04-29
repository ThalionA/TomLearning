base_dir = '/Users/theoamvr/Desktop/Experiments/TomLearning/';
data_dir = fullfile(base_dir, 'HC_V1_data');
file_list = dir(fullfile(data_dir, 'TF*_export.mat'));
n_files = length(file_list);

% Preallocate
occupancy_data = cell(1, n_files);

for i_animal = 1:n_files
    filename = file_list(i_animal).name;
    fullpath = fullfile(data_dir, filename);
    fprintf('Loading %d/%d: %-20s ... \n', i_animal, n_files, filename);
    
    % Load only the required struct
    dat = load(fullpath, 'analysis_behaviour');
    occupancy_data{i_animal} = dat.analysis_behaviour.spatialOccupation.cued.time_minspeed_gf(:);
end


% 2. Load Tom's learner index 
% (Ensure 'learnt_idx' matches the exact variable name in his .mat file)
learner_file = fullfile(data_dir, 'animal_learnt_alphabetical_order.mat');
learner_dat = load(learner_file); 
learnt_idx = learner_dat.animal_learnt_bool_sorted; 

% 2. Aggregate data and create grouping vectors for unequal trial lengths
animal_labels = cell(1, n_files);
color_labels = cell(1, n_files);

for i_animal = 1:n_files
    num_datapoints = length(occupancy_data{i_animal});
    
    % Positional grouping (Animal ID)
    animal_labels{i_animal} = repmat(i_animal, num_datapoints, 1);
    
    % Color grouping (Learner vs Non-Learner)
    if learnt_idx(i_animal) == 1
        status = 'Learner';
    else
        status = 'Non-Learner';
    end
    color_labels{i_animal} = repmat(categorical({status}), num_datapoints, 1);
end

% Flatten all arrays
all_occupancy = vertcat(occupancy_data{:});
xgroupdata = vertcat(animal_labels{:});
cgroupdata = vertcat(color_labels{:});

%% Plotting
figure('Position', [100, 100, 1200, 400]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Tile 1: Individual animals (Color-coded by learning status)
nexttile
v = violinplot(xgroupdata, all_occupancy, 'GroupByColor', cgroupdata);
axis tight

% v is now a 1x2 array (one object for Learner, one for Non-Learner).
% Categorical sorts alphabetically by default: v(1) = Learner, v(2) = Non-Learner
v(1).FaceColor = [0, 0.4470, 0.7410]; % Blue
v(2).FaceColor = [0.8500, 0.3250, 0.0980]; % Orange

ylabel('Occupancy (sec)');
xlabel('Animal Index');
title('Distribution per Animal');

% Tile 2: Aggregate distribution across all animals
nexttile
% Split the cell array data based on learning status
% (Assumes learnt_idx is a binary vector of 1s and 0s)
occ_learners = vertcat(occupancy_data{learnt_idx == 1});
occ_nonlearners = vertcat(occupancy_data{learnt_idx == 0});

hold on
histogram(occ_learners, 50, 'Normalization', 'pdf', 'FaceColor', [0, 0.4470, 0.7410], 'DisplayName', 'Learner');
histogram(occ_nonlearners, 50, 'Normalization', 'pdf', 'FaceColor', [0.8500, 0.3250, 0.0980], 'DisplayName', 'Non-Learner');
hold off

xlabel('Occupancy (sec)');
ylabel('Probability Density');
title('Aggregate Distribution');
legend('Location', 'best');

% 1. Calculate the mean occupancy within each animal (Hierarchical level 1)
animal_means = cellfun(@(x) mean(x, 'omitnan'), occupancy_data);

% 2. Separate the animal-level means by learning status
animal_means_L = animal_means(learnt_idx == 1);
animal_means_NL = animal_means(learnt_idx == 0);

% 3. Calculate group-level Mean and SEM across animals (Hierarchical level 2)
% N is now the number of animals in each group
mean_L = mean(animal_means_L, 'omitnan');
sem_L = std(animal_means_L, 'omitnan') / sqrt(sum(~isnan(animal_means_L)));

mean_NL = mean(animal_means_NL, 'omitnan');
sem_NL = std(animal_means_NL, 'omitnan') / sqrt(sum(~isnan(animal_means_NL)));

% 4. Format the text string (including N animals for clarity)
stats_str = sprintf('Learners (N=%d): %.2f \\pm %.3f s\nNon-Learners (N=%d): %.2f \\pm %.3f s', ...
    length(animal_means_L), mean_L, sem_L, length(animal_means_NL), mean_NL, sem_NL);

% 5. Add text box to the top right of the active axes
text(gca, 0.95, 0.95, stats_str, ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', ...
    'EdgeColor', 'k', ...
    'Interpreter', 'tex');