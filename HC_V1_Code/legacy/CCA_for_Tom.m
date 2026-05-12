%% Load data

load("data_for_Theo.mat")
load('period_experienced.mat')
load('licking_behaviour.mat')

%% Check areas
animal_areas_all = [];
for ianimal = 1:length(analysis_spatial)

    animal_areas_all = [animal_areas_all, categorical(units{ianimal}.regions_label)];
end

figure
histogram(animal_areas_all)


%% Parameters

n_animals = length(analysis_spatial);
area_activity = cell(1, n_animals);
n_components_reduced = 5;
n_trials_window = -3:3;
n_bins_window = -3:3;
n_shuffles = 20;
max_shift_bins = 2;

landmarks = [0, 50;...
             100, 125;...
         	 175, 200;...
             250, 275;...
         	 325, 350;...
         	 400, 425]/2.5;
learning_points = period_experienced(:, 1);

%% Run CCA
figure
t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

figure
u = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

figure
v = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

for ianimal = 1:n_animals
    animal_data = permute(analysis_spatial{ianimal}.firing.cued.freq, [1, 3, 2]);

    [~, num_bins, num_trials] = size(animal_data);

    animal_areas = units{ianimal}.regions_label;
    n_areas = numel(animal_areas);
    unit_areas = units{ianimal}.idx;

    ca1_v1_cca_all_trials = nan(n_components_reduced, num_trials);
    ca1_v1_cca_shuffled_trials = nan(n_components_reduced, num_trials, n_shuffles);
    ca1_v1_cca_shifted_trials = nan(n_components_reduced, num_trials, 2*max_shift_bins+1);
    ca1_v1_cca_shited_shuffled_trials = nan(n_components_reduced, num_trials, 2*max_shift_bins+1, n_shuffles);

    ca1_v1_cca_all_bins = nan(n_components_reduced, num_bins);
    ca1_v1_cca_shuffled_bins = nan(n_components_reduced, num_bins, n_shuffles);

    for iarea = 1:n_areas
        current_area_activity = animal_data(unit_areas(iarea, :), :, :);
        n_units_area = size(current_area_activity, 1);

        area_activity{ianimal}.(animal_areas(iarea)).spatial_binned_fr = current_area_activity;


        if n_units_area >= n_components_reduced

            [~,current_area_scores,~, ~, explained, ~] = pca(current_area_activity(:, :)', NumComponents=n_components_reduced);
            area_activity_reduced = reshape(current_area_scores', n_components_reduced, num_bins, num_trials);
            area_activity{ianimal}.(animal_areas(iarea)).area_activity_reduced = area_activity_reduced;
            fprintf('%s pca for animal %d with %d dimensions explains %.2f total variance\n', animal_areas(iarea), ianimal, n_components_reduced, sum(explained(1:n_components_reduced)))
        end
    end


    if isfield(area_activity{ianimal}, 'CA1') && isfield(area_activity{ianimal}, 'V1') && isfield(area_activity{ianimal}.CA1, 'area_activity_reduced') && isfield(area_activity{ianimal}.V1, 'area_activity_reduced') 
        for itrial = 1:num_trials
            trial_window = itrial + n_trials_window;
            trial_window = trial_window(trial_window >= 1 & trial_window <= num_trials);

            non_nan_bins = squeeze(~any(isnan(area_activity{ianimal}.CA1.area_activity_reduced(:, :, trial_window)), [1, 3]));

            ca1_reduced_pooled = area_activity{ianimal}.CA1.area_activity_reduced(:, non_nan_bins, trial_window);
            v1_reduced_pooled = area_activity{ianimal}.V1.area_activity_reduced(:, non_nan_bins, trial_window);

            try
                [~,~,r,~,~] = canoncorr(ca1_reduced_pooled(:, :)', v1_reduced_pooled(:, :)');
                n_r = min(length(r), n_components_reduced);
                ca1_v1_cca_all_trials(1:n_r, itrial) = r(1:n_r);

                for ishuffle = 1:n_shuffles
                    ca1_reduced_shuffled = ca1_reduced_pooled(:, randperm(sum(non_nan_bins), sum(non_nan_bins)), :);

                    [~,~,r,~,~] = canoncorr(ca1_reduced_shuffled(:, :)', v1_reduced_pooled(:, :)');
                    n_r = min(length(r), n_components_reduced);
                    ca1_v1_cca_shuffled_trials(1:n_r, itrial, ishuffle) = r(1:n_r);
                end



                for ishift = -max_shift_bins:max_shift_bins
                    ca1_reduced_shifted = circshift(ca1_reduced_pooled, ishift, 2);

                    [~,~,r,~,~] = canoncorr(ca1_reduced_shifted(:, :)', v1_reduced_pooled(:, :)');
                    n_r = min(length(r), n_components_reduced);
                    ca1_v1_cca_shifted_trials(1:n_r, itrial, ishift+1+max_shift_bins) = r(1:n_r);

                    for ishuffle = 1:n_shuffles
                        ca1_shifted_shuffled = ca1_reduced_shifted(:, randperm(sum(non_nan_bins), sum(non_nan_bins)), :);

                        [~,~,r,~,~] = canoncorr(ca1_shifted_shuffled(:, :)', v1_reduced_pooled(:, :)');
                        n_r = min(length(r), n_components_reduced);
                        ca1_v1_cca_shited_shuffled_trials(1:n_r, itrial, ishift+1+max_shift_bins, ishuffle) = r(1:n_r);
                    end

                end

            catch
                warning('CA1-V1 CCA failed on trial %d for animal %d', itrial, ianimal);
            end
        end

        nexttile(t)
        plot(ca1_v1_cca_all_trials(1, :)', LineWidth=1)
        hold on
        shadedErrorBar(1:num_trials, squeeze(mean(ca1_v1_cca_shuffled_trials(1, :, :), 3, "omitmissing")), squeeze(sem(ca1_v1_cca_shuffled_trials(1, :, :), 3)))
        title(sprintf('CA1-V1 - animal %d', ianimal))
        axis tight
        ylim([0, 0.7])
        xline(learning_points(ianimal))

        % figure
        % color_grad = colorGradient([0 0 1], [1, 0, 0], 2*max_shift_bins+1);
        % plot(squeeze(ca1_v1_cca_shifted_trials(1,:,:)), 'LineWidth', 1)
        % colororder(color_grad)
        % legend(strsplit(num2str(-max_shift_bins:max_shift_bins)))
        % title(sprintf('CA1-V1 - animal %d', ianimal))

        v1_to_ca1 = squeeze(mean(ca1_v1_cca_shifted_trials(1, :, 1:max_shift_bins), 3));
        neutral = squeeze(ca1_v1_cca_shifted_trials(1, :, max_shift_bins+1));
        ca1_to_v1 = squeeze(mean(ca1_v1_cca_shifted_trials(1, :, max_shift_bins+2:end), 3));
        ca1_precession_strength = (ca1_to_v1 - v1_to_ca1)./(ca1_to_v1 + v1_to_ca1);

        v1_to_ca1_shuffled = squeeze(mean(ca1_v1_cca_shited_shuffled_trials(1, :, 1:max_shift_bins, :), 3));
        ca1_to_v1_shuffled = squeeze(mean(ca1_v1_cca_shited_shuffled_trials(1, :, max_shift_bins+2:end, :), 3));
        ca1_precession_strength_shuffled = (ca1_to_v1_shuffled - v1_to_ca1_shuffled)./(ca1_to_v1_shuffled + v1_to_ca1_shuffled);

        nexttile(u)
        plot(ca1_precession_strength, LineWidth=1)
        hold on
        shadedErrorBar(1:num_trials, squeeze(mean(ca1_precession_strength_shuffled, 2, "omitmissing")), squeeze(sem(ca1_precession_strength_shuffled, 2)))
        title(sprintf('CA1-V1 - animal %d', ianimal))
        yline(0, '--')
        xline(learning_points(ianimal))
        axis tight
        ylim([-0.3, 0.3])



        non_nan_trials = squeeze(~any(isnan(area_activity{ianimal}.CA1.area_activity_reduced), [1, 2]));

        for ibin = 1:num_bins
            bin_window = ibin + n_bins_window;
            bin_window = bin_window(bin_window >= 1 & bin_window <= num_bins);

            ca1_reduced_pooled = area_activity{ianimal}.CA1.area_activity_reduced(:, bin_window, non_nan_trials);
            v1_reduced_pooled = area_activity{ianimal}.V1.area_activity_reduced(:, bin_window, non_nan_trials);

            try
                [~,~,r,~,~] = canoncorr(ca1_reduced_pooled(:, :)', v1_reduced_pooled(:, :)');
                n_r = min(length(r), n_components_reduced);
                ca1_v1_cca_all_bins(1:n_r, ibin) = r(1:n_r);

                for ishuffle = 1:n_shuffles
                    ca1_reduced_shuffled = ca1_reduced_pooled(:, :, randperm(sum(non_nan_trials), sum(non_nan_trials)));

                    [~,~,r,~,~] = canoncorr(ca1_reduced_shuffled(:, :)', v1_reduced_pooled(:, :)');
                    n_r = min(length(r), n_components_reduced);
                    ca1_v1_cca_shuffled_bins(1:n_r, ibin, ishuffle) = r(1:n_r);
                end
            catch
                warning('CA1-V1 CCA failed on bin %d for animal %d', ibin, ianimal);
            end 
        end

        nexttile(v)
        plot(ca1_v1_cca_all_bins(1, :)', LineWidth=1)
        hold on
        shadedErrorBar(1:num_bins, squeeze(mean(ca1_v1_cca_shuffled_bins(1, :, :), 3, "omitmissing")), squeeze(sem(ca1_v1_cca_shuffled_bins(1, :, :), 3)))
        title(sprintf('CA1-V1 - animal %d', ianimal))
        axis tight
        ylim([0, 0.7])
        xline(landmarks(:))
        
    end

end

xlabel(t, 'trials')
ylabel(t, 'correlation')

xlabel(u, 'trials')
ylabel(u, 'CA1 to V1 relative strength')

xlabel(v, 'spatial bins')
ylabel(v, 'correlation')


