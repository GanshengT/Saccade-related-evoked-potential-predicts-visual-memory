%% description
% This script will calculate the ERP for controls (baseline / EOG signals),
% to validate if previously identified SREPs were consistent This script is

clear;

addpath(['CircStat2012a']);
data_dir = 'data/';
sum_contact_table_filename = ['data/norm_contact_data_table.csv'];
%% parameters
% subjs = {'BJH021', 'BJH021', 'BJH024', 'BJH024', 'BJH025', 'BJH026', 'BJH026', 'BJH027', 'BJH027', 'BJH027',...
%     'BJH028', 'BJH029', 'BJH030'};
% tasks = {'BLAES_study', 'BLAES_study_2', 'BLAES_study', 'BLAES_study_2', 'BLAES_study_twosource', ...
%     'BLAES_study', 'BLAES_study_2', ...
%     'BLAES_study', 'BLAES_study_2', 'BLAES_study_3'...
%     'BLAES_study', 'BLAES_study', 'BLAES_study'};
% n_session_per_task = [2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2];
% running for one subject %
subjs = {'BJH025',};
% we can use the mean to differentiate left and right electrodes
tasks = { 'BLAES_study_twosource'};
n_session_per_task = [4];



eye_pos_in_mri_space = containers.Map();

eye_pos_in_mri_space('BJH025_left') = [-29.84, 49.21, -9.67];
eye_pos_in_mri_space('BJH025_right') = [35.51, 45.80, -9.67];
reference_subj_left_eye = [-29.84, 49.21, -9.67];
reference_subj_right_eye = [35.51, 45.80, -9.67];

BJH025_tasks = struct();
BJH025_tasks.BLAES_study_twosource = {'D8','GR4','GR5','J7'};
bad_ch_map('BJH025') = BJH025_tasks;


cluster_size_req = 0.25; % no less than 0.25 size of the other cluster
laplacian_thres = 6; % 5mm bound, switched to 6mm, same as the active processing
goodness_of_fit_thres = 0.8;
monitor_dimension = [611.2, 362.6]; % in mm
smoothing_window = 10; % 10* 1/2000s = 5ms;
% we use different smoothing_window if no change points are detected
smoothing_window_secondary = 20; % 20* 1/2000s = 10ms;
epsilon = 1e-8; % handle all 0 vector, make corr robust

trial_begin_study = 3; % 3s
sfreq = 2000; % 2000Hz
baseline_time_precedence = [-0.2, -0.09];
erp_range = [-0.01 0.04]; % 200ms
erp_pre_during_post = [-0.3 -0.04 -0.01 0.03 0.8 0.3]; % 600ms
erp_characteristic_range = [-0.05, 0.1]; % for getting change pts
erp_amp_window = 0.005; % in s, we look at max around the second change pt
alpha_test = 0.05;
low_cutoff = 4; % Low cutoff frequency (Hz) - narrow band filtering for erp
high_cutoff = 12; % High cutoff frequency (Hz)
padding_duration = 0.5; % Padding duration in seconds
color_limits = [-5, 5];

time_range_for_plotting_erp = [-0.2, 0.2];
viz_single_trial = 0; % dont plot the fitting for single trial
% load contact table
sum_contact_table = readtable(sum_contact_table_filename);
alpha_test = alpha_test / height(sum_contact_table);

window = 2048; % half a second
noverlap = 1024; 
nfft = 2 * window;
freq_range = [1 50];
% i/f fitting
freq_range_fit = [1 60]; 
narrow_bandwidth = 3; % half band width
% prepare a df, with these col names:
% subj, task, session, trial, ch_name_bci2000, blc_theta, blc_gamma


variableNames = {'subject', 'task', 'ch_name_bci2000_format', 'overall_significance', ...
    'overall_auc', 'overall_p', 'cluster_separativity', 'cluster1_corr', 'cluster2_corr',...
    'cluster_separativity_bl_control', 'cluster1_corr_bl_control', 'cluster2_corr_bl_control',...
    'n_weak_corr_trials', 'n_weak_corr_trials_bl_control',...
    'all_var_azimuth_elevation_angle', ...
    'all_cluster_4_12_hz_plv_pre_pre_pre','all_cluster_4_12_hz_plv_pre_pre', 'all_cluster_4_12_hz_plv_pre', ...
    'all_cluster_4_12_hz_plv_during', 'all_cluster_4_12_hz_plv_post',...
    'cluster1_peak_freq', 'cluster1_oscillation_amp', 'cluster1_4_12_hz_pow_pre_pre_pre','cluster1_4_12_hz_pow_pre_pre',...
    'cluster1_4_12_hz_pow_pre', 'cluster1_4_12_hz_pow_during', 'cluster1_4_12_hz_pow_post', ...
    'cluster1_no_erp_phase_pre_pre_pre','cluster1_no_erp_phase_pre_pre', 'cluster1_no_erp_phase_pre', ...
    'cluster1_4_12_hz_phase_pre_pre_pre','cluster1_4_12_hz_phase_pre_pre', 'cluster1_4_12_hz_phase_pre', ...
    'cluster1_4_12_hz_phase_during', 'cluster1_4_12_hz_phase_post', ...
    'cluster1_4_12_hz_plv_all', 'cluster1_4_12_hz_plv_pre_pre_pre', 'cluster1_4_12_hz_plv_pre_pre', ...
    'cluster1_4_12_hz_plv_pre', 'cluster1_4_12_hz_plv_during', 'cluster1_4_12_hz_plv_post', ...
    'cluster2_peak_freq', 'cluster2_oscillation_amp' 'cluster2_4_12_hz_pow_pre_pre_pre', 'cluster2_4_12_hz_pow_pre_pre', ...
    'cluster2_4_12_hz_pow_pre', 'cluster2_4_12_hz_pow_during', 'cluster2_4_12_hz_pow_post',...
    'cluster2_no_erp_phase_pre_pre_pre','cluster2_no_erp_phase_pre_pre', 'cluster2_no_erp_phase_pre', ...
    'cluster2_4_12_hz_phase_pre_pre_pre', 'cluster2_4_12_hz_phase_pre_pre','cluster2_4_12_hz_phase_pre', 'cluster2_4_12_hz_phase_during', 'cluster2_4_12_hz_phase_post', ...
    'cluster2_4_12_hz_plv_all', 'cluster2_4_12_hz_plv_pre_pre_pre',  'cluster2_4_12_hz_plv_pre_pre', ...
    'cluster2_4_12_hz_plv_pre', 'cluster2_4_12_hz_plv_during', 'cluster2_4_12_hz_plv_post', ...
    'cluster1_significance', 'cluster1_auc', 'cluster1_p', 'cluster1_auc_bl_control', 'cluster1_p_bl_control', ...
    'cluster1_latency', 'cluster1_amp', 'cluster1_amp_time', 'cluster1_polarity','cluster1_peak_to_trough', 'cluster1_trough_time', ...
    'cluster2_significance', 'cluster2_auc', 'cluster2_p', 'cluster2_auc_bl_control', 'cluster2_p_bl_control', ...
    'cluster2_latency', 'cluster2_amp', ...
    'cluster2_amp_time', 'cluster2_polarity', 'cluster2_peak_to_trough', 'cluster2_trough_time', ...
    'cluster1_mean_dir_azimuth', 'cluster1_mean_dir_elevation', ...
    'cluster1_mean_azimuth_elevation_angle', 'cluster1_var_azimuth_elevation_angle', ...
    'cluster1_mean_dir_x_screen_coord', 'cluster1_mean_dir_y_screen_coord', ...
    'cluster2_mean_dir_azimuth', 'cluster2_mean_dir_elevation', ...
    'cluster2_mean_azimuth_elevation_angle', 'cluster2_var_azimuth_elevation_angle', ...
    'cluster2_mean_dir_x_screen_coord', 'cluster2_mean_dir_y_screen_coord', ...
    're_ref_vector_x_norm', 're_ref_vector_y_norm', 're_ref_vector_z_norm', ...
    'contact_x_norm', 'contact_y_norm', 'contact_z_norm', 'left_eye_x_norm', ...
    'left_eye_y_norm', 'left_eye_z_norm', 'right_eye_x_norm', ...
    'right_eye_y_norm', 'right_eye_z_norm', ...
    'valid_side', 'contact_side', 'contact_dist_to_left_eye', 'contact_dist_to_right_eye', 'shank_ref', ...
    'weak_mean_dir_azimuth', 'weak_mean_dir_elevation'};

variableTypes = {'string', 'string', 'string', 'logical',...
    'double', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', ...
    'double', 'double', ...
    'double', ...
    'double', 'double', 'double',...
    'double', 'double', ...
    'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', ...
    'double', 'double', 'double', ...
    'double', 'double', 'double', ...
    'double', 'double', ...
    'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', ...
    'double', 'double', 'double', ...
    'double', 'double', 'double', ...
    'double', 'double', 'double', ...
    'double', 'double', 'double', 'double','double', ...
    'double', 'double', 'double',...
    'double', 'double', 'double', ...
    'logical', 'double', 'double',  'double', 'double', ...
    'double', 'double', 'double', 'string', 'double', 'double',...
    'logical', 'double', 'double',  'double', 'double', ...
    'double', 'double', ...
    'double', 'string', 'double', 'double',...
    'double', 'double', ...
    'double', 'double', ...
    'double','double', ...
    'double', 'double', ...
    'double', 'double', ...
    'double','double', ...
    'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', ...
    'double', 'double', ...
    'string', 'string', 'double', 'double', 'logical', ...
    'double', 'double'};


for i_subj = 1:length(subjs)

    responsive_erp_table = table('Size', [0, length(variableTypes)], ...
        'VariableTypes', variableTypes, ...
        'VariableNames', variableNames);
    df_evoked_potential = table('Size', [0, 44], ...
                    'VariableTypes', {'double', 'string', 'string', 'string', ...
                    'double', 'double', 'double', ...
                    'double', 'double', 'double', 'double', ...
                    'double', 'double', 'double', 'double', ...
                    'double', 'double', 'string', 'string', ...
                    'double', 'double', 'double', ...
                    'double', 'double', 'double', 'string'...
                    'double','double', 'double', 'double', 'double', ...
                    'double','double', 'double', 'double', 'double', ...
                    'double', 'double', 'double', 'double', 'double', ...
                    'double', 'double', 'double'}, ...
                    'VariableNames', {'saccade_id', 'subj', 'task', 'ch_name_bci2000_format', ...
                    'session', 'trial', 'eccentricity', ...
                    'azimuth_before', 'elevation_before', 'azimuth_after', 'elevation_after', ...
                    'gaze_x_before', 'gaze_y_before', 'gaze_x_after', 'gaze_y_after', ...
                    'duration', 'time_to_image_onset', 'behavior_condition', 'remember_condition', ...
                    'typical_erp_score', 'peak_to_trough', 'trough_time',...
                    'latency', 'amp', 'peak_time', 'polarity',...
                    '4_12_hz_pow_pre_pre_pre', '4_12_hz_pow_pre_pre', '4_12_hz_pow_pre', '4_12_hz_pow_during', '4_12_hz_pow_post',...
                    '4_12_hz_phase_pre_pre_pre', '4_12_hz_phase_pre_pre', '4_12_hz_phase_pre', '4_12_hz_phase_during', '4_12_hz_phase_post',...
                    '4_12_hz_norm_pow_pre_pre_pre', '4_12_hz_norm_pow_pre_pre', '4_12_hz_norm_pow_pre', '4_12_hz_norm_pow_during', '4_12_hz_norm_pow_post',...
                    'no_erp_phase_pre_pre_pre', 'no_erp_phase_pre_pre', 'no_erp_phase_pre'});        

    subj = subjs{i_subj};
    task = tasks{i_subj};
    n_session = n_session_per_task(i_subj);
    % get ch_label
    contact_this_subj = sum_contact_table(strcmp(sum_contact_table.SubjectID, subj), :);
    report_path = fullfile(data_dir, subj, 'report');
    coordinates = table2array(contact_this_subj(:, ...
        {'X_aligned_to_brightest_voxel', 'Y_aligned_to_brightest_voxel', 'Z_aligned_to_brightest_voxel'}));
    coordinates_norm = table2array(contact_this_subj(:, ...
        {'Norm_X', 'Norm_Y', 'Norm_Z'}));
    is_micro = contains(contact_this_subj.ShankID, '_mic');
    is_macro = ~is_micro;
    micro_indices = find(is_micro);
    macro_indices = find(is_macro);
    dist_matrix = squareform(pdist(coordinates)); % Efficient distance computation
    dist_matrix(is_micro, is_macro) = Inf;
    dist_matrix(is_macro, is_micro) = Inf;

    for idx1 = 1:length(micro_indices)
        i_micro = micro_indices(idx1);
        for idx2 = 1:length(micro_indices)
            j_micro = micro_indices(idx2);
            if i_micro ~= j_micro
                dist_matrix(i_micro, j_micro) = dist_matrix(i_micro, j_micro) + 0.1;
            end
        end
    end

    % dist_matrix = zeros(size(coordinates, 1), size(coordinates, 1));
    % 
    % for i = 1:size(coordinates, 1)
    %     for j = 1:size(coordinates, 1)
    %         dist_matrix(i, j) = sqrt(sum((coordinates(i,:) - coordinates(j,:)).^2));
    %     end
    % end
    % for i_ch = 1:height(contact_this_subj)
    % for efficiency, we run for a subset of channels
    for i_ch = 1:3
        close all;
        ch_label = contact_this_subj.labels_majority(i_ch);
        % get reref direction
        % we also need to model the eyeball
        % we use r = 1, azimuth, and elevation to model a point on the eyeball
        % then we need to transform the shank vector to this coordinates

        neighbors= find(dist_matrix(i_ch,:) < laplacian_thres & dist_matrix(i_ch,:) > 0);
        bad_ch_this_subj = bad_ch_map(subj);
        bad_ch_this_subj = bad_ch_this_subj.(task);

        original_neighbors = neighbors;
        new_neighbors = [];
        for j = 1:length(original_neighbors)
            neighbor_idx = original_neighbors(j);
            neighbor_bci2000_ch_name = contact_this_subj.ch_name_bci2000_format(neighbor_idx);
            if ismember(neighbor_bci2000_ch_name, bad_ch_this_subj)
                fprintf('Removed bad channel %d from neighbors of channel %d\n', neighbor_idx, i_ch);
            else
                new_neighbors(end+1) = neighbor_idx;
            end
        end
        neighbors = new_neighbors;

        % fprintf('Sanity Check: Channels with either 0 or >2 neighbors:\n');
        num_neighbors = length(neighbors);
        if num_neighbors == 0 || num_neighbors > 2
            fprintf('Channel %d has %d neighbor(s): %s\n', i_ch, num_neighbors, mat2str(neighbors));
        end

        current_point = coordinates(i_ch, :);
        current_point_norm = coordinates_norm(i_ch, :);
        if current_point_norm(1) > mean([reference_subj_left_eye(1), ...
                reference_subj_right_eye(1)], 'omitnan')
            contact_side = 'right';
        else
            contact_side = 'left';
        end

        neighbor_indices = neighbors;
        all_points = [current_point_norm; coordinates_norm(neighbor_indices, :)];
        % shank_ref coding; 0: no rereference, 2: small laplacian, 1:
        % bipolar; 3: not linear small laplacian
        clear shank_direction
        if size(all_points, 1) < 2 % no neighbor, seldom happen
            shank_ref = 0;
            re_ref_vector_x_norm = nan;
            re_ref_vector_y_norm = nan;
            re_ref_vector_z_norm = nan;

        else
            if size(all_points, 1) > 2 
                shank_ref = 2;
            else
                shank_ref = 1;
            end

            [coeff, score, ~, ~, explained] = pca(all_points);
            

            r_squared = explained(1) / sum(explained);
            if r_squared > goodness_of_fit_thres
                
                shank_direction = coeff(:,1)';
                re_ref_vector_x_norm = coeff(1,1);
                re_ref_vector_y_norm = coeff(2,1);
                re_ref_vector_z_norm = coeff(3,1);

            else
                shank_ref = 3;
                re_ref_vector_x_norm = nan;
                re_ref_vector_y_norm = nan;
                re_ref_vector_z_norm = nan;
            end
        end

      
        prep_signal_all_session = [];
        saccade_table_all_session = [];
        for session = 1:n_session
            filename = fullfile(data_dir, subj,  ...
                [subj, '_', task, '_session', num2str(session), ...
                '_channel', num2str(i_ch), '_saccade_onset_control_prep_signal.mat']);
            prep_signal = load(filename);
            prep_signal = prep_signal.prep_signal_saccade_onset_control;
            filename = fullfile(data_dir, subj, [subj '_' task '_session' num2str(session) ...
                    '_prep_saccade_event.mat']);
            prep_saccade_event = load(filename);
            prep_saccade_event = prep_saccade_event.data_struct_2_export;
            valid_side = prep_saccade_event.used_eye_for_saccade;
            % get eye position in MRI coordinates
            key_for_eye = [subj '_left'];
            left_eye_position_in_mri_coord = eye_pos_in_mri_space(key_for_eye);
            contact_dist_to_left_eye = norm(current_point - left_eye_position_in_mri_coord);
            key_for_eye = [subj '_right'];
            right_eye_position_in_mri_coord = eye_pos_in_mri_space(key_for_eye);
            contact_dist_to_right_eye = norm(current_point - right_eye_position_in_mri_coord);

           

            filename = fullfile(data_dir, subj, ...
                [subj, '_', task, '_session', num2str(session), '_saccade_table.mat']);
            saccade_table = load(filename);
            saccade_table = saccade_table.saccade_table;
            % verify if signal match the saccade event
            assert(size(prep_signal, 2) == sum(contains(saccade_table.label, 'control')), ...
            'saccade event and prep signal length not matched');
            prep_signal_all_session = [prep_signal_all_session, prep_signal];
            saccade_table_all_session = [saccade_table_all_session; ...
                saccade_table(contains(saccade_table.label, 'control'), :)];
        end

        % merge erp and check if responsive
        time_vector = (0:size(prep_signal_all_session, 1) - 1) / sfreq - trial_begin_study;

        task_indices = find(time_vector >= erp_range(1) & time_vector <= erp_range(2));
        task_data = prep_signal_all_session(task_indices, :);

        baseline_data = [];
        for i_trial = 1:size(task_data, 2)
            rand_start = randi([round(baseline_time_precedence(1) * sfreq),...
                round(baseline_time_precedence(2) * sfreq)], 1);
            [~, rand_start_index] = min(abs(time_vector - rand_start / sfreq));
            baseline_data = [baseline_data, prep_signal_all_session(...
                rand_start_index:rand_start_index + length(task_indices) - 1, i_trial)];
        end
        baseline_indices = find(time_vector >= baseline_time_precedence(1) & time_vector <= erp_range(1));
        baseline_data = prep_signal_all_session(baseline_indices, :);
    
        baseline_mean = mean(baseline_data, 1);
        baseline_std = std(baseline_data, [], 1);

        if any(baseline_std == 0)
            zero_std_trial_indices = find(baseline_std==0);
            saccade_table_all_session(zero_std_trial_indices, :) = [];
            prep_signal_all_session(:, zero_std_trial_indices) = [];
            task_indices = find(time_vector >= erp_range(1) & time_vector <= erp_range(2));
            task_data = prep_signal_all_session(task_indices, :);
            baseline_indices = find(time_vector >= baseline_time_precedence(1) & time_vector <= erp_range(1));
            baseline_data = prep_signal_all_session(baseline_indices, :);
        
            baseline_mean = mean(baseline_data, 1);
            baseline_std = std(baseline_data, [], 1);
        end


        zscored_baseline = (baseline_data - baseline_mean) ./ (baseline_std );
        zscored_task = (task_data - baseline_mean) ./ (baseline_std);

        baseline_corr = corrcoef(zscored_baseline);
        task_corr = corrcoef(zscored_task);
        baseline_corr_values = baseline_corr(tril(true(size(baseline_corr)), -1));
        task_corr_values = task_corr(tril(true(size(task_corr)), -1));
        
        [p_value, ~, stats] = ranksum(baseline_corr_values, task_corr_values);
        combined_values = [baseline_corr_values; task_corr_values];
        labels = [zeros(length(baseline_corr_values), 1); ones(length(task_corr_values), 1)];
        
        [sorted_values, sorted_indices] = sort(combined_values, 'ascend');
        sorted_labels = labels(sorted_indices);
        
        total_positives = sum(sorted_labels);
        total_negatives = length(sorted_labels) - total_positives;
        
        tpr = zeros(length(sorted_labels), 1);
        fpr = zeros(length(sorted_labels), 1);
        
        for i = 1:length(sorted_labels)
            tpr(i) = sum(sorted_labels(i:end)) / total_positives; % Proportion of positive labels above the current threshold
            fpr(i) = sum(~sorted_labels(i:end)) / total_negatives; % Proportion of negative labels above the current threshold
        end
        
        auc = - trapz(fpr, tpr);
        norm_concatenated_prep_signal = (prep_signal_all_session - baseline_mean) ./ (baseline_std);
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        
        time_vector = (0:length(mean_erp) - 1) / sfreq - trial_begin_study;

        

        %%  baseline control here
        weak_trials_indices = []; % To store indices of weak trials
        original_baseline_corr = baseline_corr;
        original_order = 1:size(baseline_corr, 1); % Keep track of the original indices
        dist_matrix_trial = 1 - (baseline_corr);
        Z = linkage(squareform(dist_matrix_trial), 'average');
        order = optimalleaforder(Z, dist_matrix_trial);
        clusters = cluster(Z, 'maxclust', 2);
        clusters = clusters(order);
        while true
            cluster_sizes = histcounts(clusters, 2);
            if min(cluster_sizes) >= cluster_size_req * max(cluster_sizes)
                break;
            end
            
            mean_abs_corr = mean(abs(baseline_corr), 2);
            [~, weakest_trial_idx] = min(mean_abs_corr);
            
            weak_trials_indices = [weak_trials_indices; original_order(weakest_trial_idx)];
            
            baseline_corr(weakest_trial_idx, :) = [];
            baseline_corr(:, weakest_trial_idx) = [];
            
            original_order(weakest_trial_idx) = [];
            
            dist_matrix_trial = 1 - (baseline_corr);
            Z = linkage(squareform(dist_matrix_trial), 'average');
            order = optimalleaforder(Z, dist_matrix_trial);
            clusters = cluster(Z, 'maxclust', 2);
            clusters = clusters(order);
        end
        
        final_order = [original_order(order)'; weak_trials_indices]';

        if ~isempty(weak_trials_indices)
            reordered_baseline_corr = original_baseline_corr(final_order(1:end-length(weak_trials_indices)), final_order(1:end-length(weak_trials_indices)));
            weak_baseline_corr = original_baseline_corr(final_order(1:end-length(weak_trials_indices)), weak_trials_indices);
            bottom_right_matrix = original_baseline_corr(weak_trials_indices, weak_trials_indices);
            combined_matrix = [reordered_baseline_corr, weak_baseline_corr; weak_baseline_corr', bottom_right_matrix];
        else
            combined_matrix = original_baseline_corr(final_order, final_order);
        end
        
       
        clusters_2 = cluster(Z, 'maxclust', 2);
        clusters_2 = clusters_2(order);
       
        
       

        cluster1_trials = final_order(clusters == 1);
        cluster2_trials = final_order(clusters == 2);

        within_cluster1_corr = original_baseline_corr(cluster1_trials, cluster1_trials);
        within_cluster2_corr = original_baseline_corr(cluster2_trials, cluster2_trials);
        mean_within_cluster1_corr_bl_control = mean(within_cluster1_corr(~eye(size(within_cluster1_corr))));
        mean_within_cluster2_corr_bl_control = mean(within_cluster2_corr(~eye(size(within_cluster2_corr))));
        between_cluster_corr = original_baseline_corr(cluster1_trials, cluster2_trials);
        mean_between_cluster_corr = mean(between_cluster_corr(:));
        cluster_corr_difference_bl_control = (mean_within_cluster1_corr_bl_control + mean_within_cluster2_corr_bl_control) / 2 - mean_between_cluster_corr;
        n_weak_corr_trials_bl_control = length(weak_trials_indices);

        cluster1_task_data = task_data(:, cluster1_trials);
        cluster1_baseline_mean = baseline_mean(:, cluster1_trials);
        cluster1_baseline_std = baseline_std(:, cluster1_trials);
        zscored_cluster1_baseline = (baseline_data(:, cluster1_trials) - cluster1_baseline_mean) ./ (cluster1_baseline_std);

        zscored_cluster1_task = (cluster1_task_data - cluster1_baseline_mean) ./ (cluster1_baseline_std);
        cluster1_corr = corrcoef(zscored_cluster1_task);
        cluster1_corr_values = cluster1_corr(tril(true(size(cluster1_corr)), -1));
        cluster1_baseline_corr = corrcoef(zscored_cluster1_baseline);
        cluster1_baseline_corr_values = cluster1_baseline_corr(tril(true(size(cluster1_baseline_corr)), -1));
        

        [cluster1_p_bl_control, ~, cluster1_stats] = ranksum(cluster1_corr_values, cluster1_baseline_corr_values);
        cluster1_combined_values = [cluster1_corr_values; cluster1_baseline_corr_values];
        cluster1_labels = [zeros(length(cluster1_corr_values), 1); ones(length(cluster1_baseline_corr_values), 1)];

        [sorted_values, sorted_indices] = sort(cluster1_combined_values, 'ascend');
        sorted_labels = cluster1_labels(sorted_indices);

        total_positives = sum(sorted_labels);
        total_negatives = length(sorted_labels) - total_positives;

        tpr = zeros(length(sorted_labels), 1);
        fpr = zeros(length(sorted_labels), 1);

        for i = 1:length(sorted_labels)
            tpr(i) = sum(sorted_labels(i:end)) / total_positives;
            fpr(i) = sum(~sorted_labels(i:end)) / total_negatives;
        end

        cluster1_auc_bl_control = -trapz(fpr, tpr);
        norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster1_trials) - ...
            cluster1_baseline_mean) ./ (cluster1_baseline_std);
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        
        time_vector = (0:length(mean_erp) - 1) / sfreq - trial_begin_study;
        

        cluster2_task_data = task_data(:, cluster2_trials);
        cluster2_baseline_mean = baseline_mean(:, cluster2_trials);
        cluster2_baseline_std = baseline_std(:, cluster2_trials);
        zscored_cluster2_baseline = (baseline_data(:, cluster2_trials) - cluster2_baseline_mean) ./ (cluster2_baseline_std);

        zscored_cluster2_task = (cluster2_task_data - cluster2_baseline_mean) ./ (cluster2_baseline_std);
        cluster2_corr = corrcoef(zscored_cluster2_task);
        cluster2_corr_values = cluster2_corr(tril(true(size(cluster2_corr)), -1));
        cluster2_baseline_corr = corrcoef(zscored_cluster2_baseline);
        cluster2_baseline_corr_values = cluster2_baseline_corr(tril(true(size(cluster2_baseline_corr)), -1));
        
        [cluster2_p_bl_control, ~, cluster2_stats] = ranksum(cluster2_corr_values, cluster2_baseline_corr_values);
        cluster2_combined_values = [cluster2_corr_values; cluster2_baseline_corr_values];
        cluster2_labels = [zeros(length(cluster2_corr_values), 1); ones(length(cluster2_baseline_corr_values), 1)];

        [sorted_values, sorted_indices] = sort(cluster2_combined_values, 'ascend');
        sorted_labels = cluster2_labels(sorted_indices);

        total_positives = sum(sorted_labels);
        total_negatives = length(sorted_labels) - total_positives;

        tpr = zeros(length(sorted_labels), 1);
        fpr = zeros(length(sorted_labels), 1);

        for i = 1:length(sorted_labels)
            tpr(i) = sum(sorted_labels(i:end)) / total_positives;
            fpr(i) = sum(~sorted_labels(i:end)) / total_negatives;
        end
        cluster2_auc_bl_control = -trapz(fpr, tpr);


        norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster2_trials) - ...
            cluster2_baseline_mean) ./ (cluster2_baseline_std);
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        

        baseline_corr = original_baseline_corr;

        %% baseline control done

        weak_trials_indices = []; 
        original_task_corr = task_corr;
        original_order = 1:size(task_corr, 1); 
        dist_matrix_trial = 1 - (task_corr);
        Z = linkage(squareform(dist_matrix_trial), 'average');
        order = optimalleaforder(Z, dist_matrix_trial);
        
        clusters = cluster(Z, 'maxclust', 2);
        clusters = clusters(order);
        
        while true
            cluster_sizes = histcounts(clusters, 2);
            if min(cluster_sizes) >= cluster_size_req * max(cluster_sizes)
                break;
            end
            
            mean_abs_corr = mean(abs(task_corr), 2);
            [~, weakest_trial_idx] = min(mean_abs_corr);
            
            weak_trials_indices = [weak_trials_indices; original_order(weakest_trial_idx)];
            
            task_corr(weakest_trial_idx, :) = [];
            task_corr(:, weakest_trial_idx) = [];
            
            original_order(weakest_trial_idx) = [];
            
            dist_matrix_trial = 1 - (task_corr);
            Z = linkage(squareform(dist_matrix_trial), 'average');
            order = optimalleaforder(Z, dist_matrix_trial);
            clusters = cluster(Z, 'maxclust', 2);
            clusters = clusters(order);
        end
        
        final_order = [original_order(order)'; weak_trials_indices]';
        
        if ~isempty(weak_trials_indices)
            reordered_task_corr = original_task_corr(final_order(1:end-length(weak_trials_indices)), final_order(1:end-length(weak_trials_indices)));
            weak_task_corr = original_task_corr(final_order(1:end-length(weak_trials_indices)), weak_trials_indices);
            bottom_right_matrix = original_task_corr(weak_trials_indices, weak_trials_indices);
            combined_matrix = [reordered_task_corr, weak_task_corr; weak_task_corr', bottom_right_matrix];
        else
            combined_matrix = original_task_corr(final_order, final_order);
        end
        
       
        if ~isempty(weak_trials_indices)
            weak_task_data = task_data(:, weak_trials_indices);
            
            weak_baseline_mean = baseline_mean(:, weak_trials_indices);
            weak_baseline_std = baseline_std(:, weak_trials_indices);
            zscored_weak_baseline = (baseline_data(:, weak_trials_indices) - weak_baseline_mean) ./ (weak_baseline_std);

            zscored_weak_task = (weak_task_data - weak_baseline_mean) ./ (weak_baseline_std);
            weak_corr = corrcoef(zscored_weak_task);
            weak_corr_values = weak_corr(tril(true(size(weak_corr)), -1));
            weak_baseline_corr = corrcoef(zscored_weak_baseline);
            weak_baseline_corr_values = weak_baseline_corr(tril(true(size(weak_baseline_corr)), -1));
            
            [weak_p, ~, weak_stats] = ranksum(weak_baseline_corr_values, weak_corr_values);
            weak_combined_values = [weak_baseline_corr_values; weak_corr_values];
            weak_labels = [zeros(length(weak_baseline_corr_values), 1); ones(length(weak_corr_values), 1)];
            
            [sorted_values, sorted_indices] = sort(weak_combined_values, 'ascend');
            sorted_labels = weak_labels(sorted_indices);
            
            total_positives = sum(sorted_labels);
            total_negatives = length(sorted_labels) - total_positives;
            
            tpr = zeros(length(sorted_labels), 1);
            fpr = zeros(length(sorted_labels), 1);
            
            for i = 1:length(sorted_labels)
                tpr(i) = sum(sorted_labels(i:end)) / total_positives;
                fpr(i) = sum(~sorted_labels(i:end)) / total_negatives;
            end
            
            weak_auc = -trapz(fpr, tpr);

            horizontal_angle_weak = saccade_table_all_session(weak_trials_indices, :).azimuth_after - ...
            saccade_table_all_session(weak_trials_indices, :).azimuth_before;
            verticle_angle_weak = saccade_table_all_session(weak_trials_indices, :).elevation_after - ...
            saccade_table_all_session(weak_trials_indices, :).elevation_before;
        else
            weak_auc = nan;
            horizontal_angle_weak = nan;
            verticle_angle_weak = nan;
            weak_p = nan;
        end

        cluster1_trials = final_order(clusters == 1);
        cluster2_trials = final_order(clusters == 2);

        within_cluster1_corr = original_task_corr(cluster1_trials, cluster1_trials);
        within_cluster2_corr = original_task_corr(cluster2_trials, cluster2_trials);
        mean_within_cluster1_corr = mean(within_cluster1_corr(~eye(size(within_cluster1_corr))));
        mean_within_cluster2_corr = mean(within_cluster2_corr(~eye(size(within_cluster2_corr))));
        between_cluster_corr = original_task_corr(cluster1_trials, cluster2_trials);
        mean_between_cluster_corr = mean(between_cluster_corr(:));
        cluster_corr_difference = (mean_within_cluster1_corr + mean_within_cluster2_corr) / 2 - mean_between_cluster_corr;
        n_weak_corr_trials = length(weak_trials_indices);

        cluster1_task_data = task_data(:, cluster1_trials);
        cluster1_baseline_mean = baseline_mean(:, cluster1_trials);
        cluster1_baseline_std = baseline_std(:, cluster1_trials);
        zscored_cluster1_baseline = (baseline_data(:, cluster1_trials) - cluster1_baseline_mean) ./ (cluster1_baseline_std);

        zscored_cluster1_task = (cluster1_task_data - cluster1_baseline_mean) ./ (cluster1_baseline_std);
        cluster1_corr = corrcoef(zscored_cluster1_task);
        cluster1_corr_values = cluster1_corr(tril(true(size(cluster1_corr)), -1));
        cluster1_baseline_corr = corrcoef(zscored_cluster1_baseline);
        cluster1_baseline_corr_values = cluster1_baseline_corr(tril(true(size(cluster1_baseline_corr)), -1));
        
        [cluster1_p, ~, cluster1_stats] = ranksum(cluster1_baseline_corr_values, cluster1_corr_values);
        cluster1_combined_values = [cluster1_baseline_corr_values; cluster1_corr_values];
        cluster1_labels = [zeros(length(cluster1_baseline_corr_values), 1); ones(length(cluster1_corr_values), 1)];

        [sorted_values, sorted_indices] = sort(cluster1_combined_values, 'ascend');
        sorted_labels = cluster1_labels(sorted_indices);

        total_positives = sum(sorted_labels);
        total_negatives = length(sorted_labels) - total_positives;

        tpr = zeros(length(sorted_labels), 1);
        fpr = zeros(length(sorted_labels), 1);

        for i = 1:length(sorted_labels)
            tpr(i) = sum(sorted_labels(i:end)) / total_positives;
            fpr(i) = sum(~sorted_labels(i:end)) / total_negatives;
        end

        cluster1_auc = -trapz(fpr, tpr);
        norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster1_trials) - ...
            cluster1_baseline_mean) ./ (cluster1_baseline_std);
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        
        time_vector = (0:length(mean_erp) - 1) / sfreq - trial_begin_study;

        

        horizontal_angle_cluster1 = saccade_table_all_session(cluster1_trials, :).azimuth_after - ...
            saccade_table_all_session(cluster1_trials, :).azimuth_before;
        horizontal_angle_cluster2 = saccade_table_all_session(cluster2_trials, :).azimuth_after - ...
            saccade_table_all_session(cluster2_trials, :).azimuth_before;
        vertical_angle_cluster1 = saccade_table_all_session(cluster1_trials, :).elevation_after - ...
            saccade_table_all_session(cluster1_trials, :).elevation_before;
        vertical_angle_cluster2 = saccade_table_all_session(cluster2_trials, :).elevation_after - ...
            saccade_table_all_session(cluster2_trials, :).elevation_before;
        angles1 = atan2d(vertical_angle_cluster1, horizontal_angle_cluster1);
        angles2 = atan2d(vertical_angle_cluster2, horizontal_angle_cluster2);
        
        cluster2_task_data = task_data(:, cluster2_trials);
        cluster2_baseline_mean = baseline_mean(:, cluster2_trials);
        cluster2_baseline_std = baseline_std(:, cluster2_trials);
        zscored_cluster2_baseline = (baseline_data(:, cluster2_trials) - cluster2_baseline_mean) ./ (cluster2_baseline_std);

        zscored_cluster2_task = (cluster2_task_data - cluster2_baseline_mean) ./ (cluster2_baseline_std);
        cluster2_corr = corrcoef(zscored_cluster2_task);
        cluster2_corr_values = cluster2_corr(tril(true(size(cluster2_corr)), -1));
        cluster2_baseline_corr = corrcoef(zscored_cluster2_baseline);
        cluster2_baseline_corr_values = cluster2_baseline_corr(tril(true(size(cluster2_baseline_corr)), -1));
        
        [cluster2_p, ~, cluster2_stats] = ranksum(cluster2_baseline_corr_values, cluster2_corr_values);
        cluster2_combined_values = [cluster2_baseline_corr_values; cluster2_corr_values];
        cluster2_labels = [zeros(length(cluster2_baseline_corr_values), 1); ones(length(cluster2_corr_values), 1)];

        [sorted_values, sorted_indices] = sort(cluster2_combined_values, 'ascend');
        sorted_labels = cluster2_labels(sorted_indices);

        total_positives = sum(sorted_labels);
        total_negatives = length(sorted_labels) - total_positives;

        tpr = zeros(length(sorted_labels), 1);
        fpr = zeros(length(sorted_labels), 1);

        for i = 1:length(sorted_labels)
            tpr(i) = sum(sorted_labels(i:end)) / total_positives;
            fpr(i) = sum(~sorted_labels(i:end)) / total_negatives;
        end
        cluster2_auc = -trapz(fpr, tpr);


        norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster2_trials) - ...
            cluster2_baseline_mean) ./ (cluster2_baseline_std);
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        
        gaze_x_change_cluster1 = saccade_table_all_session(cluster1_trials, :).gaze_after_x_screen_coord - ...
            saccade_table_all_session(cluster1_trials, :).gaze_before_x_screen_coord;
        gaze_x_change_cluster2 = saccade_table_all_session(cluster2_trials, :).gaze_after_x_screen_coord - ...
            saccade_table_all_session(cluster2_trials, :).gaze_before_x_screen_coord;
        gaze_y_change_cluster1 = saccade_table_all_session(cluster1_trials, :).gaze_after_y_screen_coord - ...
            saccade_table_all_session(cluster1_trials, :).gaze_before_y_screen_coord;
        gaze_y_change_cluster2 = saccade_table_all_session(cluster2_trials, :).gaze_after_y_screen_coord - ...
            saccade_table_all_session(cluster2_trials, :).gaze_before_y_screen_coord;

        gaze_x_before_cluster1 = saccade_table_all_session(cluster1_trials, :).azimuth_before;
        gaze_x_before_cluster2 = saccade_table_all_session(cluster2_trials, :).azimuth_before;
        gaze_y_before_cluster1 = saccade_table_all_session(cluster1_trials, :).elevation_before;
        gaze_y_before_cluster2 = saccade_table_all_session(cluster2_trials, :).elevation_before;
        time_vector = (0:size(prep_signal_all_session, 1) - 1) / sfreq - trial_begin_study;
        norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster1_trials) - ...
            cluster1_baseline_mean) ./ cluster1_baseline_std;
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
       
        cluster1_data_oscillation = prep_signal_all_session(:, cluster1_trials);

        psd_all_trials = [];
        for trial = 1:length(cluster1_trials)
            [psd_trial, f] = pwelch(cluster1_data_oscillation(:, trial), hamming(window), noverlap, nfft, sfreq);
            psd_all_trials = [psd_all_trials, psd_trial];
        end
        
        mean_psd = mean(psd_all_trials, 2);
        std_psd = std(psd_all_trials, [], 2);

        
        fit_range_indices = f >= freq_range_fit(1) & f <= freq_range_fit(2);
        log_f = log10(f(fit_range_indices));
        log_psd = log10(mean_psd(fit_range_indices));
        p = polyfit(log_f, log_psd, 1);
        fit_line = 10.^(polyval(p, log10(f)));

        residuals = 10*log10(mean_psd(fit_range_indices)) - 10*log10(fit_line(fit_range_indices));
        outlier_indices = residuals > 0;
        if sum(outlier_indices) > length(outlier_indices) / 2
            [~, sorted_indices] = sort(residuals, 'descend');
            outlier_indices = false(size(residuals));
            outlier_indices(sorted_indices(1:floor(length(residuals) / 2))) = true;
        end
        fit_range_indices(fit_range_indices) = ~outlier_indices;

        log_f = log10(f(fit_range_indices));
        log_psd = log10(mean_psd(fit_range_indices));
        p = polyfit(log_f, log_psd, 1);
        fit_line = 10.^(polyval(p, log10(f)));
        
        
        % time- frequency
        total_power = [];
        num_trials = size(cluster1_data_oscillation, 2);
        for trial = 1:num_trials
            [S, F, T, P] = spectrogram(cluster1_data_oscillation(:, trial), hamming(window/16), ...
                noverlap / 16, nfft / 16, sfreq, 'yaxis');
            if isempty(total_power)
                total_power = zeros(size(P));
            end
            total_power = total_power + P;
        end
        mean_power = total_power / num_trials;
        baseline_indices = find(time_vector >= baseline_time_precedence(1) & time_vector <= erp_range(1));
        baseline_time_spec = time_vector(baseline_indices) + trial_begin_study;
        baseline_indices_spec = [];
        for i = 1:length(baseline_time_spec)
            [~, idx] = min(abs(T - baseline_time_spec(i)));
            if abs(T(idx) - baseline_time_spec(i)) <= 0.01
                baseline_indices_spec = [baseline_indices_spec, idx];
            end
        end
        baseline_indices_spec = unique(baseline_indices_spec);

        baseline_power = mean(mean_power(:, baseline_indices_spec), 2);
        baseline_corrected_power = 10*log10(mean_power ./ baseline_power);
        time_vector_spec = linspace(time_vector(1), time_vector(end), length(T));

       
        valid_idx = (freq_range_fit(2) > f) & (f > high_cutoff);
        f_filtered = f(valid_idx);
        mean_psd_filtered = 10*log10(mean_psd((valid_idx)));
        fit_line_filtered = 10*log10(fit_line(valid_idx));
        diff_psd = mean_psd_filtered - fit_line_filtered;
        [cluster1_oscillation_amp, peak_idx] = max(diff_psd);
        cluster1_peak_freq = f_filtered(peak_idx);
        low_cutoff_oscillation = max(0, cluster1_peak_freq - narrow_bandwidth);
        high_cutoff_oscillation = cluster1_peak_freq + narrow_bandwidth;

        task_indices_get_phase = find(time_vector >= erp_pre_during_post(1) & time_vector <= erp_pre_during_post(end));
        cluster1_data_oscillation = prep_signal_all_session(task_indices_get_phase, cluster1_trials);
        padding_samples = padding_duration * sfreq;

        amplitude = zeros(size(cluster1_data_oscillation));
        phase_cluster1 = zeros(size(cluster1_data_oscillation));
        [b, a] = butter(2, [low_cutoff_oscillation high_cutoff_oscillation] / (sfreq / 2)); % 2nd order Butterworth filter

        for trial = 1:size(cluster1_data_oscillation, 2)
            trial_data = cluster1_data_oscillation(:, trial);
            
            padded_data = [flipud(trial_data(1:padding_samples)); trial_data; flipud(trial_data(end-padding_samples+1:end))];
            
            filtered_data = filtfilt(b, a, padded_data);
            analytic_signal = hilbert(filtered_data);
            
            analytic_signal = analytic_signal(padding_samples+1:end-padding_samples);
            amplitude(:, trial) = abs(analytic_signal);
            phase_cluster1(:, trial) = angle(analytic_signal);
        end
        
        baseline_indices = ((0:size(amplitude, 1) - 1) / sfreq + erp_pre_during_post(1) ) >= baseline_time_precedence(1) & ...
            ((0:size(amplitude, 1) - 1) / sfreq + erp_pre_during_post(1)) <= baseline_time_precedence(2);
        baseline_amp = amplitude(baseline_indices, :);
        % baseline_mean = mean([baseline_data; task_data], 1);
        baseline_amp_mean = mean(baseline_amp, 1);
        baseline_amp_std = std(baseline_amp, [], 1);
        norm_amplitude = (amplitude - ...
            baseline_amp_mean) ./ baseline_amp_std;
        mean_erp = mean(norm_amplitude, 2);
        std_error = std(norm_amplitude, [], 2) / sqrt(size(norm_amplitude, 2));
        cluster1_4_12_hz_pow = mean_erp;

       

        mean_phase = circ_mean(phase_cluster1, [], 2); % Circular mean
        std_error_phase = circ_std(phase_cluster1, [], [], 2) / sqrt(size(phase_cluster1, 2)); % Circular standard error

        % Compute the PLV
        cluster1_4_12_hz_phase = mean_phase;

        plv = abs(mean(exp(1i * phase_cluster1), 2));

        cluster1_4_12_hz_plv = plv;
        cluster1_4_12_hz_plv_all = mean(plv);

       

        time_vector = (0:size(prep_signal_all_session, 1) - 1) / sfreq - trial_begin_study;
        norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster2_trials) - ...
            cluster2_baseline_mean) ./ (cluster2_baseline_std);
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
       
        
        cluster2_data_oscillation = prep_signal_all_session(:, cluster2_trials);

        psd_all_trials = [];
        for trial = 1:length(cluster2_trials)
            [psd_trial, f] = pwelch(cluster2_data_oscillation(:, trial), hamming(window), noverlap, nfft, sfreq);
            psd_all_trials = [psd_all_trials, psd_trial];
        end
        
        mean_psd = mean(psd_all_trials, 2);
        std_psd = std(psd_all_trials, [], 2);

        
        fit_range_indices = f >= freq_range_fit(1) & f <= freq_range_fit(2);
        log_f = log10(f(fit_range_indices));
        log_psd = log10(mean_psd(fit_range_indices));
        p = polyfit(log_f, log_psd, 1);
        fit_line = 10.^(polyval(p, log10(f)));

        % Identify oscillation
        residuals = 10*log10(mean_psd(fit_range_indices)) - 10*log10(fit_line(fit_range_indices));
        outlier_indices = residuals > 0;
        if sum(outlier_indices) > length(outlier_indices) / 2
            [~, sorted_indices] = sort(residuals, 'descend');
            outlier_indices = false(size(residuals));
            outlier_indices(sorted_indices(1:floor(length(residuals) / 2))) = true;
        end
        fit_range_indices(fit_range_indices) = ~outlier_indices;

        log_f = log10(f(fit_range_indices));
        log_psd = log10(mean_psd(fit_range_indices));
        p = polyfit(log_f, log_psd, 1);
        fit_line = 10.^(polyval(p, log10(f)));
        
        
        % time- frequency
        total_power = [];
        num_trials = size(cluster2_data_oscillation, 2);
        for trial = 1:num_trials
            [S, F, T, P] = spectrogram(cluster2_data_oscillation(:, trial), hamming(window/16), ...
                noverlap / 16, nfft / 16, sfreq, 'yaxis');
            if isempty(total_power)
                total_power = zeros(size(P));
            end
            total_power = total_power + P;
        end
        mean_power = total_power / num_trials;
        baseline_indices = find(time_vector >= baseline_time_precedence(1) & time_vector <= erp_range(1));
        baseline_time_spec = time_vector(baseline_indices) + trial_begin_study;
        baseline_indices_spec = [];
        for i = 1:length(baseline_time_spec)
            [~, idx] = min(abs(T - baseline_time_spec(i)));
            if abs(T(idx) - baseline_time_spec(i)) <= 0.01
                baseline_indices_spec = [baseline_indices_spec, idx];
            end
        end
        baseline_indices_spec = unique(baseline_indices_spec);

        baseline_power = mean(mean_power(:, baseline_indices_spec), 2);
        baseline_corrected_power = 10*log10(mean_power ./ baseline_power);
        time_vector_spec = linspace(time_vector(1), time_vector(end), length(T));

       
        valid_idx = (freq_range_fit(2) > f) & (f > high_cutoff);
        f_filtered = f(valid_idx);
        mean_psd_filtered = 10*log10(mean_psd((valid_idx)));
        fit_line_filtered = 10*log10(fit_line(valid_idx));
        diff_psd = mean_psd_filtered - fit_line_filtered;
        [cluster2_oscillation_amp, peak_idx] = max(diff_psd);
        cluster2_peak_freq = f_filtered(peak_idx);
        low_cutoff_oscillation = max(0, cluster2_peak_freq - narrow_bandwidth);
        high_cutoff_oscillation = cluster2_peak_freq + narrow_bandwidth;

        task_indices_get_phase = find(time_vector >= erp_pre_during_post(1) & time_vector <= erp_pre_during_post(end));
        cluster2_data_oscillation = prep_signal_all_session(task_indices_get_phase, cluster2_trials);
        padding_samples = padding_duration * sfreq;

        amplitude = zeros(size(cluster2_data_oscillation));
        phase_cluster2 = zeros(size(cluster2_data_oscillation));
        [b, a] = butter(2, [low_cutoff_oscillation high_cutoff_oscillation] / (sfreq / 2)); % 2nd order Butterworth filter

        for trial = 1:size(cluster2_data_oscillation, 2)
            trial_data = cluster2_data_oscillation(:, trial);
            
            padded_data = [flipud(trial_data(1:padding_samples)); trial_data; flipud(trial_data(end-padding_samples+1:end))];
            
            filtered_data = filtfilt(b, a, padded_data);
            analytic_signal = hilbert(filtered_data);
            
            analytic_signal = analytic_signal(padding_samples+1:end-padding_samples);
            amplitude(:, trial) = abs(analytic_signal);
            phase_cluster2(:, trial) = angle(analytic_signal);
        end
        
        baseline_indices = find(((0:size(amplitude, 1) - 1) / sfreq + erp_pre_during_post(1) ) >= baseline_time_precedence(1) & ...
            ((0:size(amplitude, 1) - 1) / sfreq + erp_pre_during_post(1)) <= baseline_time_precedence(2));
        baseline_amp = amplitude(baseline_indices, :);
        % baseline_mean = mean([baseline_data; task_data], 1);
        baseline_amp_mean = mean(baseline_amp, 1);
        baseline_amp_std = std(baseline_amp, [], 1);
        norm_amplitude = (amplitude - ...
            baseline_amp_mean) ./ baseline_amp_std;
        mean_erp = mean(norm_amplitude, 2);
        std_error = std(norm_amplitude, [], 2) / sqrt(size(norm_amplitude, 2));
        cluster2_4_12_hz_pow = mean_erp;

       
        mean_phase = circ_mean(phase_cluster2, [], 2); % Circular mean
        std_error_phase = circ_std(phase_cluster2, [], [], 2) / sqrt(size(phase_cluster2, 2)); % Circular standard error

        % Compute the PLV
        cluster2_4_12_hz_phase = mean_phase;

        plv = abs(mean(exp(1i * phase_cluster2), 2));

        cluster2_4_12_hz_plv = plv;
        cluster2_4_12_hz_plv_all = mean(plv);

        all_phase = [phase_cluster1'; phase_cluster2']';
        all_4_12_hz_plv = abs(mean(exp(1i * all_phase), 2));

       
      

        %% characterize erp
        norm_concatenated_prep_signal = (prep_signal_all_session - baseline_mean) ./ (baseline_std);
        time_vector = (0:size(norm_concatenated_prep_signal, 1) - 1) / sfreq - trial_begin_study;
        task_indices = find(time_vector >= erp_characteristic_range(1) & time_vector <= erp_characteristic_range(2));
        task_data_for_characteristic_cluster1 = norm_concatenated_prep_signal(task_indices, cluster1_trials);
        task_data_for_characteristic_cluster2 = norm_concatenated_prep_signal(task_indices, cluster2_trials);
        % flag bad erp
        cluster1_erp_good = 1;
        cluster2_erp_good = 1;
        
        t = time_vector(find(time_vector >= erp_characteristic_range(1), 1):find(time_vector <= erp_characteristic_range(2), 1, 'last'));
        mean_erp_cluster1 = mean(task_data_for_characteristic_cluster1, 2);
        first_derivative = diff(mean(task_data_for_characteristic_cluster1, 2));
        first_derivative = smoothdata(first_derivative, 'movmean', smoothing_window);
        [change_points, ~] = findchangepts(first_derivative, 'MaxNumChanges', 2);
        if length(change_points) ~= 2
            first_derivative = smoothdata(first_derivative, 'movmean', smoothing_window_secondary);
            [change_points, ~] = findchangepts(first_derivative, 'MaxNumChanges', 2);
        end
        if length(change_points) == 2
            change_times = t(change_points);
            window_indices = find(t >= (change_times(2) - erp_amp_window) & t <= (change_times(2) + erp_amp_window));
            [erp_amplitude, max_index] = max(abs(mean_erp_cluster1(window_indices)));
            erp_polarity = 'Positive';
            if mean_erp_cluster1(window_indices(max_index)) < 0
                erp_polarity = 'Negative';
            end
            cluster1_polarity = erp_polarity;
            cluster1_latency = change_times(1);
            cluster1_amp = mean_erp_cluster1(window_indices(max_index));
            cluster1_amp_time = t(window_indices(max_index));
            if cluster1_amp_time > cluster1_latency
                trough_window_indices = find(t <= cluster1_amp_time & t >= (cluster1_latency - erp_amp_window / 2));
                if strcmp(cluster1_polarity, 'Positive')
                    [trough_amplitude, trough_index] = min(mean_erp_cluster1(trough_window_indices));
                else
                    [trough_amplitude, trough_index] = max(mean_erp_cluster1(trough_window_indices));
                end
                cluster1_peak_to_trough = cluster1_amp - trough_amplitude;
                cluster1_trough_time = t(trough_window_indices(trough_index));
            else
                cluster1_peak_to_trough = nan;
                cluster1_trough_time = nan;
            end

            
            task_indices_get_phase = find(time_vector >= -trial_begin_study & time_vector <= (cluster1_latency));
            cluster1_data_oscillation = prep_signal_all_session(task_indices_get_phase, cluster1_trials);
            padding_samples = padding_duration * sfreq;

            no_erp_phase_cluster1 = zeros(size(cluster1_data_oscillation));
            low_cutoff_oscillation = max(0, cluster1_peak_freq - narrow_bandwidth);
            high_cutoff_oscillation = cluster1_peak_freq + narrow_bandwidth;

            [b, a] = butter(2, [low_cutoff_oscillation high_cutoff_oscillation] / (sfreq / 2)); 
            for trial = 1:size(cluster1_data_oscillation, 2)
                trial_data = cluster1_data_oscillation(:, trial);
                
                padded_data = [flipud(trial_data(1:padding_samples)); trial_data; flipud(trial_data(end-padding_samples+1:end))];
                
                filtered_data = filtfilt(b, a, padded_data);
                analytic_signal = hilbert(filtered_data);
                
                analytic_signal = analytic_signal(padding_samples+1:end-padding_samples);
                no_erp_phase_cluster1(:, trial) = angle(analytic_signal);
            end
            mean_phase = circ_mean(no_erp_phase_cluster1, [], 2); 
            std_error_phase = circ_std(no_erp_phase_cluster1, [], [], 2); 
            % / sqrt(size(no_erp_phase_cluster1, 2))
            cluster1_no_erp_phase = mean_phase;
        
            plv = abs(mean(exp(1i * no_erp_phase_cluster1), 2));
        
           
            time_4_12_hz = (0:size(cluster1_no_erp_phase, 1) - 1) / sfreq - trial_begin_study;
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (- 3 / 2 / cluster1_peak_freq)));
            if any(indices_for_pre_during_post)
                cluster1_no_erp_phase_pre_pre_pre = mean(cluster1_no_erp_phase(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_no_erp_phase_pre_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (- 1 / cluster1_peak_freq)));
            if any(indices_for_pre_during_post)
                cluster1_no_erp_phase_pre_pre = mean(cluster1_no_erp_phase(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_no_erp_phase_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - ( - 1 / cluster1_peak_freq / 2)));
            if any(indices_for_pre_during_post)
                cluster1_no_erp_phase_pre = mean(cluster1_no_erp_phase(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_no_erp_phase_pre = nan;
            end

            time_4_12_hz = (0:size(cluster1_4_12_hz_pow, 1) - 1) / sfreq + erp_pre_during_post(1);

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - ( - 1 / cluster1_peak_freq)));

            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_pre_pre_pre = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_pre_pre_pre = mean(cluster1_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_plv_pre_pre_pre = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_pre_pre_pre = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_pre_pre_pre = nan;
                cluster1_4_12_hz_phase_pre_pre_pre = nan;
                cluster1_4_12_hz_plv_pre_pre_pre = nan;
                all_4_12_hz_plv_pre_pre_pre = nan;
            end
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - ( - 1 / cluster1_peak_freq / 2)));
            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_pre_pre = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_pre_pre = mean(cluster1_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_plv_pre_pre = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_pre_pre = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_pre_pre = nan;
                cluster1_4_12_hz_phase_pre_pre = nan;
                cluster1_4_12_hz_plv_pre_pre = nan;
                all_4_12_hz_plv_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz ));
            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_pre = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_pre = mean(cluster1_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_plv_pre = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_pre = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_pre = nan;
                cluster1_4_12_hz_phase_pre = nan;
                cluster1_4_12_hz_plv_pre = nan;
                all_4_12_hz_plv_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency + cluster1_amp_time) / 2));

            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_during = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_during = mean(cluster1_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_plv_during = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_during = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_during = nan;
                cluster1_4_12_hz_phase_during = nan;
                cluster1_4_12_hz_plv_during = nan;
                all_4_12_hz_plv_during = nan;
            end
            
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - ( 1 / cluster1_peak_freq)));

            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_post = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_post = mean(cluster1_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_plv_post = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_post = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_post = nan;
                cluster1_4_12_hz_phase_post = nan;
                cluster1_4_12_hz_plv_post = nan;
                all_4_12_hz_plv_post = nan;
            end
            
        else
            cluster1_no_erp_phase_pre_pre_pre = nan;
            cluster1_no_erp_phase_pre_pre = nan;
            cluster1_no_erp_phase_pre = nan;
            cluster1_polarity = nan;
            cluster1_latency = nan;
            cluster1_amp = nan;
            cluster1_amp_time = nan;
            cluster1_erp_good = 0;
            cluster1_peak_to_trough = nan;
            cluster1_trough_time = nan;
            time_4_12_hz = (0:size(cluster1_4_12_hz_pow, 1) - 1) / sfreq + erp_pre_during_post(1);

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (- 1 / cluster1_peak_freq)));


            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_pre_pre_pre = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_pre_pre_pre = mean(cluster1_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_plv_pre_pre_pre = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_pre_pre_pre = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_pre_pre_pre = nan;
                cluster1_4_12_hz_phase_pre_pre_pre = nan;
                cluster1_4_12_hz_plv_pre_pre_pre = nan;
                all_4_12_hz_plv_pre_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (- 1 / cluster1_peak_freq / 2)));

            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_pre_pre = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_pre_pre = mean(cluster1_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_plv_pre_pre = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_pre_pre = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_pre_pre = nan;
                cluster1_4_12_hz_phase_pre_pre = nan;
                cluster1_4_12_hz_plv_pre_pre = nan;
                all_4_12_hz_plv_pre_pre = nan;
            end
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz ));
            % indices_for_pre_during_post = (time_4_12_hz >= (erp_pre_during_post(3) + erp_pre_during_post(2))) & ...
                % (time_4_12_hz < (erp_pre_during_post(4) + erp_pre_during_post(2)));
            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_pre = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_pre = angle(mean(exp(1i * cluster1_4_12_hz_phase(indices_for_pre_during_post)), 'omitnan'));
                cluster1_4_12_hz_plv_pre = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_pre = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_pre = nan;
                cluster1_4_12_hz_phase_pre = nan;
                cluster1_4_12_hz_plv_pre = nan;
                all_4_12_hz_plv_pre = nan;
            end
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (erp_pre_during_post(3) + erp_pre_during_post(4))/ 2));

            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_during = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_during = angle(mean(exp(1i * cluster1_4_12_hz_phase(indices_for_pre_during_post)), 'omitnan'));
                cluster1_4_12_hz_plv_during = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_during = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_during = nan;
                cluster1_4_12_hz_phase_during = nan;
                cluster1_4_12_hz_plv_during = nan;
                all_4_12_hz_plv_during = nan;
            end
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (1 / cluster1_peak_freq)));


            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_post = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_post = angle(mean(exp(1i * cluster1_4_12_hz_phase(indices_for_pre_during_post)), 'omitnan'));
                cluster1_4_12_hz_plv_post = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_post = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_post = nan;
                cluster1_4_12_hz_phase_post = nan;
                cluster1_4_12_hz_plv_post = nan;
                all_4_12_hz_plv_post = nan;
            end
        end


        mean_erp_cluster2 = mean(task_data_for_characteristic_cluster2, 2);
       
        first_derivative = diff(mean_erp_cluster2);
        first_derivative = smoothdata(first_derivative, 'movmean', smoothing_window);
       
        [change_points, ~] = findchangepts(first_derivative, 'MaxNumChanges', 2);
        if length(change_points) ~= 2
            first_derivative = smoothdata(first_derivative, 'movmean', smoothing_window_secondary);
            [change_points, ~] = findchangepts(first_derivative, 'MaxNumChanges', 2);
        end
        if length(change_points) == 2
            change_times = t(change_points);
            window_indices = find(t >= (change_times(2) - erp_amp_window) & t <= (change_times(2) + erp_amp_window));
            [erp_amplitude, max_index] = max(abs(mean_erp_cluster2(window_indices)));
            erp_polarity = 'Positive';
            if mean_erp_cluster2(window_indices(max_index)) < 0
                erp_polarity = 'Negative';
            end
            cluster2_polarity = erp_polarity;
            cluster2_latency = change_times(1);
            cluster2_amp = mean_erp_cluster2(window_indices(max_index));
            cluster2_amp_time = t(window_indices(max_index));
            if cluster2_amp_time > cluster2_latency
                trough_window_indices = find(t <= cluster2_amp_time & t >= (cluster2_latency - erp_amp_window / 2));
                if strcmp(cluster2_polarity, 'Positive')
                    [trough_amplitude, trough_index] = min(mean_erp_cluster2(trough_window_indices));
                else
                    [trough_amplitude, trough_index] = max(mean_erp_cluster2(trough_window_indices)); 
                end
                cluster2_peak_to_trough = cluster2_amp - trough_amplitude;
                cluster2_trough_time = t(trough_window_indices(trough_index));
            else
                cluster2_peak_to_trough = nan;
                cluster2_trough_time = nan;
            end

           

            task_indices_get_phase = find(time_vector >= -trial_begin_study & time_vector <= (cluster2_latency));
            cluster2_data_oscillation = prep_signal_all_session(task_indices_get_phase, cluster2_trials);
            padding_samples = padding_duration * sfreq;

            no_erp_phase_cluster2 = zeros(size(cluster2_data_oscillation));
            low_cutoff_oscillation = max(0, cluster2_peak_freq - narrow_bandwidth);
            high_cutoff_oscillation = cluster2_peak_freq + narrow_bandwidth;

            [b, a] = butter(2, [low_cutoff_oscillation high_cutoff_oscillation] / (sfreq / 2)); 
            for trial = 1:size(cluster2_data_oscillation, 2)
                trial_data = cluster2_data_oscillation(:, trial);
                
                padded_data = [flipud(trial_data(1:padding_samples)); trial_data; flipud(trial_data(end-padding_samples+1:end))];
                
                filtered_data = filtfilt(b, a, padded_data);
                analytic_signal = hilbert(filtered_data);
                
                analytic_signal = analytic_signal(padding_samples+1:end-padding_samples);
                no_erp_phase_cluster2(:, trial) = angle(analytic_signal);
            end
            mean_phase = circ_mean(no_erp_phase_cluster2, [], 2); 
            std_error_phase = circ_std(no_erp_phase_cluster2, [], [], 2); 
            cluster2_no_erp_phase = mean_phase;
        
            plv = abs(mean(exp(1i * no_erp_phase_cluster2), 2));
        
            time_4_12_hz = (0:size(cluster2_no_erp_phase, 1) - 1) / sfreq - trial_begin_study;
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency - 1 / cluster2_peak_freq)));
            if any(indices_for_pre_during_post)
                cluster2_no_erp_phase_pre_pre_pre = mean(cluster2_no_erp_phase(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_no_erp_phase_pre_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency - 1 / cluster2_peak_freq / 2)));
            if any(indices_for_pre_during_post)
                cluster2_no_erp_phase_pre_pre = mean(cluster2_no_erp_phase(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_no_erp_phase_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency)));
            if any(indices_for_pre_during_post)
                cluster2_no_erp_phase_pre = mean(cluster2_no_erp_phase(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_no_erp_phase_pre = nan;
            end

            time_4_12_hz = (0:size(cluster2_4_12_hz_pow, 1) - 1) / sfreq + erp_pre_during_post(1);
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency - 1 / cluster2_peak_freq)));

            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_pre_pre_pre = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_pre_pre_pre = mean(cluster2_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_plv_pre_pre_pre = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_pre_pre_pre = nan;
                cluster2_4_12_hz_phase_pre_pre_pre = nan;
                cluster2_4_12_hz_plv_pre_pre_pre = nan;
            end


            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency - 1 / cluster2_peak_freq / 2)));

         
            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_pre_pre = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_pre_pre = mean(cluster2_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_plv_pre_pre = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_pre_pre = nan;
                cluster2_4_12_hz_phase_pre_pre = nan;
                cluster2_4_12_hz_plv_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency)));
         
            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_pre = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_pre = angle(mean(exp(1i * cluster2_4_12_hz_phase(indices_for_pre_during_post)), 'omitnan'));
                cluster2_4_12_hz_plv_pre = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_pre = nan;
                cluster2_4_12_hz_phase_pre = nan;
                cluster2_4_12_hz_plv_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_amp_time) ));

            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_during = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_during = angle(mean(exp(1i * cluster2_4_12_hz_phase(indices_for_pre_during_post)), 'omitnan'));
                cluster2_4_12_hz_plv_during = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_during = nan;
                cluster2_4_12_hz_phase_during = nan;
                cluster2_4_12_hz_plv_during = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_amp_time + 1 / cluster2_peak_freq)));

          
            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_post = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_post = angle(mean(exp(1i * cluster2_4_12_hz_phase(indices_for_pre_during_post)), 'omitnan'));
                cluster2_4_12_hz_plv_post = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_post = nan;
                cluster2_4_12_hz_phase_post = nan;
                cluster2_4_12_hz_plv_post = nan;
            end
            
        else
            cluster2_no_erp_phase_pre_pre_pre = nan;
            cluster2_no_erp_phase_pre_pre = nan;
            cluster2_no_erp_phase_pre = nan;
            cluster2_polarity = nan;
            cluster2_latency = nan;
            cluster2_amp = nan;
            cluster2_amp_time = nan;
            cluster2_erp_good = 0;
            cluster2_peak_to_trough = nan;
            cluster2_trough_time = nan;

            time_4_12_hz = (0:size(cluster2_4_12_hz_pow, 1) - 1) / sfreq + erp_pre_during_post(1);
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (- 1 / cluster2_peak_freq)));

            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_pre_pre_pre = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_pre_pre_pre = mean(cluster2_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_plv_pre_pre_pre = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_pre_pre_pre = nan;
                cluster2_4_12_hz_phase_pre_pre_pre = nan;
                cluster2_4_12_hz_plv_pre_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - ( - 1 / cluster2_peak_freq / 2)));

            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_pre_pre = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_pre_pre = mean(cluster2_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_plv_pre_pre = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_pre_pre = nan;
                cluster2_4_12_hz_phase_pre_pre = nan;
                cluster2_4_12_hz_plv_pre_pre = nan;
            end
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz ));


            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_pre = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_pre = mean(cluster2_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_plv_pre = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_pre = nan;
                cluster2_4_12_hz_phase_pre = nan;
                cluster2_4_12_hz_plv_pre = nan;
            end
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (erp_pre_during_post(3) + erp_pre_during_post(4)) / 2));
            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_during = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_during = mean(cluster2_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_plv_during = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_during = nan;
                cluster2_4_12_hz_phase_during = nan;
                cluster2_4_12_hz_plv_during = nan;
            end
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - ( 1 / cluster2_peak_freq)));


            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_post = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_post = mean(cluster2_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_plv_post = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_post = nan;
                cluster2_4_12_hz_phase_post = nan;
                cluster2_4_12_hz_plv_post = nan;
            end
        end




        angles_rad = deg2rad(angles1);
        mean_sin = mean(sin(angles_rad), 'omitmissing');
        mean_cos = mean(cos(angles_rad), 'omitmissing');
        mean_angle_rad = atan2(mean_sin, mean_cos);
        mean_angle_deg = rad2deg(mean_angle_rad);
        mean_angle_deg1 = mod(mean_angle_deg, 360);
        if mean_angle_deg1 > 180
            mean_angle_deg1 = mean_angle_deg1 - 360;
        elseif mean_angle_deg1 <= -180
            mean_angle_deg1 = mean_angle_deg1 + 360;
        end

        R = sqrt(mean_cos^2 + mean_sin^2);
        circular_variance1 = 1 - R;

        angles_rad = deg2rad(angles2);
        mean_sin = mean(sin(angles_rad), 'omitmissing');
        mean_cos = mean(cos(angles_rad), 'omitmissing');
        mean_angle_rad = atan2(mean_sin, mean_cos);
        mean_angle_deg = rad2deg(mean_angle_rad);
        mean_angle_deg2 = mod(mean_angle_deg, 360);
        if mean_angle_deg2 > 180
            mean_angle_deg2 = mean_angle_deg2 - 360;
        elseif mean_angle_deg2 <= -180
            mean_angle_deg2 = mean_angle_deg2 + 360;
        end
        R = sqrt(mean_cos^2 + mean_sin^2);
        circular_variance2 = 1 - R;

        angles_rad = deg2rad([angles1; angles2]);
        mean_sin = mean(sin(angles_rad), 'omitmissing');
        mean_cos = mean(cos(angles_rad), 'omitmissing');


        R = sqrt(mean_cos^2 + mean_sin^2);
        circular_variance_all = 1 - R;

        new_row = {subj, task, contact_this_subj.ch_name_bci2000_format(i_ch), p_value < alpha_test, ...
            auc, p_value, cluster_corr_difference, mean_within_cluster1_corr, mean_within_cluster2_corr,...
            cluster_corr_difference_bl_control, mean_within_cluster1_corr_bl_control, mean_within_cluster2_corr_bl_control,...
            n_weak_corr_trials, n_weak_corr_trials_bl_control,...
            circular_variance_all, ...
            all_4_12_hz_plv_pre_pre_pre, all_4_12_hz_plv_pre_pre, all_4_12_hz_plv_pre, ...
            all_4_12_hz_plv_during, all_4_12_hz_plv_post,...
            cluster1_peak_freq, cluster1_oscillation_amp, cluster1_4_12_hz_pow_pre_pre_pre, cluster1_4_12_hz_pow_pre_pre, ...
            cluster1_4_12_hz_pow_pre, cluster1_4_12_hz_pow_during, cluster1_4_12_hz_pow_post,...
            cluster1_no_erp_phase_pre_pre_pre, cluster1_no_erp_phase_pre_pre, cluster1_no_erp_phase_pre,...
            cluster1_4_12_hz_phase_pre_pre_pre, cluster1_4_12_hz_phase_pre_pre, cluster1_4_12_hz_phase_pre,...
            cluster1_4_12_hz_phase_during, cluster1_4_12_hz_phase_post,...
            cluster1_4_12_hz_plv_all, cluster1_4_12_hz_plv_pre_pre_pre, cluster1_4_12_hz_plv_pre_pre, ...
            cluster1_4_12_hz_plv_pre, cluster1_4_12_hz_plv_during, cluster1_4_12_hz_plv_post,...
            cluster2_peak_freq, cluster2_oscillation_amp, cluster2_4_12_hz_pow_pre_pre_pre, cluster2_4_12_hz_pow_pre_pre,...
            cluster2_4_12_hz_pow_pre, cluster2_4_12_hz_pow_during, cluster2_4_12_hz_pow_post,...
            cluster2_no_erp_phase_pre_pre_pre, cluster2_no_erp_phase_pre_pre, cluster2_no_erp_phase_pre,...
            cluster2_4_12_hz_phase_pre_pre_pre, cluster2_4_12_hz_phase_pre_pre, cluster2_4_12_hz_phase_pre, cluster2_4_12_hz_phase_during, cluster2_4_12_hz_phase_post,...
            cluster2_4_12_hz_plv_all, cluster2_4_12_hz_plv_pre_pre_pre, cluster2_4_12_hz_plv_pre_pre, ...
            cluster2_4_12_hz_plv_pre, cluster2_4_12_hz_plv_during, cluster2_4_12_hz_plv_post,...
            cluster1_p < alpha_test, cluster1_auc, cluster1_p, cluster1_auc_bl_control, cluster1_p_bl_control,...
            cluster1_latency, cluster1_amp, cluster1_amp_time, cluster1_polarity, cluster1_peak_to_trough, cluster1_trough_time, ...
            cluster2_p < alpha_test, cluster2_auc, cluster2_p, cluster2_auc_bl_control, cluster2_p_bl_control,...
            cluster2_latency, cluster2_amp,...
           cluster2_amp_time, cluster2_polarity, cluster2_peak_to_trough, cluster2_trough_time,...
           mean(horizontal_angle_cluster1, 'omitmissing'), mean(vertical_angle_cluster1, 'omitmissing'), ...
           mean_angle_deg1, circular_variance1, ...
           mean(gaze_x_change_cluster1, 'omitmissing'),  mean(gaze_y_change_cluster1, 'omitmissing'), ...
           mean(horizontal_angle_cluster2, 'omitmissing'), mean(vertical_angle_cluster2, 'omitmissing'), ...
           mean_angle_deg2, circular_variance2, ...
           mean(gaze_x_change_cluster2, 'omitmissing'), mean(gaze_y_change_cluster2, 'omitmissing'), ...
           re_ref_vector_x_norm, re_ref_vector_y_norm, re_ref_vector_z_norm, ...
           current_point_norm(1), current_point_norm(2), current_point_norm(3), reference_subj_left_eye(1),...
           reference_subj_left_eye(2), reference_subj_left_eye(3), reference_subj_right_eye(1),...
           reference_subj_right_eye(2), reference_subj_right_eye(3),...
           valid_side, contact_side, contact_dist_to_left_eye, contact_dist_to_right_eye, shank_ref, ...
           mean(horizontal_angle_weak, 'omitmissing'), ...
           mean(verticle_angle_weak, 'omitmissing')};
        responsive_erp_table = [responsive_erp_table; new_row];

        %% get trial property
        for i_trial = 1:size(prep_signal_all_session, 2)
            current_trial_task = norm_concatenated_prep_signal(task_indices, i_trial);
            % which cluster does it belong to
            
            if ismember(i_trial, cluster1_trials) && cluster1_erp_good
                first_derivative = diff(current_trial_task);
                first_derivative = smoothdata(first_derivative, 'movmean', smoothing_window);
                latency_window_indices = find(t >= (cluster1_latency - erp_amp_window) & ...
                    t <= (cluster1_latency + erp_amp_window));
                if latency_window_indices(end) > length(first_derivative)
                    shift_amount = latency_window_indices(end) - length(first_derivative);
                    latency_window_indices = latency_window_indices - shift_amount;
                    disp(['ch ' num2str(i_ch) 'trial ' num2str(i_trial) 'latency_window shifts'])
                end
                [~, min_index] = min(abs(t(latency_window_indices) - cluster1_latency));
                search_window_first_derivative = first_derivative(latency_window_indices);
                [change_point_index, ~] = findchangepts(search_window_first_derivative, 'MaxNumChanges', 1);
                if ~isempty(change_point_index)
                    trial_latency = t(latency_window_indices(change_point_index));
                else
                    trial_latency = cluster1_latency;
                end
                amp_window_indices = find(t >= (cluster1_amp_time - erp_amp_window) & ...
                    t <= (cluster1_amp_time + erp_amp_window));
                if amp_window_indices(end) > length(current_trial_task)
                    shift_amount = amp_window_indices(end) - length(current_trial_task);
                    amp_window_indices = amp_window_indices - shift_amount;
                    disp(['ch ' num2str(i_ch) 'trial ' num2str(i_trial) 'amp_window shifts'])
                end
                if strcmp(cluster1_polarity, 'Positive')
                    [erp_amplitude, max_index] = max(current_trial_task(amp_window_indices));
                    erp_polarity = 'Positive';
                    trial_amp_time = t(amp_window_indices(max_index));
                    erp_peak = current_trial_task(amp_window_indices(max_index));
                    if trial_amp_time > trial_latency
                        trough_window_indices = find(t <= trial_amp_time & t >= (trial_latency - erp_amp_window / 2));
                        [trough_amplitude, trough_index] = min(current_trial_task(trough_window_indices));
                        peak_to_trough = erp_peak - trough_amplitude;
                        trial_trough_time = t(trough_window_indices(trough_index));
                    else
                        peak_to_trough = nan;
                        trial_trough_time = nan;
                    end
                else
                    [erp_amplitude, max_index] = min(current_trial_task(amp_window_indices));
                    erp_polarity = 'Negative';
                    trial_amp_time = t(amp_window_indices(max_index));
                    erp_peak = current_trial_task(amp_window_indices(max_index));
                    if trial_amp_time > trial_latency
                        trough_window_indices = find(t <= trial_amp_time & t >= (trial_latency - erp_amp_window / 2));
                        [trough_amplitude, trough_index] = max(current_trial_task(trough_window_indices));
                        peak_to_trough = erp_peak - trough_amplitude;
                        trial_trough_time = t(trough_window_indices(trough_index));
                    else
                        peak_to_trough = nan;
                        trial_trough_time = nan;
                    end
                end

                time_vector = (0:size(prep_signal_all_session, 1) - 1) / sfreq - trial_begin_study;
                task_indices_get_phase = find(time_vector >= - trial_begin_study & time_vector <= erp_pre_during_post(end));
                trial_data_oscillation = prep_signal_all_session(task_indices_get_phase, i_trial);
                padding_samples = padding_duration * sfreq;
                low_cutoff_oscillation = max(0, cluster1_peak_freq - narrow_bandwidth);
                high_cutoff_oscillation = cluster1_peak_freq + narrow_bandwidth;

                [b, a] = butter(2, [low_cutoff_oscillation high_cutoff_oscillation] / (sfreq / 2));
                padded_data = [flipud(trial_data_oscillation(1:padding_samples)); trial_data_oscillation;...
                    flipud(trial_data_oscillation(end-padding_samples+1:end))];
            
                filtered_data = filtfilt(b, a, padded_data);
                analytic_signal = hilbert(filtered_data);
            
                analytic_signal = analytic_signal(padding_samples+1:end-padding_samples);
                amplitude = abs(analytic_signal);
                time_4_12_hz = (0:size(amplitude, 1) - 1) / sfreq - trial_begin_study;
                baseline_indices = find((time_4_12_hz >= baseline_time_precedence(1)) & ...
                    (time_4_12_hz <= baseline_time_precedence(2)));
                baseline_amp = amplitude(baseline_indices, :);
                baseline_amp_mean = mean(baseline_amp, 1);
                baseline_amp_std = std(baseline_amp, [], 1);
                norm_amplitude = (amplitude - ...
                    baseline_amp_mean) ./ baseline_amp_std;
                phase = angle(analytic_signal);

                task_indices_get_phase_no_erp = find(time_vector >= - trial_begin_study & time_vector <= (cluster1_latency));
                trial_data_oscillation = prep_signal_all_session(task_indices_get_phase_no_erp, i_trial);
                padding_samples = padding_duration / 2 * sfreq;

                padded_data = [flipud(trial_data_oscillation(1:padding_samples)); trial_data_oscillation; flipud(trial_data_oscillation(end-padding_samples+1:end))];
                
                filtered_data = filtfilt(b, a, padded_data);
                analytic_signal = hilbert(filtered_data);    
                analytic_signal = analytic_signal(padding_samples+1:end-padding_samples);
                no_erp_phase= angle(analytic_signal);

                time_4_12_hz = (0:size(no_erp_phase, 1) - 1) / sfreq - trial_begin_study;
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency - 1 / cluster1_peak_freq)));
                if any(indices_for_pre_during_post)
                    trial_no_erp_phase_pre_pre_pre = angle(mean(exp(1i * no_erp_phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_no_erp_phase_pre_pre_pre = nan;
                end
    
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency - 1 / cluster1_peak_freq / 2)));
                if any(indices_for_pre_during_post)
                    trial_no_erp_phase_pre_pre = angle(mean(exp(1i * no_erp_phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_no_erp_phase_pre_pre = nan;
                end
    
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency)));
                if any(indices_for_pre_during_post)
                    trial_no_erp_phase_pre = angle(mean(exp(1i * no_erp_phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_no_erp_phase_pre = nan;
                end

                time_4_12_hz = (0:size(amplitude, 1) - 1) / sfreq + erp_pre_during_post(1);
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency - 1 / cluster1_peak_freq)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_pre_pre_pre = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_pre_pre_pre = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_pre_pre_pre = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_pre_pre_pre = nan;
                    trial_4_12_hz_pow_pre_pre_pre = nan;
                    trial_4_12_hz_phase_pre_pre_pre = nan;
                end
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency - 1 / cluster1_peak_freq / 2)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_pre_pre = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_pre_pre = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_pre_pre = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_pre_pre = nan;
                    trial_4_12_hz_pow_pre_pre = nan;
                    trial_4_12_hz_phase_pre_pre = nan;
                end
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_pre = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_pre = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_pre = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_pre = nan;
                    trial_4_12_hz_pow_pre = nan;
                    trial_4_12_hz_phase_pre = nan;
                end
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_amp_time)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_during = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_during = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_during = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_during = nan;
                    trial_4_12_hz_pow_during = nan;
                    trial_4_12_hz_phase_during = nan;
                end
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_amp_time + 1 / cluster1_peak_freq / 2)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_post = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_post = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_post = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_post = nan;
                    trial_4_12_hz_pow_post = nan;
                    trial_4_12_hz_phase_post = nan;
                end

                
                new_row = {i_trial, subj, task, contact_this_subj.ch_name_bci2000_format(i_ch), ...
                    saccade_table_all_session.session(i_trial), saccade_table_all_session.trial(i_trial) ...
                    saccade_table_all_session.eccentricity(i_trial), saccade_table_all_session.azimuth_before(i_trial), ...
                           saccade_table_all_session.elevation_before(i_trial),saccade_table_all_session.azimuth_after(i_trial), ...
                           saccade_table_all_session.elevation_after(i_trial), saccade_table_all_session.gaze_before_x_screen_coord(i_trial), ...
                           saccade_table_all_session.gaze_before_y_screen_coord(i_trial),saccade_table_all_session.gaze_after_x_screen_coord(i_trial), ...
                           saccade_table_all_session.gaze_after_x_screen_coord(i_trial),saccade_table_all_session.duration(i_trial),...
                           saccade_table_all_session.time_to_image_onset(i_trial), saccade_table_all_session.behavior_condition(i_trial),...
                           saccade_table_all_session.remember_condition(i_trial), corr(mean_erp_cluster1, current_trial_task, 'Type', 'Spearman'), ...
                           peak_to_trough, trial_trough_time, ...
                           trial_latency, erp_amplitude, trial_amp_time, erp_polarity...
                           trial_4_12_hz_pow_pre_pre_pre,trial_4_12_hz_pow_pre_pre, trial_4_12_hz_pow_pre, trial_4_12_hz_pow_during, trial_4_12_hz_pow_post,...
                           trial_4_12_hz_phase_pre_pre_pre, trial_4_12_hz_phase_pre_pre, trial_4_12_hz_phase_pre, trial_4_12_hz_phase_during, trial_4_12_hz_phase_post,...
                           trial_4_12_hz_norm_pow_pre_pre_pre, trial_4_12_hz_norm_pow_pre_pre, trial_4_12_hz_norm_pow_pre, trial_4_12_hz_norm_pow_during, trial_4_12_hz_norm_pow_post,...
                           trial_no_erp_phase_pre_pre_pre, trial_no_erp_phase_pre_pre, trial_no_erp_phase_pre};

                df_evoked_potential = [df_evoked_potential; new_row];
            elseif ismember(i_trial, cluster2_trials) && cluster2_erp_good
                first_derivative = diff(current_trial_task);
                first_derivative = smoothdata(first_derivative, 'movmean', smoothing_window);

                latency_window_indices = find(t >= (cluster2_latency - erp_amp_window) & ...
                    t <= (cluster2_latency + erp_amp_window));
                if latency_window_indices(end) > length(first_derivative)
                    shift_amount = latency_window_indices(end) - length(first_derivative);
                    latency_window_indices = latency_window_indices - shift_amount;
                    disp(['ch ' num2str(i_ch) 'trial ' num2str(i_trial) 'latency_window shifts'])
                end
                [~, min_index] = min(abs(t(latency_window_indices) - cluster2_latency));
                search_window_first_derivative = first_derivative(latency_window_indices);
                [change_point_index, ~] = findchangepts(search_window_first_derivative, 'MaxNumChanges', 1);
                if ~isempty(change_point_index)
                    trial_latency = t(latency_window_indices(change_point_index));
                else
                    trial_latency = cluster2_latency;
                end
                amp_window_indices = find(t >= (cluster2_amp_time - erp_amp_window) & ...
                    t <= (cluster2_amp_time + erp_amp_window));
                if amp_window_indices(end) > length(current_trial_task)
                    shift_amount = amp_window_indices(end) - length(current_trial_task);
                    amp_window_indices = amp_window_indices - shift_amount;
                    disp(['ch ' num2str(i_ch) 'trial ' num2str(i_trial) 'amp_window shifts'])
                end
                if strcmp(cluster2_polarity, 'Positive')
                    [erp_amplitude, max_index] = max(current_trial_task(amp_window_indices));
                    erp_polarity = 'Positive';
                    trial_amp_time = t(amp_window_indices(max_index));
                    erp_peak = current_trial_task(amp_window_indices(max_index));
                    if trial_amp_time > trial_latency
                        trough_window_indices = find(t <= trial_amp_time & t >= (trial_latency - erp_amp_window / 2));
                        [trough_amplitude, trough_index] = min(current_trial_task(trough_window_indices));
                        peak_to_trough = erp_peak - trough_amplitude;
                        trial_trough_time = t(trough_window_indices(trough_index));
                    else
                        peak_to_trough = nan;
                        trial_trough_time = nan;
                    end
                else
                    [erp_amplitude, max_index] = min(current_trial_task(amp_window_indices));
                    erp_polarity = 'Negative';
                    trial_amp_time = t(amp_window_indices(max_index));
                    erp_peak = current_trial_task(amp_window_indices(max_index));
                    if trial_amp_time > trial_latency
                        trough_window_indices = find(t <= trial_amp_time & t >= (trial_latency - erp_amp_window / 2));
                        [trough_amplitude, trough_index] = max(current_trial_task(trough_window_indices));
                        peak_to_trough = erp_peak - trough_amplitude;
                        trial_trough_time = t(trough_window_indices(trough_index));
                    else
                        peak_to_trough = nan;
                        trial_trough_time = nan;
                    end
                end

                time_vector = (0:size(prep_signal_all_session, 1) - 1) / sfreq - trial_begin_study;
                task_indices_get_phase = find(time_vector >= - trial_begin_study & time_vector <= erp_pre_during_post(end));
                trial_data_oscillation = prep_signal_all_session(task_indices_get_phase, i_trial);
                padding_samples = padding_duration * sfreq;
                low_cutoff_oscillation = max(0, cluster2_peak_freq - narrow_bandwidth);
                high_cutoff_oscillation = cluster2_peak_freq + narrow_bandwidth;

                [b, a] = butter(2, [low_cutoff_oscillation high_cutoff_oscillation] / (sfreq / 2));
                padded_data = [flipud(trial_data_oscillation(1:padding_samples)); trial_data_oscillation;...
                    flipud(trial_data_oscillation(end-padding_samples+1:end))];
            
                filtered_data = filtfilt(b, a, padded_data);
                analytic_signal = hilbert(filtered_data);
            
                analytic_signal = analytic_signal(padding_samples+1:end-padding_samples);
                amplitude = abs(analytic_signal);
                time_4_12_hz = (0:size(amplitude, 1) - 1) / sfreq - trial_begin_study;
                baseline_indices = find((time_4_12_hz >= baseline_time_precedence(1)) & ...
                    (time_4_12_hz <= baseline_time_precedence(2)));
                baseline_amp = amplitude(baseline_indices, :);
                baseline_amp_mean = mean(baseline_amp, 1);
                baseline_amp_std = std(baseline_amp, [], 1);
                norm_amplitude = (amplitude - ...
                    baseline_amp_mean) ./ baseline_amp_std;
                phase = angle(analytic_signal);

                task_indices_get_phase_no_erp = find(time_vector >= - trial_begin_study & time_vector <= (cluster2_latency));
                trial_data_oscillation = prep_signal_all_session(task_indices_get_phase_no_erp, i_trial);
                padding_samples = padding_duration / 2 * sfreq;

                padded_data = [flipud(trial_data_oscillation(1:padding_samples)); trial_data_oscillation; flipud(trial_data_oscillation(end-padding_samples+1:end))];
                
                filtered_data = filtfilt(b, a, padded_data);
                analytic_signal = hilbert(filtered_data);    
                analytic_signal = analytic_signal(padding_samples+1:end-padding_samples);
                no_erp_phase= angle(analytic_signal);

                time_4_12_hz = (0:size(no_erp_phase, 1) - 1) / sfreq - trial_begin_study;
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency - 1 / cluster2_peak_freq)));
                if any(indices_for_pre_during_post)
                    trial_no_erp_phase_pre_pre_pre = angle(mean(exp(1i * no_erp_phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_no_erp_phase_pre_pre_pre = nan;
                end
    
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency - 1 / cluster2_peak_freq / 2)));
                if any(indices_for_pre_during_post)
                    trial_no_erp_phase_pre_pre = angle(mean(exp(1i * no_erp_phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_no_erp_phase_pre_pre = nan;
                end
    
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency)));
                if any(indices_for_pre_during_post)
                    trial_no_erp_phase_pre = angle(mean(exp(1i * no_erp_phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_no_erp_phase_pre = nan;
                end

                time_4_12_hz = (0:size(amplitude, 1) - 1) / sfreq + erp_pre_during_post(1);
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency - 1 / cluster2_peak_freq)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_pre_pre_pre = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_pre_pre_pre = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_pre_pre_pre = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_pre_pre_pre = nan;
                    trial_4_12_hz_pow_pre_pre_pre = nan;
                    trial_4_12_hz_phase_pre_pre_pre = nan;
                end
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency - 1 / cluster2_peak_freq / 2)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_pre_pre = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_pre_pre = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_pre_pre = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_pre_pre = nan;
                    trial_4_12_hz_pow_pre_pre = nan;
                    trial_4_12_hz_phase_pre_pre = nan;
                end
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_pre = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_pre = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_pre = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_pre = nan;
                    trial_4_12_hz_pow_pre = nan;
                    trial_4_12_hz_phase_pre = nan;
                end
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_amp_time)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_during = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_during = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_during = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_during = nan;
                    trial_4_12_hz_pow_during = nan;
                    trial_4_12_hz_phase_during = nan;
                end
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_amp_time + 1 / cluster2_peak_freq / 2)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_post = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_post = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_post = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_post = nan;
                    trial_4_12_hz_pow_post = nan;
                    trial_4_12_hz_phase_post = nan;
                end
                
                new_row = {i_trial, subj, task, contact_this_subj.ch_name_bci2000_format(i_ch), ...
                    saccade_table_all_session.session(i_trial), saccade_table_all_session.trial(i_trial) ...
                    saccade_table_all_session.eccentricity(i_trial), saccade_table_all_session.azimuth_before(i_trial), ...
                           saccade_table_all_session.elevation_before(i_trial),saccade_table_all_session.azimuth_after(i_trial), ...
                           saccade_table_all_session.elevation_after(i_trial), saccade_table_all_session.gaze_before_x_screen_coord(i_trial), ...
                           saccade_table_all_session.gaze_before_y_screen_coord(i_trial),saccade_table_all_session.gaze_after_x_screen_coord(i_trial), ...
                           saccade_table_all_session.gaze_after_x_screen_coord(i_trial),saccade_table_all_session.duration(i_trial),...
                           saccade_table_all_session.time_to_image_onset(i_trial), saccade_table_all_session.behavior_condition(i_trial),...
                           saccade_table_all_session.remember_condition(i_trial), corr(mean_erp_cluster2, current_trial_task, 'Type', 'Spearman'), ...
                           peak_to_trough, trial_trough_time,...
                           trial_latency, erp_amplitude, trial_amp_time, erp_polarity,...
                           trial_4_12_hz_pow_pre_pre_pre,trial_4_12_hz_pow_pre_pre, trial_4_12_hz_pow_pre, trial_4_12_hz_pow_during, trial_4_12_hz_pow_post,...
                           trial_4_12_hz_phase_pre_pre_pre, trial_4_12_hz_phase_pre_pre, trial_4_12_hz_phase_pre, trial_4_12_hz_phase_during, trial_4_12_hz_phase_post,...
                           trial_4_12_hz_norm_pow_pre_pre_pre, trial_4_12_hz_norm_pow_pre_pre, trial_4_12_hz_norm_pow_pre, trial_4_12_hz_norm_pow_during, trial_4_12_hz_norm_pow_post,...
                           trial_no_erp_phase_pre_pre_pre, trial_no_erp_phase_pre_pre, trial_no_erp_phase_pre};
                df_evoked_potential = [df_evoked_potential; new_row];
            else
                new_row = {i_trial, subj, task, contact_this_subj.ch_name_bci2000_format(i_ch), ...
                    saccade_table_all_session.session(i_trial), saccade_table_all_session.trial(i_trial) ...
                    saccade_table_all_session.eccentricity(i_trial), saccade_table_all_session.azimuth_before(i_trial), ...
                           saccade_table_all_session.elevation_before(i_trial),saccade_table_all_session.azimuth_after(i_trial), ...
                           saccade_table_all_session.elevation_after(i_trial), saccade_table_all_session.gaze_before_x_screen_coord(i_trial), ...
                           saccade_table_all_session.gaze_before_y_screen_coord(i_trial),saccade_table_all_session.gaze_after_x_screen_coord(i_trial), ...
                           saccade_table_all_session.gaze_after_x_screen_coord(i_trial),saccade_table_all_session.duration(i_trial),...
                           saccade_table_all_session.time_to_image_onset(i_trial), saccade_table_all_session.behavior_condition(i_trial),...
                           saccade_table_all_session.remember_condition(i_trial), corr(mean_erp_cluster2, current_trial_task, 'Type', 'Spearman'), ...
                           nan, nan, nan, nan, nan, nan, ...
                           nan, nan, nan, nan, nan,...
                           nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,...
                           nan, nan, nan};
                df_evoked_potential = [df_evoked_potential; new_row];
            end
        end
    end

    %%  add eog control
    if strcmp(subj, 'BJH025')
        erp_amp_window = 0.0025;
    end
    control_types = {'control', 'not_image_onset'};
    for i_control_type = 1: length(control_types)
        control_type = control_types{i_control_type};
        prep_signal_all_session = [];
        saccade_table_all_session = [];
        for session = 1:n_session
            filename = fullfile(data_dir, subj, ...
                [subj, '_', task, '_session', num2str(session), ...
                '_EOG_saccade_onset_', control_type, '.mat']);
            prep_signal = load(filename);
            if strcmp(control_type, 'control')
                prep_signal = prep_signal.EOG_signal_saccade_onset_control;
            else
                prep_signal = prep_signal.EOG_signal_saccade_onset_not_image_onset;
            end
            filename = fullfile(data_dir, subj,[subj '_' task '_session' num2str(session) ...
                    '_prep_saccade_event.mat']);
            prep_saccade_event = load(filename);
            prep_saccade_event = prep_saccade_event.data_struct_2_export;
            % get which eye is valid eye
            valid_side = prep_saccade_event.used_eye_for_saccade;    
    
            filename = fullfile(data_dir, subj,  ...
                [subj, '_', task, '_session', num2str(session), '_saccade_table.mat']);
            saccade_table = load(filename);
            saccade_table = saccade_table.saccade_table;
            % verify if signal match the saccade event
            if strcmp(control_type, 'control')
                assert(size(prep_signal, 2) == sum(contains(saccade_table.label, 'control')), ...
                'saccade event and prep signal length not matched');
                prep_signal_all_session = [prep_signal_all_session, prep_signal];
                saccade_table_all_session = [saccade_table_all_session; ...
                    saccade_table(contains(saccade_table.label, 'control'), :)];
            else
                assert(size(prep_signal, 2) == sum(contains(saccade_table.label, 'image_viewing_not_image_onset')), ...
                'saccade event and prep signal length not matched');
                prep_signal_all_session = [prep_signal_all_session, prep_signal];
                saccade_table_all_session = [saccade_table_all_session; ...
                    saccade_table(contains(saccade_table.label, 'image_viewing_not_image_onset'), :)];
            end
        end

        % merge erp and check if responsive
        time_vector = (0:size(prep_signal_all_session, 1) - 1) / sfreq - trial_begin_study;
        task_indices = find(time_vector >= erp_range(1) & time_vector <= erp_range(2));
        task_data = prep_signal_all_session(task_indices, :);

        baseline_data = [];
        for i_trial = 1:size(task_data, 2)
            rand_start = randi([round(baseline_time_precedence(1) * sfreq),...
                round(baseline_time_precedence(2) * sfreq)], 1);
            [~, rand_start_index] = min(abs(time_vector - rand_start / sfreq));
            baseline_data = [baseline_data, prep_signal_all_session(...
                rand_start_index:rand_start_index + length(task_indices) - 1, i_trial)];
        end
        baseline_indices = find(time_vector >= baseline_time_precedence(1) & time_vector <= erp_range(1));
        baseline_data = prep_signal_all_session(baseline_indices, :);
    
        baseline_mean = mean(baseline_data, 1);
        baseline_std = std(baseline_data, [], 1);

        if any(baseline_std == 0)
            zero_std_trial_indices = find(baseline_std==0);
            saccade_table_all_session(zero_std_trial_indices, :) = [];
            prep_signal_all_session(:, zero_std_trial_indices) = [];
            task_indices = find(time_vector >= erp_range(1) & time_vector <= erp_range(2));
            task_data = prep_signal_all_session(task_indices, :);
            baseline_indices = find(time_vector >= baseline_time_precedence(1) & time_vector <= erp_range(1));
            baseline_data = prep_signal_all_session(baseline_indices, :);
        
            baseline_mean = mean(baseline_data, 1);
            baseline_std = std(baseline_data, [], 1);
        end

        zscored_baseline = (baseline_data - baseline_mean) ./ (baseline_std);
        zscored_task = (task_data - baseline_mean) ./ (baseline_std);

        baseline_corr = corrcoef(zscored_baseline);
        task_corr = corrcoef(zscored_task);
        baseline_corr_values = baseline_corr(tril(true(size(baseline_corr)), -1));
        task_corr_values = task_corr(tril(true(size(task_corr)), -1));
        
        [p_value, ~, stats] = ranksum(baseline_corr_values, task_corr_values);
        combined_values = [baseline_corr_values; task_corr_values];
        labels = [zeros(length(baseline_corr_values), 1); ones(length(task_corr_values), 1)];
        
        [sorted_values, sorted_indices] = sort(combined_values, 'ascend');
        sorted_labels = labels(sorted_indices);
        
        total_positives = sum(sorted_labels);
        total_negatives = length(sorted_labels) - total_positives;
        
        tpr = zeros(length(sorted_labels), 1);
        fpr = zeros(length(sorted_labels), 1);
        
        for i = 1:length(sorted_labels)
            tpr(i) = sum(sorted_labels(i:end)) / total_positives; % Proportion of positive labels above the current threshold
            fpr(i) = sum(~sorted_labels(i:end)) / total_negatives; % Proportion of negative labels above the current threshold
        end
        
        auc = - trapz(fpr, tpr);
        norm_concatenated_prep_signal = (prep_signal_all_session - baseline_mean) ./ baseline_std;

   
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        
        time_vector = (0:length(mean_erp) - 1) / sfreq - trial_begin_study;

        %% baseline sub-classification control - for eog
        weak_trials_indices = []; % To store indices of weak trials
        original_baseline_corr = baseline_corr;
        original_order = 1:size(baseline_corr, 1); % Keep track of the original indices
        dist_matrix_trial = 1 - (baseline_corr);
        Z = linkage(squareform(dist_matrix_trial), 'average');
        order = optimalleaforder(Z, dist_matrix_trial);
        clusters = cluster(Z, 'maxclust', 2);
        clusters = clusters(order);
        while true
            cluster_sizes = histcounts(clusters, 2);
            if min(cluster_sizes) >= cluster_size_req * max(cluster_sizes)
                break;
            end
            
            mean_abs_corr = mean(abs(baseline_corr), 2);
            [~, weakest_trial_idx] = min(mean_abs_corr);
            
            weak_trials_indices = [weak_trials_indices; original_order(weakest_trial_idx)];
            
            baseline_corr(weakest_trial_idx, :) = [];
            baseline_corr(:, weakest_trial_idx) = [];
            
            original_order(weakest_trial_idx) = [];
            
            % re-cluster
            dist_matrix_trial = 1 - (baseline_corr);
            Z = linkage(squareform(dist_matrix_trial), 'average');
            order = optimalleaforder(Z, dist_matrix_trial);
            clusters = cluster(Z, 'maxclust', 2);
            clusters = clusters(order);
        end
        
        final_order = [original_order(order)'; weak_trials_indices]';

        if ~isempty(weak_trials_indices)
            reordered_baseline_corr = original_baseline_corr(final_order(1:end-length(weak_trials_indices)), final_order(1:end-length(weak_trials_indices)));
            weak_baseline_corr = original_baseline_corr(final_order(1:end-length(weak_trials_indices)), weak_trials_indices);
            bottom_right_matrix = original_baseline_corr(weak_trials_indices, weak_trials_indices);
            combined_matrix = [reordered_baseline_corr, weak_baseline_corr; weak_baseline_corr', bottom_right_matrix];
        else
            combined_matrix = original_baseline_corr(final_order, final_order);
        end
        
       

        cluster1_trials = final_order(clusters == 1);
        cluster2_trials = final_order(clusters == 2);

        within_cluster1_corr = original_baseline_corr(cluster1_trials, cluster1_trials);
        within_cluster2_corr = original_baseline_corr(cluster2_trials, cluster2_trials);
        mean_within_cluster1_corr_bl_control = mean(within_cluster1_corr(~eye(size(within_cluster1_corr))));
        mean_within_cluster2_corr_bl_control = mean(within_cluster2_corr(~eye(size(within_cluster2_corr))));
        between_cluster_corr = original_baseline_corr(cluster1_trials, cluster2_trials);
        mean_between_cluster_corr = mean(between_cluster_corr(:));
        cluster_corr_difference_bl_control = (mean_within_cluster1_corr_bl_control + mean_within_cluster2_corr_bl_control) / 2 - mean_between_cluster_corr;
        n_weak_corr_trials_bl_control = length(weak_trials_indices);

        cluster1_task_data = task_data(:, cluster1_trials);
        cluster1_baseline_mean = baseline_mean(:, cluster1_trials);
        cluster1_baseline_std = baseline_std(:, cluster1_trials);
        zscored_cluster1_baseline = (baseline_data(:, cluster1_trials) - cluster1_baseline_mean) ./ (cluster1_baseline_std);

        zscored_cluster1_task = (cluster1_task_data - cluster1_baseline_mean) ./ (cluster1_baseline_std);
        cluster1_corr = corrcoef(zscored_cluster1_task);
        cluster1_corr_values = cluster1_corr(tril(true(size(cluster1_corr)), -1));
        cluster1_baseline_corr = corrcoef(zscored_cluster1_baseline);
        cluster1_baseline_corr_values = cluster1_baseline_corr(tril(true(size(cluster1_baseline_corr)), -1));
        

        [cluster1_p_bl_control, ~, cluster1_stats] = ranksum(cluster1_corr_values, cluster1_baseline_corr_values);
        cluster1_combined_values = [cluster1_corr_values; cluster1_baseline_corr_values];
        cluster1_labels = [zeros(length(cluster1_corr_values), 1); ones(length(cluster1_baseline_corr_values), 1)];

        [sorted_values, sorted_indices] = sort(cluster1_combined_values, 'ascend');
        sorted_labels = cluster1_labels(sorted_indices);

        total_positives = sum(sorted_labels);
        total_negatives = length(sorted_labels) - total_positives;

        tpr = zeros(length(sorted_labels), 1);
        fpr = zeros(length(sorted_labels), 1);

        for i = 1:length(sorted_labels)
            tpr(i) = sum(sorted_labels(i:end)) / total_positives;
            fpr(i) = sum(~sorted_labels(i:end)) / total_negatives;
        end

        cluster1_auc_bl_control = -trapz(fpr, tpr);
        norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster1_trials) - ...
            cluster1_baseline_mean) ./ (cluster1_baseline_std);
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        
        time_vector = (0:length(mean_erp) - 1) / sfreq - trial_begin_study;
        
        cluster2_task_data = task_data(:, cluster2_trials);
        cluster2_baseline_mean = baseline_mean(:, cluster2_trials);
        cluster2_baseline_std = baseline_std(:, cluster2_trials);
        zscored_cluster2_baseline = (baseline_data(:, cluster2_trials) - cluster2_baseline_mean) ./ (cluster2_baseline_std);

        zscored_cluster2_task = (cluster2_task_data - cluster2_baseline_mean) ./ (cluster2_baseline_std);
        cluster2_corr = corrcoef(zscored_cluster2_task);
        cluster2_corr_values = cluster2_corr(tril(true(size(cluster2_corr)), -1));
        cluster2_baseline_corr = corrcoef(zscored_cluster2_baseline);
        cluster2_baseline_corr_values = cluster2_baseline_corr(tril(true(size(cluster2_baseline_corr)), -1));
        
        [cluster2_p_bl_control, ~, cluster2_stats] = ranksum(cluster2_corr_values, cluster2_baseline_corr_values);
        cluster2_combined_values = [cluster2_corr_values; cluster2_baseline_corr_values];
        cluster2_labels = [zeros(length(cluster2_corr_values), 1); ones(length(cluster2_baseline_corr_values), 1)];

        [sorted_values, sorted_indices] = sort(cluster2_combined_values, 'ascend');
        sorted_labels = cluster2_labels(sorted_indices);

        total_positives = sum(sorted_labels);
        total_negatives = length(sorted_labels) - total_positives;

        tpr = zeros(length(sorted_labels), 1);
        fpr = zeros(length(sorted_labels), 1);

        for i = 1:length(sorted_labels)
            tpr(i) = sum(sorted_labels(i:end)) / total_positives;
            fpr(i) = sum(~sorted_labels(i:end)) / total_negatives;
        end
        cluster2_auc_bl_control = -trapz(fpr, tpr);


        norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster2_trials) - ...
            cluster2_baseline_mean) ./ (cluster2_baseline_std);
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        

       
        baseline_corr = original_baseline_corr;
        %% baseline control done - for eog

        weak_trials_indices = []; % To store indices of weak trials
        original_task_corr = task_corr;
        original_order = 1:size(task_corr, 1); % Keep track of the original indices
        dist_matrix_trial = 1 - (task_corr);
        Z = linkage(squareform(dist_matrix_trial), 'average');
        order = optimalleaforder(Z, dist_matrix_trial);
        
        clusters = cluster(Z, 'maxclust', 2);
        clusters = clusters(order);
        
        while true
            cluster_sizes = histcounts(clusters, 2);
            if min(cluster_sizes) >= cluster_size_req * max(cluster_sizes)
                break;
            end
            
            mean_abs_corr = mean(abs(task_corr), 2);
            [~, weakest_trial_idx] = min(mean_abs_corr);

            weak_trials_indices = [weak_trials_indices; original_order(weakest_trial_idx)];
            
            task_corr(weakest_trial_idx, :) = [];
            task_corr(:, weakest_trial_idx) = [];
            
            original_order(weakest_trial_idx) = [];
            
            dist_matrix_trial = 1 - (task_corr);
            Z = linkage(squareform(dist_matrix_trial), 'average');
            order = optimalleaforder(Z, dist_matrix_trial);
            clusters = cluster(Z, 'maxclust', 2);
            clusters = clusters(order);
        end
        
        final_order = [original_order(order)'; weak_trials_indices]';
        
        if ~isempty(weak_trials_indices)
            reordered_task_corr = original_task_corr(final_order(1:end-length(weak_trials_indices)), final_order(1:end-length(weak_trials_indices)));
            weak_task_corr = original_task_corr(final_order(1:end-length(weak_trials_indices)), weak_trials_indices);
            bottom_right_matrix = original_task_corr(weak_trials_indices, weak_trials_indices);
            combined_matrix = [reordered_task_corr, weak_task_corr; weak_task_corr', bottom_right_matrix];
        else
            combined_matrix = original_task_corr(final_order, final_order);
        end
        
       
        clusters_2 = cluster(Z, 'maxclust', 2);
        clusters_2 = clusters_2(order);
     
        if ~isempty(weak_trials_indices)
            weak_task_data = task_data(:, weak_trials_indices);
            
            weak_baseline_mean = baseline_mean(:, weak_trials_indices);
            weak_baseline_std = baseline_std(:, weak_trials_indices);
            zscored_weak_baseline = (baseline_data(:, weak_trials_indices) - weak_baseline_mean) ./ (weak_baseline_std);

            zscored_weak_task = (weak_task_data - weak_baseline_mean) ./ (weak_baseline_std);
            weak_corr = corrcoef(zscored_weak_task);
            weak_corr_values = weak_corr(tril(true(size(weak_corr)), -1));
            weak_baseline_corr = corrcoef(zscored_weak_baseline);
            weak_baseline_corr_values = weak_baseline_corr(tril(true(size(weak_baseline_corr)), -1));
            
            [weak_p, ~, weak_stats] = ranksum(weak_baseline_corr_values, weak_corr_values);
            weak_combined_values = [weak_baseline_corr_values; weak_corr_values];
            weak_labels = [zeros(length(weak_baseline_corr_values), 1); ones(length(weak_corr_values), 1)];
            
            [sorted_values, sorted_indices] = sort(weak_combined_values, 'ascend');
            sorted_labels = weak_labels(sorted_indices);
            
            total_positives = sum(sorted_labels);
            total_negatives = length(sorted_labels) - total_positives;
            
            tpr = zeros(length(sorted_labels), 1);
            fpr = zeros(length(sorted_labels), 1);
            
            for i = 1:length(sorted_labels)
                tpr(i) = sum(sorted_labels(i:end)) / total_positives;
                fpr(i) = sum(~sorted_labels(i:end)) / total_negatives;
            end
            
            weak_auc = -trapz(fpr, tpr);

            horizontal_angle_weak = saccade_table_all_session(weak_trials_indices, :).horizontal_angle;
            verticle_angle_weak = saccade_table_all_session(weak_trials_indices, :).verticle_angle;
        else
            weak_auc = nan;
            horizontal_angle_weak = nan;
            verticle_angle_weak = nan;
            weak_p = nan;
        end

        cluster1_trials = final_order(clusters == 1);
        cluster2_trials = final_order(clusters == 2);

        within_cluster1_corr = original_task_corr(cluster1_trials, cluster1_trials);
        within_cluster2_corr = original_task_corr(cluster2_trials, cluster2_trials);
        mean_within_cluster1_corr = mean(within_cluster1_corr(~eye(size(within_cluster1_corr))));
        mean_within_cluster2_corr = mean(within_cluster2_corr(~eye(size(within_cluster2_corr))));
        between_cluster_corr = original_task_corr(cluster1_trials, cluster2_trials);
        mean_between_cluster_corr = mean(between_cluster_corr(:));
        cluster_corr_difference = (mean_within_cluster1_corr + mean_within_cluster2_corr) / 2 - mean_between_cluster_corr;

        cluster1_task_data = task_data(:, cluster1_trials);
        cluster1_baseline_mean = baseline_mean(:, cluster1_trials);
        cluster1_baseline_std = baseline_std(:, cluster1_trials);
        zscored_cluster1_baseline = (baseline_data(:, cluster1_trials) - cluster1_baseline_mean) ./ (cluster1_baseline_std);

        zscored_cluster1_task = (cluster1_task_data - cluster1_baseline_mean) ./ (cluster1_baseline_std);
        cluster1_corr = corrcoef(zscored_cluster1_task);
        cluster1_corr_values = cluster1_corr(tril(true(size(cluster1_corr)), -1));
        cluster1_baseline_corr = corrcoef(zscored_cluster1_baseline);
        cluster1_baseline_corr_values = cluster1_baseline_corr(tril(true(size(cluster1_baseline_corr)), -1));
        

        [cluster1_p, ~, cluster1_stats] = ranksum(cluster1_baseline_corr_values, cluster1_corr_values);
        cluster1_combined_values = [cluster1_baseline_corr_values; cluster1_corr_values];
        cluster1_labels = [zeros(length(cluster1_baseline_corr_values), 1); ones(length(cluster1_corr_values), 1)];

        [sorted_values, sorted_indices] = sort(cluster1_combined_values, 'ascend');
        sorted_labels = cluster1_labels(sorted_indices);

        total_positives = sum(sorted_labels);
        total_negatives = length(sorted_labels) - total_positives;

        tpr = zeros(length(sorted_labels), 1);
        fpr = zeros(length(sorted_labels), 1);

        for i = 1:length(sorted_labels)
            tpr(i) = sum(sorted_labels(i:end)) / total_positives;
            fpr(i) = sum(~sorted_labels(i:end)) / total_negatives;
        end

        cluster1_auc = -trapz(fpr, tpr);
        norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster1_trials) - ...
            cluster1_baseline_mean) ./ (cluster1_baseline_std);
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        
        time_vector = (0:length(mean_erp) - 1) / sfreq - trial_begin_study;

        

        horizontal_angle_cluster1 = saccade_table_all_session(cluster1_trials, :).azimuth_after - ...
            saccade_table_all_session(cluster1_trials, :).azimuth_before;
        horizontal_angle_cluster2 = saccade_table_all_session(cluster2_trials, :).azimuth_after - ...
            saccade_table_all_session(cluster2_trials, :).azimuth_before;
        vertical_angle_cluster1 = saccade_table_all_session(cluster1_trials, :).elevation_after - ...
            saccade_table_all_session(cluster1_trials, :).elevation_before;
        vertical_angle_cluster2 = saccade_table_all_session(cluster2_trials, :).elevation_after - ...
            saccade_table_all_session(cluster2_trials, :).elevation_before;
        
        angles1 = atan2d(vertical_angle_cluster1, horizontal_angle_cluster1);
        angles2 = atan2d(vertical_angle_cluster2, horizontal_angle_cluster2);
        
        

        cluster2_task_data = task_data(:, cluster2_trials);
        cluster2_baseline_mean = baseline_mean(:, cluster2_trials);
        cluster2_baseline_std = baseline_std(:, cluster2_trials);
        zscored_cluster2_baseline = (baseline_data(:, cluster2_trials) - cluster2_baseline_mean) ./ (cluster2_baseline_std);

        zscored_cluster2_task = (cluster2_task_data - cluster2_baseline_mean) ./ (cluster2_baseline_std);
        cluster2_corr = corrcoef(zscored_cluster2_task);
        cluster2_corr_values = cluster2_corr(tril(true(size(cluster2_corr)), -1));
        cluster2_baseline_corr = corrcoef(zscored_cluster2_baseline);
        cluster2_baseline_corr_values = cluster2_baseline_corr(tril(true(size(cluster2_baseline_corr)), -1));
        
        [cluster2_p, ~, cluster2_stats] = ranksum(cluster2_baseline_corr_values, cluster2_corr_values);
        cluster2_combined_values = [cluster2_baseline_corr_values; cluster2_corr_values];
        cluster2_labels = [zeros(length(cluster2_baseline_corr_values), 1); ones(length(cluster2_corr_values), 1)];

        [sorted_values, sorted_indices] = sort(cluster2_combined_values, 'ascend');
        sorted_labels = cluster2_labels(sorted_indices);

        total_positives = sum(sorted_labels);
        total_negatives = length(sorted_labels) - total_positives;

        tpr = zeros(length(sorted_labels), 1);
        fpr = zeros(length(sorted_labels), 1);

        for i = 1:length(sorted_labels)
            tpr(i) = sum(sorted_labels(i:end)) / total_positives;
            fpr(i) = sum(~sorted_labels(i:end)) / total_negatives;
        end
        cluster2_auc = -trapz(fpr, tpr);


        norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster2_trials) - ...
            cluster2_baseline_mean) ./ (cluster2_baseline_std);
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        

     

        
        gaze_x_change_cluster1 = saccade_table_all_session(cluster1_trials, :).gaze_after_x_screen_coord - ...
            saccade_table_all_session(cluster1_trials, :).gaze_before_x_screen_coord;
        gaze_x_change_cluster2 = saccade_table_all_session(cluster2_trials, :).gaze_after_x_screen_coord - ...
            saccade_table_all_session(cluster2_trials, :).gaze_before_x_screen_coord;
        gaze_y_change_cluster1 = saccade_table_all_session(cluster1_trials, :).gaze_after_y_screen_coord - ...
            saccade_table_all_session(cluster1_trials, :).gaze_before_y_screen_coord;
        gaze_y_change_cluster2 = saccade_table_all_session(cluster2_trials, :).gaze_after_y_screen_coord - ...
            saccade_table_all_session(cluster2_trials, :).gaze_before_y_screen_coord;
    
        
        gaze_x_before_cluster1 = saccade_table_all_session(cluster1_trials, :).gaze_before_x_screen_coord;
        gaze_x_before_cluster2 = saccade_table_all_session(cluster2_trials, :).gaze_before_x_screen_coord;
        gaze_y_before_cluster1 = saccade_table_all_session(cluster1_trials, :).gaze_before_y_screen_coord;
        gaze_y_before_cluster2 = saccade_table_all_session(cluster2_trials, :).gaze_before_y_screen_coord;
  
        
     
        time_vector = (0:size(prep_signal_all_session, 1) - 1) / sfreq - trial_begin_study;
        norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster1_trials) - ...
            cluster1_baseline_mean) ./ (cluster1_baseline_std);
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
       
      
        cluster1_data_oscillation = prep_signal_all_session(:, cluster1_trials);

        psd_all_trials = [];
        for trial = 1:length(cluster1_trials)
            [psd_trial, f] = pwelch(cluster1_data_oscillation(:, trial), hamming(window), noverlap, nfft, sfreq);
            psd_all_trials = [psd_all_trials, psd_trial];
        end
        
        mean_psd = mean(psd_all_trials, 2);
        std_psd = std(psd_all_trials, [], 2);

        
        fit_range_indices = f >= freq_range_fit(1) & f <= freq_range_fit(2);
        log_f = log10(f(fit_range_indices));
        log_psd = log10(mean_psd(fit_range_indices));
        p = polyfit(log_f, log_psd, 1);
        fit_line = 10.^(polyval(p, log10(f)));

        residuals = 10*log10(mean_psd(fit_range_indices)) - 10*log10(fit_line(fit_range_indices));
        outlier_indices = residuals > 0;
        if sum(outlier_indices) > length(outlier_indices) / 2
            [~, sorted_indices] = sort(residuals, 'descend');
            outlier_indices = false(size(residuals));
            outlier_indices(sorted_indices(1:floor(length(residuals) / 2))) = true;
        end
        fit_range_indices(fit_range_indices) = ~outlier_indices;

        log_f = log10(f(fit_range_indices));
        log_psd = log10(mean_psd(fit_range_indices));
        p = polyfit(log_f, log_psd, 1);
        fit_line = 10.^(polyval(p, log10(f)));
        
        
        total_power = [];
        num_trials = size(cluster1_data_oscillation, 2);
        for trial = 1:num_trials
            [S, F, T, P] = spectrogram(cluster1_data_oscillation(:, trial), hamming(window/16), ...
                noverlap / 16, nfft / 16, sfreq, 'yaxis');
            if isempty(total_power)
                total_power = zeros(size(P));
            end
            total_power = total_power + P;
        end
        mean_power = total_power / num_trials;
        baseline_indices = find(time_vector >= baseline_time_precedence(1) & time_vector <= erp_range(1));
        baseline_time_spec = time_vector(baseline_indices) + trial_begin_study;
        baseline_indices_spec = [];
        for i = 1:length(baseline_time_spec)
            [~, idx] = min(abs(T - baseline_time_spec(i)));
            if abs(T(idx) - baseline_time_spec(i)) <= 0.01
                baseline_indices_spec = [baseline_indices_spec, idx];
            end
        end
        baseline_indices_spec = unique(baseline_indices_spec);

        baseline_power = mean(mean_power(:, baseline_indices_spec), 2);
        baseline_corrected_power = 10*log10(mean_power ./ baseline_power);
        time_vector_spec = linspace(time_vector(1), time_vector(end), length(T));

        
        valid_idx = (freq_range_fit(2) > f) & (f > high_cutoff);
        f_filtered = f(valid_idx);
        mean_psd_filtered = 10*log10(mean_psd((valid_idx)));
        fit_line_filtered = 10*log10(fit_line(valid_idx));
        diff_psd = mean_psd_filtered - fit_line_filtered;
        [cluster1_oscillation_amp, peak_idx] = max(diff_psd);
        cluster1_peak_freq = f_filtered(peak_idx);
        low_cutoff_oscillation = max(0, cluster1_peak_freq - narrow_bandwidth);
        high_cutoff_oscillation = cluster1_peak_freq + narrow_bandwidth;

        task_indices_get_phase = find(time_vector >= erp_pre_during_post(1) & time_vector <= erp_pre_during_post(end));
        cluster1_data_oscillation = prep_signal_all_session(task_indices_get_phase, cluster1_trials);
        padding_samples = padding_duration * sfreq;

        amplitude = zeros(size(cluster1_data_oscillation));
        phase_cluster1 = zeros(size(cluster1_data_oscillation));
        [b, a] = butter(2, [low_cutoff_oscillation high_cutoff_oscillation] / (sfreq / 2)); % 2nd order Butterworth filter

        for trial = 1:size(cluster1_data_oscillation, 2)
            trial_data = cluster1_data_oscillation(:, trial);
            
            padded_data = [flipud(trial_data(1:padding_samples)); trial_data; flipud(trial_data(end-padding_samples+1:end))];
            
            filtered_data = filtfilt(b, a, padded_data);
            analytic_signal = hilbert(filtered_data);
            
            analytic_signal = analytic_signal(padding_samples+1:end-padding_samples);
            amplitude(:, trial) = abs(analytic_signal);
            phase_cluster1(:, trial) = angle(analytic_signal);
        end
        
        baseline_indices = find(((0:size(amplitude, 1) - 1) / sfreq + erp_pre_during_post(1) ) >= baseline_time_precedence(1) & ...
            ((0:size(amplitude, 1) - 1) / sfreq + erp_pre_during_post(1)) <= baseline_time_precedence(2));
        baseline_amp = amplitude(baseline_indices, :);
        baseline_amp_mean = mean(baseline_amp, 1);
        baseline_amp_std = std(baseline_amp, [], 1);
        norm_amplitude = (amplitude - ...
            baseline_amp_mean) ./ baseline_amp_std;
        mean_erp = mean(norm_amplitude, 2);
        std_error = std(norm_amplitude, [], 2) / sqrt(size(norm_amplitude, 2));
        cluster1_4_12_hz_pow = mean_erp;

       

        mean_phase = circ_mean(phase_cluster1, [], 2); % Circular mean
        std_error_phase = circ_std(phase_cluster1, [], [], 2) / sqrt(size(phase_cluster1, 2)); % Circular standard error

        cluster1_4_12_hz_phase = mean_phase;

        plv = abs(mean(exp(1i * phase_cluster1), 2));

        cluster1_4_12_hz_plv = plv;
        cluster1_4_12_hz_plv_all = mean(plv);

        
        time_vector = (0:size(prep_signal_all_session, 1) - 1) / sfreq - trial_begin_study;
        norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster2_trials) - ...
            cluster2_baseline_mean) ./ (cluster2_baseline_std);
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
       
       
        cluster2_data_oscillation = prep_signal_all_session(:, cluster2_trials);

        psd_all_trials = [];
        for trial = 1:length(cluster2_trials)
            [psd_trial, f] = pwelch(cluster2_data_oscillation(:, trial), hamming(window), noverlap, nfft, sfreq);
            psd_all_trials = [psd_all_trials, psd_trial];
        end
        
        mean_psd = mean(psd_all_trials, 2);
        std_psd = std(psd_all_trials, [], 2);

        
        fit_range_indices = f >= freq_range_fit(1) & f <= freq_range_fit(2);
        log_f = log10(f(fit_range_indices));
        log_psd = log10(mean_psd(fit_range_indices));
        p = polyfit(log_f, log_psd, 1);
        fit_line = 10.^(polyval(p, log10(f)));

        residuals = 10*log10(mean_psd(fit_range_indices)) - 10*log10(fit_line(fit_range_indices));
        outlier_indices = residuals > 0;
        if sum(outlier_indices) > length(outlier_indices) / 2
            [~, sorted_indices] = sort(residuals, 'descend');
            outlier_indices = false(size(residuals));
            outlier_indices(sorted_indices(1:floor(length(residuals) / 2))) = true;
        end
        fit_range_indices(fit_range_indices) = ~outlier_indices;

        log_f = log10(f(fit_range_indices));
        log_psd = log10(mean_psd(fit_range_indices));
        p = polyfit(log_f, log_psd, 1);
        fit_line = 10.^(polyval(p, log10(f)));
        
       
        total_power = [];
        num_trials = size(cluster2_data_oscillation, 2);
        for trial = 1:num_trials
            [S, F, T, P] = spectrogram(cluster2_data_oscillation(:, trial), hamming(window/16), ...
                noverlap / 16, nfft / 16, sfreq, 'yaxis');
            if isempty(total_power)
                total_power = zeros(size(P));
            end
            total_power = total_power + P;
        end
        mean_power = total_power / num_trials;
        baseline_indices = find(time_vector >= baseline_time_precedence(1) & time_vector <= erp_range(1));
        baseline_time_spec = time_vector(baseline_indices) + trial_begin_study;
        baseline_indices_spec = [];
        for i = 1:length(baseline_time_spec)
            [~, idx] = min(abs(T - baseline_time_spec(i)));
            if abs(T(idx) - baseline_time_spec(i)) <= 0.01
                baseline_indices_spec = [baseline_indices_spec, idx];
            end
        end
        baseline_indices_spec = unique(baseline_indices_spec);

        baseline_power = mean(mean_power(:, baseline_indices_spec), 2);
        baseline_corrected_power = 10*log10(mean_power ./ baseline_power);
        time_vector_spec = linspace(time_vector(1), time_vector(end), length(T));

       
        valid_idx = (freq_range_fit(2) > f) & (f > high_cutoff);
        f_filtered = f(valid_idx);
        mean_psd_filtered = 10*log10(mean_psd((valid_idx)));
        fit_line_filtered = 10*log10(fit_line(valid_idx));
        diff_psd = mean_psd_filtered - fit_line_filtered;
        [cluster2_oscillation_amp, peak_idx] = max(diff_psd);
        cluster2_peak_freq = f_filtered(peak_idx);
        low_cutoff_oscillation = max(0, cluster2_peak_freq - narrow_bandwidth);
        high_cutoff_oscillation = cluster2_peak_freq + narrow_bandwidth;

        task_indices_get_phase = find(time_vector >= erp_pre_during_post(1) & time_vector <= erp_pre_during_post(end));
        cluster2_data_oscillation = prep_signal_all_session(task_indices_get_phase, cluster2_trials);
        padding_samples = padding_duration * sfreq;

        amplitude = zeros(size(cluster2_data_oscillation));
        phase_cluster2 = zeros(size(cluster2_data_oscillation));
        [b, a] = butter(2, [low_cutoff_oscillation high_cutoff_oscillation] / (sfreq / 2)); % 2nd order Butterworth filter

        for trial = 1:size(cluster2_data_oscillation, 2)
            trial_data = cluster2_data_oscillation(:, trial);
            
            padded_data = [flipud(trial_data(1:padding_samples)); trial_data; flipud(trial_data(end-padding_samples+1:end))];
            
            filtered_data = filtfilt(b, a, padded_data);
            analytic_signal = hilbert(filtered_data);
            
            analytic_signal = analytic_signal(padding_samples+1:end-padding_samples);
            amplitude(:, trial) = abs(analytic_signal);
            phase_cluster2(:, trial) = angle(analytic_signal);
        end
        
        baseline_indices = find(((0:size(amplitude, 1) - 1) / sfreq + erp_pre_during_post(1) ) >= baseline_time_precedence(1) & ...
            ((0:size(amplitude, 1) - 1) / sfreq + erp_pre_during_post(1)) <= baseline_time_precedence(2));
        baseline_amp = amplitude(baseline_indices, :);
        baseline_amp_mean = mean(baseline_amp, 1);
        baseline_amp_std = std(baseline_amp, [], 1);
        norm_amplitude = (amplitude - ...
            baseline_amp_mean) ./ baseline_amp_std;
        mean_erp = mean(norm_amplitude, 2);
        std_error = std(norm_amplitude, [], 2) / sqrt(size(norm_amplitude, 2));
        cluster2_4_12_hz_pow = mean_erp;

        
        mean_phase = circ_mean(phase_cluster2, [], 2); % Circular mean
        std_error_phase = circ_std(phase_cluster2, [], [], 2) / sqrt(size(phase_cluster2, 2)); % Circular standard error

        % Compute the PLV
        cluster2_4_12_hz_phase = mean_phase;

        plv = abs(mean(exp(1i * phase_cluster2), 2));

        cluster2_4_12_hz_plv = plv;
        cluster2_4_12_hz_plv_all = mean(plv);

        all_phase = [phase_cluster1'; phase_cluster2']';
        all_4_12_hz_plv = abs(mean(exp(1i * all_phase), 2));

       
        %% characterize erp
        norm_concatenated_prep_signal = (prep_signal_all_session - baseline_mean) ./ (baseline_std);
        time_vector = (0:size(norm_concatenated_prep_signal, 1) - 1) / sfreq - trial_begin_study;
        task_indices = find(time_vector >= erp_characteristic_range(1) & time_vector <= erp_characteristic_range(2));
        task_data_for_characteristic_cluster1 = norm_concatenated_prep_signal(task_indices, cluster1_trials);
        task_data_for_characteristic_cluster2 = norm_concatenated_prep_signal(task_indices, cluster2_trials);
        % flag bad erp
        cluster1_erp_good = 1;
        cluster2_erp_good = 1;
        
        t = time_vector(find(time_vector >= erp_characteristic_range(1), 1):find(time_vector <= erp_characteristic_range(2), 1, 'last'));
        mean_erp_cluster1 = mean(task_data_for_characteristic_cluster1, 2);
        
        first_derivative = diff(mean(task_data_for_characteristic_cluster1, 2));
        first_derivative = smoothdata(first_derivative, 'movmean', smoothing_window);
       
        [change_points, ~] = findchangepts(first_derivative, 'MaxNumChanges', 2);
        if length(change_points) ~= 2
            first_derivative = smoothdata(first_derivative, 'movmean', smoothing_window_secondary);
            [change_points, ~] = findchangepts(first_derivative, 'MaxNumChanges', 2);
        end
        if length(change_points) == 2
            change_times = t(change_points);
            window_indices = find(t >= (change_times(2) - erp_amp_window) & t <= (change_times(2) + erp_amp_window));
            [erp_amplitude, max_index] = max(abs(mean_erp_cluster1(window_indices)));
            erp_polarity = 'Positive';
            if mean_erp_cluster1(window_indices(max_index)) < 0
                erp_polarity = 'Negative';
            end
            cluster1_polarity = erp_polarity;
            cluster1_latency = change_times(1);
            cluster1_amp = mean_erp_cluster1(window_indices(max_index));
            cluster1_amp_time = t(window_indices(max_index));
            if cluster1_amp_time > cluster1_latency
                trough_window_indices = find(t <= cluster1_amp_time & t >= (cluster1_latency - erp_amp_window / 2));
                if strcmp(cluster1_polarity, 'Positive')
                    [trough_amplitude, trough_index] = min(mean_erp_cluster1(trough_window_indices));
                else
                    [trough_amplitude, trough_index] = max(mean_erp_cluster1(trough_window_indices));
                end
                cluster1_peak_to_trough = cluster1_amp - trough_amplitude;
                cluster1_trough_time = t(trough_window_indices(trough_index));
            else
                cluster1_peak_to_trough = nan;
                cluster1_trough_time = nan;
            end

            

            % get no erp phase
            task_indices_get_phase = find(time_vector >= -trial_begin_study & time_vector <= (cluster1_latency));
            cluster1_data_oscillation = prep_signal_all_session(task_indices_get_phase, cluster1_trials);
            padding_samples = padding_duration * sfreq;

            no_erp_phase_cluster1 = zeros(size(cluster1_data_oscillation));
            low_cutoff_oscillation = max(0, cluster1_peak_freq - narrow_bandwidth);
            high_cutoff_oscillation = cluster1_peak_freq + narrow_bandwidth;

            [b, a] = butter(2, [low_cutoff_oscillation high_cutoff_oscillation] / (sfreq / 2)); 
            for trial = 1:size(cluster1_data_oscillation, 2)
                trial_data = cluster1_data_oscillation(:, trial);
                
                padded_data = [flipud(trial_data(1:padding_samples)); trial_data; flipud(trial_data(end-padding_samples+1:end))];
                
                filtered_data = filtfilt(b, a, padded_data);
                analytic_signal = hilbert(filtered_data);
                
                analytic_signal = analytic_signal(padding_samples+1:end-padding_samples);
                no_erp_phase_cluster1(:, trial) = angle(analytic_signal);
            end
            mean_phase = circ_mean(no_erp_phase_cluster1, [], 2); 
            std_error_phase = circ_std(no_erp_phase_cluster1, [], [], 2); 
            cluster1_no_erp_phase = mean_phase;
        
            plv = abs(mean(exp(1i * no_erp_phase_cluster1), 2));
        
           
            time_4_12_hz = (0:size(cluster1_no_erp_phase, 1) - 1) / sfreq - trial_begin_study;
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency - 1 / cluster1_peak_freq)));
            if any(indices_for_pre_during_post)
                cluster1_no_erp_phase_pre_pre_pre = mean(cluster1_no_erp_phase(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_no_erp_phase_pre_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency - 1 / cluster1_peak_freq / 2)));
            if any(indices_for_pre_during_post)
                cluster1_no_erp_phase_pre_pre = mean(cluster1_no_erp_phase(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_no_erp_phase_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency)));
            if any(indices_for_pre_during_post)
                cluster1_no_erp_phase_pre = mean(cluster1_no_erp_phase(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_no_erp_phase_pre = nan;
            end

            time_4_12_hz = (0:size(cluster1_4_12_hz_pow, 1) - 1) / sfreq + erp_pre_during_post(1);
           
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency - 1 / cluster1_peak_freq)));

            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_pre_pre_pre = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_pre_pre_pre = mean(cluster1_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_plv_pre_pre_pre = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_pre_pre_pre = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_pre_pre_pre = nan;
                cluster1_4_12_hz_phase_pre_pre_pre = nan;
                cluster1_4_12_hz_plv_pre_pre_pre = nan;
                all_4_12_hz_plv_pre_pre_pre = nan;
            end
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency - 1 / cluster1_peak_freq / 2)));

            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_pre_pre = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_pre_pre = mean(cluster1_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_plv_pre_pre = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_pre_pre = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_pre_pre = nan;
                cluster1_4_12_hz_phase_pre_pre = nan;
                cluster1_4_12_hz_plv_pre_pre = nan;
                all_4_12_hz_plv_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz -(cluster1_latency) ));

            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_pre = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_pre = mean(cluster1_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_plv_pre = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_pre = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_pre = nan;
                cluster1_4_12_hz_phase_pre = nan;
                cluster1_4_12_hz_plv_pre = nan;
                all_4_12_hz_plv_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_amp_time) ));

            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_during = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_during = mean(cluster1_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_plv_during = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_during = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_during = nan;
                cluster1_4_12_hz_phase_during = nan;
                cluster1_4_12_hz_plv_during = nan;
                all_4_12_hz_plv_during = nan;
            end
            
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_amp_time + 1 / cluster1_peak_freq / 2)));


            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_post = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_post = mean(cluster1_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_plv_post = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_post = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_post = nan;
                cluster1_4_12_hz_phase_post = nan;
                cluster1_4_12_hz_plv_post = nan;
                all_4_12_hz_plv_post = nan;
            end
            
        else
            cluster1_no_erp_phase_pre_pre_pre = nan;
            cluster1_no_erp_phase_pre_pre = nan;
            cluster1_no_erp_phase_pre = nan;
            cluster1_polarity = nan;
            cluster1_latency = nan;
            cluster1_amp = nan;
            cluster1_amp_time = nan;
            cluster1_erp_good = 0;
            cluster1_peak_to_trough = nan;
            cluster1_trough_time = nan;
            time_4_12_hz = (0:size(cluster1_4_12_hz_pow, 1) - 1) / sfreq + erp_pre_during_post(1);

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (erp_pre_during_post(3) - 1 / cluster1_peak_freq)));


            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_pre_pre_pre = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_pre_pre_pre = mean(cluster1_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_plv_pre_pre_pre = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_pre_pre_pre = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_pre_pre_pre = nan;
                cluster1_4_12_hz_phase_pre_pre_pre = nan;
                cluster1_4_12_hz_plv_pre_pre_pre = nan;
                all_4_12_hz_plv_pre_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (erp_pre_during_post(3) - 1 / cluster1_peak_freq / 2)));


            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_pre_pre = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_pre_pre = mean(cluster1_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_plv_pre_pre = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_pre_pre = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_pre_pre = nan;
                cluster1_4_12_hz_phase_pre_pre = nan;
                cluster1_4_12_hz_plv_pre_pre = nan;
                all_4_12_hz_plv_pre_pre = nan;
            end
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (erp_pre_during_post(3)) ));

            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_pre = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_pre = angle(mean(exp(1i * cluster1_4_12_hz_phase(indices_for_pre_during_post)), 'omitnan'));
                cluster1_4_12_hz_plv_pre = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_pre = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_pre = nan;
                cluster1_4_12_hz_phase_pre = nan;
                cluster1_4_12_hz_plv_pre = nan;
                all_4_12_hz_plv_pre = nan;
            end
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - ( erp_pre_during_post(4))));

            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_during = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_during = angle(mean(exp(1i * cluster1_4_12_hz_phase(indices_for_pre_during_post)), 'omitnan'));
                cluster1_4_12_hz_plv_during = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_during = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_during = nan;
                cluster1_4_12_hz_phase_during = nan;
                cluster1_4_12_hz_plv_during = nan;
                all_4_12_hz_plv_during = nan;
            end
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (erp_pre_during_post(3) + 1 / cluster1_peak_freq)));


            if any(indices_for_pre_during_post)
                cluster1_4_12_hz_pow_post = mean(cluster1_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster1_4_12_hz_phase_post = angle(mean(exp(1i * cluster1_4_12_hz_phase(indices_for_pre_during_post)), 'omitnan'));
                cluster1_4_12_hz_plv_post = mean(cluster1_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
                all_4_12_hz_plv_post = mean(all_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster1_4_12_hz_pow_post = nan;
                cluster1_4_12_hz_phase_post = nan;
                cluster1_4_12_hz_plv_post = nan;
                all_4_12_hz_plv_post = nan;
            end
        end


        mean_erp_cluster2 = mean(task_data_for_characteristic_cluster2, 2);

        first_derivative = diff(mean_erp_cluster2);
        first_derivative = smoothdata(first_derivative, 'movmean', smoothing_window);

        [change_points, ~] = findchangepts(first_derivative, 'MaxNumChanges', 2);
        if length(change_points) ~= 2
            first_derivative = smoothdata(first_derivative, 'movmean', smoothing_window_secondary);
            [change_points, ~] = findchangepts(first_derivative, 'MaxNumChanges', 2);
        end
        if length(change_points) == 2
            change_times = t(change_points);
            window_indices = find(t >= (change_times(2) - erp_amp_window) & t <= (change_times(2) + erp_amp_window));
            [erp_amplitude, max_index] = max(abs(mean_erp_cluster2(window_indices)));
            erp_polarity = 'Positive';
            if mean_erp_cluster2(window_indices(max_index)) < 0
                erp_polarity = 'Negative';
            end
            cluster2_polarity = erp_polarity;
            cluster2_latency = change_times(1);
            cluster2_amp = mean_erp_cluster2(window_indices(max_index));
            cluster2_amp_time = t(window_indices(max_index));
            if cluster2_amp_time > cluster2_latency
                trough_window_indices = find(t <= cluster2_amp_time & t >= (cluster2_latency - erp_amp_window / 2));
                if strcmp(cluster2_polarity, 'Positive')
                    [trough_amplitude, trough_index] = min(mean_erp_cluster2(trough_window_indices));
                else
                    [trough_amplitude, trough_index] = max(mean_erp_cluster2(trough_window_indices)); 
                end
                cluster2_peak_to_trough = cluster2_amp - trough_amplitude;
                cluster2_trough_time = t(trough_window_indices(trough_index));
            else
                cluster2_peak_to_trough = nan;
                cluster2_trough_time = nan;
            end


            task_indices_get_phase = find(time_vector >= -trial_begin_study & time_vector <= (cluster2_latency));
            cluster2_data_oscillation = prep_signal_all_session(task_indices_get_phase, cluster2_trials);
            padding_samples = padding_duration * sfreq;

            no_erp_phase_cluster2 = zeros(size(cluster2_data_oscillation));
            low_cutoff_oscillation = max(0, cluster2_peak_freq - narrow_bandwidth);
            high_cutoff_oscillation = cluster2_peak_freq + narrow_bandwidth;

            [b, a] = butter(2, [low_cutoff_oscillation high_cutoff_oscillation] / (sfreq / 2)); 
            for trial = 1:size(cluster2_data_oscillation, 2)
                trial_data = cluster2_data_oscillation(:, trial);
                
                padded_data = [flipud(trial_data(1:padding_samples)); trial_data; flipud(trial_data(end-padding_samples+1:end))];
                
                filtered_data = filtfilt(b, a, padded_data);
                analytic_signal = hilbert(filtered_data);
                
                analytic_signal = analytic_signal(padding_samples+1:end-padding_samples);
                no_erp_phase_cluster2(:, trial) = angle(analytic_signal);
            end
            mean_phase = circ_mean(no_erp_phase_cluster2, [], 2); 
            std_error_phase = circ_std(no_erp_phase_cluster2, [], [], 2); 
            cluster2_no_erp_phase = mean_phase;
        
            plv = abs(mean(exp(1i * no_erp_phase_cluster2), 2));
        
         
            time_4_12_hz = (0:size(cluster2_no_erp_phase, 1) - 1) / sfreq - trial_begin_study;
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency - 1 / cluster2_peak_freq)));
            if any(indices_for_pre_during_post)
                cluster2_no_erp_phase_pre_pre_pre = mean(cluster2_no_erp_phase(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_no_erp_phase_pre_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency - 1 / cluster2_peak_freq / 2)));
            if any(indices_for_pre_during_post)
                cluster2_no_erp_phase_pre_pre = mean(cluster2_no_erp_phase(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_no_erp_phase_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency)));
            if any(indices_for_pre_during_post)
                cluster2_no_erp_phase_pre = mean(cluster2_no_erp_phase(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_no_erp_phase_pre = nan;
            end

            time_4_12_hz = (0:size(cluster2_4_12_hz_pow, 1) - 1) / sfreq + erp_pre_during_post(1);
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency - 1 / cluster2_peak_freq)));

            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_pre_pre_pre = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_pre_pre_pre = mean(cluster2_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_plv_pre_pre_pre = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_pre_pre_pre = nan;
                cluster2_4_12_hz_phase_pre_pre_pre = nan;
                cluster2_4_12_hz_plv_pre_pre_pre = nan;
            end


            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency - 1 / cluster2_peak_freq / 2)));

         
            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_pre_pre = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_pre_pre = mean(cluster2_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_plv_pre_pre = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_pre_pre = nan;
                cluster2_4_12_hz_phase_pre_pre = nan;
                cluster2_4_12_hz_plv_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency)));
         
            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_pre = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_pre = angle(mean(exp(1i * cluster2_4_12_hz_phase(indices_for_pre_during_post)), 'omitnan'));
                cluster2_4_12_hz_plv_pre = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_pre = nan;
                cluster2_4_12_hz_phase_pre = nan;
                cluster2_4_12_hz_plv_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_amp_time) ));

            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_during = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_during = angle(mean(exp(1i * cluster2_4_12_hz_phase(indices_for_pre_during_post)), 'omitnan'));
                cluster2_4_12_hz_plv_during = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_during = nan;
                cluster2_4_12_hz_phase_during = nan;
                cluster2_4_12_hz_plv_during = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_amp_time + 1 / cluster2_peak_freq)));


          
            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_post = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_post = angle(mean(exp(1i * cluster2_4_12_hz_phase(indices_for_pre_during_post)), 'omitnan'));
                cluster2_4_12_hz_plv_post = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_post = nan;
                cluster2_4_12_hz_phase_post = nan;
                cluster2_4_12_hz_plv_post = nan;
            end
            
        else
            cluster2_no_erp_phase_pre_pre_pre = nan;
            cluster2_no_erp_phase_pre_pre = nan;
            cluster2_no_erp_phase_pre = nan;
            cluster2_polarity = nan;
            cluster2_latency = nan;
            cluster2_amp = nan;
            cluster2_amp_time = nan;
            cluster2_erp_good = 0;
            cluster2_peak_to_trough = nan;
            cluster2_trough_time = nan;

            time_4_12_hz = (0:size(cluster2_4_12_hz_pow, 1) - 1) / sfreq + erp_pre_during_post(1);
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (erp_pre_during_post(3) - 1 / cluster2_peak_freq)));

            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_pre_pre_pre = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_pre_pre_pre = mean(cluster2_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_plv_pre_pre_pre = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_pre_pre_pre = nan;
                cluster2_4_12_hz_phase_pre_pre_pre = nan;
                cluster2_4_12_hz_plv_pre_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (erp_pre_during_post(3) - 1 / cluster2_peak_freq / 2)));


            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_pre_pre = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_pre_pre = mean(cluster2_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_plv_pre_pre = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_pre_pre = nan;
                cluster2_4_12_hz_phase_pre_pre = nan;
                cluster2_4_12_hz_plv_pre_pre = nan;
            end
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (erp_pre_during_post(3) )));


            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_pre = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_pre = mean(cluster2_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_plv_pre = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_pre = nan;
                cluster2_4_12_hz_phase_pre = nan;
                cluster2_4_12_hz_plv_pre = nan;
            end
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (erp_pre_during_post(4)) ));
            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_during = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_during = mean(cluster2_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_plv_during = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_during = nan;
                cluster2_4_12_hz_phase_during = nan;
                cluster2_4_12_hz_plv_during = nan;
            end
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (erp_pre_during_post(4) + 1 / cluster2_peak_freq)));


            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_post = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_post = mean(cluster2_4_12_hz_phase(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_plv_post = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_post = nan;
                cluster2_4_12_hz_phase_post = nan;
                cluster2_4_12_hz_plv_post = nan;
            end
        end


        angles_rad = deg2rad(angles1);
        mean_sin = mean(sin(angles_rad), 'omitmissing');
        mean_cos = mean(cos(angles_rad), 'omitmissing');
        mean_angle_rad = atan2(mean_sin, mean_cos);
        mean_angle_deg = rad2deg(mean_angle_rad);
        mean_angle_deg1 = mod(mean_angle_deg, 360);
        if mean_angle_deg1 > 180
            mean_angle_deg1 = mean_angle_deg1 - 360;
        elseif mean_angle_deg1 <= -180
            mean_angle_deg1 = mean_angle_deg1 + 360;
        end

        R = sqrt(mean_cos^2 + mean_sin^2);
        circular_variance1 = 1 - R;

        angles_rad = deg2rad(angles2);
        mean_sin = mean(sin(angles_rad), 'omitmissing');
        mean_cos = mean(cos(angles_rad), 'omitmissing');
        mean_angle_rad = atan2(mean_sin, mean_cos);
        mean_angle_deg = rad2deg(mean_angle_rad);
        mean_angle_deg2 = mod(mean_angle_deg, 360);
        if mean_angle_deg2 > 180
            mean_angle_deg2 = mean_angle_deg2 - 360;
        elseif mean_angle_deg2 <= -180
            mean_angle_deg2 = mean_angle_deg2 + 360;
        end
        R = sqrt(mean_cos^2 + mean_sin^2);
        circular_variance2 = 1 - R;

        angles_rad = deg2rad([angles1; angles2]);
        mean_sin = mean(sin(angles_rad), 'omitmissing');
        mean_cos = mean(cos(angles_rad), 'omitmissing');


        R = sqrt(mean_cos^2 + mean_sin^2);
        circular_variance_all = 1 - R;

        new_row = {subj, task, ['EOG_' control_type], p_value < alpha_test, ...
            auc, p_value, cluster_corr_difference, mean_within_cluster1_corr, mean_within_cluster2_corr,...
            cluster_corr_difference_bl_control, mean_within_cluster1_corr_bl_control, mean_within_cluster2_corr_bl_control,...
            n_weak_corr_trials, n_weak_corr_trials_bl_control,...
            circular_variance_all, ...
            all_4_12_hz_plv_pre_pre_pre, all_4_12_hz_plv_pre_pre, all_4_12_hz_plv_pre, ...
            all_4_12_hz_plv_during, all_4_12_hz_plv_post,...
            cluster1_peak_freq, cluster1_oscillation_amp, cluster1_4_12_hz_pow_pre_pre_pre, cluster1_4_12_hz_pow_pre_pre, ...
            cluster1_4_12_hz_pow_pre, cluster1_4_12_hz_pow_during, cluster1_4_12_hz_pow_post,...
            cluster1_no_erp_phase_pre_pre_pre, cluster1_no_erp_phase_pre_pre, cluster1_no_erp_phase_pre,...
            cluster1_4_12_hz_phase_pre_pre_pre, cluster1_4_12_hz_phase_pre_pre, cluster1_4_12_hz_phase_pre,...
            cluster1_4_12_hz_phase_during, cluster1_4_12_hz_phase_post,...
            cluster1_4_12_hz_plv_all, cluster1_4_12_hz_plv_pre_pre_pre, cluster1_4_12_hz_plv_pre_pre, ...
            cluster1_4_12_hz_plv_pre, cluster1_4_12_hz_plv_during, cluster1_4_12_hz_plv_post,...
            cluster2_peak_freq, cluster2_oscillation_amp, cluster2_4_12_hz_pow_pre_pre_pre, cluster2_4_12_hz_pow_pre_pre,...
            cluster2_4_12_hz_pow_pre, cluster2_4_12_hz_pow_during, cluster2_4_12_hz_pow_post,...
            cluster2_no_erp_phase_pre_pre_pre, cluster2_no_erp_phase_pre_pre, cluster2_no_erp_phase_pre,...
            cluster2_4_12_hz_phase_pre_pre_pre, cluster2_4_12_hz_phase_pre_pre, cluster2_4_12_hz_phase_pre, cluster2_4_12_hz_phase_during, cluster2_4_12_hz_phase_post,...
            cluster2_4_12_hz_plv_all, cluster2_4_12_hz_plv_pre_pre_pre, cluster2_4_12_hz_plv_pre_pre, ...
            cluster2_4_12_hz_plv_pre, cluster2_4_12_hz_plv_during, cluster2_4_12_hz_plv_post,...
            cluster1_p < alpha_test, cluster1_auc, cluster1_p, cluster1_auc_bl_control, cluster1_p_bl_control,...
            cluster1_latency, cluster1_amp, cluster1_amp_time, cluster1_polarity, cluster1_peak_to_trough, cluster1_trough_time, ...
            cluster2_p < alpha_test, cluster2_auc, cluster2_p, cluster2_auc_bl_control, cluster2_p_bl_control,...
            cluster2_latency, cluster2_amp,...
           cluster2_amp_time, cluster2_polarity, cluster2_peak_to_trough, cluster2_trough_time,...
           mean(horizontal_angle_cluster1, 'omitmissing'), mean(vertical_angle_cluster1, 'omitmissing'), ...
           mean_angle_deg1, circular_variance1, ...
           mean(gaze_x_change_cluster1, 'omitmissing'),  mean(gaze_y_change_cluster1, 'omitmissing'), ...
           mean(horizontal_angle_cluster2, 'omitmissing'), mean(vertical_angle_cluster2, 'omitmissing'), ...
           mean_angle_deg2, circular_variance2, ...
           mean(gaze_x_change_cluster2, 'omitmissing'), mean(gaze_y_change_cluster2, 'omitmissing'), ...
           nan, nan, nan, ...
           (reference_subj_left_eye(1) + reference_subj_right_eye(1)) / 2, ...
           (reference_subj_left_eye(2) + reference_subj_right_eye(2)) / 2, ...
           (reference_subj_left_eye(3) + reference_subj_right_eye(3)) / 2, ...
           reference_subj_left_eye(1),...
           reference_subj_left_eye(2), reference_subj_left_eye(3), reference_subj_right_eye(1),...
           reference_subj_right_eye(2), reference_subj_right_eye(3),...
           valid_side, nan, 0, 0, 1, ...
           mean(horizontal_angle_weak, 'omitmissing'), ...
           mean(verticle_angle_weak, 'omitmissing')};
        responsive_erp_table = [responsive_erp_table; new_row];

%% get trial property        
        for i_trial = 1:size(prep_signal_all_session, 2)
            current_trial_task = norm_concatenated_prep_signal(task_indices, i_trial);

            if ismember(i_trial, cluster1_trials) && cluster1_erp_good
                first_derivative = diff(current_trial_task);
                first_derivative = smoothdata(first_derivative, 'movmean', smoothing_window);
                latency_window_indices = find(t >= (cluster1_latency - erp_amp_window) & ...
                    t <= (cluster1_latency + erp_amp_window));
                if latency_window_indices(end) > length(first_derivative)
                    shift_amount = latency_window_indices(end) - length(first_derivative);
                    latency_window_indices = latency_window_indices - shift_amount;
                    disp(['EOG trial ' num2str(i_trial) 'latency_window shifts'])
                end
                [~, min_index] = min(abs(t(latency_window_indices) - cluster1_latency));
                search_window_first_derivative = first_derivative(latency_window_indices);
                [change_point_index, ~] = findchangepts(search_window_first_derivative, 'MaxNumChanges', 1);
                if ~isempty(change_point_index)
                    trial_latency = t(latency_window_indices(change_point_index));
                else
                    trial_latency = cluster1_latency;
                end
                amp_window_indices = find(t >= (cluster1_amp_time - erp_amp_window) & ...
                    t <= (cluster1_amp_time + erp_amp_window));
                if amp_window_indices(end) > length(current_trial_task)
                    shift_amount = amp_window_indices(end) - length(current_trial_task);
                    amp_window_indices = amp_window_indices - shift_amount;
                    disp(['EOG trial ' num2str(i_trial) 'amp_window shifts'])
                end

                if strcmp(cluster1_polarity, 'Positive')
                    [erp_amplitude, max_index] = max(current_trial_task(amp_window_indices));
                    erp_polarity = 'Positive';
                    trial_amp_time = t(amp_window_indices(max_index));
                    erp_peak = current_trial_task(amp_window_indices(max_index));
                    if trial_amp_time > trial_latency
                        trough_window_indices = find(t <= trial_amp_time & t >= (trial_latency - erp_amp_window / 2));
                        [trough_amplitude, trough_index] = min(current_trial_task(trough_window_indices));
                        peak_to_trough = erp_peak - trough_amplitude;
                        trial_trough_time = t(trough_window_indices(trough_index));
                    else
                        peak_to_trough = nan;
                        trial_trough_time = nan;
                    end
                else
                    [erp_amplitude, max_index] = min(current_trial_task(amp_window_indices));
                    erp_polarity = 'Negative';
                    trial_amp_time = t(amp_window_indices(max_index));
                    erp_peak = current_trial_task(amp_window_indices(max_index));
                    if trial_amp_time > trial_latency
                        trough_window_indices = find(t <= trial_amp_time & t >= (trial_latency - erp_amp_window / 2));
                        [trough_amplitude, trough_index] = max(current_trial_task(trough_window_indices));
                        peak_to_trough = erp_peak - trough_amplitude;
                        trial_trough_time = t(trough_window_indices(trough_index));
                    else
                        peak_to_trough = nan;
                        trial_trough_time = nan;
                    end
                end
                time_vector = (0:size(prep_signal_all_session, 1) - 1) / sfreq - trial_begin_study;
                task_indices_get_phase = find(time_vector >= - trial_begin_study & time_vector <= erp_pre_during_post(end));
                trial_data_oscillation = prep_signal_all_session(task_indices_get_phase, i_trial);
                padding_samples = padding_duration * sfreq;
                low_cutoff_oscillation = max(0, cluster1_peak_freq - narrow_bandwidth);
                high_cutoff_oscillation = cluster1_peak_freq + narrow_bandwidth;

                [b, a] = butter(2, [low_cutoff_oscillation high_cutoff_oscillation] / (sfreq / 2));
                padded_data = [flipud(trial_data_oscillation(1:padding_samples)); trial_data_oscillation;...
                    flipud(trial_data_oscillation(end-padding_samples+1:end))];
            
                filtered_data = filtfilt(b, a, padded_data);
                analytic_signal = hilbert(filtered_data);
            
                analytic_signal = analytic_signal(padding_samples+1:end-padding_samples);
                amplitude = abs(analytic_signal);
                time_4_12_hz = (0:size(amplitude, 1) - 1) / sfreq - trial_begin_study;
                baseline_indices = find((time_4_12_hz >= baseline_time_precedence(1)) & ...
                    (time_4_12_hz <= baseline_time_precedence(2)));
                baseline_amp = amplitude(baseline_indices, :);
                baseline_amp_mean = mean(baseline_amp, 1);
                baseline_amp_std = std(baseline_amp, [], 1);
                norm_amplitude = (amplitude - ...
                    baseline_amp_mean) ./ baseline_amp_std;
                phase = angle(analytic_signal);

                task_indices_get_phase_no_erp = find(time_vector >= - trial_begin_study & time_vector <= (cluster1_latency));
                trial_data_oscillation = prep_signal_all_session(task_indices_get_phase_no_erp, i_trial);
                padding_samples = padding_duration / 2 * sfreq;

                padded_data = [flipud(trial_data_oscillation(1:padding_samples)); trial_data_oscillation; flipud(trial_data_oscillation(end-padding_samples+1:end))];
                
                filtered_data = filtfilt(b, a, padded_data);
                analytic_signal = hilbert(filtered_data);    
                analytic_signal = analytic_signal(padding_samples+1:end-padding_samples);
                no_erp_phase= angle(analytic_signal);

                time_4_12_hz = (0:size(no_erp_phase, 1) - 1) / sfreq - trial_begin_study;
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency - 1 / cluster1_peak_freq)));
                if any(indices_for_pre_during_post)
                    trial_no_erp_phase_pre_pre_pre = angle(mean(exp(1i * no_erp_phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_no_erp_phase_pre_pre_pre = nan;
                end
    
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency - 1 / cluster1_peak_freq / 2)));
                if any(indices_for_pre_during_post)
                    trial_no_erp_phase_pre_pre = angle(mean(exp(1i * no_erp_phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_no_erp_phase_pre_pre = nan;
                end
    
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency) ));
                if any(indices_for_pre_during_post)
                    trial_no_erp_phase_pre = angle(mean(exp(1i * no_erp_phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_no_erp_phase_pre = nan;
                end

                time_4_12_hz = (0:size(amplitude, 1) - 1) / sfreq + erp_pre_during_post(1);
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency - 1 / cluster1_peak_freq)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_pre_pre_pre = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_pre_pre_pre = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_pre_pre_pre = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_pre_pre_pre = nan;
                    trial_4_12_hz_pow_pre_pre_pre = nan;
                    trial_4_12_hz_phase_pre_pre_pre = nan;
                end
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency - 1 / cluster1_peak_freq / 2)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_pre_pre = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_pre_pre = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_pre_pre = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_pre_pre = nan;
                    trial_4_12_hz_pow_pre_pre = nan;
                    trial_4_12_hz_phase_pre_pre = nan;
                end
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency) ));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_pre = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_pre = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_pre = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_pre = nan;
                    trial_4_12_hz_pow_pre = nan;
                    trial_4_12_hz_phase_pre = nan;
                end
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_amp_time)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_during = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_during = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_during = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_during = nan;
                    trial_4_12_hz_pow_during = nan;
                    trial_4_12_hz_phase_during = nan;
                end
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_amp_time + 1 / cluster1_peak_freq / 2)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_post = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_post = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_post = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_post = nan;
                    trial_4_12_hz_pow_post = nan;
                    trial_4_12_hz_phase_post = nan;
                end



                new_row = {i_trial, subj, task, ['EOG_' control_type], ...
                    saccade_table_all_session.session(i_trial), saccade_table_all_session.trial(i_trial) ...
                    saccade_table_all_session.eccentricity(i_trial), saccade_table_all_session.azimuth_before(i_trial), ...
                           saccade_table_all_session.elevation_before(i_trial),saccade_table_all_session.azimuth_after(i_trial), ...
                           saccade_table_all_session.elevation_after(i_trial), saccade_table_all_session.gaze_before_x_screen_coord(i_trial), ...
                           saccade_table_all_session.gaze_before_y_screen_coord(i_trial),saccade_table_all_session.gaze_after_x_screen_coord(i_trial), ...
                           saccade_table_all_session.gaze_after_x_screen_coord(i_trial),saccade_table_all_session.duration(i_trial),...
                           saccade_table_all_session.time_to_image_onset(i_trial), saccade_table_all_session.behavior_condition(i_trial),...
                           saccade_table_all_session.remember_condition(i_trial), corr(mean_erp_cluster1, current_trial_task, 'Type', 'Spearman'), ...
                           peak_to_trough, trial_trough_time, ...
                           trial_latency, erp_amplitude, trial_amp_time, erp_polarity...
                           trial_4_12_hz_pow_pre_pre_pre, trial_4_12_hz_pow_pre_pre, trial_4_12_hz_pow_pre, trial_4_12_hz_pow_during, trial_4_12_hz_pow_post,...
                           trial_4_12_hz_phase_pre_pre_pre, trial_4_12_hz_phase_pre_pre, trial_4_12_hz_phase_pre, trial_4_12_hz_phase_during, trial_4_12_hz_phase_post,...
                           trial_4_12_hz_norm_pow_pre_pre_pre, trial_4_12_hz_norm_pow_pre_pre, trial_4_12_hz_norm_pow_pre, trial_4_12_hz_norm_pow_during,... 
                           trial_4_12_hz_norm_pow_post,...   
                           trial_no_erp_phase_pre_pre_pre, trial_no_erp_phase_pre_pre, trial_no_erp_phase_pre};

                df_evoked_potential = [df_evoked_potential; new_row];

            elseif ismember(i_trial, cluster2_trials) && cluster2_erp_good
                first_derivative = diff(current_trial_task);
                first_derivative = smoothdata(first_derivative, 'movmean', smoothing_window);

                latency_window_indices = find(t >= (cluster2_latency - erp_amp_window) & ...
                    t <= (cluster2_latency + erp_amp_window));
                if latency_window_indices(end) > length(first_derivative)
                    shift_amount = latency_window_indices(end) - length(first_derivative);
                    latency_window_indices = latency_window_indices - shift_amount;
                    disp(['EOG trial ' num2str(i_trial) 'latency_window shifts'])
                end
                [~, min_index] = min(abs(t(latency_window_indices) - cluster2_latency));
                search_window_first_derivative = first_derivative(latency_window_indices);
                [change_point_index, ~] = findchangepts(search_window_first_derivative, 'MaxNumChanges', 1);
                if ~isempty(change_point_index)
                    trial_latency = t(latency_window_indices(change_point_index));
                else
                    trial_latency = cluster2_latency;
                end
                amp_window_indices = find(t >= (cluster2_amp_time - erp_amp_window) & ...
                    t <= (cluster2_amp_time + erp_amp_window));
                if amp_window_indices(end) > length(current_trial_task)
                    shift_amount = amp_window_indices(end) - length(current_trial_task);
                    amp_window_indices = amp_window_indices - shift_amount;
                    disp(['EOG trial ' num2str(i_trial) 'amp_window shifts'])
                end
                if strcmp(cluster2_polarity, 'Positive')
                    [erp_amplitude, max_index] = max(current_trial_task(amp_window_indices));
                    erp_polarity = 'Positive';
                    trial_amp_time = t(amp_window_indices(max_index));
                    erp_peak = current_trial_task(amp_window_indices(max_index));
                    if trial_amp_time > trial_latency
                        trough_window_indices = find(t <= trial_amp_time & t >= (trial_latency - erp_amp_window / 2));
                        [trough_amplitude, trough_index] = min(current_trial_task(trough_window_indices));
                        peak_to_trough = erp_peak - trough_amplitude;
                        trial_trough_time = t(trough_window_indices(trough_index));
                    else
                        peak_to_trough = nan;
                        trial_trough_time = nan;
                    end
                else
                    [erp_amplitude, max_index] = min(current_trial_task(amp_window_indices));
                    erp_polarity = 'Negative';
                    trial_amp_time = t(amp_window_indices(max_index));
                    erp_peak = current_trial_task(amp_window_indices(max_index));
                    if trial_amp_time > trial_latency
                        trough_window_indices = find(t <= trial_amp_time & t >= (trial_latency - erp_amp_window / 2));
                        [trough_amplitude, trough_index] = max(current_trial_task(trough_window_indices));
                        peak_to_trough = erp_peak - trough_amplitude;
                        trial_trough_time = t(trough_window_indices(trough_index));
                    else
                        peak_to_trough = nan;
                        trial_trough_time = nan;
                    end
                end

                time_vector = (0:size(prep_signal_all_session, 1) - 1) / sfreq - trial_begin_study;
                task_indices_get_phase = find(time_vector >= - trial_begin_study & time_vector <= erp_pre_during_post(end));
                trial_data_oscillation = prep_signal_all_session(task_indices_get_phase, i_trial);
                padding_samples = padding_duration * sfreq;
                low_cutoff_oscillation = max(0, cluster2_peak_freq - narrow_bandwidth);
                high_cutoff_oscillation = cluster2_peak_freq + narrow_bandwidth;

                [b, a] = butter(2, [low_cutoff_oscillation high_cutoff_oscillation] / (sfreq / 2));
                padded_data = [flipud(trial_data_oscillation(1:padding_samples)); trial_data_oscillation;...
                    flipud(trial_data_oscillation(end-padding_samples+1:end))];
            
                filtered_data = filtfilt(b, a, padded_data);
                analytic_signal = hilbert(filtered_data);
            
                analytic_signal = analytic_signal(padding_samples+1:end-padding_samples);
                amplitude = abs(analytic_signal);
                time_4_12_hz = (0:size(amplitude, 1) - 1) / sfreq - trial_begin_study;
                baseline_indices = find((time_4_12_hz >= baseline_time_precedence(1)) & ...
                    (time_4_12_hz <= baseline_time_precedence(2)));
                baseline_amp = amplitude(baseline_indices, :);
                baseline_amp_mean = mean(baseline_amp, 1);
                baseline_amp_std = std(baseline_amp, [], 1);
                norm_amplitude = (amplitude - ...
                    baseline_amp_mean) ./ baseline_amp_std;
                phase = angle(analytic_signal);

                task_indices_get_phase_no_erp = find(time_vector >= - trial_begin_study & time_vector <= (cluster2_latency));
                trial_data_oscillation = prep_signal_all_session(task_indices_get_phase_no_erp, i_trial);
                padding_samples = padding_duration / 2 * sfreq;

                padded_data = [flipud(trial_data_oscillation(1:padding_samples)); trial_data_oscillation; flipud(trial_data_oscillation(end-padding_samples+1:end))];
                
                filtered_data = filtfilt(b, a, padded_data);
                analytic_signal = hilbert(filtered_data);    
                analytic_signal = analytic_signal(padding_samples+1:end-padding_samples);
                no_erp_phase= angle(analytic_signal);

                time_4_12_hz = (0:size(no_erp_phase, 1) - 1) / sfreq - trial_begin_study;
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency - 1 / cluster2_peak_freq)));
                if any(indices_for_pre_during_post)
                    trial_no_erp_phase_pre_pre_pre = angle(mean(exp(1i * no_erp_phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_no_erp_phase_pre_pre_pre = nan;
                end
    
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency - 1 / cluster2_peak_freq / 2)));
                if any(indices_for_pre_during_post)
                    trial_no_erp_phase_pre_pre = angle(mean(exp(1i * no_erp_phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_no_erp_phase_pre_pre = nan;
                end
    
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - cluster2_latency));
                if any(indices_for_pre_during_post)
                    trial_no_erp_phase_pre = angle(mean(exp(1i * no_erp_phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_no_erp_phase_pre = nan;
                end

                time_4_12_hz = (0:size(amplitude, 1) - 1) / sfreq + erp_pre_during_post(1);
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency - 1 / cluster2_peak_freq)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_pre_pre_pre = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_pre_pre_pre = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_pre_pre_pre = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_pre_pre_pre = nan;
                    trial_4_12_hz_pow_pre_pre_pre = nan;
                    trial_4_12_hz_phase_pre_pre_pre = nan;
                end
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency - 1 / cluster2_peak_freq / 2)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_pre_pre = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_pre_pre = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_pre_pre = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_pre_pre = nan;
                    trial_4_12_hz_pow_pre_pre = nan;
                    trial_4_12_hz_phase_pre_pre = nan;
                end
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_pre = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_pre = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_pre = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_pre = nan;
                    trial_4_12_hz_pow_pre = nan;
                    trial_4_12_hz_phase_pre = nan;
                end
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_amp_time)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_during = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_during = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_during = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_during = nan;
                    trial_4_12_hz_pow_during = nan;
                    trial_4_12_hz_phase_during = nan;
                end
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_amp_time + 1 / cluster2_peak_freq / 2)));

                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_post = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_post = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_post = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_post = nan;
                    trial_4_12_hz_pow_post = nan;
                    trial_4_12_hz_phase_post = nan;
                end
              
                new_row = {i_trial, subj, task, ['EOG_' control_type], ...
                    saccade_table_all_session.session(i_trial), saccade_table_all_session.trial(i_trial) ...
                    saccade_table_all_session.eccentricity(i_trial), saccade_table_all_session.azimuth_before(i_trial), ...
                           saccade_table_all_session.elevation_before(i_trial),saccade_table_all_session.azimuth_after(i_trial), ...
                           saccade_table_all_session.elevation_after(i_trial), saccade_table_all_session.gaze_before_x_screen_coord(i_trial), ...
                           saccade_table_all_session.gaze_before_y_screen_coord(i_trial),saccade_table_all_session.gaze_after_x_screen_coord(i_trial), ...
                           saccade_table_all_session.gaze_after_x_screen_coord(i_trial),saccade_table_all_session.duration(i_trial),...
                           saccade_table_all_session.time_to_image_onset(i_trial), saccade_table_all_session.behavior_condition(i_trial),...
                           saccade_table_all_session.remember_condition(i_trial), corr(mean_erp_cluster2, current_trial_task, 'Type', 'Spearman'), ...
                           peak_to_trough, trial_trough_time,...
                           trial_latency, erp_amplitude, trial_amp_time, erp_polarity,...
                           trial_4_12_hz_pow_pre_pre_pre, trial_4_12_hz_pow_pre_pre, trial_4_12_hz_pow_pre, trial_4_12_hz_pow_during, trial_4_12_hz_pow_post,...
                           trial_4_12_hz_phase_pre_pre_pre, trial_4_12_hz_phase_pre_pre, trial_4_12_hz_phase_pre, trial_4_12_hz_phase_during, trial_4_12_hz_phase_post,...
                           trial_4_12_hz_norm_pow_pre_pre_pre, trial_4_12_hz_norm_pow_pre_pre, trial_4_12_hz_norm_pow_pre, trial_4_12_hz_norm_pow_during, trial_4_12_hz_norm_pow_post,...
                           trial_no_erp_phase_pre_pre_pre, trial_no_erp_phase_pre_pre, trial_no_erp_phase_pre};


                df_evoked_potential = [df_evoked_potential; new_row];
            else
                new_row = {i_trial, subj, task, ['EOG_' control_type], ...
                    saccade_table_all_session.session(i_trial), saccade_table_all_session.trial(i_trial) ...
                    saccade_table_all_session.eccentricity(i_trial), saccade_table_all_session.azimuth_before(i_trial), ...
                           saccade_table_all_session.elevation_before(i_trial),saccade_table_all_session.azimuth_after(i_trial), ...
                           saccade_table_all_session.elevation_after(i_trial), saccade_table_all_session.gaze_before_x_screen_coord(i_trial), ...
                           saccade_table_all_session.gaze_before_y_screen_coord(i_trial),saccade_table_all_session.gaze_after_x_screen_coord(i_trial), ...
                           saccade_table_all_session.gaze_after_x_screen_coord(i_trial),saccade_table_all_session.duration(i_trial),...
                           saccade_table_all_session.time_to_image_onset(i_trial), saccade_table_all_session.behavior_condition(i_trial),...
                           saccade_table_all_session.remember_condition(i_trial), corr(mean_erp_cluster2, current_trial_task, 'Type', 'Spearman'), ...
                           nan, nan, nan, nan, nan, nan,...
                           nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,...
                           nan, nan, nan};
                df_evoked_potential = [df_evoked_potential; new_row];
            end

        end
    end
    % change it back
    if strcmp(subj, 'BJH025')
        erp_amp_window = 0.005;
    end

    filename = fullfile('result', ['control_responsive_erp_' subj '_' task '.csv']);
    writetable(responsive_erp_table, filename);
    
    filename = fullfile('result', ['control_df_evoked_potential_' subj '_' task '.csv']);
    writetable(df_evoked_potential, filename);
end
