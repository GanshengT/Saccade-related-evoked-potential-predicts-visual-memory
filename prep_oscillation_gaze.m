%% description
% This script will calculate the activation (ERP) for each contact

clear;
% addpath(['BLEAS/updated_mex_files']);
% addpath(['BLEAS/theta_saccade/Cluster Fix Package_GT_version']);
addpath(['CircStat2012a']);
data_dir = 'data';
sum_contact_table_filename = 'data/norm_contact_data_table.csv';
%% parameters
% subjs = {'BJH028', 'BJH029', 'BJH021', 'BJH021', 'BJH024', 'BJH024', 'BJH025', 'BJH026', 'BJH026', 'BJH027', 'BJH027', 'BJH027',...
%     'BJH030'};
% % we can use the mean to differentiate left and right electrodes
% tasks = {'BLAES_study', 'BLAES_study', 'BLAES_study', 'BLAES_study_2', 'BLAES_study', 'BLAES_study_2', 'BLAES_study_twosource', ...
%     'BLAES_study', 'BLAES_study_2', ...
%     'BLAES_study', 'BLAES_study_2', 'BLAES_study_3'...
%     'BLAES_study'};
% n_session_per_task = [2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2];

% running for one subject %
subjs = {'BJH025',};
% we can use the mean to differentiate left and right electrodes
tasks = { 'BLAES_study_twosource'};
n_session_per_task = [4];
%%%%%%
eye_pos_in_mri_space = containers.Map();
eye_pos_in_mri_space('BJH025_left') = [-29.84, 49.21, -9.67];
eye_pos_in_mri_space('BJH025_right') = [35.51, 45.80, -9.67];

reference_subj_left_eye = [-29.84, 49.21, -9.67];
reference_subj_right_eye = [35.51, 45.80, -9.67];

bad_ch_map = containers.Map();

BJH025_tasks = struct();
BJH025_tasks.BLAES_study_twosource = {'D8','GR4','GR5','J7'};
bad_ch_map('BJH025') = BJH025_tasks;

cluster_size_req = 0.25; % no less than 0.25 size of the other cluster
laplacian_thres = 6; % 5mm bound, we are changing laplacian threshold to 6 for 28 and 29
goodness_of_fit_thres = 0.8;
monitor_dimension = [611.2, 362.6]; % in mm
smoothing_window = 10; % 10* 1/2000s = 5ms;
epsilon = 1e-8; % handle all 0 vector, make corr robust

trial_begin_study = 3; % 3s
sfreq = 2000; % 2000Hz
baseline_time_precedence = [-0.2, -0.09];
erp_range = [-0.01 0.04]; % 200ms
erp_pre_during_post = [-0.3 -0.04 -0.01 0.03 0.8 0.3]; % 600ms
erp_characteristic_range = [-0.05, 0.1]; % for getting change pts
erp_amp_window = 0.01; % in s, we look at max around the second change pt
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
% prepare a df, with these col names:
% subj, task, session, trial, ch_name_bci2000, blc_theta, blc_gamma
% psd verification
window = 2048; % half a second
noverlap = 1024; 
nfft = 2 * window;
freq_range = [1 50];
% i/f fitting
freq_range_fit = [1 60]; 
narrow_bandwidth = 3; % half band width
smoothing_window_secondary = 20; % 20* 1/2000s = 10ms;
plot_individual_trial = 0;

%% gaze properties
df_gaze_properties = table('Size', [0, 21], ...
    'VariableTypes', {'string', 'string', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'string', 'double', 'double', ...
    'string', 'string'}, ...
    'VariableNames', {'subj', 'task', 'session', 'trial', ...
    'id', 'gaze_before_x_screen_coord', 'gaze_before_y_screen_coord', ...
    'gaze_after_x_screen_coord', 'gaze_after_y_screen_coord',...
    'azimuth_before', 'elevation_before', 'azimuth_after', 'elevation_after',...
    'eccentricity', 'duration', 'time_to_image_onset', 'event', 'amgd_stim_flag', 'pre_amgd_stim_flag', ...
    'behavior_condition', 'remember_condition'});


for i_subj = 1:length(subjs)
    subj = subjs{i_subj};
    task = tasks{i_subj};
    n_session = n_session_per_task(i_subj);
    % get ch_label
    contact_this_subj = sum_contact_table(strcmp(sum_contact_table.SubjectID, subj), :);
    report_path = fullfile(data_dir, subj, 'report');
    id_gaze_event = 1;
    for session = 1:n_session
        % load event table and saccade table
        filename = fullfile(data_dir, subj, ...
                [subj, '_', task, '_session', num2str(session), '_saccade_table.mat']);
        saccade_table = load(filename);
        saccade_table = saccade_table.saccade_table;
        filename = fullfile(data_dir, subj, [subj '_' task '_session' num2str(session) ...
                    '_prep_saccade_event.mat']);
        prep_saccade_event = load(filename);
        prep_saccade_event = prep_saccade_event.data_struct_2_export;
        % get fixation duration iterate over fixationstat
        for i_trial = 1:length(prep_saccade_event.trial_event_struct)
            if prep_saccade_event.trial_event_struct(i_trial).stimulation_condition == '1'
                amgd_stim_flag = 1;
            else
                amgd_stim_flag = 0;
            end
            if (i_trial > 1) && ...
                    (prep_saccade_event.trial_event_struct(i_trial - 1).stimulation_condition == '1')
                pre_amgd_stim_flag = 1;
            else
                pre_amgd_stim_flag = 0;
            end
            for i_eye = 1:size(prep_saccade_event.fixationstats{i_trial}.fixationtimes, 2)
                if all(prep_saccade_event.fixationstats{i_trial}.valid_fixation_label(:, i_eye))
                    new_row_gaze = {subj, task, session, i_trial, ...
                        id_gaze_event, prep_saccade_event.fixationstats{i_trial}.fixations(1, i_eye),...
                        prep_saccade_event.fixationstats{i_trial}.fixations(2, i_eye),...
                        prep_saccade_event.fixationstats{i_trial}.fixations(1, i_eye),...
                        prep_saccade_event.fixationstats{i_trial}.fixations(2, i_eye),...
                        nan, nan, nan, nan, nan...
                        (prep_saccade_event.fixationstats{i_trial}.fixationtimes(2, i_eye) - ...
                        prep_saccade_event.fixationstats{i_trial}.fixationtimes(1, i_eye)) / sfreq,...
                        prep_saccade_event.fixationstats{i_trial}.fixationtimes(1, i_eye) / prep_saccade_event.sfreq + ...
                            prep_saccade_event.trial_event_begin,...
                        'fixation', amgd_stim_flag, pre_amgd_stim_flag,...
                        prep_saccade_event.trial_event_struct(i_trial).behavior_condition,...
                        prep_saccade_event.trial_event_struct(i_trial).remember_condition};
                    id_gaze_event = id_gaze_event + 1;
                    df_gaze_properties = [df_gaze_properties; new_row_gaze];
                end

            end
            % get sub-saccade table
            subset_index = (strcmp(saccade_table.subj, subj)) & ...
                (saccade_table.session == session) & strcmp(saccade_table.task, task) & ...
                (saccade_table.trial == i_trial);
            sub_saccade_table = saccade_table(subset_index, :);
            for i_eye = 1:height(sub_saccade_table)
                new_row_gaze = {subj, task, session, i_trial, ...
                    id_gaze_event, sub_saccade_table.gaze_before_x_screen_coord(i_eye),...
                    sub_saccade_table.gaze_before_y_screen_coord(i_eye),...
                    sub_saccade_table.gaze_after_x_screen_coord(i_eye),...
                    sub_saccade_table.gaze_after_y_screen_coord(i_eye),...
                    sub_saccade_table.azimuth_before(i_eye), sub_saccade_table.elevation_before(i_eye),...
                    sub_saccade_table.azimuth_after(i_eye), sub_saccade_table.elevation_after(i_eye), ...
                    sub_saccade_table.eccentricity(i_eye), sub_saccade_table.duration(i_eye),...
                    sub_saccade_table.time_to_image_onset(i_eye),'saccade', amgd_stim_flag, pre_amgd_stim_flag,...
                    prep_saccade_event.trial_event_struct(i_trial).behavior_condition,...
                    prep_saccade_event.trial_event_struct(i_trial).remember_condition};   
                id_gaze_event = id_gaze_event + 1;
                df_gaze_properties = [df_gaze_properties; new_row_gaze];
            end
        end
    end
% save df as csv and mat
    filename = fullfile('result', [subj '_' task '_gaze_properties.csv']);
    writetable(df_gaze_properties, filename);
end



%% identify responsive channel
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

 
% ref_dir calculation : get neighbors (small laplacian), fit a line, if a
% line can not be fitted (use contact in other shanks, flag shank_ref as 0)

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

    % for i_ch = 1:height(contact_this_subj)
    % uncomment for running for all channel, the provided script is being
    % tested in a subset of channel for efficiency
    for i_ch = 1:6
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

        % fit a line for neighbors including this ch
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
        % shank_ref coding; 0 no rereference, 2: small laplacian, 1,
        % bipolar; 3: not linear small laplacian
        clear shank_direction
        if size(all_points, 1) < 2 % no neighbor, seldom happen
            shank_ref = 0;
            re_ref_vector_x_norm = nan;
            re_ref_vector_y_norm = nan;
            re_ref_vector_z_norm = nan;
        else
            if size(all_points, 1) > 2 % meaning more than 2 neighbors, small laplacian
                shank_ref = 2;
            else
                shank_ref = 1;
            end
            % Perform linear regression
            [coeff, score, ~, ~, explained] = pca(all_points);
            
            % Check the goodness of fit
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
            filename = fullfile(data_dir, subj, ...
                [subj, '_', task, '_session', num2str(session), ...
                '_channel', num2str(i_ch), '_saccade_onset_not_image_onset_prep_signal.mat']);
            prep_signal = load(filename);
            prep_signal = prep_signal.prep_signal_saccade_onset_not_image_onset;
            filename = fullfile(data_dir, subj, [subj '_' task '_session' num2str(session) ...
                    '_prep_saccade_event.mat']);
            prep_saccade_event = load(filename);
            prep_saccade_event = prep_saccade_event.data_struct_2_export;
            % get which eye is valid eye
            valid_side = prep_saccade_event.used_eye_for_saccade;
            % get eye position in MRI coordinates
            key_for_eye = [subj '_left'];
            left_eye_position_in_mri_coord = eye_pos_in_mri_space(key_for_eye);
            % % get distance vector between eye and contact in mm
            contact_dist_to_left_eye = norm(current_point - left_eye_position_in_mri_coord);
            key_for_eye = [subj '_right'];
            right_eye_position_in_mri_coord = eye_pos_in_mri_space(key_for_eye);
            % % get distance vector between eye and contact in mm
            contact_dist_to_right_eye = norm(current_point - right_eye_position_in_mri_coord);

            filename = fullfile(data_dir, subj, ...
                [subj, '_', task, '_session', num2str(session), '_saccade_table.mat']);
            saccade_table = load(filename);
            saccade_table = saccade_table.saccade_table;
            % verify if signal match the saccade event
            assert(size(prep_signal, 2) == sum(contains(saccade_table.label, 'image_viewing_not_image_onset')), ...
            'saccade event and prep signal length not matched');
            prep_signal_all_session = [prep_signal_all_session, prep_signal];
            saccade_table_all_session = [saccade_table_all_session; ...
                saccade_table(contains(saccade_table.label, 'image_viewing_not_image_onset'), :)];
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
        if any(baseline_std==0)
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

        zscored_baseline = (baseline_data - baseline_mean) ./ baseline_std;
        zscored_task = (task_data - baseline_mean) ./ baseline_std;

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
        
        %  TPR and FPR as the number of labels larger than the current value
        for i = 1:length(sorted_labels)
            tpr(i) = sum(sorted_labels(i:end)) / total_positives; % Proportion of positive labels above the current threshold
            fpr(i) = sum(~sorted_labels(i:end)) / total_negatives; % Proportion of negative labels above the current threshold
        end
        
        auc = - trapz(fpr, tpr);
        norm_concatenated_prep_signal = (prep_signal_all_session - baseline_mean) ./ baseline_std;
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        
        time_vector = (0:length(mean_erp) - 1) / sfreq - trial_begin_study;

        % verification plot

        % fig = figure('Visible', 'off');
        % set(fig, 'position', [100, 100, 1200, 800])
        % subplot(3, 2, [1, 2]);
        % hold on;
        % 
        % plot(time_vector, mean_erp, 'b', 'LineWidth', 2);
        % 
        % fill([time_vector, fliplr(time_vector)], ...
        %  [mean_erp' + std_error', fliplr(mean_erp' - std_error')], ...
        %  'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % xline(0, 'r--', 'LineWidth', 1.5);
        % 
        % xlabel('Time (s)');
        % ylabel('Amplitude (\muV)');
        % title('Event-Related Potential (ERP) - Full Signal');
        % legend({'ERP', 'Standard Error', 'saccade onset'});
        % 
        % set(gca, 'FontSize', 12);
        % grid on;
        % hold off;
        % 
        % subplot(3, 2, [3, 4]);
        % hold on;
        % 
        % range_start = find(time_vector >= -0.2, 1);
        % range_end = find(time_vector <= 0.2, 1, 'last');
        % 
        % plot(time_vector(range_start:range_end), mean_erp(range_start:range_end), 'b', 'LineWidth', 2);
        % 
        % fill([time_vector(range_start:range_end), fliplr(time_vector(range_start:range_end))], ...
        %  [mean_erp(range_start:range_end)' + std_error(range_start:range_end)', ...
        %  fliplr((mean_erp(range_start:range_end))' - (std_error(range_start:range_end))')], ...
        %  'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % xline(0, 'r--', 'LineWidth', 1.5);
        % 
        % xlabel('Time (s)');
        % ylabel('Amplitude (\muV)');
        % title('Event-Related Potential (ERP) - Around Saccade Onset');
        % legend({'ERP', 'Standard Error', 'saccade onset'});
        % 
        % 
        % set(gca, 'FontSize', 12);
        % grid on;
        % hold off;
        % 
        % subplot(3, 2, 5);
        % hold on;
        % num_bins = 20; % Specify the number of bins
        % histogram(baseline_corr_values, num_bins, 'Normalization', 'probability', 'FaceAlpha', 0.5, 'DisplayName', 'Baseline');
        % histogram(task_corr_values, num_bins, 'Normalization', 'probability', 'FaceAlpha', 0.5, 'DisplayName', 'Task');
        % xlabel('Correlation Values');
        % ylabel('Probability');
        % title('Distribution of Baseline and Task Correlation Values');
        % legend('show');
        % hold off;
        % subplot(3, 2, 6);
        % hold on;
        % plot(fpr, tpr, 'b', 'LineWidth', 2);
        % plot([0 1], [0 1], 'r--'); % Diagonal line for reference
        % xlabel('False Positive Rate');
        % ylabel('True Positive Rate');
        % title(sprintf('ROC Curve (AUC = %.2f, p = %.3f)', auc, p_value));
        % hold off;

        %% baseline sub-classification control
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
            
            % identify trial with weakest connection
            mean_abs_corr = mean(abs(baseline_corr), 2);
            [~, weakest_trial_idx] = min(mean_abs_corr);
            
            %  the index of the weakest trial
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
        
        % verification plot
        % fig = figure('Visible','off');
        % set(fig, 'position', [100, 100, 900, 900])
        % 
        % subplot(2, 2, 1);
        % imagesc(original_baseline_corr);
        % colormap('coolwarm');
        % clim([-1, 1]);
        % axis square;
        % title('Original Baseline Correlation Matrix');
        % xlabel('Trials');
        % ylabel('Trials');
        % colorbar;
        % 
        % subplot(2, 2, 2);
        % dendrogram(Z, 0);
        % title('Dendrogram');
        % xlabel('Trials');
        % ylabel('Distance');
        % 
        % subplot(2, 2, 3);
        % imagesc(combined_matrix);
        % colormap('coolwarm');
        % clim([-1, 1]);
        % axis square;
        % title('Reordered Baseline Correlation Matrix with Clusters');
        % xlabel('Trials');
        % ylabel('Trials');
        % 
        % hold on;
        % clusters_2 = cluster(Z, 'maxclust', 2);
        % clusters_2 = clusters_2(order);
        % for i = 1:max(clusters_2)
        %     idx = find(clusters_2 == i);
        %     rectangle('Position', [idx(1)-0.5, idx(1)-0.5, length(idx), length(idx)], ...
        %               'EdgeColor', 'w', 'LineWidth', 2, 'LineStyle', '--');
        %     text(mean(idx), mean(idx), num2str(i), 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        % end
        % 
        % if ~isempty(weak_trials_indices)
        %     weak_idx_start = size(reordered_baseline_corr, 1) + 1;
        %     weak_idx_end = size(combined_matrix, 1);
        %     rectangle('Position', [weak_idx_start - 0.5, weak_idx_start - 0.5, weak_idx_end - weak_idx_start + 1, weak_idx_end - weak_idx_start + 1], ...
        %               'EdgeColor', 'k', 'LineWidth', 2, 'LineStyle', '--');
        %     text((weak_idx_start + weak_idx_end) / 2, (weak_idx_start + weak_idx_end) / 2, ...
        %          'Weak', 'Color', 'k', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        % end
        % 
        % subplot(2, 2, 4);
        % axis off; % Turn off the axis
        % c = colorbar('southoutside');
        % c.Label.String = 'Correlation Coefficient';
        % c.Label.FontSize = 12;
        % colormap('coolwarm');
        % clim([-1, 1]);

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
        zscored_cluster1_baseline = (baseline_data(:, cluster1_trials) - cluster1_baseline_mean) ./ cluster1_baseline_std;

        zscored_cluster1_task = (cluster1_task_data - cluster1_baseline_mean) ./ cluster1_baseline_std;
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
            cluster1_baseline_mean) ./ cluster1_baseline_std;
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        
        time_vector = (0:length(mean_erp) - 1) / sfreq - trial_begin_study;
        % verification
        % fig = figure('visible', 'off');
        % set(fig, 'position', [100, 100, 1500, 700])
        % range_start = find(time_vector >= time_range_for_plotting_erp(1), 1);
        % range_end = find(time_vector <= time_range_for_plotting_erp(2), 1, 'last');
        % subplot(2, 5, [1 2]);
        % plot(time_vector(range_start:range_end), mean_erp(range_start:range_end), 'b', 'LineWidth', 2);
        % hold on;
        % fill([time_vector(range_start:range_end), fliplr(time_vector(range_start:range_end))], ...
        %  [mean_erp(range_start:range_end)' + std_error(range_start:range_end)', ...
        %  fliplr((mean_erp(range_start:range_end))' - (std_error(range_start:range_end))')], ...
        %  'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % xline(0, 'r--', 'LineWidth', 1.5);
        % 
        % xlabel('Time (s)');
        % ylabel('Amplitude (\muV)');
        % title(['Cluster 1 ERP n_trial: ' num2str(size(norm_concatenated_prep_signal, 2))], 'Interpreter', 'none');
        % legend({'ERP', 'Standard Error', 'saccade onset'});
        % grid on;
        % hold off;
        % 
        % subplot(2, 5, 3)
        % hold on;
        % num_bins = 20; % Specify the number of bins
        % histogram(cluster1_baseline_corr_values, num_bins, 'Normalization', 'probability', 'FaceAlpha', 0.5, 'DisplayName', 'Baseline');
        % histogram(cluster1_corr_values, num_bins, 'Normalization', 'probability', 'FaceAlpha', 0.5, 'DisplayName', 'Task');
        % xlabel('Correlation Values');
        % ylabel('Probability');
        % title('Cluster 1 Distribution of Correlation Values');
        % legend('show');
        % hold off;
        % 
        % subplot(2, 5, 4);
        % hold on;
        % plot(fpr, tpr, 'b', 'LineWidth', 2);
        % plot([0 1], [0 1], 'r--'); % Diagonal line for reference
        % xlabel('False Positive Rate');
        % ylabel('True Positive Rate');
        % title(sprintf('Cluster 1 ROC Curve (AUC = %.2f, p = %.3f)', cluster1_auc_bl_control, cluster1_p_bl_control));
        % hold off;

        cluster2_task_data = task_data(:, cluster2_trials);
        cluster2_baseline_mean = baseline_mean(:, cluster2_trials);
        cluster2_baseline_std = baseline_std(:, cluster2_trials);
        zscored_cluster2_baseline = (baseline_data(:, cluster2_trials) - cluster2_baseline_mean) ./ cluster2_baseline_std;

        zscored_cluster2_task = (cluster2_task_data - cluster2_baseline_mean) ./ cluster2_baseline_std;
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
            cluster2_baseline_mean) ./ cluster2_baseline_std;
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        

        % subplot(2, 5, [6 7]);
        % range_start = find(time_vector >= time_range_for_plotting_erp(1), 1);
        % range_end = find(time_vector <= time_range_for_plotting_erp(2), 1, 'last');
        % plot(time_vector(range_start:range_end), mean_erp(range_start:range_end), 'b', 'LineWidth', 2);
        % hold on;
        % fill([time_vector(range_start:range_end), fliplr(time_vector(range_start:range_end))], ...
        %  [mean_erp(range_start:range_end)' + std_error(range_start:range_end)', ...
        %  fliplr((mean_erp(range_start:range_end))' - (std_error(range_start:range_end))')], ...
        %  'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % xline(0, 'r--', 'LineWidth', 1.5);
        % xlabel('Time (s)');
        % ylabel('Amplitude (\muV)');
        % title(['Cluster 2 ERP n_trial: ' num2str(size(norm_concatenated_prep_signal, 2))], 'Interpreter', 'none');
        % legend({'ERP', 'Standard Error', 'saccade onset'});
        % grid on;
        % hold off;
        % 
        % subplot(2, 5, 8)
        % hold on;
        % num_bins = 20; % Specify the number of bins
        % histogram(cluster2_baseline_corr_values, num_bins, 'Normalization', 'probability', 'FaceAlpha', 0.5, 'DisplayName', 'Baseline');
        % histogram(cluster2_corr_values, num_bins, 'Normalization', 'probability', 'FaceAlpha', 0.5, 'DisplayName', 'Task');
        % xlabel('Correlation Values');
        % ylabel('Probability');
        % title('Cluster 2 Distribution of Correlation Values');
        % legend('show');
        % hold off;
        % 
        % subplot(2, 5, 9);
        % hold on;
        % plot(fpr, tpr, 'b', 'LineWidth', 2);
        % plot([0 1], [0 1], 'r--'); % Diagonal line for reference
        % xlabel('False Positive Rate');
        % ylabel('True Positive Rate');
        % title(sprintf('Cluster 2 ROC Curve (AUC = %.2f, p = %.3f)', cluster2_auc_bl_control, cluster2_p_bl_control));
        % hold off;
        
        % sgtitle(['Channel ', ch_label, ' Subclassification of baseline-control analysis']);

        % we will use baseline_corr later (we need to get the order back)
        baseline_corr = original_baseline_corr;
        %% baseline control done
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
            
            % identify trial with weakest connection
            mean_abs_corr = mean(abs(task_corr), 2);
            [~, weakest_trial_idx] = min(mean_abs_corr);       
            weak_trials_indices = [weak_trials_indices; original_order(weakest_trial_idx)];
            
            % Remove the weakest trial 
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
        
        % fig = figure('Visible','off');
        % set(fig, 'position', [100, 100, 900, 900])
        % 
        % subplot(2, 2, 1);
        % imagesc(original_task_corr);
        % colormap('coolwarm');
        % clim([-1, 1]);
        % axis square;
        % title('Original Task Correlation Matrix');
        % xlabel('Trials');
        % ylabel('Trials');
        % colorbar;
        % 
        % subplot(2, 2, 2);
        % dendrogram(Z, 0);
        % title('Dendrogram');
        % xlabel('Trials');
        % ylabel('Distance');
        % 
        % subplot(2, 2, 3);
        % imagesc(combined_matrix);
        % colormap('coolwarm');
        % clim([-1, 1]);
        % axis square;
        % title('Reordered Task Correlation Matrix with Clusters');
        % xlabel('Trials');
        % ylabel('Trials');
        % 
        % hold on;
        % clusters_2 = cluster(Z, 'maxclust', 2);
        % clusters_2 = clusters_2(order);
        % for i = 1:max(clusters_2)
        %     idx = find(clusters_2 == i);
        %     rectangle('Position', [idx(1)-0.5, idx(1)-0.5, length(idx), length(idx)], ...
        %               'EdgeColor', 'w', 'LineWidth', 2, 'LineStyle', '--');
        %     text(mean(idx), mean(idx), num2str(i), 'Color', 'w', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        % end
        % 
        % if ~isempty(weak_trials_indices)
        %     weak_idx_start = size(reordered_task_corr, 1) + 1;
        %     weak_idx_end = size(combined_matrix, 1);
        %     rectangle('Position', [weak_idx_start - 0.5, weak_idx_start - 0.5, weak_idx_end - weak_idx_start + 1, weak_idx_end - weak_idx_start + 1], ...
        %               'EdgeColor', 'k', 'LineWidth', 2, 'LineStyle', '--');
        %     text((weak_idx_start + weak_idx_end) / 2, (weak_idx_start + weak_idx_end) / 2, ...
        %          'Weak', 'Color', 'k', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        % end
        % 
        % subplot(2, 2, 4);
        % axis off; % Turn off the axis
        % c = colorbar('southoutside');
        % c.Label.String = 'Correlation Coefficient';
        % c.Label.FontSize = 12;
        % colormap('coolwarm');
        % clim([-1, 1]);

        % perform erp analysis for sub-classification
        if length(weak_trials_indices) >1
            weak_task_data = task_data(:, weak_trials_indices);
            
            weak_baseline_mean = baseline_mean(:, weak_trials_indices);
            weak_baseline_std = baseline_std(:, weak_trials_indices);
            zscored_weak_baseline = (baseline_data(:, weak_trials_indices) - weak_baseline_mean) ./ weak_baseline_std;

            zscored_weak_task = (weak_task_data - weak_baseline_mean) ./ weak_baseline_std;
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
        zscored_cluster1_baseline = (baseline_data(:, cluster1_trials) - cluster1_baseline_mean) ./ cluster1_baseline_std;

        zscored_cluster1_task = (cluster1_task_data - cluster1_baseline_mean) ./ cluster1_baseline_std;
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
            cluster1_baseline_mean) ./ cluster1_baseline_std;
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        
        time_vector = (0:length(mean_erp) - 1) / sfreq - trial_begin_study;

        % enable below for figures
        % fig = figure('visible', 'off');
        % set(fig, 'position', [100, 100, 1500, 700])
        % range_start = find(time_vector >= time_range_for_plotting_erp(1), 1);
        % range_end = find(time_vector <= time_range_for_plotting_erp(2), 1, 'last');
        % subplot(2, 5, [1 2]);
        % plot(time_vector(range_start:range_end), mean_erp(range_start:range_end), 'b', 'LineWidth', 2);
        % hold on;
        % fill([time_vector(range_start:range_end), fliplr(time_vector(range_start:range_end))], ...
        %  [mean_erp(range_start:range_end)' + std_error(range_start:range_end)', ...
        %  fliplr((mean_erp(range_start:range_end))' - (std_error(range_start:range_end))')], ...
        %  'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % xline(0, 'r--', 'LineWidth', 1.5);
        % 
        % xlabel('Time (s)');
        % ylabel('Amplitude (\muV)');
        % title(['Cluster 1 ERP n_trial: ' num2str(size(norm_concatenated_prep_signal, 2))], 'Interpreter', 'none');
        % legend({'ERP', 'Standard Error', 'saccade onset'});
        % grid on;
        % hold off;
        % 
        % subplot(2, 5, 3)
        % hold on;
        % num_bins = 20; % Specify the number of bins
        % histogram(cluster1_baseline_corr_values, num_bins, 'Normalization', 'probability', 'FaceAlpha', 0.5, 'DisplayName', 'Baseline');
        % histogram(cluster1_corr_values, num_bins, 'Normalization', 'probability', 'FaceAlpha', 0.5, 'DisplayName', 'Task');
        % xlabel('Correlation Values');
        % ylabel('Probability');
        % title('Cluster 1 Distribution of Correlation Values');
        % legend('show');
        % hold off;
        % 
        % subplot(2, 5, 4);
        % hold on;
        % plot(fpr, tpr, 'b', 'LineWidth', 2);
        % plot([0 1], [0 1], 'r--'); % Diagonal line for reference
        % xlabel('False Positive Rate');
        % ylabel('True Positive Rate');
        % title(sprintf('Cluster 1 ROC Curve (AUC = %.2f, p = %.3f)', cluster1_auc, cluster1_p));
        % hold off;

        horizontal_angle_cluster1 = saccade_table_all_session(cluster1_trials, :).azimuth_after - ...
            saccade_table_all_session(cluster1_trials, :).azimuth_before;
        horizontal_angle_cluster2 = saccade_table_all_session(cluster2_trials, :).azimuth_after - ...
            saccade_table_all_session(cluster2_trials, :).azimuth_before;
        vertical_angle_cluster1 = saccade_table_all_session(cluster1_trials, :).elevation_after - ...
            saccade_table_all_session(cluster1_trials, :).elevation_before;
        vertical_angle_cluster2 = saccade_table_all_session(cluster2_trials, :).elevation_after - ...
            saccade_table_all_session(cluster2_trials, :).elevation_before;
        % subplot(2, 5, 5);
        % % 3d surf plot
        % [counts, centers1] = hist3([horizontal_angle_cluster1, ...
        %     vertical_angle_cluster1], 'Nbins', [10, 10]);

        % [counts2, centers2] = hist3([horizontal_angle_cluster2, ...
        %     vertical_angle_cluster2], 'Nbins', [10, 10]);
        % 
        % scatter(horizontal_angle_cluster1, vertical_angle_cluster1, ...
        %     100, 'r', 'filled', 'MarkerFaceAlpha', 0.1);
        % hold on;
        % % Cluster 2
        % scatter(horizontal_angle_cluster2, vertical_angle_cluster2, ...
        %     100, 'b', 'filled', 'MarkerFaceAlpha', 0.1);
        % title('C1/2 saccade angle');
        % xlabel('Azimuth (degrees)');
        % ylabel('Elevation (degrees)');
        % % zlabel('Density');
        % legend({'C1', 'C2'});
        % hold off;
        % 
        % subplot(2, 5, 10);
        % in this subplot, we will do polarhist for the angle between
        % horizontal and verticle, use different colors to mark cluster 1
        % or 2, angles1 is in degree between -180 and 180, 0 means right,
        % -90 means down
        angles1 = atan2d(vertical_angle_cluster1, horizontal_angle_cluster1);
        angles2 = atan2d(vertical_angle_cluster2, horizontal_angle_cluster2);
        
        % % Polar histogram for Cluster 1
        % polarhistogram(deg2rad(angles1), 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        % hold on;
        % polarhistogram(deg2rad(angles2), 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        % 
        % ax = gca;
        % ax.ThetaTick = [0 90 180 270];
        % ax.ThetaTickLabel = {'Right', 'Up', 'Left', 'Down'};
        % if verticle angle = 1, horizongtal angle = 0, angles1 = 90, which
        % is up
        % show the distance vector and shank direction vector
        % Calculate azimuth (in degrees)
        % shank_angle = atan2d(atan2d(shank_direction(3), sqrt(shank_direction(1)^2 + shank_direction(2)^2)),...
        %     atan2d(shank_direction(2), shank_direction(1)));
        
        % Calculate elevation (in degrees)
        % phi = atan2d(shank_direction(3), sqrt(shank_direction(1)^2 + shank_direction(2)^2));
        % shank_angle = atan2d(shank_elevation, shank_azimuth);
        % polarplot([0, deg2rad(shank_angle)], [0, max(max(sum(counts,1)), max(sum(counts2,1)))], 'k-', 'LineWidth', 2);
        % if exist('shank_direction', 'var') && ~isempty(shank_direction)
        %     % Calculate the shank angle (in degrees)
        %     shank_angle = atan2d(...
        %         atan2d(shank_direction(3), sqrt(shank_direction(1)^2 + shank_direction(2)^2)), ...
        %         atan2d(shank_direction(2), shank_direction(1)) ...
        %     );
        % 
        %     % Plot the distance vector and shank direction vector
        %     polarplot([0, deg2rad(shank_angle)], ...
        %               [0, max(max(sum(counts, 1)), max(sum(counts2, 1)))], ...
        %               'k-', 'LineWidth', 2);
        %     title('Polar Histogram of atan(v,h)');
        % else
        %     title('Polar Hist of atan(v,h) - not linear reref');
        % end
        % legend('C1', 'C2');

       

        cluster2_task_data = task_data(:, cluster2_trials);
        cluster2_baseline_mean = baseline_mean(:, cluster2_trials);
        cluster2_baseline_std = baseline_std(:, cluster2_trials);
        zscored_cluster2_baseline = (baseline_data(:, cluster2_trials) - cluster2_baseline_mean) ./ cluster2_baseline_std;

        zscored_cluster2_task = (cluster2_task_data - cluster2_baseline_mean) ./ cluster2_baseline_std;
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
            cluster2_baseline_mean) ./ cluster2_baseline_std;
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        

        % subplot(2, 5, [6 7]);
        % range_start = find(time_vector >= time_range_for_plotting_erp(1), 1);
        % range_end = find(time_vector <= time_range_for_plotting_erp(2), 1, 'last');
        % plot(time_vector(range_start:range_end), mean_erp(range_start:range_end), 'b', 'LineWidth', 2);
        % hold on;
        % fill([time_vector(range_start:range_end), fliplr(time_vector(range_start:range_end))], ...
        %  [mean_erp(range_start:range_end)' + std_error(range_start:range_end)', ...
        %  fliplr((mean_erp(range_start:range_end))' - (std_error(range_start:range_end))')], ...
        %  'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % xline(0, 'r--', 'LineWidth', 1.5);
        % xlabel('Time (s)');
        % ylabel('Amplitude (\muV)');
        % title(['Cluster 2 ERP n_trial: ' num2str(size(norm_concatenated_prep_signal, 2))], 'Interpreter', 'none');
        % legend({'ERP', 'Standard Error', 'saccade onset'});
        % grid on;
        % hold off;
        % 
        % subplot(2, 5, 8)
        % hold on;
        % num_bins = 20; % Specify the number of bins
        % histogram(cluster2_baseline_corr_values, num_bins, 'Normalization', 'probability', 'FaceAlpha', 0.5, 'DisplayName', 'Baseline');
        % histogram(cluster2_corr_values, num_bins, 'Normalization', 'probability', 'FaceAlpha', 0.5, 'DisplayName', 'Task');
        % xlabel('Correlation Values');
        % ylabel('Probability');
        % title('Cluster 2 Distribution of Correlation Values');
        % legend('show');
        % hold off;
        % 
        % subplot(2, 5, 9);
        % hold on;
        % plot(fpr, tpr, 'b', 'LineWidth', 2);
        % plot([0 1], [0 1], 'r--'); % Diagonal line for reference
        % xlabel('False Positive Rate');
        % ylabel('True Positive Rate');
        % title(sprintf('Cluster 2 ROC Curve (AUC = %.2f, p = %.3f)', cluster2_auc, cluster2_p));
        % hold off;
        % 
        % sgtitle(['Channel ', ch_label, ' Saccade Directions and ERP Analysis']);


        % plot individual trials
        % if plot_individual_trial
        %     fig = figure('visible', 'off');
        %     set(fig, 'position', [100, 100, 1000, 800])
        %     norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster1_trials) - ...
        %         cluster1_baseline_mean) ./ cluster1_baseline_std;
        %     mean_erp = mean(norm_concatenated_prep_signal, 2);
        %     std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        % 
        %     time_vector = (0:length(mean_erp) - 1) / sfreq - trial_begin_study;
        %     range_start = find(time_vector >= time_range_for_plotting_erp(1), 1);
        %     range_end = find(time_vector <= time_range_for_plotting_erp(2), 1, 'last');
        %     subplot(2, 1, 1);
        %     hold on;
        %     for trial = 1:size(norm_concatenated_prep_signal, 2)
        %         plot(time_vector(range_start:range_end), norm_concatenated_prep_signal(range_start:range_end, trial), ...
        %             'Color', [0.6, 0.6, 1, 0.05], 'LineWidth', 1);  % Semi-transparent blue lines
        %     end
        %     plot(time_vector(range_start:range_end), mean_erp(range_start:range_end), 'b', 'LineWidth', 4);
        %     % hold on;
        %     % fill([time_vector(range_start:range_end), fliplr(time_vector(range_start:range_end))], ...
        %     %  [mean_erp(range_start:range_end)' + std_error(range_start:range_end)', ...
        %     %  fliplr((mean_erp(range_start:range_end))' - (std_error(range_start:range_end))')], ...
        %     %  'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        %     xline(0, 'r--', 'LineWidth', 1.5);
        %     xlabel('Time (s)');
        %     ylabel('Amplitude (\muV)');
        %     title(['Cluster 1 ERP n_trial: ' num2str(size(norm_concatenated_prep_signal, 2))], 'Interpreter', 'none');
        %     legend({'ERP', 'Standard Error', 'saccade onset'});
        %     grid on;
        %     hold off;
        %     ylim([-3.5, 3.5]);
        % 
        %     norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster2_trials) - ...
        %         cluster2_baseline_mean) ./ cluster2_baseline_std;
        %     mean_erp = mean(norm_concatenated_prep_signal, 2);
        %     std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
        % 
        %     subplot(2, 1, 2);
        %     range_start = find(time_vector >= time_range_for_plotting_erp(1), 1);
        %     range_end = find(time_vector <= time_range_for_plotting_erp(2), 1, 'last');
        %     hold on;
        %     for trial = 1:size(norm_concatenated_prep_signal, 2)
        %         plot(time_vector(range_start:range_end), norm_concatenated_prep_signal(range_start:range_end, trial), ...
        %             'Color', [0.6, 0.6, 1, 0.05], 'LineWidth', 1);  % Semi-transparent blue lines
        %     end
        %     plot(time_vector(range_start:range_end), mean_erp(range_start:range_end), 'b', 'LineWidth', 4);
        % 
        %     % fill([time_vector(range_start:range_end), fliplr(time_vector(range_start:range_end))], ...
        %     %  [mean_erp(range_start:range_end)' + std_error(range_start:range_end)', ...
        %     %  fliplr((mean_erp(range_start:range_end))' - (std_error(range_start:range_end))')], ...
        %     %  'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        %     xline(0, 'r--', 'LineWidth', 1.5);
        %     xlabel('Time (s)');
        %     ylabel('Amplitude (\muV)');
        %     title(['Cluster 2 ERP n_trial: ' num2str(size(norm_concatenated_prep_signal, 2))], 'Interpreter', 'none');
        %     legend({'ERP', 'Standard Error', 'saccade onset'});
        %     grid on;
        %     hold off;
        %     ylim([-3.5, 3.5]);
        %     sgtitle(['Channel ', ch_label, ' trial - ERP Analysis']);
           
            % we need to create a new figure, use two color to represent the two clusters
            % fig = figure('visible', 'off');
            % set(fig, 'position', [100, 100, 1000, 600])
            % norm_concatenated_prep_signal = (prep_signal_all_session - ...
            %     baseline_mean) ./ baseline_std;
            % mean_erp = mean(norm_concatenated_prep_signal, 2);
            % % std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
            % 
            % time_vector = (0:length(mean_erp) - 1) / sfreq - trial_begin_study;
            % range_start = find(time_vector >= time_range_for_plotting_erp(1), 1);
            % range_end = find(time_vector <= (time_range_for_plotting_erp(2) + 0.1), 1, 'last');
            % hold on;
            % plot(time_vector(range_start:range_end), mean_erp(range_start:range_end), 'k', 'LineWidth', 4);
            % 
            % % plot cluster 1
            % norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster1_trials) - ...
            %     cluster1_baseline_mean) ./ cluster1_baseline_std;
            % mean_erp = mean(norm_concatenated_prep_signal, 2);
            % std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
            % 
            % time_vector = (0:length(mean_erp) - 1) / sfreq - trial_begin_study;
            % % range_start = find(time_vector >= time_range_for_plotting_erp(1), 1);
            % % range_end = find(time_vector <= (time_range_for_plotting_erp(2) + 0.2), 1, 'last');
            % hold on;
            % for trial = 1:size(norm_concatenated_prep_signal, 2)
            %     plot(time_vector(range_start:range_end), norm_concatenated_prep_signal(range_start:range_end, trial), ...
            %         'Color', 'r', 'LineWidth', 0.5);  % Semi-transparent blue lines
            % end
            % plot(time_vector(range_start:range_end), mean_erp(range_start:range_end), 'r', 'LineWidth', 4);
            % 
            % % plot cluster 2
            % norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster2_trials) - ...
            %     cluster2_baseline_mean) ./ cluster2_baseline_std;
            % mean_erp = mean(norm_concatenated_prep_signal, 2);
            % std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
            % 
            % % range_start = find(time_vector >= time_range_for_plotting_erp(1), 1);
            % % range_end = find(time_vector <= time_range_for_plotting_erp(2), 1, 'last');
            % hold on;
            % for trial = 1:size(norm_concatenated_prep_signal, 2)
            %     plot(time_vector(range_start:range_end), norm_concatenated_prep_signal(range_start:range_end, trial), ...
            %         'Color', 'b', 'LineWidth', 0.5);  % Semi-transparent blue lines
            % end
            % plot(time_vector(range_start:range_end), mean_erp(range_start:range_end), 'b', 'LineWidth', 4);
            % xline(0, 'r--', 'LineWidth', 1.5);
            % xlabel('Time (s)');
            % ylabel('Amplitude (\muV)');
            % sgtitle(['Channel ', ch_label, ' trial trace- ERP Analysis r - cluster 2']);
            % grid on;
            % hold off;
            % ylim([-3.5, 3.5]);
        % end


        % see what determines the polarity
        % fig = figure('visible', 'off');
        % set(fig, 'position', [100, 100, 1500, 700])
        % subplot(2,3,1)
        % remember_cluster1 = sum(strcmp(saccade_table_all_session(cluster1_trials, :).remember_condition, 'remember'));
        % remember_cluster2 = sum(strcmp(saccade_table_all_session(cluster2_trials, :).remember_condition, 'remember'));
        % forgotten_cluster1 = sum(strcmp(saccade_table_all_session(cluster1_trials, :).remember_condition, 'forgotten'));
        % forgotten_cluster2 = sum(strcmp(saccade_table_all_session(cluster2_trials, :).remember_condition, 'forgotten'));
        % total_cluster1 = remember_cluster1 + forgotten_cluster1;
        % total_cluster2 = remember_cluster2 + forgotten_cluster2;
        % 
        % percent_remember_cluster1 = (remember_cluster1 / total_cluster1) * 100;
        % percent_forgotten_cluster1 = (forgotten_cluster1 / total_cluster1) * 100;
        % percent_remember_cluster2 = (remember_cluster2 / total_cluster2) * 100;
        % percent_forgotten_cluster2 = (forgotten_cluster2 / total_cluster2) * 100;
        % 
        % conf_matrix = [remember_cluster1, remember_cluster2; forgotten_cluster1, forgotten_cluster2];
        % percent_matrix = [percent_remember_cluster1, percent_remember_cluster2; percent_forgotten_cluster1, percent_forgotten_cluster2];
        % 
        % h = heatmap({'Cluster 1', 'Cluster 2'}, {'Remember', 'Forgotten'}, conf_matrix, ...
        %   'Colormap', sky, 'ColorbarVisible', 'on');
        % 
        % h.CellLabelColor = 'none';
        % hs = struct(h);
        % hs.Colorbar.Label.String = 'Number of Trials';
        % title('Cluster and Remember');
        % ax = gca;
        % secondaryAx = axes('Position', ax.Position, 'Color', 'none', 'XAxisLocation', ...
        %     'top', 'YAxisLocation', 'right', 'XTick', [], 'YTick', []);
        % textStrings = strcat(num2str(conf_matrix(:)), ' (', num2str(percent_matrix(:), '%.1f'), '%)');
        % textStrings = strtrim(cellstr(textStrings));
        % [x, y] = meshgrid(1:2);
        % x = x(:);
        % y = y(:);
        % y = 3 - y;
        % for k = 1:length(textStrings)
        %     text(x(k), y(k), textStrings{k}, 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', 'Parent', secondaryAx);
        % end
        % secondaryAx.XLim = [0.7, 2.7];
        % secondaryAx.YLim = [0.5, 2.5];
        % 
        % subplot(2,3,2)
        % like_cluster1 = sum(strcmp(saccade_table_all_session(cluster1_trials, :).behavior_condition, 'like'));
        % like_cluster2 = sum(strcmp(saccade_table_all_session(cluster2_trials, :).behavior_condition, 'like'));
        % dislike_cluster1 = sum(strcmp(saccade_table_all_session(cluster1_trials, :).behavior_condition, 'dislike'));
        % dislike_cluster2 = sum(strcmp(saccade_table_all_session(cluster2_trials, :).behavior_condition, 'dislike'));
        % total_cluster1 = like_cluster1 + dislike_cluster1;
        % total_cluster2 = like_cluster2 + dislike_cluster2;
        % 
        % percent_remember_cluster1 = (like_cluster1 / total_cluster1) * 100;
        % percent_forgotten_cluster1 = (dislike_cluster1 / total_cluster1) * 100;
        % percent_remember_cluster2 = (like_cluster2 / total_cluster2) * 100;
        % percent_forgotten_cluster2 = (dislike_cluster2 / total_cluster2) * 100;
        % 
        % conf_matrix = [like_cluster1, like_cluster2; dislike_cluster1, dislike_cluster2];
        % percent_matrix = [percent_remember_cluster1, percent_remember_cluster2; percent_forgotten_cluster1, percent_forgotten_cluster2];
        % 
        % h = heatmap({'Cluster 1', 'Cluster 2'}, {'Like', 'Dislike'}, conf_matrix, ...
        %   'Colormap', sky, 'ColorbarVisible', 'on');
        % 
        % h.CellLabelColor = 'none';
        % hs = struct(h);
        % hs.Colorbar.Label.String = 'Number of Trials';
        % title('Cluster and Like');
        % ax = gca;
        % secondaryAx = axes('Position', ax.Position, 'Color', 'none', 'XAxisLocation', ...
        %     'top', 'YAxisLocation', 'right', 'XTick', [], 'YTick', []);
        % textStrings = strcat(num2str(conf_matrix(:)), ' (', num2str(percent_matrix(:), '%.1f'), '%)');
        % textStrings = strtrim(cellstr(textStrings));
        % [x, y] = meshgrid(1:2);
        % x = x(:);
        % y = y(:);
        % y = 3 - y;
        % for k = 1:length(textStrings)
        %     text(x(k), y(k), textStrings{k}, 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', 'Parent', secondaryAx);
        % end
        % secondaryAx.XLim = [0.7, 2.7];
        % secondaryAx.YLim = [0.5, 2.5];
        % 
        % subplot(2,3,3)
        gaze_x_change_cluster1 = saccade_table_all_session(cluster1_trials, :).gaze_after_x_screen_coord - ...
            saccade_table_all_session(cluster1_trials, :).gaze_before_x_screen_coord;
        gaze_x_change_cluster2 = saccade_table_all_session(cluster2_trials, :).gaze_after_x_screen_coord - ...
            saccade_table_all_session(cluster2_trials, :).gaze_before_x_screen_coord;
        gaze_y_change_cluster1 = saccade_table_all_session(cluster1_trials, :).gaze_after_y_screen_coord - ...
            saccade_table_all_session(cluster1_trials, :).gaze_before_y_screen_coord;
        gaze_y_change_cluster2 = saccade_table_all_session(cluster2_trials, :).gaze_after_y_screen_coord - ...
            saccade_table_all_session(cluster2_trials, :).gaze_before_y_screen_coord;

        % scatter(gaze_x_change_cluster1, gaze_y_change_cluster1, ...
        %     100, 'r', 'filled', 'MarkerFaceAlpha', 0.1);
        % hold on;
        % % Cluster 2
        % scatter(gaze_x_change_cluster2, gaze_y_change_cluster2, ...
        %     100, 'b', 'filled', 'MarkerFaceAlpha', 0.1);
        % title('C1/2 saccade angle');
        % xlabel('delta x (screen)');
        % ylabel('delta y (screen)');
        % % zlabel('Density');
        % legend({'C1', 'C2'});
        % hold off;
        % 
        % subplot(2,3,4)
        % time_to_image_onset_cluster1 = saccade_table_all_session(cluster1_trials, :).time_to_image_onset;
        % time_to_image_onset_cluster2 = saccade_table_all_session(cluster2_trials, :).time_to_image_onset;
        % hold on;
        % num_bins = 20; 
        % histogram(time_to_image_onset_cluster1, num_bins, 'Normalization', ...
        %     'probability', 'FaceAlpha', 0.5, 'DisplayName', 'C1');
        % histogram(time_to_image_onset_cluster2, num_bins, 'Normalization', 'probability', ...
        %     'FaceAlpha', 0.5, 'DisplayName', 'C2');
        % xlabel('time to image onset');
        % ylabel('Probability');
        % title('Distribution of time');
        % legend('show');
        % hold off;
        % 
        % subplot(2,3,5)
        % gaze_x_before_cluster1 = saccade_table_all_session(cluster1_trials, :).gaze_before_x_screen_coord;
        % gaze_x_before_cluster2 = saccade_table_all_session(cluster2_trials, :).gaze_before_x_screen_coord;
        % gaze_y_before_cluster1 = saccade_table_all_session(cluster1_trials, :).gaze_before_y_screen_coord;
        % gaze_y_before_cluster2 = saccade_table_all_session(cluster2_trials, :).gaze_before_y_screen_coord;
        % 
        % scatter(gaze_x_before_cluster1, gaze_y_before_cluster1, ...
        %     100, 'r', 'filled', 'MarkerFaceAlpha', 0.1);
        % hold on;
        % % Cluster 2
        % scatter(gaze_x_before_cluster2, gaze_y_before_cluster2, ...
        %     100, 'b', 'filled', 'MarkerFaceAlpha', 0.1);
        % title('C1/2 initial saccade location');
        % xlabel('x (screen)');
        % ylabel('y (screen)');
        % % zlabel('Density');
        % legend({'C1', 'C2'});
        % hold off;
        % set(gca, 'YDir', 'reverse');
        % 
        % subplot(2,3,6)
        gaze_x_before_cluster1 = saccade_table_all_session(cluster1_trials, :).azimuth_before;
        gaze_x_before_cluster2 = saccade_table_all_session(cluster2_trials, :).azimuth_before;
        gaze_y_before_cluster1 = saccade_table_all_session(cluster1_trials, :).elevation_before;
        gaze_y_before_cluster2 = saccade_table_all_session(cluster2_trials, :).elevation_before;
        % scatter(gaze_x_before_cluster1, gaze_y_before_cluster1, ...
        %     100, 'r', 'filled', 'MarkerFaceAlpha', 0.1);
        % hold on;
        % % Cluster 2
        % scatter(gaze_x_before_cluster2, gaze_y_before_cluster2, ...
        %     100, 'b', 'filled', 'MarkerFaceAlpha', 0.1);
        % title('C1/2 initial saccade location');
        % xlabel('azimuth (eye coord)');
        % ylabel('elevation (eye coord)');
        % % zlabel('Density');
        % legend({'C1', 'C2'});
        % hold off;
        % set(gca, 'XDir', 'reverse');
        % sgtitle(['Channel ', ch_label, ' factor for polarity']);


        time_vector = (0:size(prep_signal_all_session, 1) - 1) / sfreq - trial_begin_study;
        norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster1_trials) - ...
            cluster1_baseline_mean) ./ cluster1_baseline_std;
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
       
        % fig = figure('Visible', 'off');
        % set(fig, 'position', [100, 100, 1500, 700])
        % range_start = find(time_vector >= erp_pre_during_post(1), 1);
        % range_end = find(time_vector <= erp_pre_during_post(end), 1, 'last');
        % subplot(2, 4, 1);
        % plot(time_vector(range_start:range_end), mean_erp(range_start:range_end), 'b', 'LineWidth', 2);
        % hold on;
        % fill([time_vector(range_start:range_end), fliplr(time_vector(range_start:range_end))], ...
        %  [mean_erp(range_start:range_end)' + std_error(range_start:range_end)', ...
        %  fliplr((mean_erp(range_start:range_end))' - (std_error(range_start:range_end))')], ...
        %  'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % xline(0, 'r--', 'LineWidth', 1.5);
        % 
        % xlabel('Time (s)');
        % ylabel('Amplitude (\muV)');
        % title(['Cluster 1 ERP n_trial: ' num2str(size(norm_concatenated_prep_signal, 2))], 'Interpreter', 'none');
        % legend({'ERP', 'Standard Error', 'saccade onset'});
        % grid on;
        % hold off;
        % 
        % subplot(2, 4, 5);
        % imagesc(time_vector(range_start:range_end), 1:size(norm_concatenated_prep_signal, 2), ...
        %     norm_concatenated_prep_signal(range_start: range_end, :)');
        % colorbar;
        % title('Baseline Corrected Raw Signal');
        % ylabel('Trial');
        % xlabel('Time (s)');
        % clim(color_limits);

        % compute overall psd
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
        
        % subplot(2,4,2)
        % idx = f >= freq_range(1) & f <= freq_range(2);
        % plot(f(idx), 10*log10(mean_psd(idx)), 'b', 'LineWidth', 2);
        % hold on;
        % std_psd_db = 10*log10(mean_psd + std_psd) - 10*log10(mean_psd); 
        % 
        % upper_bound_db = 10*log10(mean_psd(idx)) + std_psd_db(idx);
        % lower_bound_db = 10*log10(mean_psd(idx)) - std_psd_db(idx);
        % 
        % fill([f(idx); flipud(f(idx))], [upper_bound_db; flipud(lower_bound_db)], ...
        % 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % 
        % plot(f(idx), 10*log10(fit_line(idx)), 'r--', 'LineWidth', 1.5);
        % legend({'Mean PSD', 'Standard Error', '1/f Fit'});
        % 
        % xlabel('Frequency (Hz)');
        % ylabel('Power/Frequency (dB/Hz)');
        % title('Power Spectral Density (PSD)');
        % grid on;
        % hold off;
        % legend({'Mean PSD', 'Standard Error', '1/f Fit'});
        % text_x = f(find(idx, 1, 'last')); 
        % text_y = 10*log10(fit_line(find(idx, 1, 'last')));
        % fit_text = sprintf('psd = %.2f*log(f) + %.2f', p(1), p(2));
        % text(text_x, text_y, fit_text, 'Color', 'r', 'FontSize', 12, 'HorizontalAlignment', 'right');

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

        % subplot(2, 4, 6);
        % imagesc(time_vector_spec(find(time_vector_spec >= erp_pre_during_post(1), 1):find(time_vector_spec <= erp_pre_during_post(end), 1, 'last')), ...
        %     F, baseline_corrected_power(:, find(time_vector_spec >= erp_pre_during_post(1), 1):find(time_vector_spec <= erp_pre_during_post(end), 1, 'last')));
        % axis xy;
        % ylim(freq_range);
        % xlabel('Time (s)');
        % ylabel('Frequency (Hz)');
        % title('Baseline Corrected Time-Frequency Plot');
        % cb = colorbar;
        % clim([-2, 2]);
        % ylabel(cb, 'Power (dB)'); 

      
        valid_idx = (freq_range_fit(2) > f) & (f > low_cutoff);
        f_filtered = f(valid_idx);
        mean_psd_filtered = 10*log10(mean_psd((valid_idx)));
        fit_line_filtered = 10*log10(fit_line(valid_idx));
        diff_psd = mean_psd_filtered - fit_line_filtered;
        [cluster1_oscillation_amp, peak_idx] = max(diff_psd);
        cluster1_peak_freq = f_filtered(peak_idx);
        low_cutoff_oscillation = max(0, cluster1_peak_freq - narrow_bandwidth);
        high_cutoff_oscillation = cluster1_peak_freq + narrow_bandwidth;

        % get power and phase
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
        % baseline_mean = mean([baseline_data; task_data], 1);
        baseline_amp_mean = mean(baseline_amp, 1);
        baseline_amp_std = std(baseline_amp, [], 1);
        norm_amplitude = (amplitude - ...
            baseline_amp_mean) ./ baseline_amp_std;
        mean_erp = mean(norm_amplitude, 2);
        std_error = std(norm_amplitude, [], 2) / sqrt(size(norm_amplitude, 2));
        cluster1_4_12_hz_pow = mean_erp;

        % subplot(2, 4, 3);
        % plot(time_vector(range_start:range_end), mean_erp, 'b', 'LineWidth', 2);
        % hold on;
        % fill([time_vector(range_start:range_end), fliplr(time_vector(range_start:range_end))], ...
        %  [mean_erp' + std_error', ...
        %  fliplr((mean_erp)' - (std_error)')], ...
        %  'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % xline(0, 'r--', 'LineWidth', 1.5);
        % 
        % xlabel('Time (s)');
        % ylabel('Baseline corrected Power Amplitude (a.u.)');
        % title([sprintf('%.1f', low_cutoff_oscillation) '-' sprintf('%.1f', high_cutoff_oscillation) 'Hz power'], 'Interpreter', 'none');
        % legend({'ERP', 'Standard Error', 'saccade onset'});
        % grid on;
        % hold off;

        % subplot(2, 4, 7);
        % imagesc(time_vector(range_start:range_end), 1:size(norm_amplitude, 2), ...
        %     norm_amplitude');
        % colorbar;
        % title('Baseline Corrected Power');
        % ylabel('Trial');
        % xlabel('Time (s)');
        % clim(color_limits);
        % 
        mean_phase = circ_mean(phase_cluster1, [], 2); % Circular mean
        std_error_phase = circ_std(phase_cluster1, [], [], 2) / sqrt(size(phase_cluster1, 2)); % Circular standard error

        % Compute the PLV
        cluster1_4_12_hz_phase = mean_phase;
        plv = abs(mean(exp(1i * phase_cluster1), 2));

        cluster1_4_12_hz_plv = plv;
        cluster1_4_12_hz_plv_all = mean(plv);

        % subplot(2, 4, 4);
        % plot(time_vector(range_start:range_end), mean_phase, 'b', 'LineWidth', 2);
        % hold on;
        % plot(time_vector(range_start:range_end), plv, 'r', 'LineWidth', 2);
        % plot(time_vector(range_start:range_end), std_error_phase * 10, 'green', 'LineWidth', 2);
        % fill([time_vector(range_start:range_end), fliplr(time_vector(range_start:range_end))], ...
        %  [mean_phase' + std_error_phase', ...
        %  fliplr((mean_phase)' - (std_error_phase)')], ...
        %  'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % xline(0, 'r--', 'LineWidth', 1.5);
        % 
        % xlabel('Time (s)');
        % ylabel('Phase (rad)');
        % title([sprintf('%.1f', low_cutoff_oscillation) '-' sprintf('%.1f', high_cutoff_oscillation) 'Hz phase'], 'Interpreter', 'none');
        % legend({'ERP', 'PLV', ' var', 'Standard Error', 'saccade onset'});
        % grid on;
        % hold off;
        % xlim([-0.25, 0.25]);
        % 
        % subplot(2, 4, 8);
        % imagesc(time_vector(range_start:range_end), 1:size(phase_cluster1, 2), ...
        %     phase_cluster1');
        % colorbar;
        % title('Phase');
        % ylabel('Trial');
        % xlabel('Time (s)');
        % xlim([-0.25, 0.25]);
        % 
        % sgtitle(['Channel ', ch_label, ' oscillation Analysis']);


  % cluster 2
        time_vector = (0:size(prep_signal_all_session, 1) - 1) / sfreq - trial_begin_study;
        norm_concatenated_prep_signal = (prep_signal_all_session(:, cluster2_trials) - ...
            cluster2_baseline_mean) ./ cluster2_baseline_std;
        mean_erp = mean(norm_concatenated_prep_signal, 2);
        std_error = std(norm_concatenated_prep_signal, [], 2) / sqrt(size(norm_concatenated_prep_signal, 2));
       
        % fig = figure('Visible', 'off');
        % set(fig, 'position', [100, 100, 1500, 700])
        % range_start = find(time_vector >= erp_pre_during_post(1), 1);
        % range_end = find(time_vector <= erp_pre_during_post(end), 1, 'last');
        % subplot(2, 4, 1);
        % plot(time_vector(range_start:range_end), mean_erp(range_start:range_end), 'b', 'LineWidth', 2);
        % hold on;
        % fill([time_vector(range_start:range_end), fliplr(time_vector(range_start:range_end))], ...
        %  [mean_erp(range_start:range_end)' + std_error(range_start:range_end)', ...
        %  fliplr((mean_erp(range_start:range_end))' - (std_error(range_start:range_end))')], ...
        %  'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % xline(0, 'r--', 'LineWidth', 1.5);
        % 
        % xlabel('Time (s)');
        % ylabel('Amplitude (\muV)');
        % title(['Cluster 2 ERP n_trial: ' num2str(size(norm_concatenated_prep_signal, 2))], 'Interpreter', 'none');
        % legend({'ERP', 'Standard Error', 'saccade onset'});
        % grid on;
        % hold off;
        % 
        % subplot(2, 4, 5);
        % imagesc(time_vector(range_start:range_end), 1:size(norm_concatenated_prep_signal, 2), ...
        %     norm_concatenated_prep_signal(range_start: range_end, :)');
        % colorbar;
        % title('Baseline Corrected Raw Signal');
        % ylabel('Trial');
        % xlabel('Time (s)');
        % clim(color_limits);

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
        
        % subplot(2,4,2)
        % idx = f >= freq_range(1) & f <= freq_range(2);
        % plot(f(idx), 10*log10(mean_psd(idx)), 'b', 'LineWidth', 2);
        % hold on;
        % std_psd_db = 10*log10(mean_psd + std_psd) - 10*log10(mean_psd); 
        % 
        % upper_bound_db = 10*log10(mean_psd(idx)) + std_psd_db(idx);
        % lower_bound_db = 10*log10(mean_psd(idx)) - std_psd_db(idx);
        % 
        % fill([f(idx); flipud(f(idx))], [upper_bound_db; flipud(lower_bound_db)], ...
        % 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % 
        % plot(f(idx), 10*log10(fit_line(idx)), 'r--', 'LineWidth', 1.5);
        % legend({'Mean PSD', 'Standard Error', '1/f Fit'});
        % 
        % xlabel('Frequency (Hz)');
        % ylabel('Power/Frequency (dB/Hz)');
        % title('Power Spectral Density (PSD)');
        % grid on;
        % hold off;
        % legend({'Mean PSD', 'Standard Error', '1/f Fit'});
        % text_x = f(find(idx, 1, 'last')); 
        % text_y = 10*log10(fit_line(find(idx, 1, 'last')));
        % fit_text = sprintf('psd = %.2f*log(f) + %.2f', p(1), p(2));
        % text(text_x, text_y, fit_text, 'Color', 'r', 'FontSize', 12, 'HorizontalAlignment', 'right');

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

        % subplot(2, 4, 6);
        % imagesc(time_vector_spec(find(time_vector_spec >= erp_pre_during_post(1), 1):find(time_vector_spec <= erp_pre_during_post(end), 1, 'last')), ...
        %     F, baseline_corrected_power(:, find(time_vector_spec >= erp_pre_during_post(1), 1):find(time_vector_spec <= erp_pre_during_post(end), 1, 'last')));
        % axis xy;
        % ylim(freq_range);
        % xlabel('Time (s)');
        % ylabel('Frequency (Hz)');
        % title('Baseline Corrected Time-Frequency Plot');
        % cb = colorbar;
        % clim([-2, 2]);
        % ylabel(cb, 'Power (dB)'); 

        valid_idx = (freq_range_fit(2) > f) & (f > low_cutoff);
        f_filtered = f(valid_idx);
        mean_psd_filtered = 10*log10(mean_psd((valid_idx)));
        fit_line_filtered = 10*log10(fit_line(valid_idx));
        diff_psd = mean_psd_filtered - fit_line_filtered;
        [cluster2_oscillation_amp, peak_idx] = max(diff_psd);
        cluster2_peak_freq = f_filtered(peak_idx);
        low_cutoff_oscillation = max(0, cluster2_peak_freq - narrow_bandwidth);
        high_cutoff_oscillation = cluster2_peak_freq + narrow_bandwidth;

        % get power and phrase
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

        % subplot(2, 4, 3);
        % plot(time_vector(range_start:range_end), mean_erp, 'b', 'LineWidth', 2);
        % hold on;
        % fill([time_vector(range_start:range_end), fliplr(time_vector(range_start:range_end))], ...
        %  [mean_erp' + std_error', ...
        %  fliplr((mean_erp)' - (std_error)')], ...
        %  'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % xline(0, 'r--', 'LineWidth', 1.5);
        % 
        % xlabel('Time (s)');
        % ylabel('Baseline corrected Power Amplitude (a.u.)');
        % title([sprintf('%.1f', low_cutoff_oscillation) '-' sprintf('%.1f', high_cutoff_oscillation) 'Hz power'], 'Interpreter', 'none');
        % legend({'ERP', 'Standard Error', 'saccade onset'});
        % grid on;
        % hold off;
        % 
        % subplot(2, 4, 7);
        % imagesc(time_vector(range_start:range_end), 1:size(norm_amplitude, 2), ...
        %     norm_amplitude');
        % colorbar;
        % title('Baseline Corrected Power');
        % ylabel('Trial');
        % xlabel('Time (s)');
        % clim(color_limits);

        mean_phase = circ_mean(phase_cluster2, [], 2); % Circular mean
        std_error_phase = circ_std(phase_cluster2, [], [], 2) / sqrt(size(phase_cluster2, 2)); % Circular standard error

        % Compute the PLV
        cluster2_4_12_hz_phase = mean_phase;

        plv = abs(mean(exp(1i * phase_cluster2), 2));

        cluster2_4_12_hz_plv = plv;
        cluster2_4_12_hz_plv_all = mean(plv);

        all_phase = [phase_cluster1'; phase_cluster2']';
        all_4_12_hz_plv = abs(mean(exp(1i * all_phase), 2));

        % subplot(2, 4, 4);
        % plot(time_vector(range_start:range_end), mean_phase, 'b', 'LineWidth', 2);
        % hold on;
        % plot(time_vector(range_start:range_end), plv, 'r', 'LineWidth', 2);
        % plot(time_vector(range_start:range_end), std_error_phase * 10, 'green', 'LineWidth', 2);
        % fill([time_vector(range_start:range_end), fliplr(time_vector(range_start:range_end))], ...
        %  [mean_phase' + std_error_phase', ...
        %  fliplr((mean_phase)' - (std_error_phase)')], ...
        %  'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % xline(0, 'r--', 'LineWidth', 1.5);
        % 
        % xlabel('Time (s)');
        % ylabel('Phase (rad)');
        % title([sprintf('%.1f', low_cutoff_oscillation) '-' sprintf('%.1f', high_cutoff_oscillation) 'Hz phase'], 'Interpreter', 'none');
        % legend({'ERP', 'PLV', 'var', 'Standard Error', 'saccade onset'});
        % grid on;
        % hold off;
        % xlim([-0.25, 0.25]);
        % 
        % subplot(2, 4, 8);
        % imagesc(time_vector(range_start:range_end), 1:size(phase_cluster2, 2), ...
        %     phase_cluster2');
        % colorbar;
        % title('Phase');
        % ylabel('Trial');
        % xlabel('Time (s)');
        % xlim([-0.25, 0.25]);
        % 
        % sgtitle(['Channel ', ch_label, ' oscillation Analysis']);


        %% characterize erp
        norm_concatenated_prep_signal = (prep_signal_all_session - baseline_mean) ./ baseline_std;
        time_vector = (0:size(norm_concatenated_prep_signal, 1) - 1) / sfreq - trial_begin_study;
        task_indices = find(time_vector >= erp_characteristic_range(1) & time_vector <= erp_characteristic_range(2));
        task_data_for_characteristic_cluster1 = norm_concatenated_prep_signal(task_indices, cluster1_trials);
        task_data_for_characteristic_cluster2 = norm_concatenated_prep_signal(task_indices, cluster2_trials);
        % flag bad erp
        cluster1_erp_good = 1;
        cluster2_erp_good = 1;

        % fig = figure('visible', 'off');
        % set(fig, 'position', [100, 100, 1500, 700]);
        % subplot(2, 3, 1)
        t = time_vector(find(time_vector >= erp_characteristic_range(1), 1):find(time_vector <= erp_characteristic_range(2), 1, 'last'));
        mean_erp_cluster1 = mean(task_data_for_characteristic_cluster1, 2);
        % plot(t, mean_erp_cluster1, 'b', 'LineWidth', 2);
        first_derivative = diff(mean(task_data_for_characteristic_cluster1, 2));
        first_derivative = smoothdata(first_derivative, 'movmean', smoothing_window);
        % hold on;
        % plot(t(1:end-1), 10 * first_derivative);
        [change_points, ~] = findchangepts(first_derivative, 'MaxNumChanges', 2);
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

            % for i = 1:length(change_times)
            %     xline(change_times(i), 'r--', 'LineWidth', 1.5, 'Label', sprintf('Change at %.1f %s', change_times(i) * 1e3), ...
            %         'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right');
            % end
            % plot(t(window_indices(max_index)), mean_erp_cluster1(window_indices(max_index)), 'go', 'MarkerFaceColor', 'g');
            % text(t(window_indices(max_index)), mean_erp_cluster1(window_indices(max_index)), ...
            %     sprintf('  ERP Amp: %.2f %s', mean_erp_cluster1(window_indices(max_index)), erp_polarity(1)), ...
            %     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
            % if ~isnan(cluster1_peak_to_trough)
            %     trough_time = t(trough_window_indices(trough_index));
            %     plot(trough_time, trough_amplitude, 'ro', 'MarkerFaceColor', 'r');
            %     text(trough_time, trough_amplitude, sprintf('  Trough: %.2f \n  Time: %.2f ms', trough_amplitude, trough_time * 1e3), ...
            %          'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
            % end
            % legend({'Mean ERP', 'First Derivative', 'Change Points'}, 'Location', 'Best');
            % xlabel('Time (s)');
            % ylabel('Amplitude (\muV)');
            % title('Cluster 1: Mean ERP charateristics');
            % grid on;
            % hold off;

            % define pre, during, and post
            % during would be between latency and peak_time, pre would be
            % between erp_pre_during_post(2) and latency, post would be
            % between peak_time and erp_pre_during_post(5)

            % get no erp phase
            task_indices_get_phase = find(time_vector >= -trial_begin_study & time_vector <= cluster1_latency);
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
        
            % subplot(2, 3, 2);
            % range_start = find(time_vector >= -trial_begin_study, 1);
            % range_end = find(time_vector <= (cluster1_latency), 1, 'last');
            % plot(time_vector(range_start:range_end), mean_phase, 'b', 'LineWidth', 2);
            % hold on;
            % plot(time_vector(range_start:range_end), plv, 'r', 'LineWidth', 2);
            % plot(time_vector(range_start:range_end), mean(no_erp_phase_cluster1, 2), 'green', 'LineWidth', 2);
            % fill([time_vector(range_start:range_end), fliplr(time_vector(range_start:range_end))], ...
            %  [mean_phase' + std_error_phase', ...
            %  fliplr((mean_phase)' - (std_error_phase)')], ...
            %  'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            % xline(0, 'r--', 'LineWidth', 1.5);
            % 
            % xlabel('Time (s)');
            % ylabel('Phase (rad)');
            % title(['C1 ' sprintf('%.1f', low_cutoff_oscillation) '-' sprintf('%.1f', high_cutoff_oscillation) 'Hz phase'], 'Interpreter', 'none');
            % legend({'ERP', 'PLV', 'mean', 'Standard Error', 'saccade onset'});
            % grid on;
            % hold off;
            % xlim([-0.2, 0]);
            % 
            % subplot(2, 3, 5);
            % imagesc(time_vector(range_start:range_end), 1:size(no_erp_phase_cluster1, 2), ...
            %     no_erp_phase_cluster1');
            % colorbar;
            % title('Phase');
            % ylabel('Trial');
            % xlabel('Time (s)');
            % xlim([-0.2, 0]);

            time_4_12_hz = (0:size(cluster1_no_erp_phase, 1) - 1) / sfreq - trial_begin_study;
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency - 1 / cluster1_peak_freq)));
            if any(indices_for_pre_during_post)
                cluster1_no_erp_phase_pre_pre_pre = mean(cluster1_no_erp_phase(indices_for_pre_during_post),'omitnan');
            else
                cluster1_no_erp_phase_pre_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency - 1 / cluster1_peak_freq / 2)));
            if any(indices_for_pre_during_post)
                cluster1_no_erp_phase_pre_pre = mean(cluster1_no_erp_phase(indices_for_pre_during_post),'omitnan');
            else
                cluster1_no_erp_phase_pre_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency)));
            if any(indices_for_pre_during_post)
                cluster1_no_erp_phase_pre = mean(cluster1_no_erp_phase(indices_for_pre_during_post),'omitnan');
            else
                cluster1_no_erp_phase_pre = nan;
            end

            time_4_12_hz = (0:size(cluster1_4_12_hz_pow, 1) - 1) / sfreq + erp_pre_during_post(1);
            % during is
            % at half between latency and amp_time, post is at end time,
            % prepre is latency - half cycle, preprepre is - cycle
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

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency)));
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

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_amp_time)));
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
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (erp_pre_during_post(3))));
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
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (erp_pre_during_post(4)) ));

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
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (erp_pre_during_post(3) + 1 / cluster1_peak_freq / 2)));

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
        % subplot(2, 3, 4)
        % plot(t, mean_erp_cluster2, 'b', 'LineWidth', 2);
        first_derivative = diff(mean_erp_cluster2);
        first_derivative = smoothdata(first_derivative, 'movmean', smoothing_window);
        % hold on;
        % plot(t(1:end-1), 10 * first_derivative);
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

            % for i = 1:length(change_times)
            %     xline(change_times(i), 'r--', 'LineWidth', 1.5, 'Label', sprintf('Change at %.1f %s', change_times(i) * 1e3), ...
            %         'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right');
            % end
            % plot(t(window_indices(max_index)), mean_erp_cluster2(window_indices(max_index)), 'go', 'MarkerFaceColor', 'g');
            % text(t(window_indices(max_index)), mean_erp_cluster2(window_indices(max_index)), ...
            %     sprintf('  ERP Amp: %.2f %s', mean_erp_cluster2(window_indices(max_index)), erp_polarity(1)), ...
            %     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
            % if ~isnan(cluster2_peak_to_trough)
            %     trough_time = t(trough_window_indices(trough_index));
            %     plot(trough_time, trough_amplitude, 'ro', 'MarkerFaceColor', 'r');
            %     text(trough_time, trough_amplitude, sprintf('  Trough: %.2f \n  Time: %.2f ms', trough_amplitude, trough_time * 1e3), ...
            %          'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
            % end
            % legend({'Mean ERP', 'First Derivative', 'Change Points'}, 'Location', 'Best');
            % xlabel('Time (s)');
            % ylabel('Amplitude (\muV)');
            % title('Cluster 2: Mean ERP charateristics');
            % grid on;
            % hold off;

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
        
            % subplot(2, 3, 3);
            % range_start = find(time_vector >= -trial_begin_study, 1);
            % range_end = find(time_vector <= (cluster2_latency), 1, 'last');
            % plot(time_vector(range_start:range_end), mean_phase, 'b', 'LineWidth', 2);
            % hold on;
            % plot(time_vector(range_start:range_end), plv, 'r', 'LineWidth', 2);
            % plot(time_vector(range_start:range_end), mean(no_erp_phase_cluster2, 2), 'green', 'LineWidth', 2);
            % fill([time_vector(range_start:range_end), fliplr(time_vector(range_start:range_end))], ...
            %  [mean_phase' + std_error_phase', ...
            %  fliplr((mean_phase)' - (std_error_phase)')], ...
            %  'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            % xline(0, 'r--', 'LineWidth', 1.5);
            % 
            % xlabel('Time (s)');
            % ylabel('Phase (rad)');
            % title(['C2 ' sprintf('%.1f', low_cutoff_oscillation) '-' sprintf('%.1f', high_cutoff_oscillation) 'Hz phase'], 'Interpreter', 'none');
            % legend({'ERP', 'PLV', 'mean', 'Standard Error', 'saccade onset'});
            % grid on;
            % hold off;
            % xlim([-0.2, 0]);

        
            % subplot(2, 3, 6);
            % imagesc(time_vector(range_start:range_end), 1:size(no_erp_phase_cluster2, 2), ...
            %     no_erp_phase_cluster2');
            % colorbar;
            % title('Phase');
            % ylabel('Trial');
            % xlabel('Time (s)');
            % xlim([-0.2, 0]);

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

            time_4_12_hz = (0:size(cluster2_4_12_hz_pow, 1) - 1) / sfreq - trial_begin_study;
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
            
            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_latency) ));
         
            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_pre = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_pre = angle(mean(exp(1i * cluster2_4_12_hz_phase(indices_for_pre_during_post)), 'omitnan'));
                cluster2_4_12_hz_plv_pre = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_pre = nan;
                cluster2_4_12_hz_phase_pre = nan;
                cluster2_4_12_hz_plv_pre = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_amp_time)));

            if any(indices_for_pre_during_post)
                cluster2_4_12_hz_pow_during = mean(cluster2_4_12_hz_pow(indices_for_pre_during_post), 'omitnan');
                cluster2_4_12_hz_phase_during = angle(mean(exp(1i * cluster2_4_12_hz_phase(indices_for_pre_during_post)), 'omitnan'));
                cluster2_4_12_hz_plv_during = mean(cluster2_4_12_hz_plv(indices_for_pre_during_post), 'omitnan');
            else
                cluster2_4_12_hz_pow_during = nan;
                cluster2_4_12_hz_phase_during = nan;
                cluster2_4_12_hz_plv_during = nan;
            end

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster2_amp_time + 1 / cluster2_peak_freq / 2)));

          
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

            time_4_12_hz = (0:size(cluster2_4_12_hz_pow, 1) - 1) / sfreq - trial_begin_study;
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

            [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - ( erp_pre_during_post(3) - 1 / cluster2_peak_freq / 2)));

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
  

        % so first change pt represents the latency, second change pt
        % represents the erp_location, the max(abs(around second pt)) would
        % be the erp amp, and depending on the polarity, we say it's pos
        % erp or neg erp.

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
            % alternative is to compare corr to decide when cluster it belongs to.
            if ismember(i_trial, cluster1_trials) && cluster1_erp_good
                % if cluster1, use the latency and amp_time to guide 
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

        %         task_indices_get_phase = find(time_vector >= erp_pre_during_post(1) & time_vector <= erp_pre_during_post(end));
        % cluster1_data_oscillation = prep_signal_all_session(task_indices_get_phase, cluster1_trials);

                % phase1 = phase_cluster1(:, (cluster1_trials == i_trial));

                % get no erp phase
                task_indices_get_phase_no_erp = find(time_vector >= - trial_begin_study & time_vector <= cluster1_latency);
                trial_data_oscillation = prep_signal_all_session(task_indices_get_phase_no_erp, i_trial);
                padding_samples = padding_duration * sfreq;

                padded_data = [flipud(trial_data_oscillation(1:padding_samples)); trial_data_oscillation; flipud(trial_data_oscillation(end-padding_samples+1:end))];
                
                filtered_data_no_erp = filtfilt(b, a, padded_data);
                analytic_signal = hilbert(filtered_data_no_erp);    
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

                time_4_12_hz = (0:size(amplitude, 1) - 1) / sfreq - trial_begin_study;
                [~, indices_for_pre_during_post] = min(abs(time_4_12_hz - (cluster1_latency - 1 / cluster1_peak_freq)));
                % indices_for_pre_during_post = (time_4_12_hz >= (trial_latency + 3 * erp_pre_during_post(2))) & ...
                % (time_4_12_hz < (trial_amp_time + 3 * erp_pre_during_post(2)));
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
                % indices_for_pre_during_post = (time_4_12_hz >= (trial_latency + 2 * erp_pre_during_post(2))) & ...
                % (time_4_12_hz < (trial_amp_time + 2 * erp_pre_during_post(2)));
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
                % indices_for_pre_during_post = (time_4_12_hz >= (trial_latency + 1 * erp_pre_during_post(2))) & ...
                % (time_4_12_hz < (trial_amp_time + 1 * erp_pre_during_post(2)));
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

                % indices_for_pre_during_post = (time_4_12_hz >= trial_latency) & (time_4_12_hz < trial_amp_time);
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
                % indices_for_pre_during_post = (time_4_12_hz >= (trial_latency - erp_pre_during_post(2))) & ...
                %     (time_4_12_hz < (trial_amp_time - erp_pre_during_post(2)));
                if any(indices_for_pre_during_post)
                    trial_4_12_hz_norm_pow_post = mean(norm_amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_pow_post = mean(amplitude(indices_for_pre_during_post), 'omitnan');
                    trial_4_12_hz_phase_post = angle(mean(exp(1i * phase(indices_for_pre_during_post)), 'omitnan'));
                else
                    trial_4_12_hz_norm_pow_post = nan;
                    trial_4_12_hz_pow_post = nan;
                    trial_4_12_hz_phase_post = nan;
                end

                % if viz_single_trial
                %     fig = figure('visible', 'off');
                %     set(fig, 'position', [100, 100, 1500, 700]);
                %     subplot(3, 1, 1);
                %     plot(t, current_trial_task, 'b', 'LineWidth', 1.5);
                %     hold on;
                %     plot(time_vector(task_indices_get_phase), filtered_data(padding_samples+1:end-padding_samples), 'b', 'LineWidth', 1.5);
                %     % Annotate the latency with xline
                %     xline(trial_latency, 'r--', 'LineWidth', 1.5, 'Label', sprintf('Latency: %.2f ms', trial_latency * 1e3), ...
                %           'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right');
                %     % Annotate the peak
                %     plot(trial_amp_time, erp_peak, 'go', 'MarkerFaceColor', 'g');
                %     text(trial_amp_time, erp_peak, sprintf('  Peak: %.2f \n  Time: %.2f ms', erp_peak, trial_amp_time * 1e3), ...
                %          'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
                %     if ~isnan(peak_to_trough)
                %         trough_time = t(trough_window_indices(trough_index));
                %         plot(trough_time, trough_amplitude, 'ro', 'MarkerFaceColor', 'r');
                %         text(trough_time, trough_amplitude, sprintf('  Trough: %.2f \n  Time: %.2f ms', trough_amplitude, trough_time * 1e3), ...
                %              'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
                %     end
                %     xlabel('Time (s)');
                %     ylabel('z-score (a.u.)');
                %     title(sprintf('C1 Trial %d - ERP Characterization', i_trial));
                %     legend({'Current Trial', 'Latency', 'ERP Peak'}, 'Location', 'Best');
                %     grid on;
                %     hold off;
                % 
                %     subplot(3, 1, 2);
                %     plot(time_vector(task_indices_get_phase), norm_amplitude, 'b', 'LineWidth', 1.5);
                %     hold on;
                %     plot(time_vector(task_indices_get_phase), phase, 'r', 'LineWidth', 1.5);
                %     xline(trial_latency - 1 / cluster1_peak_freq / 2, 'k--', 'LineWidth', 1.5, 'Label', 'Prepre', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right');
                %     xline(trial_latency, 'r--', 'LineWidth', 1.5, 'Label', sprintf('Latency: %.2f ms', trial_latency * 1e3), ...
                %           'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right');
                %     xline(trial_amp_time, 'g--', 'LineWidth', 1.5, 'Label', sprintf('Peak: %.2f ms', trial_amp_time * 1e3), ...
                %           'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right');
                %     xline(trial_latency - 1 / cluster1_peak_freq, 'k--', 'LineWidth', 1.5, 'Label', 'Prerere', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
                %     xlabel('Time (s)');
                %     ylabel('z-score (a.u.)');
                %     title([cluster1_polarity sprintf(' Trial %d - Oscillation Amplitude', i_trial)]);
                %     legend({'Amplitude', 'Phase', 'Pre', 'Latency', 'Peak', 'Post'}, 'Location', 'Best');
                %     grid on;
                %     hold off;
                %     xlim([-0.1, 0.15])
                % 
                %     subplot(3, 1, 3);
                %     plot(time_vector(task_indices_get_phase_no_erp), trial_data_oscillation, 'b', 'LineWidth', 1.5);
                %     hold on;
                %     plot(time_vector(task_indices_get_phase_no_erp), filtered_data_no_erp(padding_samples+1:end-padding_samples) * 2, 'g', 'LineWidth', 1.5);
                %     plot(time_vector(task_indices_get_phase_no_erp), no_erp_phase * 10, 'r', 'LineWidth', 1.5);
                %     xlabel('Time (s)');
                %     ylabel('Amplitude (a.u.)');
                %     title([sprintf(' Trial-%d ', i_trial), sprintf('%.1f', low_cutoff_oscillation), ' to ', ...
                %         sprintf('%.1f', high_cutoff_oscillation)]);
                %     legend({'Raw Signal', 'Filtered Signal', 'Phase'}, 'Location', 'Best');
                %     grid on;
                %     hold off;
                %     xlim([-0.5, 0])
                % end

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
                           trial_latency, erp_amplitude, trial_amp_time, erp_polarity,...
                           trial_4_12_hz_pow_pre_pre_pre, trial_4_12_hz_pow_pre_pre, trial_4_12_hz_pow_pre, trial_4_12_hz_pow_during, trial_4_12_hz_pow_post,...
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

                % get no erp phase
                task_indices_get_phase_no_erp = find(time_vector >= - trial_begin_study & time_vector <= (cluster2_latency));
                trial_data_oscillation = prep_signal_all_session(task_indices_get_phase_no_erp, i_trial);
                padding_samples = padding_duration * sfreq;

                padded_data = [flipud(trial_data_oscillation(1:padding_samples)); trial_data_oscillation; flipud(trial_data_oscillation(end-padding_samples+1:end))];
                
                filtered_data_no_erp = filtfilt(b, a, padded_data);
                analytic_signal = hilbert(filtered_data_no_erp);    
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

                time_4_12_hz = (0:size(amplitude, 1) - 1) / sfreq - trial_begin_study;
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
                % if viz_single_trial
                %     fig = figure('visible', 'off');
                %     set(fig, 'position', [100, 100, 1500, 700]);
                %     subplot(3, 1, 1);
                %     plot(t, current_trial_task, 'b', 'LineWidth', 1.5);
                %     hold on;
                %     % Annotate the latency with xline
                %     xline(trial_latency, 'r--', 'LineWidth', 1.5, 'Label', sprintf('Latency: %.2f ms', trial_latency * 1e3), ...
                %           'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right');
                %     % Annotate the peak
                %     plot(trial_amp_time, erp_peak, 'go', 'MarkerFaceColor', 'g');
                %     text(trial_amp_time, erp_peak, sprintf('  Peak: %.2f \n  Time: %.2f ms', erp_peak, trial_amp_time * 1e3), ...
                %          'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
                %     if ~isnan(peak_to_trough)
                %         trough_time = t(trough_window_indices(trough_index));
                %         plot(trough_time, trough_amplitude, 'ro', 'MarkerFaceColor', 'r');
                %         text(trough_time, trough_amplitude, sprintf('  Trough: %.2f \n  Time: %.2f ms', trough_amplitude, trough_time * 1e3), ...
                %              'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
                %     end
                %     xlabel('Time (s)');
                %     ylabel('z-score (a.u.)');
                %     title(sprintf('C2 Trial %d - ERP Characterization', i_trial));
                %     legend({'Current Trial', 'Latency', 'ERP Peak'}, 'Location', 'Best');
                %     grid on;
                %     hold off;
                % 
                %     subplot(3, 1, 2);
                %     plot(time_vector(task_indices_get_phase), norm_amplitude, 'b', 'LineWidth', 1.5);
                %     hold on;
                %     plot(time_vector(task_indices_get_phase), phase, 'r', 'LineWidth', 1.5);
                %     xline(trial_latency - 1 / cluster1_peak_freq / 2, 'k--', 'LineWidth', 1.5, 'Label', 'Prepre', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right');
                %     xline(trial_latency, 'r--', 'LineWidth', 1.5, 'Label', sprintf('Latency: %.2f ms', trial_latency * 1e3), ...
                %           'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right');
                %     xline(trial_amp_time, 'g--', 'LineWidth', 1.5, 'Label', sprintf('Peak: %.2f ms', trial_amp_time * 1e3), ...
                %           'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right');
                %     xline(trial_latency - 1 / cluster1_peak_freq, 'k--', 'LineWidth', 1.5, 'Label', 'Prerere', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
                %     xlabel('Time (s)');
                %     ylabel('z-score (a.u.)');
                %     title([cluster1_polarity sprintf(' Trial %d - Oscillation Amplitude', i_trial)]);
                %     legend({'Amplitude', 'Phase', 'Pre', 'Latency', 'Peak', 'Post'}, 'Location', 'Best');
                %     grid on;
                %     hold off;
                %     xlim([-0.1, 0.15]);
                % 
                %     subplot(3, 1, 3);
                %     plot(time_vector(task_indices_get_phase_no_erp), trial_data_oscillation, 'b', 'LineWidth', 1.5);
                %     hold on;
                %     plot(time_vector(task_indices_get_phase_no_erp), filtered_data_no_erp(padding_samples+1:end-padding_samples) * 2, 'g', 'LineWidth', 1.5);
                %     plot(time_vector(task_indices_get_phase_no_erp), no_erp_phase * 10, 'r', 'LineWidth', 1.5);
                %     xlabel('Time (s)');
                %     ylabel('Amplitude (a.u.)');
                %     title([sprintf(' Trial-%d ', i_trial), sprintf('%.1f', low_cutoff_oscillation), ' to ', ...
                %         sprintf('%.1f', high_cutoff_oscillation)]);
                %     legend({'Raw Signal', 'Filtered Signal', 'Phase'}, 'Location', 'Best');
                %     grid on;
                %     hold off;
                %     xlim([-0.5, 0]);
                % end
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
                           saccade_table_all_session.remember_condition(i_trial), nan, ...
                           nan, nan, nan, nan, nan, nan, ...
                           nan, nan, nan, nan, nan,...
                           nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan};
                df_evoked_potential = [df_evoked_potential; new_row];
            end
        end
    end
    filename = fullfile('result', ['responsive_erp_' subj '_' task '.mat']);
    save(filename, 'responsive_erp_table');

    filename = fullfile('result', ['df_evoked_potential_' subj '_' task '.mat']);
    save(filename, 'df_evoked_potential');

end

