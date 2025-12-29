%% description
% identify channels with epileptic discharge
% identifying events where the envelope of the high-pass filtered (>250 Hz) time-course exceeds a +5 SD threshold or 
% where both the amplitude and gradient surpass a +5 SD threshold
%% parameters
clear;
trial_begin_study = 3; % 3s
sfreq = 2000; % 2000Hz
viz_epilepsy_detection = 0;

addpath(['CircStat2012a']);
data_dir = 'data/';
% subjs = {'BJH021', 'BJH021', 'BJH028', 'BJH029', 'BJH030', 'BJH025', 'BJH026', 'BJH026', 'BJH027', 'BJH027', 'BJH027', 'BJH024', 'BJH024'};
% tasks = {'BLAES_study', 'BLAES_study_2', ...
%     'BLAES_study', 'BLAES_study', 'BLAES_study', 'BLAES_study_twosource', 'BLAES_study', 'BLAES_study_2', ...
%     'BLAES_study', 'BLAES_study_2', 'BLAES_study_3', 'BLAES_study', 'BLAES_study_2'};
% n_session_per_task = [2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4];
subjs = {'BJH025'};
tasks = {'BLAES_study_twosource'};
n_session_per_task = [4];

sum_contact_table_filename = ['data/norm_contact_data_table.csv'];
sum_contact_table = readtable(sum_contact_table_filename);


csv_dir = 'result';



all_summary_tables = table();
for i_subj = 1:length(subjs)
    trial_results = struct('subj', {}, 'ch_id', {}, 'ch_name', {}, ...
    'task', {}, 'session', {}, 'encoding_trial_id', {}, 'saccade_id_saccade_table', {},...
    'saccade_trial_id', {}, 'perc_epileptic', {}, 'epileptic_events', {});
    subj = subjs{i_subj};
    task = tasks{i_subj};
    n_session = n_session_per_task(i_subj);

    % get ch_label
    contact_this_subj = sum_contact_table(strcmp(sum_contact_table.SubjectID, subj), :);
    % similarly, we tested this examplary script with a subset of data for
    % efficiency
    for i_ch = 1:6
    % for i_ch = 1:height(contact_this_subj) 
        ch_label = contact_this_subj.labels_majority(i_ch);
        ch_name_bci2000_format_this_ch = contact_this_subj.ch_name_bci2000_format(i_ch);
        prep_signal_all_session = [];
        saccade_table_all_session = [];
        for session =1:n_session
        % get data
            prep_signal_filename = fullfile(data_dir, subj, [subj, '_', task, '_session', ...
                num2str(session), '_channel', num2str(i_ch), '_saccade_onset_not_image_onset_prep_signal.mat']);
            prep_signal = load(prep_signal_filename);
            prep_signal = prep_signal.prep_signal_saccade_onset_not_image_onset;
            filename = fullfile(data_dir, subj,  ...
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
        time_vector = (0:size(prep_signal_all_session, 1) - 1) / sfreq - trial_begin_study;
        n_trials = size(prep_signal_all_session, 2);
        for i_trial = 1:n_trials
            signal = prep_signal_all_session(:, i_trial);
            
            % High-pass filter the signal (>250 Hz) using a 4th order Butterworth filter
            [b, a] = butter(4, 250/(sfreq/2), 'high');
            hpSignal = filtfilt(b, a, signal);
            
            % Criterion 1: Compute envelope using Hilbert transform
            analytic_signal = hilbert(hpSignal);
            envelope_signal = abs(analytic_signal);
            threshold_env = mean(envelope_signal) + 5*std(envelope_signal);
            idx_env = find(abs(envelope_signal) > threshold_env);
            
            % Criterion 2: Amplitude and gradient thresholds
            threshold_amp = mean(signal) + 5*std(signal);
            idx_amp = find(abs(signal) > threshold_amp);
            gradient_signal = [0; diff(hpSignal)];
            threshold_grad = mean(gradient_signal) + 5*std(gradient_signal);
            idx_grad = find(gradient_signal > threshold_grad);
            idx_combined = intersect(idx_amp, idx_grad);
            
            idx_epileptic = union(idx_env, idx_combined);
            
            %  continuous epileptic events using a binary mask and connectivity analysis
            binary_signal = false(size(signal));
            binary_signal(idx_epileptic) = true;
            % requires Image Processing Toolbox)
            CC = bwconncomp(binary_signal);
            epileptic_events = zeros(CC.NumObjects, 2); % each row: [start_time, end_time]
            for i_event = 1:CC.NumObjects
                event_idx = CC.PixelIdxList{i_event};
                epileptic_events(i_event, 1) = time_vector(event_idx(1));
                epileptic_events(i_event, 2) = time_vector(event_idx(end));
            end
            
            % percentage of time classified as epileptic
            % discharge (already percentage)
            perc_epileptic = (length(idx_epileptic) / length(signal)) * 100;
            
            %% Save 
            trial_entry.subj = subj;
            trial_entry.ch_id = i_ch;
            trial_entry.ch_name = ch_name_bci2000_format_this_ch{1};
            trial_entry.task = task;
            trial_entry.session = saccade_table_all_session.session(i_trial);
            trial_entry.encoding_trial_id = saccade_table_all_session.trial(i_trial);
            trial_entry.saccade_id_saccade_table = saccade_table_all_session.id(i_trial);
            trial_entry.saccade_trial_id = i_trial;
            trial_entry.perc_epileptic = perc_epileptic;
            trial_entry.epileptic_events = epileptic_events;
            trial_results(end+1) = trial_entry;

            % diagnostic plot
            
            % if viz_epilepsy_detection
            %     fig = figure('Visible','off');
            %     subplot(3,1,1);
            %     plot(time_vector, signal, 'k', 'LineWidth', 1); hold on;
            %     yline(threshold_amp, '--g', 'LineWidth', 1);
            %     yline(threshold_amp - 10*std(signal), '--g', 'LineWidth', 1);
            %     ylims = ylim;
            %     for j = 1:size(epileptic_events,1)
            %         rect_start = epileptic_events(j, 1);
            %         rect_width = epileptic_events(j, 2) - epileptic_events(j, 1);
            %         rectangle('Position', [rect_start, ylims(1), rect_width, ylims(2)-ylims(1)], ...
            %             'FaceColor', [0.5, 0, 0.5, 0.3], 'EdgeColor', 'none');
            %     end
            %     title(sprintf('Subj: %s | Ch: %s | Trial: %d\nOriginal Signal with Amplitude Thresholds', ...
            %         subj, ch_name_bci2000_format_this_ch{1}, i_trial));
            %     xlabel('Time (s)'); ylabel('Amplitude');
            %     hold off;
            % 
            % 
            %     subplot(3,1,2);
            %     plot(time_vector, gradient_signal, 'b', 'LineWidth', 1);
            % 
            %     yline(threshold_grad, '--g', 'LineWidth', 1);
            %     yline(threshold_grad - 10*std(gradient_signal), '--g', 'LineWidth', 1);
            % 
            %     ylims = ylim;
            %     for j = 1:size(epileptic_events,1)
            %         rect_start = epileptic_events(j, 1);
            %         rect_width = epileptic_events(j, 2) - epileptic_events(j, 1);
            %         rectangle('Position', [rect_start, ylims(1), rect_width, ylims(2)-ylims(1)], ...
            %             'FaceColor', [0.5, 0, 0.5, 0.3], 'EdgeColor', 'none');
            %     end
            %     title(sprintf('Subj: %s | Ch: %s | Trial: %d\n Gradient Signal with Gradient Thresholds', ...
            %         subj, ch_name_bci2000_format_this_ch{1}, i_trial));
            %     xlabel('Time (s)'); ylabel('Amplitude');
            %     hold off;
            % 
            %     subplot(3,1,3);
            %     plot(time_vector, hpSignal, 'b', 'LineWidth', 1);
            %     yline(threshold_env, '--r', 'LineWidth', 1);
            %     yline(threshold_env - 10*std(hpSignal), '--r', 'LineWidth', 1);
            %     ylims = ylim;
            %     for j = 1:size(epileptic_events,1)
            %         rect_start = epileptic_events(j, 1);
            %         rect_width = epileptic_events(j, 2) - epileptic_events(j, 1);
            %         rectangle('Position', [rect_start, ylims(1), rect_width, ylims(2)-ylims(1)], ...
            %             'FaceColor', [0.5, 0, 0.5, 0.3], 'EdgeColor', 'none');
            %     end
            %     title(sprintf('Subj: %s | Ch: %s | Trial: %d\n High-Passed Signal with Envelope Threshold', ...
            %         subj, ch_name_bci2000_format_this_ch{1}, i_trial));
            %     xlabel('Time (s)'); ylabel('Amplitude');
            %     hold off;

            % end
        end
    end
    trial_table = struct2table(trial_results);

    summary_table = groupsummary(trial_table, {'subj','task','ch_id','ch_name'}, 'mean', 'perc_epileptic');

    csv_filename = fullfile(csv_dir, sprintf('%s_epileptic_activity_results_all_trials.csv', subj));
    writetable(trial_table, csv_filename);
   
end


