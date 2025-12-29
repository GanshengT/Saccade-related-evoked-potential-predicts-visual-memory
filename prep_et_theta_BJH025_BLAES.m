%% Description
% This script preprocess SEEG, and visualize time-frequency representation at image onset, 
% we will also visualize trial coherence, then we will find correlation
% between theta and saccade
% output, preprocessed signals and theta amplitude, phase over time

clear;

addpath(['/Users/ganshengtan/Library/CloudStorage/Box-Box/Washu/' ...
    'projects/BLEAS/updated_mex_files']);
data_dir = '/Users/ganshengtan/Library/CloudStorage/Box-Box/Washu/projects/BLEAS/theta_saccade/data/';
task_dir = ['/Users/ganshengtan/Library/CloudStorage/Box-Box/' ...
    'Washu/projects/BLEAS/all_tasks/BLAES'];
%% change here for new subj
subj_id = 'BJH025';
tasks_session_dict = containers.Map();
tasks_session_dict('BLAES_study_twosource') = [1, 2, 3, 4];
tasks_session_dict('BLAES_test_twosource') = [1, 2, 3, 4];
% get seeg contact indices
mapping_imaging_datfile_ch = containers.Map();
mapping_imaging_datfile_ch('A_IFG_OFG_SG') = 'A';
mapping_imaging_datfile_ch('B_SFG_ACC_SG') = 'B';
mapping_imaging_datfile_ch('C_MFG_ALNS_PIRIFC') = 'C';
mapping_imaging_datfile_ch('D_SUBCE_DLNS') = 'D';
mapping_imaging_datfile_ch('E_MTG_SUPTPOLE') = 'E';
mapping_imaging_datfile_ch('F_ITG_INFTPOLE') = 'F';
mapping_imaging_datfile_ch('GL_STG_VLNS_AMY') = 'GL';
mapping_imaging_datfile_ch('GR_STG_VLNS_AMY') = 'GR';
mapping_imaging_datfile_ch('HL_MTG_AHC_BEFR') = 'HL';
mapping_imaging_datfile_ch('HR_MTG_AHC_BEFR') = 'HR';
mapping_imaging_datfile_ch('I_MTG_PHC') = 'I';
mapping_imaging_datfile_ch('J_ITG_FUG_EC') = 'J';
mapping_imaging_datfile_ch('K_ITG_FUG_PHG') = 'K';
mapping_imaging_datfile_ch('L_MTG_EC_LESION') = 'L';

bad_ch_map = containers.Map();
BJH025_tasks = struct();
BJH025_tasks.BLAES_study_twosource = {'D8','GR4','GR5','J7'};
bad_ch_map('BJH025') = BJH025_tasks;

eeg_ch_names = {'FP1', 'F3', 'C3', 'P3', 'O1', 'FP2', 'F4', 'C4', ...
    'P4', 'O2', 'F7', 'T7', 'P7', 'F8', 'T8', 'P8', 'F9', 'F10', ...
    'FPZ', 'FZ', 'CZ', 'PZ', 'OZ'};
EOG_ch = {'FP1', 'FP2'};

imaging_processing_temp_folder = [data_dir, '/', subj_id, '/imaging_process'];
photodiode_name = 'DC03';
photodiode_ch_id_mic = 19;
n_contact_per_micro_shank = 8;
% check with openvar('parameters.ChannelNames')
% seeg_contact_indices = [7, 176];
% reorder coordinate based on ch_name in seeg recording
dat_file2match_seeg_indices = [data_dir, subj_id, '/BLAES_study_twosource/' ...
    'ECOG001/ECOGS001R01.dat'];
[signal, states, parameters] = load_bcidat(dat_file2match_seeg_indices);
ch_name_bci2000 = parameters.ChannelNames.Value;
contact_info_filename = fullfile(data_dir, subj_id, 'imaging_process/contact_data_table.mat');
contact_info = load(contact_info_filename);
contact_info = contact_info.contact_data_table;

contact_info_seeg = contact_info(contact_info.ContactID ~= 99, :);
contact_info_seeg.ch_name_bci2000_format = cellfun(@(x, y) [x{1} num2str(y)], ...
                                  regexp(contact_info_seeg.ShankID, '^[^_]*', 'match'), ...
                                  num2cell(contact_info_seeg.ContactID), ...
                                  'UniformOutput', false);

seeg_contact_indices = [];
rows_to_delete = [];
for i = 1:height(contact_info_seeg)
    shank_in_contact_info = contact_info_seeg.ShankID{i};
    ch_name_prefix = mapping_imaging_datfile_ch(shank_in_contact_info);
    ch_name = [ch_name_prefix num2str(contact_info_seeg.ContactID(i))];
    idx = find(strcmp(ch_name_bci2000, ch_name));
    % If found, store the index, otherwise store NaN
    if ~isempty(idx)
        seeg_contact_indices(end + 1) = idx;
    else
        rows_to_delete(end + 1) = i;
        disp([ch_name 'not found in bci dat file'])
    end
end
contact_info_seeg(rows_to_delete, :) = [];


% identify photodiode for sync
if ismember(99, contact_info.ContactID)
    filename = [data_dir, subj_id, '/BLAES_study_twosource/' ...
    'ECOG001/ECOGS001R01_1.dat'];
    [signal_mic, states_mic, parameters_mic] = load_bcidat(filename);
    figure,set(gcf,'position',[0 0 1500 700]);
    plot(double(states.(photodiode_name) / 10),  'DisplayName', photodiode_name);
    title(['is ' photodiode_name ' photodiode?']);
    set(gca,'FontSize',12);
    hold on;
    plot(states.StimulusCode, 'DisplayName', 'StimulusCode');
    legend('show');
    hold off;
    potential_photodiode_mic_id = find(~strncmp(parameters_mic.ChannelNames.Value,'chan',4));
    figure,set(gcf,'position',[0 0 1000 500])
    for ii = 1:length(potential_photodiode_mic_id)
        tmpD = single(signal_mic(:,potential_photodiode_mic_id(ii)));
        subplot(3,3,ii),
        plot(tmpD)
        title(num2str(potential_photodiode_mic_id(ii)))
        set(gca,'FontSize',12)
    end
    sgtitle(['is ' num2str(photodiode_ch_id_mic) ' photodiode?']);
end

% small laplacian referencing threshold
laplacian_thres = 6; % 5mm bound, 
% DIXI microdeep electrode contact distance is 1.5mm, macro contact length is 2mm
num_cycle_wavelet = 7;
freqs = 2:70; % 1Hz oscillation wavelet is longer than the epoch, beyond 70Hz is maybe broadband activity
prep_signal_image_onset_length = 5; % 5s after image onset

saccade_onset_control_image_onset = 0.2; % get saccade not within 200ms around the image onset

temp_save_path = fullfile(data_dir, subj_id, 'prep_data');
report_path = fullfile(data_dir, subj_id, 'report');
if ~exist(temp_save_path, 'dir')
    mkdir(temp_save_path);
end

if ~exist(report_path, 'dir')
    mkdir(report_path);
end

temp_save_tf_path = fullfile(data_dir, subj_id, 'prep_data', 'time_frequency');
if ~exist(temp_save_tf_path, 'dir')
    mkdir(temp_save_tf_path);
end


unique_tasks = unique(tasks_session_dict.keys);
for i_task = 1:length(unique_tasks)
    current_task = unique_tasks{i_task};
    sessions_during_this_task = tasks_session_dict(current_task);
    for i_session = 1:length(sessions_during_this_task)
        current_session = sessions_during_this_task(i_session);
        filename = fullfile([data_dir, subj_id, '/', current_task, ...
            '/ECOG001/ECOGS001R0', num2str(current_session), '.dat']);
        [signal, states, parameters] = load_bcidat(filename);
        if ismember(99, contact_info.ContactID)
            filename = fullfile([data_dir, subj_id, '/', current_task, ...
                '/ECOG001/ECOGS001R0', num2str(current_session), '_1.dat']);
            [signal_mic, states_mic, parameters_mic] = load_bcidat(filename);
            PD_Mic = single(signal_mic(:,photodiode_ch_id_mic));
            PD_UpSam = upsample(double(states.(photodiode_name)),...
                parameters_mic.SamplingRate.NumericValue/parameters.SamplingRate.NumericValue);
            for i = 2:length(PD_UpSam)
                if PD_UpSam(i) == 0
                    PD_UpSam(i) = PD_UpSam(i-1);
                end
            end
            [r,lags] = xcorr(diff(PD_Mic),diff(PD_UpSam));
            % figure,plot(r)
            [~,lagInd] = max(r);
            lagSample = lags(lagInd); %micro lags behind of macro
            lagLatency = lagSample/parameters_mic.SamplingRate.NumericValue;
            disp(['Micro recorded ',num2str(lagLatency),'s earlier' ]);
            StimCode_UpSample = upsample(double(states.StimulusCode),...
                parameters_mic.SamplingRate.NumericValue/parameters.SamplingRate.NumericValue);
            for i = 2:length(StimCode_UpSample)
                if StimCode_UpSample(i) == 0
                    StimCode_UpSample(i) = StimCode_UpSample(i-1);
                end
            end
            if lagSample>0 %micro earlier
                StimCode_Mic = [zeros(lagSample,1);StimCode_UpSample]; 
            else %marco earlier - that should be the case
                StimCode_Mic = StimCode_UpSample(lagSample*-1+1:end);
            end
            [~,stimOnset] = findpeaks(StimCode_Mic,'MinPeakDistance',30000*1.25); %sample index for each event
        
            %plot PD and the digital triggers to see alignments
            fig = figure;
            set(fig,'position',[0 0 1200 400]);
            viz_time = 30;
            subplot(3, 1, 1)
            plot(-PD_Mic(1: parameters_mic.SamplingRate.NumericValue*viz_time) / max(PD_Mic));
            hold on;
            plot(PD_UpSam(1: parameters_mic.SamplingRate.NumericValue*viz_time) / max(PD_UpSam));
            legend({'PD micro','PD macro'})
            set(gca,'FontSize',12,'LineWidth',2);
            title('prior to alignment');
            subplot(3, 1, 2)
            hold on
            if lagSample>0
                plot(-PD_Mic(lagSample: lagSample + parameters_mic.SamplingRate.NumericValue*viz_time) / max(PD_Mic));
                plot(PD_UpSam(1: parameters_mic.SamplingRate.NumericValue*viz_time) / max(PD_UpSam));
            else
                plot(-PD_Mic(1: 1 + parameters_mic.SamplingRate.NumericValue*viz_time) / max(PD_Mic));
                plot(PD_UpSam(-lagSample: -lagSample+parameters_mic.SamplingRate.NumericValue*viz_time) / max(PD_UpSam));
            end
            legend({'PD micro','PD macro'})
            set(gca,'FontSize',12,'LineWidth',2);
            title('after alignment');
            subplot(3, 1, 3)
            hold on
            if lagSample>0
                plot(-StimCode_UpSample(lagSample:lagSample+parameters_mic.SamplingRate.NumericValue*viz_time)/max(StimCode_UpSample),'r-')
                plot(PD_Mic(lagSample:lagSample+parameters_mic.SamplingRate.NumericValue*viz_time)/max(PD_Mic),'b-')
            else
                plot(-StimCode_Mic(1:viz_time*parameters_mic.SamplingRate.NumericValue)/max(StimCode_Mic),'r-') %normalize to -1 to 0
                plot(PD_Mic(1:viz_time*parameters_mic.SamplingRate.NumericValue)/max(PD_Mic),'b-') %normalize to 0-1
            end
            legend({'StimOn','PD'})
            set(gca,'FontSize',12,'LineWidth',2);
            title('mic stim code and PD');
            box off
            filename = [report_path '/mic_macro_alignment_' ...
                current_task '_' num2str(current_session) '.svg'];
            saveas(fig, filename);
            % do delay correction, then downsample signal_mic, and make it
            % same length as signal
            corrected_signal_mic = circshift(double(signal_mic), -lagSample, 1); 
            downsampled_signal_mic = resample(corrected_signal_mic, ...
                parameters.SamplingRate.NumericValue, parameters_mic.SamplingRate.NumericValue);
        
            if size(downsampled_signal_mic, 1) > size(signal, 1)
                downsampled_signal_mic = downsampled_signal_mic(1:size(signal, 1), :);
            elseif size(downsampled_signal_mic, 1) < size(signal, 1)
                downsampled_signal_mic = [downsampled_signal_mic; zeros(size(signal, 1) - size(downsampled_signal_mic, 1), size(downsampled_signal_mic, 2))];
            end
        end

        filename = fullfile(temp_save_path, [subj_id '_' current_task '_session' num2str(current_session) ...
            '_prep_saccade_event.mat']);
        event_info = load(filename);
        event_info = event_info.data_struct_2_export;

        % find EOG from frontal electrode
        % one channel can get blink eog
        ch_EOG_one = double(signal(:,strcmp(parameters.ChannelNames.Value, EOG_ch{1})));
        ch_EOG_two = double(signal(:,strcmp(parameters.ChannelNames.Value, EOG_ch{2})));

        % all_ch_names = parameters.ChannelNames.Value;
        % eeg_indices = find(ismember(all_ch_names, eeg_ch_names));
        % EOG_ref = mean(double(signal(:, eeg_indices)), 2);
        EOG_signal = ch_EOG_one - ch_EOG_two;

        % get contacts that have recordings, and get seeg signals ordered as the contact table
        variableTypes = varfun(@class, contact_info_seeg, 'OutputFormat', 'cell');
        contact_info_seeg_matched_bci2000_signal =contact_info_seeg(1, :);
        contact_info_seeg_matched_bci2000_signal(1, :) = [];

        seeg_signal = zeros(size(signal, 1), length(seeg_contact_indices));
        for i_ch = 1:size(seeg_signal, 2)
            if isnan(seeg_contact_indices(i_ch))
                continue
            else
                seeg_signal(:, i_ch) = signal(:, seeg_contact_indices(i_ch));
                contact_info_seeg_matched_bci2000_signal = [...
                    contact_info_seeg_matched_bci2000_signal; contact_info_seeg(i_ch, :)];
                shank_in_contact_info = contact_info_seeg.ShankID{i_ch};
                ch_name_prefix = mapping_imaging_datfile_ch(shank_in_contact_info);
                contact_info_seeg_matched_bci2000_signal.ch_name_bci2000_format{i_ch} = [ch_name_prefix num2str(contact_info_seeg.ContactID(i_ch))];
            end
        end

        % check
        assert(sum(contact_info_seeg_matched_bci2000_signal.ContactID ~= 99) == size(seeg_signal, 2), "check seeg channel");
        % we will add lfp from micro contact
        signal_mic_row_id = 1;
        for i_micro_shank = find(contact_info.ContactID == 99)'
            % populate the contact_info_seeg_matched_bci2000_signal
            for i_micro_contact = 1:n_contact_per_micro_shank
                new_row = contact_info(i_micro_shank, :);
                % signal_mic_row_id = (i_micro_shank - 1) * n_contact_per_micro_shank + i_micro_contact;
                new_row.ch_name_bci2000_format = ['chan' num2str(signal_mic_row_id)];
                new_row.ShankID = [new_row.ShankID{1} '_mic'];
                contact_info_seeg_matched_bci2000_signal = [contact_info_seeg_matched_bci2000_signal; new_row];
                contact_info_seeg_matched_bci2000_signal.ContactID(end) = i_micro_contact;
                id_channel_name_in_mic = find(cellfun(@(x) strcmp(x, ['chan' num2str(signal_mic_row_id)]), ...
                    parameters_mic.ChannelNames.Value));
                seeg_signal = [seeg_signal, downsampled_signal_mic(:, id_channel_name_in_mic)];
                signal_mic_row_id = signal_mic_row_id + 1;
            end
        end

        % save if the file does not exist
        table_filename = [imaging_processing_temp_folder, ...
            '/contact_data_table_matched_bci2000_signal.mat'];
        if exist(table_filename, 'file') ~= 2
            save(table_filename, 'contact_info_seeg_matched_bci2000_signal', '-v7.3');
            contact_data_struct = table2struct(contact_info_seeg_matched_bci2000_signal);
            table_filename = [imaging_processing_temp_folder, ...
                '/contact_data_table_struct_matched_bci2000_signal.mat'];
            save(table_filename, 'contact_data_struct', '-v7.3');
        else
            disp('contact info is matched to signal');
        end

        % % reorder coordinate based on ch_name in seeg recording
        % ch_name_bci2000 = parameters.ChannelNames.Value(seeg_contact_indices(1): seeg_contact_indices(2));
        % contact_info_seeg = contact_info(contact_info.ContactID ~= 99, :);
        % contact_info_seeg.ch_name_bci2000_format = cellfun(@(x, y) [x{1} num2str(y)], ...
        %                                   regexp(contact_info_seeg.ShankID, '^[^_]*', 'match'), ...
        %                                   num2cell(contact_info_seeg.ContactID), ...
        %                                   'UniformOutput', false);
        % 
        % variableTypes = varfun(@class, contact_info_seeg, 'OutputFormat', 'cell');
        % 
        % contact_info_seeg_matched_bci2000_signal =contact_info_seeg(1, :);
        % contact_info_seeg_matched_bci2000_signal(1, :) = [];
        % for i = 1:length(ch_name_bci2000)
        %     current_name = ch_name_bci2000{i};
        %     match_idx = strcmp(contact_info_seeg.ch_name_bci2000_format, current_name);
        % 
        %     if sum(match_idx) == 1
        %         matched_rows = contact_info_seeg(match_idx, :);
        %         contact_info_seeg_matched_bci2000_signal = [contact_info_seeg_matched_bci2000_signal; matched_rows];
        %     elseif sum(match_idx) > 1
        %         fprintf('Multiple matches found for %s\n', current_name);
        %     else
        %         fprintf('No match found for %s\n', current_name);
        %     end
        % end
        % table_filename = [imaging_processing_temp_folder, ...
        %     '/contact_data_table_matched_bci2000_signal.mat'];
        % if exist(table_filename, 'file') ~= 2
        %     save(table_filename, 'contact_info_seeg_matched_bci2000_signal', '-v7.3');
        %     contact_data_struct = table2struct(contact_info_seeg_matched_bci2000_signal);
        %     table_filename = [imaging_processing_temp_folder, ...
        %         '/contact_data_table_struct_matched_bci2000_signal.mat'];
        %     save(table_filename, 'contact_data_struct', '-v7.3');
        % else
        %     disp('contact info is matched to signal');
        % end

        % bipolar rereferncing
        % coordinates = table2array(contact_info(contact_info.ContactID ~= 99, ...
        %     {'ShankID', 'ContactID', 'X_aligned_to_brightest_voxel', 'Y_aligned_to_brightest_voxel', 'Z_aligned_to_brightest_voxel'}));
        % % Calculate the Euclidean distances between each pair of electrodes
        % coordinates = table2array(contact_info_seeg_matched_bci2000_signal(:, ...
        %     {'X_aligned_to_brightest_voxel', 'Y_aligned_to_brightest_voxel', 'Z_aligned_to_brightest_voxel'}));
        % 
        % dist_matrix = zeros(size(seeg_signal, 2), size(seeg_signal, 2));
        % 
        % for i = 1:size(seeg_signal, 2)
        %     for j = 1:size(seeg_signal, 2)
        %         dist_matrix(i, j) = sqrt(sum((coordinates(i,:) - coordinates(j,:)).^2));
        %     end
        % end
        % % add neighbors for micro wires
        % for i = (height(contact_info_seeg) + 1):size(seeg_signal, 2)
        %     for j = (height(contact_info_seeg) + 1):size(seeg_signal, 2)
        %         if i ~= j
        %             dist_matrix(i, j) = dist_matrix(i, j) + 0.1;
        %         end
        %     end
        % end
        % 
        % neighbors = cell(size(seeg_signal, 2), 1);
        % for i = 1:size(seeg_signal, 2)
        %     neighbors{i} = find(dist_matrix(i,:) < laplacian_thres & dist_matrix(i,:) > 0);
        % end
        is_micro = contains(contact_info_seeg_matched_bci2000_signal.ShankID, '_mic');
        is_macro = ~is_micro;
        
        micro_indices = find(is_micro);
        macro_indices = find(is_macro);
        dc_offsets_micro = mean(seeg_signal(:, micro_indices), 1); % Mean across baseline for each microwire
        seeg_signal(:, micro_indices) = seeg_signal(:, micro_indices) - dc_offsets_micro;

        coordinates = table2array(contact_info_seeg_matched_bci2000_signal(:, ...
            {'X_aligned_to_brightest_voxel', 'Y_aligned_to_brightest_voxel', 'Z_aligned_to_brightest_voxel'}));

        dist_matrix = squareform(pdist(coordinates)); % Efficient distance computation
        dist_matrix(is_micro, is_macro) = Inf;
        dist_matrix(is_macro, is_micro) = Inf;

        for i = (height(contact_info_seeg) + 1):size(seeg_signal, 2)
            for j = (height(contact_info_seeg) + 1):size(seeg_signal, 2)
                if i ~= j
                    dist_matrix(i, j) = dist_matrix(i, j) + 0.1;
                end
            end
        end

        neighbors = cell(size(seeg_signal, 2), 1);
        for i = 1:size(seeg_signal, 2)
            neighbors{i} = find(dist_matrix(i,:) < laplacian_thres & dist_matrix(i,:) > 0);
        end

        bad_ch_this_subj = bad_ch_map(subj_id);
        bad_ch_this_subj = bad_ch_this_subj.(current_task);

        for i = 1:length(neighbors)
            original_neighbors = neighbors{i};
            new_neighbors = [];
            for j = 1:length(original_neighbors)
                neighbor_idx = original_neighbors(j);
                neighbor_bci2000_ch_name = contact_info_seeg_matched_bci2000_signal.ch_name_bci2000_format(neighbor_idx);
                if ismember(neighbor_bci2000_ch_name, bad_ch_this_subj)
                    fprintf('Removed bad channel %d from neighbors of channel %d\n', neighbor_idx, i);
                else
                    new_neighbors(end+1) = neighbor_idx;
                end
            end
            neighbors{i} = new_neighbors;
        end

        fprintf('Sanity Check: Channels with either 0 or >2 neighbors:\n');
        for i = 1:length(neighbors)
            num_neighbors = length(neighbors{i});
            if num_neighbors == 0 || num_neighbors > 2
                fprintf('Channel %d has %d neighbor(s): %s\n', i, num_neighbors, mat2str(neighbors{i}));
            end
        end

        laplacian_referenced_signal = zeros(size(seeg_signal));
        for i = 1:size(seeg_signal, 2)
            if isempty(neighbors{i})
                laplacian_referenced_signal(:,i) = seeg_signal(:,i); % No neighbors found, use original signal
            else
                laplacian_referenced_signal(:,i) = seeg_signal(:,i) - mean(seeg_signal(:,neighbors{i}), 2);
            end
        end
        % figure; plot(laplacian_referenced_signal(:, 57)); ylim([-200,200])
        clear seeg_signal

        % epoch the data, and viz time-frequency, it's good have epoch of
        % same size
        event_table = event_info.trial_event_struct;

        if contains(current_task, 'BLAES_test') 
            max_length = 0;  
            for i = 1:size(event_table, 1)
                onset_sample = event_table(i).fixed_image_onset_photodiode_sample + ...
                               event_info.trial_event_begin * event_info.sfreq;
                offset_sample = event_table(i).response_image_offset_photodiode_sample + ...
                                event_info.trial_event_end * event_info.sfreq;
                epoch_length = offset_sample - onset_sample + 1;
                if epoch_length > max_length
                    max_length = epoch_length;  
                end
            end
    
            epochs = cell(size(event_table, 1), 1);
            for i = 1:size(event_table, 1)
                onset_sample = event_table(i).fixed_image_onset_photodiode_sample + ...
                               event_info.trial_event_begin * event_info.sfreq;
                offset_sample = onset_sample + max_length - 1;  
            
                if offset_sample > size(laplacian_referenced_signal, 1)
                    offset_sample = size(laplacian_referenced_signal, 1); 
                end
            
                epochs{i} = laplacian_referenced_signal(onset_sample:offset_sample, :);
                
                if length(epochs{i}) < max_length
                    padding = repmat(epochs{i}(end, :), max_length - length(epochs{i}), 1);  
                    epochs{i} = [epochs{i}; padding]; 
                end
            end
            
        elseif contains(current_task, 'BLAES_study') 
            max_length = 0;  
            for i = 1:size(event_table, 1)
                onset_sample = event_table(i).image_onset_photodiode_sample + ...
                               event_info.trial_event_begin * event_info.sfreq;
                offset_sample = event_table(i).image_offset_photodiode_sample + ...
                                event_info.trial_event_end * event_info.sfreq;
                epoch_length = offset_sample - onset_sample + 1;
                if epoch_length > max_length
                    max_length = epoch_length;  
                end
            end
    
            epochs = cell(size(event_table, 1), 1);
            for i = 1:size(event_table, 1)
                onset_sample = event_table(i).image_onset_photodiode_sample + ...
                               event_info.trial_event_begin * event_info.sfreq;
                offset_sample = onset_sample + max_length - 1;  
            
                if offset_sample > length(laplacian_referenced_signal)
                    offset_sample = length(laplacian_referenced_signal); 
                end
            
                epochs{i} = laplacian_referenced_signal(onset_sample:offset_sample, :);
                
                if length(epochs{i}) < max_length
                    padding = repmat(epochs{i}(end, :), max_length - length(epochs{i}), 1);  
                    epochs{i} = [epochs{i}; padding]; 
                end
            end
        end
        
        num_trials = size(epochs, 1);
        if size(epochs, 1) > 0
            [uniform_length, num_channels] = size(epochs{1});
            all_data = zeros(num_trials, uniform_length, num_channels);
        
            for i = 1:num_trials
                all_data(i, :, :) = epochs{i};
            end
        
            for ch = 1:num_channels
                prep_signal_image_onset = squeeze(all_data(:, :, ch));
        
                filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_channel%d_image_onset_prep_signal.mat', ...
                            subj_id, current_task, current_session, ch));
                if ~exist(filename, 'file')
                        save(filename, 'prep_signal_image_onset', '-v7.3');
                else
                    fprintf('File %s already exists. Skip\n', filename);
                end
       
            end
        end
        clear all_data
      
        %% for full processing, uncomment the following
        % % epoch * ch * freq * time
        % % It is painful to have too many files, I will save one file per
        % % channel including all epochs
        % % Iterate over each channel
        % for ch = 1:size(epochs{1}, 2)
        %     ifft_wavelet_filtered = zeros(numel(epochs), length(freqs), size(epochs{1}, 1));
        %     filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_channel%d_image_onset.mat', ...
        %             subj_id, current_task, current_session, ch));
        %     if ~exist(filename, 'file')
        %         for i_epoch = 1:numel(epochs)
        %             reconstructed_signal = zeros(size(epochs{i_epoch}(:, ch)));
        %             % we should not have phase shift after correction
        %             for f_idx = 1:length(freqs)
        %                 f = freqs(f_idx);
        %                 [wavelet, wavelet_length] = morlet_wavelet_cohen(f, event_info.sfreq, num_cycle_wavelet);
        %                 wavelet = wavelet / sum(abs(wavelet));
        %                 padded_signal = [zeros(wavelet_length, 1); epochs{i_epoch}(:, ch)];
        %                 lenConv = length(padded_signal) + wavelet_length;
        %                 fft_signal = fft(padded_signal, lenConv);
        %                 fft_wavelet = fft(wavelet, lenConv);
        %                 conv_result = ifft(fft_signal .* fft_wavelet');
        % 
        %                 start_idx = floor(wavelet_length / 2) + 1;
        %                 end_idx = start_idx + length(epochs{i_epoch}(:, ch)) - 1;
        %                 ifft_wavelet_filtered(i_epoch, f_idx, :) = conv_result(start_idx:end_idx);
        % 
        %                 reconstructed_signal = reconstructed_signal + squeeze(real(ifft_wavelet_filtered(i_epoch, f_idx, :)));
        %             end
        % 
        %             energy_original = sum(epochs{i_epoch}(:, ch).^2);
        %             energy_reconstructed = sum(reconstructed_signal.^2);
        %             scale_factor = sqrt(energy_original / energy_reconstructed);
        %             ifft_wavelet_filtered(i_epoch, :, :) = ifft_wavelet_filtered(i_epoch, :, :) * scale_factor;
        %         end
        %         save(filename, 'ifft_wavelet_filtered', '-v7.3');
        %     else
        %         fprintf('File %s already exists. Skipping.\n', filename);
        % 
        %     end
        % end

        %%%%%%%%%%%%%%%%
        
        % % Loop through each epoch and frequency
        % for i_epoch = 1:numel(epochs)
        %     for ch = 1:size(epochs{1}, 2)
        %         % Initialize the filtered data matrix for current channel
        %         filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_epoch%d_channel%d_image_onset.mat', ...
        %             subj_id, current_task, current_session, i_epoch, ch));
        %         if ~exist(filename, 'file')
        %             ifft_wavelet_filtered = zeros(length(freqs), size(epochs{1}, 1));
        %             % Initialize a variable for the reconstructed signal
        %             reconstructed_signal = zeros(size(epochs{i_epoch}(:, ch)));
        %             for f_idx = 1:length(freqs)
        %                 f = freqs(f_idx);
        %                 % d = designfilt('bandpassfir', 'FilterOrder', 20, ...
        %                 %                'CutoffFrequency1', f-0.5, 'CutoffFrequency2', f+0.5, ...
        %                 %                'SampleRate', event_info.sfreq);
        %                 % filtered_signal = filtfilt(d, epochs{i_epoch}(:, ch));  
        %                 % analytic_signal = hilbert(filtered_signal);
        %                 % amplitude_envelope = abs(analytic_signal);  % Amplitude envelope from Hilbert transform
        %                 % narrow_filtered(i, :, f) = amplitude_envelope;
        % 
        %                 % % Wavelet filter
        %                 % [wavelet, wavelet_length] = morlet_wavelet_cohen(f, event_info.sfreq, num_cycle_wavelet);
        %                 % wavelet_center = ceil(wavelet_length / 2);
        %                 % half_wavelet = floor(wavelet_length / 2); 
        %                 % filtered_wavelet = conv(epochs{i}(:, 1), wavelet, 'same');
        %                 % wavelet_filtered(i, :, f) = circshift(abs(filtered_wavelet), wavelet_center);
        % 
        %                 % Frequency domain filtering
        %                 [wavelet, wavelet_length] = morlet_wavelet_cohen(f, event_info.sfreq, num_cycle_wavelet);
        %                 wavelet = wavelet / sum(abs(wavelet));
        %                 padded_signal = [zeros(wavelet_length, 1); epochs{i_epoch}(:, ch)]; 
        %                 lenConv = length(padded_signal) + wavelet_length;
        %                 fft_signal = fft(padded_signal, lenConv);
        %                 % lenConv = length(epochs{i_epoch}(:, ch)) + wavelet_length;
        %                 % fft_signal = fft(epochs{i_epoch}(:, ch), lenConv);
        %                 fft_wavelet = fft(wavelet, lenConv);
        %                 conv_result = ifft(fft_signal .* fft_wavelet');
        %                 % Extract the relevant portion of the convolution result
        %                 start_idx = floor(wavelet_length / 2) + 1;
        %                 % start_idx = 1;
        %                 end_idx = start_idx + length(epochs{i_epoch}(:, ch)) - 1;
        %                 % if end_idx > length(conv_result)
        %                 %     end_idx = length(conv_result);
        %                 % end
        %                 ifft_wavelet_filtered(f_idx, :) = conv_result(start_idx:end_idx);
        %                 % validate
        %                 % figure; plot(epochs{i_epoch}(:, ch)); hold on; plot(real(ifft_wavelet_filtered(f_idx, :)))
        % 
        %                 reconstructed_signal = reconstructed_signal + real(ifft_wavelet_filtered(f_idx, :))';
        % 
        %                 % to align the phase ifft_wavelet_filtered with the
        %                 % original signal, one can use plot(abs(filtered_wavelet).*cos(angle(filtered_wavelet)));
        %                 % or take the real part, these two are equivalent
        %             end
        %             energy_original = sum(epochs{i_epoch}(:, ch).^2);
        %             energy_reconstructed = sum(reconstructed_signal.^2);
        %             scale_factor = sqrt(energy_original / energy_reconstructed);
        %             ifft_wavelet_filtered = ifft_wavelet_filtered * scale_factor;
        % 
        %             % check signal
        %             % figure;
        %             % plot(epochs{i_epoch}(:, ch));
        %             % hold on 
        %             % plot(real(ifft_wavelet_filtered(60,:)));
        % 
        %             % Save filtered data for the current channel of the current epoch
        % 
        %             save(filename, 'ifft_wavelet_filtered', '-v7.3');
        %         else
        %             fprintf('File %s already exists. Skip\n', filename);
        %         end
        % 
        %     end
        % end

        % process saccade onset epoch
        if contains(current_task, 'BLAES_test') 
            saccade_onsets_image_viewing = [];
            saccade_onsets_image_viewing_not_image_onset = [];
            saccade_onsets_control = [];

            saccade_dist_image_viewing = [];
            saccade_dist_image_viewing_not_image_onset = [];
            saccade_dist_control = [];
            saccade_table = table([], [], [], [], [], [], [], [], [], [], [], ...
                [], [], [], [], [], [], [], [], 'VariableNames', ...
                 {'id', 'subj', 'task', 'session', 'trial', 'label', ...
                 'eccentricity', 'gaze_before_x_screen_coord', 'gaze_before_y_screen_coord', ...
                 'gaze_after_x_screen_coord', 'gaze_after_y_screen_coord','azimuth_before',...
                 'azimuth_after','elevation_before', 'elevation_after', ...
                 'duration', 'time_to_image_onset', 'behavior_condition', 'remember_condition'});

            % n_epoch_saccade_onset_image_viewing = 0; 
            % n_epoch_saccade_onset_control = 0; 
            saccade_id = 1;
            for i = 1:size(event_table, 1)
                fixationstats = event_info.fixationstats{i};
                for i_saccade = 1:size(fixationstats.saccadetimes, 2)

                    trial_onset_sample = event_table(i).fixed_image_onset_photodiode_sample + ...
                                   event_info.trial_event_begin * event_info.sfreq;
                    trial_offset_sample = event_table(i).response_image_offset_photodiode_sample + ...
                                    event_info.trial_event_end * event_info.sfreq;
                    if (fixationstats.valid_saccade_label(1, i_saccade) == 0) ||...
                            (fixationstats.valid_saccade_label(2, i_saccade) == 0) 
                        continue
                    end

                    saccade_onset_sample = trial_onset_sample + fixationstats.saccadetimes(1, i_saccade) + ...
                                   event_info.trial_event_begin * event_info.sfreq;
                    

                    % saccade during image viewing
                    if (fixationstats.saccadetimes(1, i_saccade) > ...
                            -event_info.trial_event_begin * event_info.sfreq) && ( ...
                            fixationstats.saccadetimes(1, i_saccade) < ...
                            (trial_offset_sample - trial_onset_sample - event_info.trial_event_end * event_info.sfreq))
                        saccade_onsets_image_viewing = [saccade_onsets_image_viewing saccade_onset_sample];
                        saccade_dist_image_viewing = [saccade_dist_image_viewing fixationstats.SaaccadeClusterValues(i_saccade, 3)];
                        
                        if (fixationstats.saccadetimes(1, i_saccade) > ...
                            -event_info.trial_event_begin * event_info.sfreq + saccade_onset_control_image_onset * event_info.sfreq) 
                            saccade_onsets_image_viewing_not_image_onset = [saccade_onsets_image_viewing_not_image_onset saccade_onset_sample];
                            saccade_dist_image_viewing_not_image_onset = [saccade_dist_image_viewing_not_image_onset fixationstats.SaaccadeClusterValues(i_saccade, 3)];
                            newRow = {saccade_id, subj_id, current_task, current_session, i, ...
                                'image_viewing_not_image_onset', ...
                                fixationstats.eccentricity(i_saccade),  ...
                                fixationstats.gaze_pos_before(1, i_saccade),  ...
                                fixationstats.gaze_pos_before(2, i_saccade),...
                                fixationstats.gaze_pos_after(1, i_saccade),  ...
                                fixationstats.gaze_pos_after(2, i_saccade),...
                                fixationstats.azimuth_before(i_saccade),...
                                fixationstats.azimuth_after(i_saccade),...
                                fixationstats.elevation_before(i_saccade),...
                                fixationstats.elevation_after(i_saccade),...
                                (fixationstats.saccadetimes(2, i_saccade) - ...
                                    fixationstats.saccadetimes(1, i_saccade)) / event_info.sfreq,...
                                fixationstats.saccadetimes(1, i_saccade) / event_info.sfreq + ...
                                event_info.trial_event_begin,...
                                event_table(i).behavior_condition,...
                                event_table(i).remember_condition};
                        else
                            newRow = {saccade_id, subj_id, current_task, current_session, i, ...
                            'image_viewing_around_image_onset', ...
                            fixationstats.eccentricity(i_saccade),  ...
                            fixationstats.gaze_pos_before(1, i_saccade),  ...
                            fixationstats.gaze_pos_before(2, i_saccade),...
                            fixationstats.gaze_pos_after(1, i_saccade),  ...
                            fixationstats.gaze_pos_after(2, i_saccade),...
                            fixationstats.azimuth_before(i_saccade),...
                            fixationstats.azimuth_after(i_saccade),...
                            fixationstats.elevation_before(i_saccade),...
                            fixationstats.elevation_after(i_saccade),...
                            (fixationstats.saccadetimes(2, i_saccade) - ...
                                fixationstats.saccadetimes(1, i_saccade)) / event_info.sfreq,...
                            fixationstats.saccadetimes(1, i_saccade) / event_info.sfreq + ...
                                event_info.trial_event_begin,...
                            event_table(i).behavior_condition,...
                            event_table(i).remember_condition};

                        end
                        % n_epoch_saccade_onset_image_viewing = n_epoch_saccade_onset_image_viewing + 1;
                    else
                        saccade_onsets_control = [saccade_onsets_control saccade_onset_sample];
                        saccade_dist_control = [saccade_dist_control fixationstats.SaaccadeClusterValues(i_saccade, 3)];
                        % n_epoch_saccade_onset_control = n_epoch_saccade_onset_control + 1;
                        newRow = {saccade_id, subj_id, current_task, current_session, i, ...
                                'control', ...
                                fixationstats.eccentricity(i_saccade),  ...
                                fixationstats.gaze_pos_before(1, i_saccade),  ...
                                fixationstats.gaze_pos_before(2, i_saccade),...
                                fixationstats.gaze_pos_after(1, i_saccade),  ...
                                fixationstats.gaze_pos_after(2, i_saccade),...
                                fixationstats.azimuth_before(i_saccade),...
                                fixationstats.azimuth_after(i_saccade),...
                                fixationstats.elevation_before(i_saccade),...
                                fixationstats.elevation_after(i_saccade),...
                                (fixationstats.saccadetimes(2, i_saccade) - ...
                                    fixationstats.saccadetimes(1, i_saccade)) / event_info.sfreq,...
                                fixationstats.saccadetimes(1, i_saccade) / event_info.sfreq + ...
                                    event_info.trial_event_begin,...
                                event_table(i).behavior_condition,...
                                event_table(i).remember_condition};

                    end

                    saccade_table = [saccade_table; newRow];
                    saccade_id = saccade_id + 1;
                end
            end
            saccade_distance.saccade_dist_image_viewing = saccade_dist_image_viewing;
            saccade_distance.saccade_dist_image_viewing_not_image_onset = saccade_dist_image_viewing_not_image_onset;
            saccade_distance.saccade_dist_control = saccade_dist_control;

            filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_saccade_distance.mat', ...
                            subj_id, current_task, current_session));

            save(filename, 'saccade_distance', '-v7.3');
            filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_saccade_table.mat', ...
                subj_id, current_task, current_session));
            save(filename, 'saccade_table', '-v7.3');

            % epoch EOG
            EOG_signal_saccade_onset_around_image_onset = zeros((1 - event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ...
                    (size(saccade_onsets_image_viewing, 2) - size(saccade_onsets_image_viewing_not_image_onset, 2)));
                
            EOG_signal_saccade_onset_not_image_onset = zeros((1 - event_info.trial_event_begin * event_info.sfreq + ...
                               event_info.trial_event_end * event_info.sfreq), size(saccade_onsets_image_viewing_not_image_onset, 2));
           
            index_around_image_onset = 1;
            index_control = 1;
            index_not_image_onset = 1;
            for i_epoch = 1:length(saccade_onsets_image_viewing)
                saccade_onset_sample = saccade_onsets_image_viewing(i_epoch);

                if ismember(saccade_onset_sample, saccade_onsets_image_viewing_not_image_onset)

                    EOG_signal_saccade_onset_not_image_onset(:, index_not_image_onset) = EOG_signal(...
                        saccade_onset_sample:(saccade_onset_sample - ...
                               event_info.trial_event_begin * event_info.sfreq + ...
                               event_info.trial_event_end * event_info.sfreq));
                    index_not_image_onset = index_not_image_onset + 1;
               else

                    EOG_signal_saccade_onset_around_image_onset(:, index_around_image_onset) = EOG_signal(...
                        saccade_onset_sample:(saccade_onset_sample - ...
                               event_info.trial_event_begin * event_info.sfreq + ...
                               event_info.trial_event_end * event_info.sfreq));
                    index_around_image_onset = index_around_image_onset + 1;
                    % filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_epoch%d_channel%d_saccade_onset.mat', ...
                    %     subj_id, current_task, current_session, i_epoch, ch));
                end
            end

            filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_EOG_saccade_onset_not_image_onset.mat', ...
                        subj_id, current_task, current_session));   
            if ~exist(filename, 'file')
                    save(filename, 'EOG_signal_saccade_onset_not_image_onset', '-v7.3');
            else
                fprintf('File %s already exists. Skip\n', filename);
            end

            filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_EOG_saccade_onset_around_image_onset.mat', ...
                        subj_id, current_task, current_session));
            if ~exist(filename, 'file')
                    save(filename, 'EOG_signal_saccade_onset_around_image_onset', '-v7.3');
            else
                fprintf('File %s already exists. Skip\n', filename);
            end
            EOG_signal_saccade_onset_control = zeros((1 - event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ...
                                   size(saccade_onsets_control, 2));
            for i_epoch = 1:length(saccade_onsets_control)

                saccade_onset_sample = saccade_onsets_control(i_epoch);
                EOG_signal_saccade_onset_control(:, index_control) = EOG_signal(...
                        saccade_onset_sample:(saccade_onset_sample - ...
                               event_info.trial_event_begin * event_info.sfreq + ...
                               event_info.trial_event_end * event_info.sfreq));
                index_control = index_control + 1;

            end

            filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_EOG_saccade_onset_control.mat', ...
                        subj_id, current_task, current_session));
            if ~exist(filename, 'file')
                    save(filename, 'EOG_signal_saccade_onset_control', '-v7.3');
            else
                fprintf('File %s already exists. Skip\n', filename);
            end

            for ch = 1:size(epochs{1}, 2)
                reconstructed_signal = zeros(size(laplacian_referenced_signal, 1), 1);

                ifft_wavelet_filtered_ch = zeros(length(freqs), size(laplacian_referenced_signal, 1));
                for f_idx = 1:length(freqs)
                    f = freqs(f_idx);
                    [wavelet, wavelet_length] = morlet_wavelet_cohen(f, event_info.sfreq, num_cycle_wavelet);
                    wavelet = wavelet / sum(abs(wavelet));
                    padded_signal = [zeros(wavelet_length, 1); laplacian_referenced_signal(:, ch)]; 
                    lenConv = length(padded_signal) + wavelet_length;
                    fft_signal = fft(padded_signal, lenConv);
                    fft_wavelet = fft(wavelet, lenConv);
                    conv_result = ifft(fft_signal .* fft_wavelet');
                    % Extract the relevant portion of the convolution result
                    start_idx = floor(wavelet_length / 2) + 1;
                    end_idx = start_idx + size(laplacian_referenced_signal, 1) - 1;

                    ifft_wavelet_filtered_ch(f_idx, :) = conv_result(start_idx:end_idx);
                    
                    reconstructed_signal = reconstructed_signal + real(ifft_wavelet_filtered_ch(f_idx, :))';
                end
                energy_original = sum(laplacian_referenced_signal(:, ch).^2);
                energy_reconstructed = sum(reconstructed_signal.^2);
                scale_factor = sqrt(energy_original / energy_reconstructed);
                ifft_wavelet_filtered_ch = ifft_wavelet_filtered_ch * scale_factor;
                % validate
                % figure; plot(laplacian_referenced_signal(:, ch)); hold on; plot(real(ifft_wavelet_filtered_ch(f_idx, :)))
                % control - during fixation cross stimulus

                ifft_wavelet_filtered_saccade_onset_around_image_onset = zeros(length(freqs), ...
                    (1 - event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ...
                    (size(saccade_onsets_image_viewing, 2) - size(saccade_onsets_image_viewing_not_image_onset, 2)));
                
                ifft_wavelet_filtered_saccade_onset_not_image_onset = zeros(length(freqs), ...
                    (1 - event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), size(saccade_onsets_image_viewing_not_image_onset, 2));
                
                prep_signal_saccade_onset_around_image_onset = zeros((1 - event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ...
                    (size(saccade_onsets_image_viewing, 2) - size(saccade_onsets_image_viewing_not_image_onset, 2)));
                
                prep_signal_saccade_onset_not_image_onset = zeros((1 - event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), size(saccade_onsets_image_viewing_not_image_onset, 2));
              

                index_around_image_onset = 1;
                index_control = 1;
                index_not_image_onset = 1;
                for i_epoch = 1:length(saccade_onsets_image_viewing)
                    saccade_onset_sample = saccade_onsets_image_viewing(i_epoch);
                    ifft_wavelet_filtered = ifft_wavelet_filtered_ch(:, saccade_onset_sample:(saccade_onset_sample - ...
                                   event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq));
                    if ismember(saccade_onset_sample, saccade_onsets_image_viewing_not_image_onset)
                        ifft_wavelet_filtered_saccade_onset_not_image_onset(:, :, index_not_image_onset) = ...
                            ifft_wavelet_filtered;
                        prep_signal_saccade_onset_not_image_onset(:, index_not_image_onset) = laplacian_referenced_signal(...
                            saccade_onset_sample:(saccade_onset_sample - ...
                                   event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ch);
                        index_not_image_onset = index_not_image_onset + 1;

                        % filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_epoch%d_channel%d_saccade_onset_not_image_onset.mat', ...
                        %     subj_id, current_task, current_session, i_epoch, ch));
                    else
                        ifft_wavelet_filtered_saccade_onset_around_image_onset(:, :, index_around_image_onset) = ...
                            ifft_wavelet_filtered;
                        prep_signal_saccade_onset_around_image_onset(:, index_around_image_onset) = laplacian_referenced_signal(...
                            saccade_onset_sample:(saccade_onset_sample - ...
                                   event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ch);
                        index_around_image_onset = index_around_image_onset + 1;
                        % filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_epoch%d_channel%d_saccade_onset.mat', ...
                        %     subj_id, current_task, current_session, i_epoch, ch));
                    end
                end
                filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_channel%d_saccade_onset_not_image_onset.mat', ...
                            subj_id, current_task, current_session, ch));   
                if ~exist(filename, 'file')
                        save(filename, 'ifft_wavelet_filtered_saccade_onset_not_image_onset', '-v7.3');
                else
                    fprintf('File %s already exists. Skip\n', filename);
                end

                filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_channel%d_saccade_onset_around_image_onset.mat', ...
                            subj_id, current_task, current_session, ch));
                if ~exist(filename, 'file')
                        save(filename, 'ifft_wavelet_filtered_saccade_onset_around_image_onset', '-v7.3');
                else
                    fprintf('File %s already exists. Skip\n', filename);
                end
                filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_channel%d_saccade_onset_around_image_onset_prep_signal.mat', ...
                            subj_id, current_task, current_session, ch));
                if ~exist(filename, 'file')
                        save(filename, 'prep_signal_saccade_onset_around_image_onset', '-v7.3');
                else
                    fprintf('File %s already exists. Skip\n', filename);
                end

                filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_channel%d_saccade_onset_not_image_onset_prep_signal.mat', ...
                            subj_id, current_task, current_session, ch));
                if ~exist(filename, 'file')
                        save(filename, 'prep_signal_saccade_onset_not_image_onset', '-v7.3');
                else
                    fprintf('File %s already exists. Skip\n', filename);
                end
                clear ifft_wavelet_filtered_saccade_onset_not_image_onset
                clear ifft_wavelet_filtered_saccade_onset_around_image_onset
                clear prep_signal_saccade_onset_not_image_onset
                clear prep_signal_saccade_onset_around_image_onset

                ifft_wavelet_filtered_saccade_onset_control = zeros(length(freqs), (1 - event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ...
                     size(saccade_onsets_control, 2));
                prep_signal_saccade_onset_control = zeros((1 - event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ...
                                   size(saccade_onsets_control, 2));

                for i_epoch = 1:length(saccade_onsets_control)
                    % filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_epoch%d_channel%d_saccade_onset_control.mat', ...
                    %     subj_id, current_task, current_session, i_epoch, ch));
                    % if ~exist(filename, 'file')
                    saccade_onset_sample = saccade_onsets_control(i_epoch);
                    ifft_wavelet_filtered = ifft_wavelet_filtered_ch(:, saccade_onset_sample:(saccade_onset_sample - ...
                                   event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq));
                    ifft_wavelet_filtered_saccade_onset_control(:, :, index_control) = ifft_wavelet_filtered;
                    prep_signal_saccade_onset_control(:, index_control) = laplacian_referenced_signal(...
                            saccade_onset_sample:(saccade_onset_sample - ...
                                   event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ch);
                    index_control = index_control + 1;

                end

                filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_channel%d_saccade_onset_control.mat', ...
                            subj_id, current_task, current_session, ch));
                if ~exist(filename, 'file')
                        save(filename, 'ifft_wavelet_filtered_saccade_onset_control', '-v7.3');
                else
                    fprintf('File %s already exists. Skip\n', filename);
                end
                filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_channel%d_saccade_onset_control_prep_signal.mat', ...
                            subj_id, current_task, current_session, ch));
                if ~exist(filename, 'file')
                        save(filename, 'prep_signal_saccade_onset_control', '-v7.3');
                else
                    fprintf('File %s already exists. Skip\n', filename);
                end
                clear ifft_wavelet_filtered_saccade_onset_control
                clear prep_signal_saccade_onset_control
            end

        elseif contains(current_task, 'BLAES_study') 
            saccade_onsets_image_viewing = [];
            saccade_onsets_image_viewing_not_image_onset = [];
            saccade_onsets_control = [];

            saccade_dist_image_viewing = [];
            saccade_dist_image_viewing_not_image_onset = [];
            saccade_dist_control = [];

            saccade_table = table([], [], [], [], [], [], [], [], [], [], [], ...
                [], [], [], [], [], [], [], [], 'VariableNames', ...
                 {'id', 'subj', 'task', 'session', 'trial', 'label', ...
                 'eccentricity', 'gaze_before_x_screen_coord', 'gaze_before_y_screen_coord', ...
                 'gaze_after_x_screen_coord', 'gaze_after_y_screen_coord','azimuth_before',...
                 'azimuth_after','elevation_before', 'elevation_after', ...
                 'duration','time_to_image_onset', 'behavior_condition', 'remember_condition'});

            saccade_id = 1;
            for i = 1:size(event_table, 1)
                % we exclude stimulation trials
                if strcmp(event_table(i).stimulation_condition, '1')
                    continue
                end
                fixationstats = event_info.fixationstats{i};
                for i_saccade = 1:size(fixationstats.saccadetimes, 2)

                    trial_onset_sample = event_table(i).image_onset_photodiode_sample + ...
                                   event_info.trial_event_begin * event_info.sfreq;
                    trial_offset_sample = event_table(i).image_offset_photodiode_sample + ...
                                    event_info.trial_event_end * event_info.sfreq;
                    if (fixationstats.valid_saccade_label(1, i_saccade) == 0) ||...
                            (fixationstats.valid_saccade_label(2, i_saccade) == 0) 
                        continue
                    end
                    % The onset of saccade trial, not the onset of saccade
                    saccade_onset_sample = trial_onset_sample + fixationstats.saccadetimes(1, i_saccade) + ...
                                   event_info.trial_event_begin * event_info.sfreq;

                    % saccade during image viewing
                    if (fixationstats.saccadetimes(1, i_saccade) > ...
                            -event_info.trial_event_begin * event_info.sfreq) && ( ...
                            fixationstats.saccadetimes(1, i_saccade) < ...
                            (trial_offset_sample - trial_onset_sample - event_info.trial_event_end * event_info.sfreq))
                        saccade_onsets_image_viewing = [saccade_onsets_image_viewing saccade_onset_sample];
                        saccade_dist_image_viewing = [saccade_dist_image_viewing fixationstats.SaaccadeClusterValues(i_saccade, 3)];

                        % n_epoch_saccade_onset_image_viewing = n_epoch_saccade_onset_image_viewing + 1;
                        if (fixationstats.saccadetimes(1, i_saccade) > ...
                            -event_info.trial_event_begin * event_info.sfreq + saccade_onset_control_image_onset * event_info.sfreq) 
                            saccade_onsets_image_viewing_not_image_onset = [saccade_onsets_image_viewing_not_image_onset saccade_onset_sample];
                            saccade_dist_image_viewing_not_image_onset = [saccade_dist_image_viewing_not_image_onset fixationstats.SaaccadeClusterValues(i_saccade, 3)];
                            newRow = {saccade_id, subj_id, current_task, current_session, i, ...
                                'image_viewing_not_image_onset', ...
                                fixationstats.eccentricity(i_saccade),  ...
                                fixationstats.gaze_pos_before(1, i_saccade),  ...
                                fixationstats.gaze_pos_before(2, i_saccade),...
                                fixationstats.gaze_pos_after(1, i_saccade),  ...
                                fixationstats.gaze_pos_after(2, i_saccade),...
                                fixationstats.azimuth_before(i_saccade),...
                                fixationstats.azimuth_after(i_saccade),...
                                fixationstats.elevation_before(i_saccade),...
                                fixationstats.elevation_after(i_saccade),...
                                (fixationstats.saccadetimes(2, i_saccade) - ...
                                    fixationstats.saccadetimes(1, i_saccade)) / event_info.sfreq,...
                                fixationstats.saccadetimes(1, i_saccade) / event_info.sfreq + ...
                                event_info.trial_event_begin,...
                                event_table(i).behavior_condition,...
                                event_table(i).remember_condition};
                        else
                            newRow = {saccade_id, subj_id, current_task, current_session, i, ...
                                'image_viewing_around_image_onset', ...
                                fixationstats.eccentricity(i_saccade),  ...
                                fixationstats.gaze_pos_before(1, i_saccade),  ...
                                fixationstats.gaze_pos_before(2, i_saccade),...
                                fixationstats.gaze_pos_after(1, i_saccade),  ...
                                fixationstats.gaze_pos_after(2, i_saccade),...
                                fixationstats.azimuth_before(i_saccade),...
                                fixationstats.azimuth_after(i_saccade),...
                                fixationstats.elevation_before(i_saccade),...
                                fixationstats.elevation_after(i_saccade),...
                                (fixationstats.saccadetimes(2, i_saccade) - ...
                                    fixationstats.saccadetimes(1, i_saccade)) / event_info.sfreq,...
                                fixationstats.saccadetimes(1, i_saccade) / event_info.sfreq + ...
                                event_info.trial_event_begin,...
                                event_table(i).behavior_condition,...
                                event_table(i).remember_condition};
                        end
                    else
                        saccade_onsets_control = [saccade_onsets_control saccade_onset_sample];
                        saccade_dist_control = [saccade_dist_control fixationstats.SaaccadeClusterValues(i_saccade, 3)];
                        % n_epoch_saccade_onset_control = n_epoch_saccade_onset_control + 1;
                        newRow = {saccade_id, subj_id, current_task, current_session, i, ...
                                'control', ...
                                fixationstats.eccentricity(i_saccade),  ...
                                fixationstats.gaze_pos_before(1, i_saccade),  ...
                                fixationstats.gaze_pos_before(2, i_saccade),...
                                fixationstats.gaze_pos_after(1, i_saccade),  ...
                                fixationstats.gaze_pos_after(2, i_saccade),...
                                fixationstats.azimuth_before(i_saccade),...
                                fixationstats.azimuth_after(i_saccade),...
                                fixationstats.elevation_before(i_saccade),...
                                fixationstats.elevation_after(i_saccade),...
                                (fixationstats.saccadetimes(2, i_saccade) - ...
                                    fixationstats.saccadetimes(1, i_saccade)) / event_info.sfreq,...
                                fixationstats.saccadetimes(1, i_saccade) / event_info.sfreq + ...
                                event_info.trial_event_begin,...
                                event_table(i).behavior_condition,...
                                event_table(i).remember_condition};
                    end

                    saccade_table = [saccade_table; newRow];
                    saccade_id = saccade_id + 1;
                end
            end
            saccade_distance.saccade_dist_image_viewing = saccade_dist_image_viewing;
            saccade_distance.saccade_dist_image_viewing_not_image_onset = saccade_dist_image_viewing_not_image_onset;
            saccade_distance.saccade_dist_control = saccade_dist_control;

            filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_saccade_distance.mat', ...
                            subj_id, current_task, current_session));

            save(filename, 'saccade_distance', '-v7.3');

            filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_saccade_table.mat', ...
                subj_id, current_task, current_session));
            save(filename, 'saccade_table', '-v7.3');

            % epoch EOG
            EOG_signal_saccade_onset_around_image_onset = zeros((1 - event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ...
                    (size(saccade_onsets_image_viewing, 2) - size(saccade_onsets_image_viewing_not_image_onset, 2)));
                
            EOG_signal_saccade_onset_not_image_onset = zeros((1 - event_info.trial_event_begin * event_info.sfreq + ...
                               event_info.trial_event_end * event_info.sfreq), size(saccade_onsets_image_viewing_not_image_onset, 2));
           
            index_around_image_onset = 1;
            index_control = 1;
            index_not_image_onset = 1;
            for i_epoch = 1:length(saccade_onsets_image_viewing)
                saccade_onset_sample = saccade_onsets_image_viewing(i_epoch);

                if ismember(saccade_onset_sample, saccade_onsets_image_viewing_not_image_onset)

                    EOG_signal_saccade_onset_not_image_onset(:, index_not_image_onset) = EOG_signal(...
                        saccade_onset_sample:(saccade_onset_sample - ...
                               event_info.trial_event_begin * event_info.sfreq + ...
                               event_info.trial_event_end * event_info.sfreq));
                    index_not_image_onset = index_not_image_onset + 1;
               else

                    EOG_signal_saccade_onset_around_image_onset(:, index_around_image_onset) = EOG_signal(...
                        saccade_onset_sample:(saccade_onset_sample - ...
                               event_info.trial_event_begin * event_info.sfreq + ...
                               event_info.trial_event_end * event_info.sfreq));
                    index_around_image_onset = index_around_image_onset + 1;
                    % filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_epoch%d_channel%d_saccade_onset.mat', ...
                    %     subj_id, current_task, current_session, i_epoch, ch));
                end
            end

            filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_EOG_saccade_onset_not_image_onset.mat', ...
                        subj_id, current_task, current_session));   
            if ~exist(filename, 'file')
                    save(filename, 'EOG_signal_saccade_onset_not_image_onset', '-v7.3');
            else
                fprintf('File %s already exists. Skip\n', filename);
            end

            filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_EOG_saccade_onset_around_image_onset.mat', ...
                        subj_id, current_task, current_session));
            if ~exist(filename, 'file')
                    save(filename, 'EOG_signal_saccade_onset_around_image_onset', '-v7.3');
            else
                fprintf('File %s already exists. Skip\n', filename);
            end
            EOG_signal_saccade_onset_control = zeros((1 - event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ...
                                   size(saccade_onsets_control, 2));
            for i_epoch = 1:length(saccade_onsets_control)

                saccade_onset_sample = saccade_onsets_control(i_epoch);
                EOG_signal_saccade_onset_control(:, index_control) = EOG_signal(...
                        saccade_onset_sample:(saccade_onset_sample - ...
                               event_info.trial_event_begin * event_info.sfreq + ...
                               event_info.trial_event_end * event_info.sfreq));
                index_control = index_control + 1;

            end

            filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_EOG_saccade_onset_control.mat', ...
                        subj_id, current_task, current_session));
            if ~exist(filename, 'file')
                    save(filename, 'EOG_signal_saccade_onset_control', '-v7.3');
            else
                fprintf('File %s already exists. Skip\n', filename);
            end


            for ch = 1:size(epochs{1}, 2)
                reconstructed_signal = zeros(size(laplacian_referenced_signal, 1), 1);

                ifft_wavelet_filtered_ch = zeros(length(freqs), size(laplacian_referenced_signal, 1));
                for f_idx = 1:length(freqs)
                    f = freqs(f_idx);
                    [wavelet, wavelet_length] = morlet_wavelet_cohen(f, event_info.sfreq, num_cycle_wavelet);
                    wavelet = wavelet / sum(abs(wavelet));
                    padded_signal = [zeros(wavelet_length, 1); laplacian_referenced_signal(:, ch)]; 
                    lenConv = length(padded_signal) + wavelet_length;
                    fft_signal = fft(padded_signal, lenConv);
                    fft_wavelet = fft(wavelet, lenConv);
                    conv_result = ifft(fft_signal .* fft_wavelet');
                    % Extract the relevant portion of the convolution result
                    start_idx = floor(wavelet_length / 2) + 1;
                    end_idx = start_idx + size(laplacian_referenced_signal, 1) - 1;

                    ifft_wavelet_filtered_ch(f_idx, :) = conv_result(start_idx:end_idx);
                    
                    reconstructed_signal = reconstructed_signal + real(ifft_wavelet_filtered_ch(f_idx, :))';
                end
                energy_original = sum(laplacian_referenced_signal(:, ch).^2);
                energy_reconstructed = sum(reconstructed_signal.^2);
                scale_factor = sqrt(energy_original / energy_reconstructed);
                ifft_wavelet_filtered_ch = ifft_wavelet_filtered_ch * scale_factor;
                % validate
                % figure; plot(laplacian_referenced_signal(:, ch)); hold on; plot(real(ifft_wavelet_filtered_ch(f_idx, :)))
                
                % control - during fixation cross stimulus
                ifft_wavelet_filtered_saccade_onset_around_image_onset = zeros(length(freqs), ...
                    (1 - event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ...
                    (size(saccade_onsets_image_viewing, 2) - size(saccade_onsets_image_viewing_not_image_onset, 2)));
                
                ifft_wavelet_filtered_saccade_onset_not_image_onset = zeros(length(freqs), ...
                    (1 - event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), size(saccade_onsets_image_viewing_not_image_onset, 2));
                
                prep_signal_saccade_onset_around_image_onset = zeros((1 - event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ...
                    (size(saccade_onsets_image_viewing, 2) - size(saccade_onsets_image_viewing_not_image_onset, 2)));
                
                prep_signal_saccade_onset_not_image_onset = zeros((1 - event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), size(saccade_onsets_image_viewing_not_image_onset, 2));
               
                
                index_around_image_onset = 1;
                index_control = 1;
                index_not_image_onset = 1;
                for i_epoch = 1:length(saccade_onsets_image_viewing)
                    saccade_onset_sample = saccade_onsets_image_viewing(i_epoch);
                    ifft_wavelet_filtered = ifft_wavelet_filtered_ch(:, saccade_onset_sample:(saccade_onset_sample - ...
                                   event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq));
                    if ismember(saccade_onset_sample, saccade_onsets_image_viewing_not_image_onset)
                        ifft_wavelet_filtered_saccade_onset_not_image_onset(:, :, index_not_image_onset) = ...
                            ifft_wavelet_filtered;
                        prep_signal_saccade_onset_not_image_onset(:, index_not_image_onset) = laplacian_referenced_signal(...
                            saccade_onset_sample:(saccade_onset_sample - ...
                                   event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ch);
                        index_not_image_onset = index_not_image_onset + 1;

                        % filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_epoch%d_channel%d_saccade_onset_not_image_onset.mat', ...
                        %     subj_id, current_task, current_session, i_epoch, ch));
                    else
                        ifft_wavelet_filtered_saccade_onset_around_image_onset(:, :, index_around_image_onset) = ...
                            ifft_wavelet_filtered;
                        prep_signal_saccade_onset_around_image_onset(:, index_around_image_onset) = laplacian_referenced_signal(...
                            saccade_onset_sample:(saccade_onset_sample - ...
                                   event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ch);
                        index_around_image_onset = index_around_image_onset + 1;
                        % filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_epoch%d_channel%d_saccade_onset.mat', ...
                        %     subj_id, current_task, current_session, i_epoch, ch));
                    end
                end
                filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_channel%d_saccade_onset_not_image_onset.mat', ...
                            subj_id, current_task, current_session, ch));   
                if ~exist(filename, 'file')
                        save(filename, 'ifft_wavelet_filtered_saccade_onset_not_image_onset', '-v7.3');
                else
                    fprintf('File %s already exists. Skip\n', filename);
                end

                filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_channel%d_saccade_onset_around_image_onset.mat', ...
                            subj_id, current_task, current_session, ch));
                if ~exist(filename, 'file')
                        save(filename, 'ifft_wavelet_filtered_saccade_onset_around_image_onset', '-v7.3');
                else
                    fprintf('File %s already exists. Skip\n', filename);
                end

                filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_channel%d_saccade_onset_around_image_onset_prep_signal.mat', ...
                            subj_id, current_task, current_session, ch));
                if ~exist(filename, 'file')
                        save(filename, 'prep_signal_saccade_onset_around_image_onset', '-v7.3');
                else
                    fprintf('File %s already exists. Skip\n', filename);
                end

                filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_channel%d_saccade_onset_not_image_onset_prep_signal.mat', ...
                            subj_id, current_task, current_session, ch));
                if ~exist(filename, 'file')
                        save(filename, 'prep_signal_saccade_onset_not_image_onset', '-v7.3');
                else
                    fprintf('File %s already exists. Skip\n', filename);
                end

                clear ifft_wavelet_filtered_saccade_onset_not_image_onset
                clear ifft_wavelet_filtered_saccade_onset_around_image_onset
                clear prep_signal_saccade_onset_around_image_onset
                clear prep_signal_saccade_onset_not_image_onset

                ifft_wavelet_filtered_saccade_onset_control = zeros(length(freqs), (1 - event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ...
                     size(saccade_onsets_control, 2));
                prep_signal_saccade_onset_control = zeros((1 - event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ...
                                   size(saccade_onsets_control, 2));
  

                for i_epoch = 1:length(saccade_onsets_control)
                    % filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_epoch%d_channel%d_saccade_onset_control.mat', ...
                    %     subj_id, current_task, current_session, i_epoch, ch));
                    % if ~exist(filename, 'file')
                    saccade_onset_sample = saccade_onsets_control(i_epoch);
                    ifft_wavelet_filtered = ifft_wavelet_filtered_ch(:, saccade_onset_sample:(saccade_onset_sample - ...
                                   event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq));
                    ifft_wavelet_filtered_saccade_onset_control(:, :, index_control) = ifft_wavelet_filtered;
                    prep_signal_saccade_onset_control(:, index_control) = laplacian_referenced_signal(...
                            saccade_onset_sample:(saccade_onset_sample - ...
                                   event_info.trial_event_begin * event_info.sfreq + ...
                                   event_info.trial_event_end * event_info.sfreq), ch);
                    index_control = index_control + 1;
                    %     save(filename, 'ifft_wavelet_filtered', '-v7.3');
                    % else
                    %     fprintf('File %s already exists. Skip\n', filename);
                    % end
                end

                filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_channel%d_saccade_onset_control.mat', ...
                            subj_id, current_task, current_session, ch));
                if ~exist(filename, 'file')
                        save(filename, 'ifft_wavelet_filtered_saccade_onset_control', '-v7.3');
                else
                    fprintf('File %s already exists. Skip\n', filename);
                end

                filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_channel%d_saccade_onset_control_prep_signal.mat', ...
                            subj_id, current_task, current_session, ch));
                if ~exist(filename, 'file')
                        save(filename, 'prep_signal_saccade_onset_control', '-v7.3');
                else
                    fprintf('File %s already exists. Skip\n', filename);
                end

                clear ifft_wavelet_filtered_saccade_onset_control
                clear prep_signal_saccade_onset_control
            end
        end
        close all;
    end
end

%% discard
%         % for visualization purpose
%         figure;
%         colormap jet;  % Set the color map for better visual distinction in imagesc
%         channel_indices = [10 ,20];
%         % Total number of subplots
%         total_plots = length(channel_indices) * 2;
% 
%         % Loop through selected channels
%         for idx = 1:length(channel_indices)
%             % Raw Signal Plot
%             subplot(2, length(channel_indices), idx);
%             plot(epochs{i_epoch}(:, channel_indices(idx)));  
%             title(sprintf('Raw Signal - Channel %d', channel_indices(idx)));
%             xlabel('Time (samples)');
%             ylabel('Amplitude');
% 
%             % Time-Frequency Representation
%             subplot(2, length(channel_indices), idx + length(channel_indices));
%             channel_data = squeeze(abs(ifft_wavelet_filtered(channel_indices(idx), :, :)));
%             channel_data = (channel_data - mean(channel_data(1:2000))) / std(channel_data(1:2000));
%             imagesc('XData', 1:size(channel_data, 2), 'YData', freqs, 'CData', channel_data);
%             axis xy;  % Flip the axis so lower frequencies are at the bottom
%             colorbar;  % Add a color bar to indicate scale
%             title(sprintf('Time-Frequency - Channel %d', channel_indices(idx)));
%             xlabel('Time (samples)');
%             ylabel('Frequency (Hz)');
%         end
% 
%         sgtitle('Raw Signal and Time-Frequency Representation of Selected Channels');  
% 
% %% test wavelet methods - not relevant
% % % Signal parameters
% fs = 1000; % Sampling frequency in Hz
% T = 5; % Total duration in seconds
% t = 0:1/fs:T-1/fs; % Time vector
% % 
% % % Amplitude modulation for 4 Hz and 10 Hz components
% amp4 = linspace(1, 10, length(t)/2);
% amp4 = [amp4, fliplr(amp4)]; % Varying amplitude for 4 Hz
% % amp10 = circshift(amp4, 2000);
% % % Create the signal
% signal4Hz = amp4 .* sin(2 * pi * 4 * t); % 4 Hz component
% % signal10Hz = amp10 .* sin(2 * pi * 10 * t); % 10 Hz component
% % signal = signal4Hz + signal10Hz; % Combined signal
% % 
% % % Plot the original combined signal
% % figure;
% % plot(t, signal);
% % xlabel('Time (s)');
% % ylabel('Amplitude');
% % title('Combined 4 Hz and 10 Hz Signal with Varying Amplitude');
% % 
% % figure;
% % plot(t, amp4);
% % hold on;% Number of cycles in the wavelet (can adjust for specificity)
% % plot(t,amp10);
% % 
% % % Define the Morlet wavelet for 4 Hz
% [wavelet, wavelet_length] = morlet_wavelet_cohen(4, fs, num_cycle_wavelet);
% wavelet = wavelet / sum(abs(wavelet));
% lenConv = length(signal4Hz) + wavelet_length - 1;
% fft_signal = fft(signal4Hz, lenConv);
% fft_wavelet = fft(wavelet, lenConv);
% conv_result = ifft(fft_signal .* fft_wavelet);
% % Extract the relevant portion of the convolution result
% start_idx = floor(wavelet_length / 2) + 1;
% end_idx = start_idx + length(signal4Hz) - 1;
% figure;
% plot(squeeze(real(conv_result(start_idx:end_idx)))*10);
% hold on; plot(signal4Hz)
% 
% 
% 
% % Perform convolution with the signal
% filtered_wavelet = conv(signal, wavelet, 'same');
% wavelet_amplitude = circshift(abs(filtered_wavelet), - ceil(wavelet_length / 2 ));
% wavelet_amplitude = abs(filtered_wavelet);  % Extract amplitude (envelope)
% 
% % Plot the wavelet-filtered amplitude (4 Hz component)
% plot(t, wavelet_amplitude / 100);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Filtered Amplitude Isolating 10 Hz Component');
% 

% 
% wavelet_amplitude = abs(filtered_wavelet) * (padded_length / length(signal)); 
% figure;
% plot(t, amp4);
% hold on;% Number of cycles in the wavelet (can adjust for specificity)
% plot(t,amp10);
% 
% % Plot the wavelet-filtered amplitude (4 Hz component)
% plot(t, wavelet_amplitude / 100);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Filtered Amplitude Isolating 10 Hz Component');
% 
% % narrow band - not good
% d = designfilt('bandpassfir', 'FilterOrder', 50, ...
%                    'CutoffFrequency1', 10-0.5, 'CutoffFrequency2', 10+0.5, ...
%                    'SampleRate', fs);
% filtered_signal = filtfilt(d, signal);  
% analytic_signal = hilbert(filtered_signal);
% amplitude = abs(analytic_signal);
% % Plot the wavelet-filtered amplitude (4 Hz component)
% figure;
% plot(t, amp4);
% hold on;% Number of cycles in the wavelet (can adjust for specificity)
% plot(t,amp10);
% plot(t, amplitude);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Filtered Amplitude Isolating 10 Hz Component');
% 
% figure;
% plot(signal);
% hold on;
% plot(filtered_signal);
% 
% 
% amplitude = abs(filtered_wavelet);  % Amplitude envelope
% phase = angle(filtered_wavelet);    
% 
% reconstructed_10Hz = amplitude .* cos(2 * pi * 10 * t + phase);
% 
% figure;
% plot(t, reconstructed_10Hz, 'b-', 'DisplayName', 'Reconstructed 10 Hz Component');
% hold on;
% plot(t, signal10Hz, 'r-', 'DisplayName', 'Original 10 Hz Component');
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Comparison of Reconstructed and Original 4 Hz Components');
% legend show;
% 
% 
% % wavelet_amplitude = circshift(wavelet_amplitude, -wavelet_center);
% 
% % Plot the wavelet-filtered amplitude
% % figure;
% plot(t, wavelet_amplitude);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Amplitude from Wavelet Filtering');
% 
% 
% 
% % filter_methods = {narrow_filtered, wavelet_filtered, ifft_wavelet_filtered};
% % titles = {'Average TF Representation (Narrowband Filter)', ...
% %           'Average TF Representation (Wavelet Filter)', ...
% %           'Average TF Representation (IFFT Wavelet Filter)'};
% % 
% % figure;
% % num_methods = numel(filter_methods);
% % 
% % for i = 1:num_methods
% %     avg_data = mean(filter_methods{i}, 1);  
% %     avg_data = squeeze(avg_data); 
% % 
% %     subplot(1, num_methods, i);
% %     imagesc(1:size(avg_data, 2), freqs, 10*log10(avg_data)); % Convert to dB
% %     axis xy;  
% %     xlabel('Time (samples)');
% %     ylabel('Frequency (Hz)');
% %     title(titles{i});
% %     colorbar;  
% % end
% 
% 
% method_titles = {'Narrowband 4 Hz', 'Wavelet 4 Hz', 'IFFT Wavelet 4 Hz'};
% epochs_to_plot = [1, 2];  
% num_methods = numel(filter_methods) + 1; 
% 
% 
% for e = 1:numel(epochs_to_plot)
% epoch = epochs_to_plot(e);
% 
% figure;  
% subplot(2, num_methods, 1 + (e-1)*num_methods);
% plot(t, squeeze(epochs{epoch}(:, f_index)));
% title(sprintf('Epoch %d: Raw Signal', e));
% xlabel('Time (seconds)');
% ylabel('Amplitude');
% 
% for i = 1:numel(filter_methods)
%     subplot(2, num_methods, i+1 + (e-1)*num_methods);  % Positioning each subplot
%     plot(t, squeeze(filter_methods{i}(epoch, :, f_index)));
%     title(sprintf('Epoch %d: %s', e, method_titles{i}));
%     xlabel('Time (seconds)');
%     ylabel('Amplitude');
% end
% end
% 



            % for saccade onset
            % epochs = cell(n_epoch_saccade_onset_image_viewing, 1);
            % epoch_counter = 1;
            % for i = 1:size(event_table, 1)
            %     fixationstats = event_info.fixationstats{i};
            %     for i_saccade = 1:size(fixationstats.saccadetimes, 2)
            % 
            %         trial_onset_sample = event_table(i).fixed_image_onset_photodiode_sample + ...
            %                        event_info.trial_event_begin * event_info.sfreq;
            %         trial_offset_sample = event_table(i).response_image_offset_photodiode_sample + ...
            %                         event_info.trial_event_end * event_info.sfreq;
            %         if fixationstats.valid_saccade_label(1, i_saccade) == 0
            %             continue
            %         end
            % 
            %         saccade_onset_sample = trial_onset_sample + fixationstats.saccadetimes(1, i_saccade) + ...
            %                        event_info.trial_event_begin * event_info.sfreq;
            % 
            % 
            % 
            %         % saccade during image viewing
            %         if (fixationstats.saccadetimes(1, i_saccade) > ...
            %                 -event_info.trial_event_begin * event_info.sfreq) && ( ...
            %                 fixationstats.saccadetimes(1, i_saccade) < ...
            %                 trial_offset_sample - event_info.trial_event_end * event_info.sfreq)
            %             epochs{epoch_counter} = laplacian_referenced_signal(saccade_onset_sample:(saccade_onset_sample - ...
            %                        event_info.trial_event_begin * event_info.sfreq + ...
            %                        event_info.trial_event_end * event_info.sfreq), :);
            %             epoch_counter = epoch_counter+ 1;
            %         end
            %     end
            % end
            % 
            % for i_epoch = 1:numel(epochs)
            %     for ch = 1:size(epochs{1}, 2)
            %         ifft_wavelet_filtered = zeros(length(freqs), size(epochs{1}, 1));
            % 
            %         reconstructed_signal = zeros(size(epochs{i_epoch}(:, ch)));
            %         for f_idx = 1:length(freqs)
            %             f = freqs(f_idx);
            % 
            %             [wavelet, wavelet_length] = morlet_wavelet_cohen(f, event_info.sfreq, num_cycle_wavelet);
            %             wavelet = wavelet / sum(abs(wavelet));
            %             padded_signal = [zeros(wavelet_length, 1); epochs{i_epoch}(:, ch)]; 
            %             lenConv = length(padded_signal) + wavelet_length;
            %             fft_signal = fft(padded_signal, lenConv);
            %             fft_wavelet = fft(wavelet, lenConv);
            %             conv_result = ifft(fft_signal .* fft_wavelet');
            %             % Extract the relevant portion of the convolution result
            %             start_idx = floor(wavelet_length / 2) + 1;
            %             end_idx = start_idx + length(epochs{i_epoch}(:, ch)) - 1;
            % 
            %             ifft_wavelet_filtered(f_idx, :) = conv_result(start_idx:end_idx);
            %             % validate
            %             % figure; plot(epochs{i_epoch}(:, ch)); hold on; plot(real(ifft_wavelet_filtered(f_idx, :)))
            % 
            %             reconstructed_signal = reconstructed_signal + real(ifft_wavelet_filtered(f_idx, :))';
            % 
            %         end
            %         energy_original = sum(epochs{i_epoch}(:, ch).^2);
            %         energy_reconstructed = sum(reconstructed_signal.^2);
            %         scale_factor = sqrt(energy_original / energy_reconstructed);
            %         ifft_wavelet_filtered = ifft_wavelet_filtered * scale_factor;
            % 
            %         filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_epoch%d_channel%d_saccade_onset.mat', ...
            %             subj_id, current_task, current_session, i_epoch, ch));
            %         save(filename, 'ifft_wavelet_filtered', '-v7.3');
            %     end
            % end
            % 
            % epochs = cell(n_epoch_saccade_onset_control, 1);
            % epoch_counter = 1;
            % for i = 1:size(event_table, 1)
            %     fixationstats = event_info.fixationstats{i};
            %     for i_saccade = 1:size(fixationstats.saccadetimes, 2)
            % 
            %         trial_onset_sample = event_table(i).fixed_image_onset_photodiode_sample + ...
            %                        event_info.trial_event_begin * event_info.sfreq;
            %         trial_offset_sample = event_table(i).response_image_offset_photodiode_sample + ...
            %                         event_info.trial_event_end * event_info.sfreq;
            %         if fixationstats.valid_saccade_label(1, i_saccade) == 0
            %             continue
            %         end
            % 
            %         saccade_onset_sample = trial_onset_sample + fixationstats.saccadetimes(1, i_saccade) + ...
            %                        event_info.trial_event_begin * event_info.sfreq;
            % 
            % 
            %         % saccade during image viewing
            %         if ~((fixationstats.saccadetimes(1, i_saccade) > ...
            %                 -event_info.trial_event_begin * event_info.sfreq) && ( ...
            %                 fixationstats.saccadetimes(1, i_saccade) < ...
            %                 trial_offset_sample - event_info.trial_event_end * event_info.sfreq))
            %             epochs{epoch_counter} = laplacian_referenced_signal(saccade_onset_sample:(saccade_onset_sample - ...
            %                        event_info.trial_event_begin * event_info.sfreq + ...
            %                        event_info.trial_event_end * event_info.sfreq), :);
            %             epoch_counter = epoch_counter+ 1;
            %         end
            %     end
            % end
            % 
            % for i_epoch = 1:numel(epochs)
            %     for ch = 1:size(epochs{1}, 2)
            %         ifft_wavelet_filtered = zeros(length(freqs), size(epochs{1}, 1));
            % 
            %         reconstructed_signal = zeros(size(epochs{i_epoch}(:, ch)));
            %         for f_idx = 1:length(freqs)
            %             f = freqs(f_idx);
            % 
            %             [wavelet, wavelet_length] = morlet_wavelet_cohen(f, event_info.sfreq, num_cycle_wavelet);
            %             wavelet = wavelet / sum(abs(wavelet));
            %             padded_signal = [zeros(wavelet_length, 1); epochs{i_epoch}(:, ch)]; 
            %             lenConv = length(padded_signal) + wavelet_length;
            %             fft_signal = fft(padded_signal, lenConv);
            %             fft_wavelet = fft(wavelet, lenConv);
            %             conv_result = ifft(fft_signal .* fft_wavelet');
            %             % Extract the relevant portion of the convolution result
            %             start_idx = floor(wavelet_length / 2) + 1;
            %             end_idx = start_idx + length(epochs{i_epoch}(:, ch)) - 1;
            % 
            %             ifft_wavelet_filtered(f_idx, :) = conv_result(start_idx:end_idx);
            %             % validate
            %             % figure; plot(epochs{i_epoch}(:, ch)); hold on; plot(real(ifft_wavelet_filtered(f_idx, :)))
            % 
            %             reconstructed_signal = reconstructed_signal + real(ifft_wavelet_filtered(f_idx, :))';
            % 
            %         end
            %         energy_original = sum(epochs{i_epoch}(:, ch).^2);
            %         energy_reconstructed = sum(reconstructed_signal.^2);
            %         scale_factor = sqrt(energy_original / energy_reconstructed);
            %         ifft_wavelet_filtered = ifft_wavelet_filtered * scale_factor;
            % 
            %         filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_epoch%d_channel%d_saccade_onset_control.mat', ...
            %             subj_id, current_task, current_session, i_epoch, ch));
            %         save(filename, 'ifft_wavelet_filtered', '-v7.3');
            %     end
            % end

% 

 % we could do better, do wavelet transform first, then use
            % % saccade onset to get ifft for each ch, each epoch
            % n_epoch_saccade_onset_image_viewing = 0; 
            % n_epoch_saccade_onset_control = 0; 
            % for i = 1:size(event_table, 1)
            %     fixationstats = event_info.fixationstats{i};
            %     for i_saccade = 1:size(fixationstats.saccadetimes, 2)
            % 
            %         trial_onset_sample = event_table(i).image_onset_photodiode_sample + ...
            %                        event_info.trial_event_begin * event_info.sfreq;
            %         trial_offset_sample = event_table(i).image_offset_photodiode_sample + ...
            %                         event_info.trial_event_end * event_info.sfreq;
            %         if fixationstats.valid_saccade_label(1, i_saccade) == 0
            %             continue
            %         end
            % 
            %         saccade_onset_sample = trial_onset_sample + fixationstats.saccadetimes(1, i_saccade) + ...
            %                        event_info.trial_event_begin * event_info.sfreq;
            % 
            %         % saccade during image viewing
            %         if (fixationstats.saccadetimes(1, i_saccade) > ...
            %                 -event_info.trial_event_begin * event_info.sfreq) && ( ...
            %                 fixationstats.saccadetimes(1, i_saccade) < ...
            %                 trial_offset_sample - event_info.trial_event_end * event_info.sfreq)
            %             n_epoch_saccade_onset_image_viewing = n_epoch_saccade_onset_image_viewing + 1;
            %         else
            %             n_epoch_saccade_onset_control = n_epoch_saccade_onset_control + 1;
            %         end
            %     end
            % end
            % epochs = cell(n_epoch_saccade_onset_image_viewing, 1);
            % epoch_counter = 1;
            % for i = 1:size(event_table, 1)
            %     fixationstats = event_info.fixationstats{i};
            %     for i_saccade = 1:size(fixationstats.saccadetimes, 2)
            % 
            %         trial_onset_sample = event_table(i).image_onset_photodiode_sample + ...
            %                        event_info.trial_event_begin * event_info.sfreq;
            %         trial_offset_sample = event_table(i).image_offset_photodiode_sample + ...
            %                         event_info.trial_event_end * event_info.sfreq;
            %         if fixationstats.valid_saccade_label(1, i_saccade) == 0
            %             continue
            %         end
            % 
            %         saccade_onset_sample = trial_onset_sample + fixationstats.saccadetimes(1, i_saccade) + ...
            %                        event_info.trial_event_begin * event_info.sfreq;
            % 
            %         % saccade during image viewing
            %         if (fixationstats.saccadetimes(1, i_saccade) > ...
            %                 -event_info.trial_event_begin * event_info.sfreq) && ( ...
            %                 fixationstats.saccadetimes(1, i_saccade) < ...
            %                 trial_offset_sample - event_info.trial_event_end * event_info.sfreq)
            %             epochs{epoch_counter} = laplacian_referenced_signal(saccade_onset_sample:(saccade_onset_sample - ...
            %                        event_info.trial_event_begin * event_info.sfreq + ...
            %                        event_info.trial_event_end * event_info.sfreq), :);
            %             epoch_counter = epoch_counter+ 1;
            %         end
            %     end
            % end
            % 
            % for i_epoch = 1:numel(epochs)
            %     for ch = 1:size(epochs{1}, 2)
            %         ifft_wavelet_filtered = zeros(length(freqs), size(epochs{1}, 1));
            % 
            %         reconstructed_signal = zeros(size(epochs{i_epoch}(:, ch)));
            %         for f_idx = 1:length(freqs)
            %             f = freqs(f_idx);
            % 
            %             [wavelet, wavelet_length] = morlet_wavelet_cohen(f, event_info.sfreq, num_cycle_wavelet);
            %             wavelet = wavelet / sum(abs(wavelet));
            %             padded_signal = [zeros(wavelet_length, 1); epochs{i_epoch}(:, ch)]; 
            %             lenConv = length(padded_signal) + wavelet_length;
            %             fft_signal = fft(padded_signal, lenConv);
            %             fft_wavelet = fft(wavelet, lenConv);
            %             conv_result = ifft(fft_signal .* fft_wavelet');
            %             % Extract the relevant portion of the convolution result
            %             start_idx = floor(wavelet_length / 2) + 1;
            %             end_idx = start_idx + length(epochs{i_epoch}(:, ch)) - 1;
            % 
            %             ifft_wavelet_filtered(f_idx, :) = conv_result(start_idx:end_idx);
            %             % validate
            %             % figure; plot(epochs{i_epoch}(:, ch)); hold on; plot(real(ifft_wavelet_filtered(f_idx, :)))
            % 
            %             reconstructed_signal = reconstructed_signal + real(ifft_wavelet_filtered(f_idx, :))';
            % 
            %         end
            %         energy_original = sum(epochs{i_epoch}(:, ch).^2);
            %         energy_reconstructed = sum(reconstructed_signal.^2);
            %         scale_factor = sqrt(energy_original / energy_reconstructed);
            %         ifft_wavelet_filtered = ifft_wavelet_filtered * scale_factor;
            % 
            %         filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_epoch%d_channel%d_saccade_onset.mat', ...
            %             subj_id, current_task, current_session, i_epoch, ch));
            %         save(filename, 'ifft_wavelet_filtered', '-v7.3');
            %     end
            % end
            % 
            % epochs = cell(n_epoch_saccade_onset_control, 1);
            % epoch_counter = 1;
            % for i = 1:size(event_table, 1)
            %     fixationstats = event_info.fixationstats{i};
            %     for i_saccade = 1:size(fixationstats.saccadetimes, 2)
            % 
            %         trial_onset_sample = event_table(i).image_onset_photodiode_sample + ...
            %                        event_info.trial_event_begin * event_info.sfreq;
            %         trial_offset_sample = event_table(i).image_offset_photodiode_sample + ...
            %                         event_info.trial_event_end * event_info.sfreq;
            %         if fixationstats.valid_saccade_label(1, i_saccade) == 0
            %             continue
            %         end
            % 
            %         saccade_onset_sample = trial_onset_sample + fixationstats.saccadetimes(1, i_saccade) + ...
            %                        event_info.trial_event_begin * event_info.sfreq;
            % 
            % 
            %         % saccade during image viewing
            %         if ~((fixationstats.saccadetimes(1, i_saccade) > ...
            %                 -event_info.trial_event_begin * event_info.sfreq) && ( ...
            %                 fixationstats.saccadetimes(1, i_saccade) < ...
            %                 trial_offset_sample - event_info.trial_event_end * event_info.sfreq))
            %             epochs{epoch_counter} = laplacian_referenced_signal(saccade_onset_sample:(saccade_onset_sample - ...
            %                        event_info.trial_event_begin * event_info.sfreq + ...
            %                        event_info.trial_event_end * event_info.sfreq), :);
            %             epoch_counter = epoch_counter+ 1;
            %         end
            %     end
            % end
            % 
            % for i_epoch = 1:numel(epochs)
            %     for ch = 1:size(epochs{1}, 2)
            %         ifft_wavelet_filtered = zeros(length(freqs), size(epochs{1}, 1));
            % 
            %         reconstructed_signal = zeros(size(epochs{i_epoch}(:, ch)));
            %         for f_idx = 1:length(freqs)
            %             f = freqs(f_idx);
            % 
            %             [wavelet, wavelet_length] = morlet_wavelet_cohen(f, event_info.sfreq, num_cycle_wavelet);
            %             wavelet = wavelet / sum(abs(wavelet));
            %             padded_signal = [zeros(wavelet_length, 1); epochs{i_epoch}(:, ch)]; 
            %             lenConv = length(padded_signal) + wavelet_length;
            %             fft_signal = fft(padded_signal, lenConv);
            %             fft_wavelet = fft(wavelet, lenConv);
            %             conv_result = ifft(fft_signal .* fft_wavelet');
            %             % Extract the relevant portion of the convolution result
            %             start_idx = floor(wavelet_length / 2) + 1;
            %             end_idx = start_idx + length(epochs{i_epoch}(:, ch)) - 1;
            % 
            %             ifft_wavelet_filtered(f_idx, :) = conv_result(start_idx:end_idx);
            %             % validate
            %             % figure; plot(epochs{i_epoch}(:, ch)); hold on; plot(real(ifft_wavelet_filtered(f_idx, :)))
            % 
            %             reconstructed_signal = reconstructed_signal + real(ifft_wavelet_filtered(f_idx, :))';
            % 
            %         end
            %         energy_original = sum(epochs{i_epoch}(:, ch).^2);
            %         energy_reconstructed = sum(reconstructed_signal.^2);
            %         scale_factor = sqrt(energy_original / energy_reconstructed);
            %         ifft_wavelet_filtered = ifft_wavelet_filtered * scale_factor;
            % 
            %         filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_epoch%d_channel%d_saccade_onset_control.mat', ...
            %             subj_id, current_task, current_session, i_epoch, ch));
            %         save(filename, 'ifft_wavelet_filtered', '-v7.3');
            %     end
            % end

            % alternative, save one file for each channel:
            % It is painful to have too many files, I will save one file per
        % channel including all epochs
        % Iterate over each channel
        % for ch = 1:size(epochs{1}, 2)
        %     ifft_wavelet_filtered = zeros(numel(epochs), length(freqs), size(epochs{1}, 1));
        %     filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_channel%d_image_onset.mat', ...
        %             subj_id, current_task, current_session, ch));
        %     if ~exist(filename, 'file')
        %         for i_epoch = 1:numel(epochs)
        %             reconstructed_signal = zeros(size(epochs{i_epoch}(:, ch)));
        %             for f_idx = 1:length(freqs)
        %                 f = freqs(f_idx);
        %                 [wavelet, wavelet_length] = morlet_wavelet_cohen(f, event_info.sfreq, num_cycle_wavelet);
        %                 wavelet = wavelet / sum(abs(wavelet));
        %                 padded_signal = [zeros(wavelet_length, 1); epochs{i_epoch}(:, ch)];
        %                 lenConv = length(padded_signal) + wavelet_length;
        %                 fft_signal = fft(padded_signal, lenConv);
        %                 fft_wavelet = fft(wavelet, lenConv);
        %                 conv_result = ifft(fft_signal .* fft_wavelet');
        % 
        %                 start_idx = floor(wavelet_length / 2) + 1;
        %                 end_idx = start_idx + length(epochs{i_epoch}(:, ch)) - 1;
        %                 ifft_wavelet_filtered(i_epoch, f_idx, :) = conv_result(start_idx:end_idx);
        % 
        %                 reconstructed_signal = reconstructed_signal + real(ifft_wavelet_filtered(i_epoch, f_idx, :));
        %             end
        % 
        %             energy_original = sum(epochs{i_epoch}(:, ch).^2);
        %             energy_reconstructed = sum(reconstructed_signal.^2);
        %             scale_factor = sqrt(energy_original / energy_reconstructed);
        %             ifft_wavelet_filtered(i_epoch, :, :) = ifft_wavelet_filtered(i_epoch, :, :) * scale_factor;
        %         end
        %             save(filename, 'ifft_wavelet_filtered', '-v7.3');
        %     else
        %             fprintf('File %s already exists. Skipping.\n', filename);
        % 
        %     end
        % end


        % for ch = 1:size(epochs{1}, 2)
        %         reconstructed_signal = zeros(size(laplacian_referenced_signal, 1), 1);
        %         ifft_wavelet_filtered_ch = zeros(length(freqs), size(laplacian_referenced_signal, 1));
        %         for f_idx = 1:length(freqs)
        %             f = freqs(f_idx);
        %             [wavelet, wavelet_length] = morlet_wavelet_cohen(f, event_info.sfreq, num_cycle_wavelet);
        %             wavelet = wavelet / sum(abs(wavelet));
        %             padded_signal = [zeros(wavelet_length, 1); laplacian_referenced_signal(:, ch)]; 
        %             lenConv = length(padded_signal) + wavelet_length;
        %             fft_signal = fft(padded_signal, lenConv);
        %             fft_wavelet = fft(wavelet, lenConv);
        %             conv_result = ifft(fft_signal .* fft_wavelet');
        %             % Extract the relevant portion of the convolution result
        %             start_idx = floor(wavelet_length / 2) + 1;
        %             end_idx = start_idx + size(laplacian_referenced_signal, 1) - 1;
        % 
        %             ifft_wavelet_filtered_ch(f_idx, :) = conv_result(start_idx:end_idx);
        %             % validate
        %             % figure; plot(epochs{i_epoch}(:, ch)); hold on; plot(real(ifft_wavelet_filtered(f_idx, :)))
        %             reconstructed_signal = reconstructed_signal + real(ifft_wavelet_filtered_ch(f_idx, :))';
        %         end
        %         energy_original = sum(laplacian_referenced_signal(:, ch).^2);
        %         energy_reconstructed = sum(reconstructed_signal.^2);
        %         scale_factor = sqrt(energy_original / energy_reconstructed);
        %         ifft_wavelet_filtered_ch = ifft_wavelet_filtered_ch * scale_factor;
        % 
        %         for i_epoch = 1:length(saccade_onsets_image_viewing)
        %             saccade_onset_sample = saccade_onsets_image_viewing(i_epoch);
        %             ifft_wavelet_filtered = ifft_wavelet_filtered_ch(:, saccade_onset_sample:(saccade_onset_sample - ...
        %                            event_info.trial_event_begin * event_info.sfreq + ...
        %                            event_info.trial_event_end * event_info.sfreq));
        %             if ismember(saccade_onset_sample, saccade_onsets_image_viewing_not_image_onset)
        %                 filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_epoch%d_channel%d_saccade_onset_not_image_onset.mat', ...
        %                     subj_id, current_task, current_session, i_epoch, ch));
        %             else
        %                 filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_epoch%d_channel%d_saccade_onset.mat', ...
        %                     subj_id, current_task, current_session, i_epoch, ch));
        %             end
        %             if ~exist(filename, 'file')
        %                 save(filename, 'ifft_wavelet_filtered', '-v7.3');
        %             else
        %                 fprintf('File %s already exists. Skip\n', filename);
        %             end
        %         end
        % 
        %         for i_epoch = 1:length(saccade_onsets_control)
        %             saccade_onset_sample = saccade_onsets_control(i_epoch);
        %             ifft_wavelet_filtered = ifft_wavelet_filtered_ch(:, saccade_onset_sample:(saccade_onset_sample - ...
        %                            event_info.trial_event_begin * event_info.sfreq + ...
        %                            event_info.trial_event_end * event_info.sfreq));
        %             filename = fullfile(temp_save_tf_path, sprintf('%s_%s_session%d_epoch%d_channel%d_saccade_onset_control.mat', ...
        %                 subj_id, current_task, current_session, i_epoch, ch));
        %             if ~exist(filename, 'file')
        %                 save(filename, 'ifft_wavelet_filtered', '-v7.3');
        %             else
        %                 fprintf('File %s already exists. Skip\n', filename);
        %             end
        %         end