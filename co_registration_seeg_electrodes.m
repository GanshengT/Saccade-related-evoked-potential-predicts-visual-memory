%% Description
% This script project patient MRI to the reference space, require
% freesurfer, imaging processing script from https://github.com/GanshengT/intracranial_contact_loc

% addpath(genpath(['imaging_processing_scripts']));
data_path = '/Users/ganshengtan/Library/CloudStorage/Box-Box/Washu/projects/BLEAS/theta_saccade/data';

% subj info
subj_ids = {'BJH025', 'BJH026', 'BJH027', 'BJH024', 'BJH021', 'BJH028', 'BJH029', 'BJH030'};
MRI_names = {'SAG_T1_MPRAGE', 'SAG_T1_MPRAGE', 'SAG_T1_MPRAGE','x3D_T1_MS_P',...
    'SAG_T1_MPRAGE', 'SAG_T1_MPR', 'x3D_T1_MS_P', 'x3D_T1_MS_P'};

reference_subj = 1;
% we will use this subject as fixed image
for i_subj_id = 1:length(subj_ids)
    subj_id = subj_ids{i_subj_id};
    imaging_processing_temp_folder = [data_path, '/', subj_id, '/imaging_process'];
    % we will only need MRI and electrode coordinates
    % we will output a 3D MNI brain, with RAS axis
    % Then we will project patient MRI and contact coordinates to MNI
    % space.
    % we will output three projection plot with dots representing coverage,
    % further coverage is smaller
    
    % Ensure subfolder for each subject exists
    subj_folder = [coverage_viz_path, '/', subj_id];
    if ~exist(subj_folder, 'dir')
        mkdir(subj_folder);
    end

    mri_file_name = fullfile(imaging_processing_temp_folder, [MRI_names{i_subj_id} '_converted_from_orig_mgz.nii']);
    if i_subj_id == reference_subj
        mri_target_file = fullfile(subj_folder, 'mri_reference.nii');  
    else
        mri_target_file = fullfile(subj_folder, 'mri_patient_space.nii');
    end
    copyfile(mri_file_name, mri_target_file);
end

    % copy the MRI to coverage folder for processing
    % if the subj is reference subject, we rename it the ref_mri

    %% processing
    % set up fixed_image
fixed_image = fullfile(coverage_viz_path, subj_ids{reference_subj}, ...
    'mri_reference.nii');
% need to rename path
norm_script_path = ['imaging_processing_scripts/GT_helper/normalize_to_reference_mri.sh'];
ants_path = '/Users/ganshengtan/Library/CloudStorage/Box-Box/Washu/projects/BLEAS/theta_saccade/install/bin';
setenv('PATH', [ants_path, ':', getenv('PATH')]);
get_norm_image_script_path = ['imaging_processing_scripts/GT_helper/get_moving_image_at_reference_space.sh'];
chmod_cmd = sprintf('chmod +x %s', norm_script_path);
system(chmod_cmd);
chmod_cmd = sprintf('chmod +x %s', get_norm_image_script_path);
% require ANT
antsRegistrationSyN_path = ['bin/antsRegistration'];
chmod_cmd = sprintf('chmod +x %s', antsRegistrationSyN_path);
system(chmod_cmd);
antsRegistrationSyN_path = ['bin/PrintHeader'];
chmod_cmd = sprintf('chmod +x %s', antsRegistrationSyN_path);
system(chmod_cmd);
antsRegistrationSyN_path = ['bin/antsRegistrationSyN.sh'];
chmod_cmd = sprintf('chmod +x %s', antsRegistrationSyN_path);
system(chmod_cmd);
antsApplyTransforms_path = ['bin/antsApplyTransforms'];
chmod_cmd = sprintf('chmod +x %s', antsApplyTransforms_path);
system(chmod_cmd);

norm_contact_data_table = [];


for i_subj_id = 1:length(subj_ids)
    imaging_processing_temp_folder = [data_path, '/', subj_ids{i_subj_id}, '/imaging_process'];

    if i_subj_id == reference_subj
        % get contact coordinates
        if exist(fullfile( ...
            imaging_processing_temp_folder,...
            'contact_data_table_matched_bci2000_signal.mat'), 'file') == 2
            contact_data_table = load(fullfile( ...
                imaging_processing_temp_folder,...
                'contact_data_table_matched_bci2000_signal.mat'));
            contact_data_table = contact_data_table.contact_info_seeg_matched_bci2000_signal;
        else
            contact_data_table = load(fullfile( ...
                imaging_processing_temp_folder,...
                'contact_data_table.mat'));
            contact_data_table.ch_name_bci2000_format = repmat({'none'}, height(contact_data_table), 1);

        end
        contact_data_table.Norm_X = contact_data_table.X_aligned_to_brightest_voxel;
        contact_data_table.Norm_Y = contact_data_table.Y_aligned_to_brightest_voxel;
        contact_data_table.Norm_Z = contact_data_table.Z_aligned_to_brightest_voxel;
        contact_data_table.SubjectID = repmat(subj_ids{i_subj_id}, size(contact_data_table, 1), 1);
        norm_contact_data_table = [norm_contact_data_table; contact_data_table];

 
    else
        if exist(fullfile( ...
            imaging_processing_temp_folder,...
            'contact_data_table_matched_bci2000_signal.mat'), 'file') == 2
            contact_data_table = load(fullfile( ...
                imaging_processing_temp_folder,...
                'contact_data_table_matched_bci2000_signal.mat'));
            contact_data_table = contact_data_table.contact_info_seeg_matched_bci2000_signal;
        else
            contact_data_table = load(fullfile( ...
                imaging_processing_temp_folder,...
                'contact_data_table.mat'));
            contact_data_table = contact_data_table.contact_data_table;
            % add col corresponding to 
            contact_data_table.ch_name_bci2000_format = repmat({'none'}, height(contact_data_table), 1);

        end
        moving_image = fullfile(coverage_viz_path, subj_ids{i_subj_id}, 'mri_patient_space.nii');
        output_prefix = fullfile(coverage_viz_path, subj_ids{i_subj_id}, 'norm_output_');
    
        if exist([coverage_viz_path, '/', subj_ids{i_subj_id}, ...
            '/norm_output_0GenericAffine.mat'], 'file')
            disp('norm has been done previously, affine matrix found')
            % affine_matrix_struct = load([coverage_viz_path, '/', subj_ids{i_subj_id}, ...
            % '/norm_output_0GenericAffine.mat']);
        else
            cmd = sprintf('%s %s %s %s', norm_script_path, moving_image, fixed_image, output_prefix);
            % comment the follow line if normalization has been run
            [status, result] = system(cmd);
            assert(status == 0, 'Normalization is not successful')
            disp(result);
        end
        % require coverage_viz_path
        affine_matrix_struct = load([coverage_viz_path, '/', subj_ids{i_subj_id}, ...
            '/norm_output_0GenericAffine.mat']);

        % transform contacts
        X = contact_data_table{:, 6};
        Y = contact_data_table{:, 7};
        Z = contact_data_table{:, 8};
        coords = [X, Y, Z, ones(length(X), 1)];  
        
        % doc for understanding 0GenericAffine output, page 6 from http://stnava.github.io/ANTsDoc/
         % ANTS will give us warp and inverse warp, we can use warp and
        % norm_output_0GenericAffine to map the patient space to reference
        % space. inverse warp is used to map from the reference space to
        % patient space
        
        affine_matrix = ea_antsmat2mat(affine_matrix_struct.AffineTransform_double_3_3, ...
            affine_matrix_struct.fixed);
        % 

        transformed_coords = coords * affine_matrix';

        contact_data_table.Norm_X = transformed_coords(:, 1);
        contact_data_table.Norm_Y = transformed_coords(:, 2);
        contact_data_table.Norm_Z = transformed_coords(:, 3);
        contact_data_table.SubjectID = repmat(subj_ids{i_subj_id}, size(contact_data_table, 1), 1);

        % export the norm electrode location and check in freeview
        norm_contacts_coordinates_folder =  fullfile(coverage_viz_path, ...
            subj_ids{i_subj_id}, 'norm_contact_coordinate');
        if ~exist(norm_contacts_coordinates_folder, 'dir')
            mkdir(norm_contacts_coordinates_folder);
        end
        unique_shankIDs = unique(contact_data_table.ShankID);
        for i_shankID = 1:length(unique_shankIDs)
            shankID = unique_shankIDs{i_shankID};
            filtered_table = contact_data_table(strcmp(contact_data_table.ShankID, shankID), :);
            X = filtered_table.Norm_X;
            Y = filtered_table.Norm_Y;
            Z = filtered_table.Norm_Z;
            file_name = sprintf('%s_norm.dat', shankID);
        
            fid = fopen([norm_contacts_coordinates_folder, '/', file_name], 'w');
            fprintf(fid, '\n');
            for i_contact = 1:height(filtered_table)
                fprintf(fid, [num2str([X(i_contact), Y(i_contact), Z(i_contact)]) '\n']);
            end
            fprintf(fid, 'info\n');
            fprintf(fid, 'numpoints %d\n', height(filtered_table));
            fprintf(fid, 'useRealRAS 1\n');
            fclose(fid);
        end

        cmd = sprintf('%s %s %s %s', get_norm_image_script_path, moving_image, ...
            fixed_image, fullfile(coverage_viz_path, subj_ids{i_subj_id}));
        [status, result] = system(cmd);
        disp(result);
        norm_contact_data_table = [norm_contact_data_table; contact_data_table];
    end 
end

% save the table


