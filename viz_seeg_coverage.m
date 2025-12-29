
output_csv_file = fullfile('data', 'norm_contact_data_table.csv');
norm_contact_data_table = readtable(output_csv_file);
%% exclude bad chs, and micro contacts
bad_ch_map = containers.Map();
BJH021_tasks = struct();
BJH021_tasks.BLAES_study = {'A3', 'A4', 'E8', 'O1'};
BJH021_tasks.BLAES_study_2 = {'A4', 'A5', 'A6', 'E8', 'O1'};
bad_ch_map('BJH021') = BJH021_tasks;

BJH024_tasks = struct();
BJH024_tasks.BLAES_study = {'FL8', 'ER2', 'GR1', 'GR2', 'HL14', 'HR3', 'chan1', 'chan2', 'chan3', 'chan4', ...
    'chan5', 'chan6', 'chan7', 'chan8', 'chan9', 'chan10', 'chan11', 'chan12', ...
    'chan13', 'chan14', 'chan15', 'chan16'};
BJH024_tasks.BLAES_study_2 = {'FL8', 'GL3', 'GL4'};
bad_ch_map('BJH024') = BJH024_tasks;

BJH025_tasks = struct();
BJH025_tasks.BLAES_study_twosource = {'D8','GR4','GR5','J7'};
bad_ch_map('BJH025') = BJH025_tasks;

BJH026_tasks = struct();
BJH026_tasks.BLAES_study = {'DL16', 'GL1', 'GL2', 'GR14', 'chan16'};
BJH026_tasks.BLAES_study_2 = {'DL16', 'GL3', 'GL4', 'GR5', 'GR4', 'chan16', 'chan15', 'chan14'};
bad_ch_map('BJH026') = BJH026_tasks;

BJH027_tasks = struct();
BJH027_tasks.BLAES_study = {'GL3', 'GL4', 'KL1'};
BJH027_tasks.BLAES_study_2 = {'GL10', 'GR4', 'GR5', 'KL1'};
BJH027_tasks.BLAES_study_3 = {'GL4', 'GL5', 'GL10','KL1', 'chan11'};
bad_ch_map('BJH027') = BJH027_tasks;

BJH028_tasks = struct();
BJH028_tasks.BLAES_study = {'A_left1', 'A_left2', 'B_left1', 'B_left2', 'B_left3', 
                 'B_left4', 'B_left5', 'B_left6',  'I_left7', 'K_left1'};
bad_ch_map('BJH028') = BJH028_tasks;

BJH029_tasks = struct();
BJH029_tasks.BLAES_study = {'A_left4', 'A_left5'};
bad_ch_map('BJH029') = BJH029_tasks;

BJH030_tasks = struct();
BJH030_tasks.BLAES_study = {'DR13', 'KBL3', 'KBL4'};
bad_ch_map('BJH030') = BJH030_tasks;

subj_keys = keys(bad_ch_map);          
common_bad_map = containers.Map();     
for k = 1:length(subj_keys)
    subj = subj_keys{k};
    tasks = bad_ch_map(subj);
    task_fields = fieldnames(tasks);
    
    common_bad = tasks.(task_fields{1});
    
    for j = 2:length(task_fields)
        current_task = tasks.(task_fields{j});
        common_bad = intersect(common_bad, current_task);
    end
    common_bad_map(subj) = common_bad;

end

rows_to_exclude = false(height(norm_contact_data_table), 1);

for i = 1:height(norm_contact_data_table)
    curr_subj = norm_contact_data_table.SubjectID{i};
    curr_channel = norm_contact_data_table.ch_name_bci2000_format{i};
    
    if contains(curr_channel, 'chan')
        rows_to_exclude(i) = true;
        continue; 
    end
    if isKey(bad_ch_map, curr_subj)
        if any(strcmp(curr_channel, common_bad_map(curr_subj)))
            rows_to_exclude(i) = true;
        end
    end
end

norm_contact_data_table_bad_ch_excluded = norm_contact_data_table(~rows_to_exclude, :);

% parameters
isovalue_to_get_surface = 0.5;
% first 3D glass brain with axis and origin
% require freesurfer
freesurfer_path = '/Applications/freesurfer/7.4.1';
lut_path = [freesurfer_path, '/luts/FreeSurferColorLUT.txt'];
color_table = load_freesurfer_color_table(lut_path);
% require subject-specific segmentation
% get segmentation vol
segmentation_folder = [data_path, '/', subj_ids{reference_subj}, '/IMAGING/segmentation'];
aparc_aseg_mgz_path = [segmentation_folder, '/mri/aparc+aseg.mgz'];
aseg_mgz_path = [segmentation_folder, '/mri/aseg.mgz'];
% we can also use the mri_read function with freesurfer
segmentation_volume = MRIread(aseg_mgz_path);

% predefine surface that we want to plot
% brain_region2plot = {'Left-Cerebral-Cortex', 'Right-Cerebral-Cortex',
%     'Left-Cerebral-White-Matter', 'Right-Cerebral-White-Matter'};
brain_region2plot = {'Cerebral-Cortex', 'Cerebral-White-Matter'};
% , 'Cerebral-White-Matter'
sides = {'Left', 'Right'};
fig = figure;
% get white matter color
region_color = color_table.RGB(3, :) / 255;
for side_id = 1:length(sides)
    binary_vol = zeros(size(segmentation_volume.vol));
    side = sides{side_id};
    for i_region = 1:length(brain_region2plot)
        region_label = [side '-' brain_region2plot{i_region}];
        region_label_id = find(strcmp(color_table.Label, region_label));
        region_id = color_table.ID(region_label_id);

        binary_vol(segmentation_volume.vol == region_id)=true;
    end
    if sum(binary_vol(:)) > 0
        [col_voxel_space, row_voxel_space, slice_voxel_space] = ...
            meshgrid(1:size(binary_vol, 2), 1:size(binary_vol, 1), 1:size(binary_vol, 3));
        [faces, vert_voxel_space]=isosurface(col_voxel_space, row_voxel_space, ...
            slice_voxel_space, binary_vol, isovalue_to_get_surface); 
        vert_voxel_space_homogeneous=[vert_voxel_space ones(size(vert_voxel_space, 1), 1)]; 
        vert_ras_space_homogeneous = segmentation_volume.vox2ras1 * ...
        vert_voxel_space_homogeneous';
        vert_ras_space = vert_ras_space_homogeneous(1:3, :)';
        patch('Faces', faces, 'Vertices', vert_ras_space, ...
            'FaceColor', region_color, 'FaceAlpha', 0.4, ...
             'EdgeAlpha', 0.01);

    end
end
axis equal;
xlabel('X (left <--> right) [mm]');
ylabel('Y (posterior <--> anterior) [mm]');
zlabel('Z (inferior <--> superior) [mm]');
title('3D Brain Visualization with Contacts');
hold on

% Plot the contact locations (the )
scatter3(norm_contact_data_table_bad_ch_excluded.Norm_X, norm_contact_data_table_bad_ch_excluded.Norm_Y, norm_contact_data_table_bad_ch_excluded.Norm_Z, ...
    50, 'r', 'filled');
% create 3D video

%% 2d 
fig = figure;
set(fig, 'Position', [100, 100, 1000, 800]);  % [left, bottom, width, height]
set(fig, 'Renderer', 'opengl');
% get white matter color
region_color = color_table.RGB(3, :) / 255;

for side_id = 1:length(sides)
    binary_vol = zeros(size(segmentation_volume.vol));
    side = sides{side_id};
    for i_region = 1:length(brain_region2plot)
        region_label = [side '-' brain_region2plot{i_region}];
        region_label_id = find(strcmp(color_table.Label, region_label));
        region_id = color_table.ID(region_label_id);
        binary_vol(segmentation_volume.vol == region_id)=true;
    end
    
    if sum(binary_vol(:)) > 0
        binary_vol = smooth3(binary_vol, 'gaussian', 5); % reduce imaging noise
        binary_vol = smooth3(binary_vol, 'gaussian', 5);

        
        [col_voxel_space, row_voxel_space, slice_voxel_space] = ...
            meshgrid(1:size(binary_vol, 2), 1:size(binary_vol, 1), 1:size(binary_vol, 3));
        [faces, vert_voxel_space] = isosurface(col_voxel_space, row_voxel_space, ...
            slice_voxel_space, binary_vol, isovalue_to_get_surface); 
        
        vert_voxel_space_homogeneous = [vert_voxel_space ones(size(vert_voxel_space, 1), 1)];
        vert_ras_space_homogeneous = segmentation_volume.vox2ras1 * vert_voxel_space_homogeneous';
        vert_ras_space = vert_ras_space_homogeneous(1:3, :)';
        
        z_values = vert_ras_space(:, 3);
        min_z = min(z_values);
        max_z = max(z_values);
        normalized_z = (z_values - min_z) / (max_z  - min_z);  % Normalize Z between 0 and 1
        alpha_vals = exp(2 * (1 - normalized_z));       
        alpha_vals = 0.5 * (alpha_vals - min(alpha_vals)) / (max(alpha_vals) - min(alpha_vals));
  
        p = patch('Faces', faces, 'Vertices', vert_ras_space, ...
                  'FaceColor', region_color, 'FaceVertexAlphaData', alpha_vals, 'FaceAlpha', 'interp','AlphaDataMapping', 'none',...
                  'EdgeColor', 'none');  % No edges ('EdgeColor', 'none')


        lighting gouraud;  % Gouraud lighting for smooth lighting effects
        % material([ka kd ks n sc]): Diffuse Reflection (kd): Provides
        % directional lighting, which highlights surface details and
        % contours based on the direction of the light source.
        % 
        material([0.2 0.2 0.3]);  
        light('Position', [1, 20, -150], 'Style', 'infinite'); 
    end
end

axis equal;
xlabel('X (left <--> right) [mm]');
ylabel('Y (posterior <--> anterior) [mm]');
zlabel('Z (inferior <--> superior) [mm]');
title('3D Brain Visualization with Contacts');
hold on
min_z = min(norm_contact_data_table_bad_ch_excluded.Norm_Z);
max_z = max(norm_contact_data_table_bad_ch_excluded.Norm_Z);
normalized_z = (norm_contact_data_table_bad_ch_excluded.Norm_Z - min_z) / (max_z - min_z);  

alpha_vals = 0.5 * exp(1 - normalized_z);   

scatter3(norm_contact_data_table_bad_ch_excluded.Norm_X, norm_contact_data_table_bad_ch_excluded.Norm_Y, norm_contact_data_table_bad_ch_excluded.Norm_Z, ...
    30, 'r', 'filled', 'MarkerFaceAlpha', 'flat', 'MarkerEdgeAlpha', 'flat', ...
    'AlphaData', alpha_vals);  % Use 'flat' to apply different alpha values

view([0 0.2 -1]);  % Bottom-up view

