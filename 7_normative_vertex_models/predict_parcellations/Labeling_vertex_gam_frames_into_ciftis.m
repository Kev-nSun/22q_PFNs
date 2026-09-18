% based on Zaixu's code "Step_8th_Visualize_Workbench_Atlas"
%% Initialize SPM and load in files
clear

% Add paths to Connectome Workbench and SPM (if needed)
addpath('C:\workbench\bin_windows64'); % Connectome Workbench binaries
addpath('C:\workbench\spm12');         % SPM (if needed)

% Initialize SPM (if needed)
spm('Defaults','fmri');
spm_jobman('initcfg');

% ---- Base paths ----
MatFolder = 'C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/TCA_PFN_loadings_gam_predictions/091426_normed_gam_hard_parcellation_frames';
VisualizeRoot = 'C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/TCA_PFN_loadings_gam_predictions/091426_VisualizeFolder';

% ---- Network naming + colors (same each frame) ----
SystemName = {'Net 1: DM', 'Net 2: SM', 'Net 3: FP', 'Net 4: SM', 'Net 5: DA', ...
              'Net 6: VS', 'Net 7: VA', 'Net 8: DM', 'Net 9: VA', 'Net 10: VS', 'Net 11: SM', ...
              'Net 12: DM', 'Net 13: SM', 'Net 14: DA', 'Net 15: FP', 'Net 16: AU', 'Net 17: FP'};
ColorPlate = {'242 139 168', '173 216 230', '244 197 115', '73 143 191', ...
              '65 171 93', '137 63 153', '217 117 242', '226 57 93', ...
              '206 28 249', '102 5 122', '33 113 181', '170 12 61', ...
              '7 69 132', '0 109 44', '216 144 72', '78 49 168', '204 109 14'};
          
%% ---- Loop frames ----
    for frame = 1:100

        Map = sprintf('Frame_%03d', frame);
        VisualizeFolder = fullfile(VisualizeRoot, Map);
        if ~exist(VisualizeFolder, 'dir')
            mkdir(VisualizeFolder);
        end

        % Load L/R mats for this frame
        PFN_L = load(fullfile(MatFolder, sprintf('Frame_%03d_PFN_GAM_hard_parcellation_L.mat', frame)));
        PFN_R = load(fullfile(MatFolder, sprintf('Frame_%03d_PFN_GAM_hard_parcellation_R.mat', frame)));

        sbj_AtlasLabel_lh = PFN_L.hard_parcellation;
        sbj_AtlasLabel_rh = PFN_R.hard_parcellation;

        %% PFN ID labeling
        ColorInfo_Atlas = fullfile(VisualizeFolder, 'name_Atlas.txt');
        if exist(ColorInfo_Atlas, 'file'); delete(ColorInfo_Atlas); end

        for net = 1:17
            system(['echo ' SystemName{net} ' >> "' ColorInfo_Atlas '"']);
            system(['echo ' num2str(net) ' ' ColorPlate{net} ' 1 >> "' ColorInfo_Atlas '"']);
        end

        % ---- Left hemi ----
        V_lh = gifti;
        V_lh.cdata = sbj_AtlasLabel_lh;
        V_lh_File = fullfile(VisualizeFolder, [Map '_lh.func.gii']);
        save(V_lh, V_lh_File);
        pause(1);

        V_lh_Label_File = fullfile(VisualizeFolder, [Map '_lh_AtlasLabel.label.gii']);
        cmd = ['wb_command -metric-label-import "' V_lh_File '" "' ColorInfo_Atlas '" "' V_lh_Label_File '"'];
        system(cmd);

        % ---- Right hemi ----
        V_rh = gifti;
        V_rh.cdata = sbj_AtlasLabel_rh;
        V_rh_File = fullfile(VisualizeFolder, [Map '_rh.func.gii']);
        save(V_rh, V_rh_File);
        pause(1);

        V_rh_Label_File = fullfile(VisualizeFolder, [Map '_rh_AtlasLabel.label.gii']);
        cmd = ['wb_command -metric-label-import "' V_rh_File '" "' ColorInfo_Atlas '" "' V_rh_Label_File '"'];
        system(cmd);

        % ---- Convert into cifti file ----
        DlabelOut = fullfile(VisualizeFolder, [Map '_AtlasLabel.dlabel.nii']);
        cmd = ['wb_command -cifti-create-label "' DlabelOut '" -left-label "' V_lh_Label_File '" -right-label "' V_rh_Label_File '"'];
        system(cmd);
        pause(1);

        % Optional cleanup (uncomment if you want to remove intermediates):
        % delete(V_lh_File); delete(V_rh_File); delete(V_lh_Label_File); delete(V_rh_Label_File);

    end