%% Final processing step
% Based on Zaixu's code "Step_8th_Visualize_Workbench_Atlas"

clear
clc

%% ================= USER SETTINGS =================

% Software paths
WorkbenchDir = 'C:\workbench\bin_windows64';
SPMDir       = 'C:\workbench\spm12';

% Input/output folder
VisualizeFolder = ['C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/Analyses/Results/TCA_PFN_loadings_gam_predictions/Normed_delta_maps/'];

% Map naming components
Frames     = {'stable','delta_1st','delta_99th'};
MapSuffix = 'scaling_pred_topo';

% Number of PFNs/networks
NumNetworks = 17;

% Delete intermediate GIFTI files after CIFTI creation?
DeleteIntermediateFiles = false;

%% =================================================


%% Initialize software

addpath(WorkbenchDir);
addpath(SPMDir);

spm('Defaults','fmri');
spm_jobman('initcfg');

% Full path to wb_command
WB_Command = fullfile(WorkbenchDir, 'wb_command.exe');


%% Create atlas label/color file

ColorInfo_Atlas = fullfile(VisualizeFolder, 'name_Atlas.txt');

if exist(ColorInfo_Atlas, 'file')
    delete(ColorInfo_Atlas);
end

SystemName = { ...
    'Net 1: DM', ...
    'Net 2: SM', ...
    'Net 3: FP', ...
    'Net 4: SM', ...
    'Net 5: DA', ...
    'Net 6: VS', ...
    'Net 7: VA', ...
    'Net 8: DM', ...
    'Net 9: VA', ...
    'Net 10: VS', ...
    'Net 11: SM', ...
    'Net 12: DM', ...
    'Net 13: SM', ...
    'Net 14: DA', ...
    'Net 15: FP', ...
    'Net 16: AU', ...
    'Net 17: FP'};

ColorPlate = { ...
    '242 139 168', ...
    '173 216 230', ...
    '244 197 115', ...
    '73 143 191', ...
    '65 171 93', ...
    '137 63 153', ...
    '217 117 242', ...
    '226 57 93', ...
    '206 28 249', ...
    '102 5 122', ...
    '33 113 181', ...
    '170 12 61', ...
    '7 69 132', ...
    '0 109 44', ...
    '216 144 72', ...
    '78 49 168', ...
    '204 109 14'};

if numel(SystemName) ~= NumNetworks || numel(ColorPlate) ~= NumNetworks
    error('SystemName and ColorPlate must each contain %d entries.', ...
        NumNetworks);
end

fid = fopen(ColorInfo_Atlas, 'w');

if fid == -1
    error('Could not create label file: %s', ColorInfo_Atlas);
end

for i = 1:NumNetworks

    fprintf(fid, '%s\n', SystemName{i});
    fprintf(fid, '%d %s 1\n', i, ColorPlate{i});

end

fclose(fid);


%% Process each map

for m = 1:numel(Frames)

    frame = Frames{m};

    fprintf('\n========================================\n');
    fprintf('Processing %s map\n', frame);
    fprintf('========================================\n');

    %% Construct map name automatically

    Map = sprintf('%s_%s', ...
        frame, MapSuffix);


    %% Construct input filenames

    L_mat = fullfile(VisualizeFolder, ...
        [Map '_L.mat']);

    R_mat = fullfile(VisualizeFolder, ...
        [Map '_R.mat']);


    %% Check input files

    if ~exist(L_mat, 'file')
        error('Left hemisphere file not found:\n%s', L_mat);
    end

    if ~exist(R_mat, 'file')
        error('Right hemisphere file not found:\n%s', R_mat);
    end


    %% Load hemisphere data

    PFN_L = load(L_mat);
    PFN_R = load(R_mat);


    %% Construct MAT variable names automatically

    lh_field = [frame '_lh'];
    rh_field = [frame '_rh'];


    %% Check that expected variables exist

    if ~isfield(PFN_L, 'x')
        error('Variable "x" not found in:\n%s', L_mat);
    end

    if ~isfield(PFN_R, 'x')
        error('Variable "x" not found in:\n%s', R_mat);
    end

    if ~isfield(PFN_L.x, lh_field)
        error('Field "%s" not found in:\n%s', ...
            lh_field, L_mat);
    end

    if ~isfield(PFN_R.x, rh_field)
        error('Field "%s" not found in:\n%s', ...
            rh_field, R_mat);
    end


    %% Extract atlas labels

    sbj_AtlasLabel_lh = PFN_L.x.(lh_field);
    sbj_AtlasLabel_rh = PFN_R.x.(rh_field);


    %% Left hemisphere

    V_lh = gifti;
    V_lh.cdata = sbj_AtlasLabel_lh;

    V_lh_File = fullfile(VisualizeFolder, ...
        [Map '_lh.func.gii']);

    save(V_lh, V_lh_File);

    V_lh_Label_File = fullfile(VisualizeFolder, ...
        [Map '_lh_AtlasLabel.label.gii']);

    cmd = sprintf( ...
        '"%s" -metric-label-import "%s" "%s" "%s"', ...
        WB_Command, ...
        V_lh_File, ...
        ColorInfo_Atlas, ...
        V_lh_Label_File);

    status = system(cmd);

    if status ~= 0
        error('Workbench failed for left hemisphere: %s', Map);
    end


    %% Right hemisphere

    V_rh = gifti;
    V_rh.cdata = sbj_AtlasLabel_rh;

    V_rh_File = fullfile(VisualizeFolder, ...
        [Map '_rh.func.gii']);

    save(V_rh, V_rh_File);

    V_rh_Label_File = fullfile(VisualizeFolder, ...
        [Map '_rh_AtlasLabel.label.gii']);

    cmd = sprintf( ...
        '"%s" -metric-label-import "%s" "%s" "%s"', ...
        WB_Command, ...
        V_rh_File, ...
        ColorInfo_Atlas, ...
        V_rh_Label_File);

    status = system(cmd);

    if status ~= 0
        error('Workbench failed for right hemisphere: %s', Map);
    end


    %% Create final CIFTI label

    CIFTI_File = fullfile(VisualizeFolder, ...
        [Map '_AtlasLabel.dlabel.nii']);

    cmd = sprintf( ...
        '"%s" -cifti-create-label "%s" -left-label "%s" -right-label "%s"', ...
        WB_Command, ...
        CIFTI_File, ...
        V_lh_Label_File, ...
        V_rh_Label_File);

    status = system(cmd);

    if status ~= 0
        error('Workbench failed while creating CIFTI: %s', Map);
    end


    %% Optional cleanup

    if DeleteIntermediateFiles

        delete(V_lh_File);
        delete(V_rh_File);
        delete(V_lh_Label_File);
        delete(V_rh_Label_File);

    end


    fprintf('Created:\n%s\n', CIFTI_File);

end


fprintf('\nAll requested maps finished successfully.\n');