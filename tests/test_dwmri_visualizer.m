%% Set environment
addpath(genpath('nifti_utils'));
addpath(genpath('dwmri_visualizer'));
data_path = 'data';

%% Volume test

% Load data
fa_vol = nifti_utils.load_untouch_nii_vol_scaled(fullfile(data_path,'fa.nii.gz'),'double');
mask_vol = nifti_utils.load_untouch_nii_vol_scaled(fullfile(data_path,'mask.nii.gz'),'logical');
xform_RAS = nifti_utils.get_voxel_RAS_xform(fullfile(data_path,'dwmri.nii.gz'));

% Get visualizer
dv = dwmri_visualizer([],...
                      fa_vol, ...
                      mask_vol, ...
                      xform_RAS, ...
                      'vol');
                  
% Make a plot
figure;
dv.plot_slice(26,'axial','slice',[],subplot(1,3,1));
axis image
dv.plot_slice(47,'coronal','slice',[],subplot(1,3,2));
axis image
dv.plot_slice(47,'sagittal','slice',[],subplot(1,3,3));
axis image

%% Outline test

fa_vol = nifti_utils.load_untouch_nii_vol_scaled(fullfile(data_path,'fa.nii.gz'),'double');
mask_vol = nifti_utils.load_untouch_nii_vol_scaled(fullfile(data_path,'mask.nii.gz'),'logical');
mask2_vol = nifti_utils.load_untouch_nii_vol_scaled(fullfile(data_path,'mask2.nii.gz'),'logical');
xform_RAS = nifti_utils.get_voxel_RAS_xform(fullfile(data_path,'dwmri.nii.gz'));

% Get visualizer
dv = dwmri_visualizer({mask_vol,mask2_vol}, ...
                      fa_vol, ...
                      mask_vol, ...
                      xform_RAS, ...
                      'outlines', ...
                      {'r','g';...
                       1,1});

% Make a plot
figure;
dv.plot_slice(26,'axial','boundingbox',[],subplot(1,3,1));
axis image
dv.plot_slice(47,'coronal','boundingbox',[],subplot(1,3,2));
axis image
dv.plot_slice(47,'sagittal','boundingbox',[],subplot(1,3,3));
axis image

%% Colorized FA test

% Load data
v1_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(data_path,'v1.nii.gz'),'double');
fa_vol = nifti_utils.load_untouch_nii_vol_scaled(fullfile(data_path,'fa.nii.gz'),'double');
xform_RAS = nifti_utils.get_voxel_RAS_xform(fullfile(data_path,'dwmri.nii.gz'));

% Get visualizer
dv = dwmri_visualizer({v1_vol,fa_vol}, ...
                      [], ...
                      fa_vol, ...
                      xform_RAS, ...
                      'colorized_FA', ...
                      {'w',1});

% Make a plot
dv.plot_slice(26,'axial','bottom-centroid',1);
axis image

%% Directions test - bedpostX

% Load data
dyads1_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(data_path,'dyads1.nii.gz'),'double');
dyads2_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(data_path,'dyads2.nii.gz'),'double');
mean_f1samples_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(data_path,'mean_f1samples.nii.gz'),'double');
mean_f2samples_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(data_path,'mean_f2samples.nii.gz'),'double');
mean_fsumsamples_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(data_path,'mean_fsumsamples.nii.gz'),'double');
mask_vol = nifti_utils.load_untouch_nii_vol_scaled(fullfile(data_path,'mask.nii.gz'),'logical');
xform_RAS = nifti_utils.get_voxel_RAS_xform(fullfile(data_path,'dwmri.nii.gz'));

% Get visualizer
dv = dwmri_visualizer({dyads1_vol,dyads2_vol; ...
                       ones(size(mask_vol)),mean_f2samples_vol}, ...
                      mean_fsumsamples_vol, ... 
                      mask_vol, ...
                      xform_RAS, ...
                      'directions', ...
                      {'r','b'; ...
                        0.25,0.25; ...
                        1,0.05});

% Make a plot
dv.plot_slice(47,'coronal','boundingbox');
axis image

%% sh_coefs test - qball

% Load data
sh_coefs_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(data_path,'qball_sh_coefs.nii.gz'),'double');
fa_vol = nifti_utils.load_untouch_nii_vol_scaled(fullfile(data_path,'fa.nii.gz'),'double');
xform_RAS = nifti_utils.get_voxel_RAS_xform(fullfile(data_path,'dwmri.nii.gz'));

% Get visualizer
dv = dwmri_visualizer(sh_coefs_vol, ...
                      fa_vol, ...
                      fa_vol*1.5, ...
                      xform_RAS, ...
                      'sh_coefs', ...
                      {6,20,true});

% Make a plot
dv.plot_slice(47,'coronal','top-centroid',1);
axis image
light('Position', [5, 5, 5], 'Style', 'infinite')

%% sh_coefs test - MT_SCSD

% Load data
dwmri_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(data_path,'dwmri.nii.gz'),'double');
bvals = dlmread(fullfile(data_path,'dwmri.bval'));
sh_coefs_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(data_path,'mt_scsd_sh_coefs.nii.gz'),'double');
fa_vol = nifti_utils.load_untouch_nii_vol_scaled(fullfile(data_path,'fa.nii.gz'),'double');
xform_RAS = nifti_utils.get_voxel_RAS_xform(fullfile(data_path,'dwmri.nii.gz'));

% Get visualizer
dv = dwmri_visualizer(sh_coefs_vol, ...
                      mean(dwmri_vol(:,:,:,bvals == max(bvals)),4)-0.25, ...
                      fa_vol*2, ...
                      xform_RAS, ...
                      'sh_coefs', ...
                      {8,60,true});

% Make a plot
dv.plot_slice(47,'coronal','top-centroid',1);
axis image
light('Position', [5, 5, 5], 'Style', 'infinite')

%% DTI 

% Load data
dwmri_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(data_path,'dwmri.nii.gz'),'double');
bvals = dlmread(fullfile(data_path,'dwmri.bval'));
dt_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(data_path,'dt.nii.gz'),'double');
fa_vol = nifti_utils.load_untouch_nii_vol_scaled(fullfile(data_path,'fa.nii.gz'),'double');
xform_RAS = nifti_utils.get_voxel_RAS_xform(fullfile(data_path,'dwmri.nii.gz'));

% Get visualizer
dv = dwmri_visualizer({dt_vol(:,:,:,3:end)}, ...    
                      mean(dwmri_vol(:,:,:,bvals == max(bvals)),4)-0.25, ...
                      fa_vol*1.75, ...
                      xform_RAS, ...
                      'DT', ...
                      {40,1});

% Make a plot
dv.plot_slice(47,'coronal','top-centroid',1);
axis image

%% DTI  - colorized by v1 direction

% Load data
dwmri_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(data_path,'dwmri.nii.gz'),'double');
bvals = dlmread(fullfile(data_path,'dwmri.bval'));
dt_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(data_path,'dt.nii.gz'),'double');
v1_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(data_path,'v1.nii.gz'),'double');
fa_vol = nifti_utils.load_untouch_nii_vol_scaled(fullfile(data_path,'fa.nii.gz'),'double');
xform_RAS = nifti_utils.get_voxel_RAS_xform(fullfile(data_path,'dwmri.nii.gz'));

% Get visualizer
dv = dwmri_visualizer({dt_vol(:,:,:,3:end),v1_vol}, ...    
                      mean(dwmri_vol(:,:,:,bvals == max(bvals)),4)-0.25, ...
                      fa_vol*1.75, ...
                      xform_RAS, ...
                      'DT', ...
                      {40,1});

% Make a plot
dv.plot_slice(26,'axial','boundingbox');
axis image
light('Position', [5, 5, 5], 'Style', 'infinite')

%% FOD test - SD

% Load data
fod_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(data_path,'sd_fod.nii.gz'),'double');
fa_vol = nifti_utils.load_untouch_nii_vol_scaled(fullfile(data_path,'fa.nii.gz'),'double');
xform_RAS = nifti_utils.get_voxel_RAS_xform(fullfile(data_path,'dwmri.nii.gz'));

% Get visualizer
dv = dwmri_visualizer(fod_vol, ...
                      fa_vol, ...
                      fa_vol*1.5, ...
                      xform_RAS, ...
                      'OD', ...
                      false);

% Make a plot
dv.plot_slice(47,'coronal','boundingbox',1);
axis image

%% PASMRI test - PASMRI

% Load data
pas_vol = nifti_utils.load_untouch_nii4D_vol_scaled(fullfile(data_path,'pas.nii.gz'),'double');
fa_vol = nifti_utils.load_untouch_nii_vol_scaled(fullfile(data_path,'fa.nii.gz'),'double');
xform_RAS = nifti_utils.get_voxel_RAS_xform(fullfile(data_path,'dwmri.nii.gz'));

% Get visualizer
dv = dwmri_visualizer(pas_vol, ...
                      fa_vol, ...
                      fa_vol*1.5, ...
                      xform_RAS, ...
                      'PAS', ...
                      {1.4,16,'~/camino',120,false});

% Make a plot
dv.plot_slice(47,'coronal','top-centroid',1);
axis image
