%% Set environment
addpath(genpath('nifti_utils'));
addpath(genpath('dwmri_visualizer'));
data_path = 'data';
       
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

% Make the plot
dwmri_hdr = load_untouch_header_only(fullfile(data_path,'dwmri.nii.gz'));
f = dwmri_axial_coronal_plot(dv, ...
                             {'MT_SCSD', ...
                             ['   B-values: ' num2str(unique(bvals)) '; lmax: 8; mrtrix path: ~/mrtrix3'], ...
                             ['   Gradient Directions: ' num2str(length(find(bvals ~= 0)))], ...
                             ['   Slice Dimensions: ' num2str(dwmri_hdr.dime.dim(2)) ',' num2str(dwmri_hdr.dime.dim(3))], ...
                             ['   Slices: ' num2str(dwmri_hdr.dime.dim(4))], ...
                             ['   Voxel Resolution: ' num2str(round(dwmri_hdr.dime.pixdim(2),1)) ' x ' num2str(round(dwmri_hdr.dime.pixdim(3),1)) ' x ' num2str(round(dwmri_hdr.dime.pixdim(4),1))], ...
                             'Justin Blaber', ...
                             'justin.blaber@vanderbilt.edu', ...
                             'Multi-tissue constrained spherical deconvolution for improved analysis of multi-shell diffusion MRI data (Tournier et al 2014)', ...
                             'N/A', ...
                             'MT_SCSD_v1_1_0'});
                     
% Apply light to axes
all_axes = findall(f,'type','axes');
for i = 1:length(all_axes)
    light(all_axes(i), 'Position', [5, 5, 5], 'Style', 'infinite');
end
