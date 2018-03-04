# dwmri_visualizer
Visualization toolbox for DWMRI related processing, including: DTI, spherical harmonics, etc... 

# Installation instructions:
1) Install [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
2) (optional - only if plotting PASMRI) Install [camino](http://camino.cs.ucl.ac.uk/)
3) Download repos and (optional) example data:
```
git clone https://github.com/justinblaber/nifti_utils.git
git clone https://github.com/justinblaber/dwmri_visualizer.git

# Optionally download example data (note this is quite large... sorry)
wget https://justinblaber.org/downloads/github/dwmri_visualizer/data.zip
unzip data.zip
```
3) In MATLAB:
```
>> addpath(genpath('nifti_utils'));
>> addpath(genpath('dwmri_visualizer'));
>> data_path = 'data';
```
Then, run specific cells in each script:
```
>> edit test_dwmri_visualizer.m
>> edit test_dwmri_plot.m
```
Select cells from `test_dwmri_visualizer.m` and `test_dwmri_plot.m`  are explained below:
```
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
```
<p align="center">
  <img width="708" height="611"  src="https://i.imgur.com/i1c0Uwq.png">
</p>

---

```
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
```
<p align="center">
  <img width="686" height="916" src="https://i.imgur.com/HAEeMGR.png">
</p>

---

```
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
```
<p align="center">
  <img width="763" height="761" src="https://i.imgur.com/x6dPglj.png">
</p>

---
```
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
```

<p align="center">
  <a href="https://justinblaber.org/downloads/dwmri_visualizer/MT_SCSD.pdf"><img width="611" height="791" src="https://i.imgur.com/rhafPYz.png"></a>
</p>
