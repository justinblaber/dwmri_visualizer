classdef dwmri_visualizer < handle %#ok<*PROP,*INUSL,*PROPLC>
    properties (Access = public)
        % For all plots with lines, glyphs, etc... assumes all data 
        % processing was done with respect to "voxel convention". If 
        % processing was done with dwmri_modelfit library, then everything 
        % should work properly.
        
        % Inputs
        data        % data corresponding to type
        bg_vol      % volume which is displayed as the background - for display purposes this should be scaled between 0 and 1
        mask_vol    % mask - can be used to modulate size of lines/glyphs as well, if input is double
        xform_RAS   % 3x3 transform (no translation) to "voxel" RAS, usually obtained from the nifti header
        type        % string which stores the "type"
        info        % info depends on the type
        
        % Supported types:
        %   1) 'vol':                                                      used for plotting volumes
        %       data - none
        %       bg_vol - vol
        %       info - none
        %
        %   2) 'outlines':                                                 used for plotting outlines of logical volume inputs
        %       data - {outline_vol1, ..., outline_volN}
        %       bg_vol - vol
        %       info - {line_color1, ..., line_colorN; 
        %               line_width1, ..., line_widthN}                     widths are optional
        %
        %   3) 'colorized_FA':
        %       data - {V1_vol, FA_vol}
        %       bg_vol - none
        %       info - {line_color, line_width}                            line_color and line_width are optional
        %
        %   4) 'directions':                                               for plotting peaks
        %       data - {direction_vol1, ..., direction_volN; 
        %               threshold_vol1, ..., threshold_volN}               thresholds are optional
        %       bg_vol - vol
        %       info - {line_color1, ..., line_colorN; 
        %               line_width1, ..., line_widthN;                     widths are optional
        %               threshold1, ..., thresholdN;                       thresholds are optional
        %               scale_factor1, ..., scale_factorN}                 optional input for scaling the lines as they are displayed
        %
        %   5) 'sh_coefs':
        %       data - sh_coefs_vol
        %       bg_vol - vol
        %       info - {lmax, sphere_num, min_max_normalization}
        %
        %   6) 'OD':                                                       This is orientational distribution, basically for the SD method, which does not return spherical harmonics
        %       data - OD_vol
        %       bg_vol - vol
        %       info - min_max_normalization
        %
        %   7) 'PAS':                                                      Assumes reduced encoding was used
        %       data - PAS_vol
        %       bg_vol - vol
        %       info - {pas_radius, reduced_encoding, camino_path, sphere_num, min_max_normalization}
        %
        %   8) 'DT':                                                       diffusion tensor - stored values must be 4D with 4th dim of: [D_xx D_xy D_xz D_yy D_yz D_zz]
        %       data - {DT_vol, V1_vol}                                    V1_vol is optional; if provided, the diffusion tensor will be colorized by v1 orientation 
        %       bg_vol - vol
        %       info - {sphere_num,exponent}                               exponent is optional
        %
        % Updates:
        %       v1.1.0 (8  Mar 2017): Added diffusion tensor glyphs, more 
        %           positions ('slice', 'boundingbox', etc...) to plot, 
        %           sagittal view, and "vol" and "outlines" input types
        %       v1.2.0 (24 May 2017): Put an explicit background volume
        %           option. Also added the option to colorize DT glyphs by 
        %           V1 and also modulate size through the mask.        
    end
        
    methods (Static) 
        function xform_plot = get_xform_plot(orientation)            
            % xform_plot is the transform applied AFTER RAS transform has 
            % been applied. This will get the plots to appear as specified 
            % in the plot_slice() function header.
            switch orientation
                case 'axial'                 
                    xform_plot = [0  -1   0; ...
                                  1   0   0; ...
                                  0   0   1];
                case 'coronal'                   
                    xform_plot = [0   0  -1; ... 
                                  1   0   0; ...
                                  0  -1   0];
                case 'sagittal'                   
                    xform_plot = [0   0  -1; ... 
                                  0   1   0; ...
                                  1   0   0];
            end   
        end
                
        function slice = apply_xfm_3D_slice(slice,xform)     
            % Save size
            size_slice = size(slice);
            
            % Reshape and apply xform
            slice = xform * reshape(slice,[],3)';
            
            % Reshape back to slice
            slice = reshape(slice',size_slice(1),size_slice(2),[]);
        end 
        
        function plot_direction(i,j,direction,line_size,line_color,line_width,parent_axes)
            % Input to line(): (i,j) => (y,x)
            p1x = j - (line_size/2)*direction(2);
            p1y = i - (line_size/2)*direction(1);
            p1z = 1 - (line_size/2)*direction(3); % Add 1 so line is "above" image and can be seen fully

            p2x = j + (line_size/2)*direction(2);
            p2y = i + (line_size/2)*direction(1);
            p2z = 1 + (line_size/2)*direction(3); % Add 1 so line is "above" image and can be seen fully

            line([p1x p2x],[p1y p2y],[p1z p2z],'color',line_color,'LineWidth',line_width,'parent',parent_axes);
        end    
        
        function plot_OD(i,j,OD,sphere_coords,orientation,min_max_normalized,OD_size,OD_color,parent_axes)
            % Remove negative values
            OD(OD < 0) = 0; 

            % Normalize OD
            if min_max_normalized % This applies min-max normalization, usually for Q-ball
                OD_normalized = (OD-min(OD(:)))./(max(OD(:))-min(OD(:)));
            else
                OD_normalized = OD./max(OD(:));
            end
            
            % Set the size
            OD_normalized = OD_normalized.*OD_size;

            % Create surface
            OD_surf_x = sphere_coords(:,:,2).*OD_normalized;  % x and y are flipped
            OD_surf_y = sphere_coords(:,:,1).*OD_normalized;
            OD_surf_z = sphere_coords(:,:,3).*OD_normalized;

            % Set color array  
            if isempty(OD_color)
                switch orientation
                    case 'axial'
                        color_idx = [1 2 3];
                    case 'coronal'
                        color_idx = [1 3 2];
                    case 'sagittal'
                        color_idx = [3 1 2];
                end
                OD_surf_norm = sqrt(OD_surf_x.^2 + OD_surf_y.^2 + OD_surf_z.^2);
                OD_surf_color = cat(3,abs(OD_surf_x./OD_surf_norm),abs(OD_surf_y./OD_surf_norm),abs(OD_surf_z./OD_surf_norm));
                OD_surf_color = OD_surf_color(:,:,color_idx);
            else
                % Typically for colorized diffusion tensor by V1
                OD_surf_color = repmat(permute(OD_color,[1 3 2]),size(OD_normalized,1),size(OD_normalized,2));
            end

            % Plot as surface
            scale_factor = 2.5;
            surf(parent_axes, ...
                 OD_surf_x./scale_factor+j, ...
                 OD_surf_y./scale_factor+i, ...
                 OD_surf_z./scale_factor+1, ... % Add one to raise it above imshow
                 OD_surf_color); 
        end
    end
    
    methods (Access = public)
        function obj = dwmri_visualizer(data, bg_vol, mask_vol, xform_RAS, type, info)
            obj.data = data;
            obj.mask_vol = mask_vol;
            obj.xform_RAS = xform_RAS;
            obj.type = type;
            if exist('info','var')
                obj.info = info;
            end
            
            % Set bg image
            if strcmp(type,'colorized_FA')
                % bg is colorized FA - this gets set to RAS orientation later
                obj.bg_vol = abs(bsxfun(@times,data{1},data{2})); 
            else
                % Every other type has a manual volume set; store as rgb
                if size(bg_vol,4) == 1
                    obj.bg_vol = repmat(bg_vol,[1 1 1 3]);
                elseif size(bg_vol,4) == 3
                    obj.bg_vol = bg_vol;
                else
                    error(['4th dimension of bg_vol must have size of 1 or 3; size of: ' num2str(size(bg_vol,4)) ' was found.' ]);
                end
            end
        end       
        
        function data = get_data(obj)
            data = obj.data;
        end
        
        function mask_vol = get_mask_vol(obj)
            mask_vol = obj.mask_vol;
        end
        
        function xform_RAS = get_xform_RAS(obj)
            xform_RAS = obj.xform_RAS;
        end
        
        function type = get_type(obj)    
            type = obj.type;
        end   
        
        function info = get_info(obj)              
            info = obj.info;
        end            
                
        function bg_vol = get_bg_vol(obj)
            bg_vol = obj.bg_vol;
        end       
                        
        function slice = get_slice_in_plot_orientation(obj, slice_num, vol, orientation)
            % This will grab the slice and orient it correctly into "plot"
            % orientation. slice_num is WRT RAS orientation.
            xform_RAS = obj.get_xform_RAS(); 
            
            % Get slice dimension 
            switch orientation
                case 'axial'
                    slice_dim = find(xform_RAS(3,:));  
                    slice_dim = xform_RAS(3,slice_dim)*slice_dim;
                case 'coronal'
                    slice_dim = find(xform_RAS(2,:)); 
                    slice_dim = xform_RAS(2,slice_dim)*slice_dim;
                case 'sagittal'
                    slice_dim = find(xform_RAS(1,:)); 
                    slice_dim = xform_RAS(1,slice_dim)*slice_dim;
            end   
            
            % Grab slice - make sure not to squeeze() here since this will
            % get reoriented.
            switch slice_dim
                case 1
                    slice = vol(slice_num,:,:,:);
                case 2
                    slice = vol(:,slice_num,:,:);
                case 3
                    slice = vol(:,:,slice_num,:);
                case -1
                    slice = vol(size(vol,1)-slice_num+1,:,:,:);
                case -2
                    slice = vol(:,size(vol,2)-slice_num+1,:,:);
                case -3
                    slice = vol(:,:,size(vol,3)-slice_num+1,:);
            end
                        
            % Get "total" xform to get from storage orientation to plot orientation.
            xform_total = dwmri_visualizer.get_xform_plot(orientation) * xform_RAS;              
                                
            % Convert slice to plot orientation - still may need to apply 
            % xform to components, but that is left to the caller.       
            slice = vol_utils.convert_to_xform(slice,xform_total);
            
            % You can't squeeze() here since one of the dimensions of the 
            % slice could be 1. Just remove 3rd dimension.
            slice = reshape(slice,[size(slice,1) size(slice,2) size(slice,4)]);
        end
                
        function plot_slice(obj,slice_num,orientation,position,window_scale,parent_axes)
            % Inputs:
            %   slice_num is WRT orientation
            %   orientation can either be: 'axial', 'coronal', or 'sagittal'
            %   position: (NOTE: this is relative to plot position)
            %       centroid
            %       top-centroid
            %       bottom-centroid
            %       left-centroid
            %       right-centroid
            %       slice
            %       boundingbox
            %       TODO: Add more positions
            %   window_scale determines how "zoomed-in" the figure is. Only
            %       applicable to *centroid input
            %   parent_axes is optional; if omitted, a new figure will open.
            %
            % Plots are:
            %
            %   Axial:          Coronal:        Sagittal:
            %       A               S               S
            %       |               |               |
            %   L ----- R       L ----- R       P ----- A          
            %       |               |               |
            %       P               I               I
            %
            
            if ~exist('parent_axes','var')
                f = figure();
                parent_axes = axes('parent',f);
            end    
            
            % ------------------------------------------------------------%
            % Get info ---------------------------------------------------%
            % ------------------------------------------------------------%
            
            data = obj.get_data();
            bg_vol = obj.get_bg_vol();
            mask_vol = obj.get_mask_vol();
            xform_RAS = obj.get_xform_RAS();
            type = obj.get_type();
            info = obj.get_info();
            
            % Get "total" xform to get from storage orientation to plot orientation.
            xform_total = dwmri_visualizer.get_xform_plot(orientation) * xform_RAS; 
            
            % ------------------------------------------------------------%
            % Get slice  -------------------------------------------------%
            % ------------------------------------------------------------%
            disp('Getting slice...');
                                                
            % Grab background and mask slice - handle the components on a per-type basis.
            bg_slice = obj.get_slice_in_plot_orientation(slice_num,bg_vol,orientation);
            mask_slice = obj.get_slice_in_plot_orientation(slice_num,mask_vol,orientation);
                        
            % Finish slices
            if strcmp(type,'outlines')                                       
                % Handle data slice      
                data_slice = {};
                for i = 1:length(data)
                    data_slice{i} = obj.get_slice_in_plot_orientation(slice_num,data{i},orientation); %#ok<AGROW>
                end                   
            elseif strcmp(type,'colorized_FA')
                % Handle data slice 
                data_slice = obj.get_slice_in_plot_orientation(slice_num,data{1},orientation);
                data_slice = dwmri_visualizer.apply_xfm_3D_slice(data_slice,xform_total);

                % Handle background slice - NOTE: this uses xform_RAS, NOT xform_total
                bg_slice = abs(dwmri_visualizer.apply_xfm_3D_slice(bg_slice,xform_RAS));                                                              
            elseif strcmp(type,'directions')                                          
                % Handle data slice      
                data_slice = {};
                for i = 1:size(data,2)           
                    data_slice{1,i} = obj.get_slice_in_plot_orientation(slice_num,data{1,i},orientation); %#ok<AGROW>
                    data_slice{1,i} = dwmri_visualizer.apply_xfm_3D_slice(data_slice{1,i},xform_total); %#ok<AGROW>
                    if size(data,1) >= 2
                        % Threshold is scalar so no reorientation of
                        % components is needed.
                        data_slice{2,i} = obj.get_slice_in_plot_orientation(slice_num,data{2,i},orientation); %#ok<AGROW>
                    end
                end   
            elseif strcmp(type,'sh_coefs') % For ODs just apply xform to sphere coords - this is the simplest/fastest way to reorient in matlab
                data_slice{1} = obj.get_slice_in_plot_orientation(slice_num,data,orientation);
                % Get sphere_coords
                [sphere_x,sphere_y,sphere_z] = sphere(info{2});
                data_slice{2} = dwmri_visualizer.apply_xfm_3D_slice(cat(3,sphere_x,sphere_y,sphere_z),xform_total);
            elseif strcmp(type,'OD')
                data_slice{1} = obj.get_slice_in_plot_orientation(slice_num,data,orientation);
                % Get sphere_coords
                [sphere_x,sphere_y,sphere_z] = sphere(sqrt(length(data_slice{1}))-1);
                data_slice{2} = dwmri_visualizer.apply_xfm_3D_slice(cat(3,sphere_x,sphere_y,sphere_z),xform_total);
            elseif strcmp(type,'PAS')
                data_slice{1} = obj.get_slice_in_plot_orientation(slice_num,data,orientation);
                % Get sphere_coords
                [sphere_x,sphere_y,sphere_z] = sphere(info{4});
                data_slice{2} = dwmri_visualizer.apply_xfm_3D_slice(cat(3,sphere_x,sphere_y,sphere_z),xform_total);      
            elseif strcmp(type,'DT')
                data_slice{1} = obj.get_slice_in_plot_orientation(slice_num,data{1},orientation);  
                % Get sphere_coords
                [sphere_x,sphere_y,sphere_z] = sphere(info{1});
                data_slice{2} = dwmri_visualizer.apply_xfm_3D_slice(cat(3,sphere_x,sphere_y,sphere_z),xform_total);      
                
                % Get V1 in RAS configuration if V1 is provided
                if length(data) >= 2
                    % NOTE: this uses xform_RAS, NOT xform_total
                    data_slice{3} = obj.get_slice_in_plot_orientation(slice_num,data{2},orientation); 
                    data_slice{3} = abs(dwmri_visualizer.apply_xfm_3D_slice(data_slice{3},xform_RAS)); 
                end
            end
                        
            % ------------------------------------------------------------%
            % Plot slice -------------------------------------------------%
            % ------------------------------------------------------------%
            disp('Plotting slice...');
            
            position_split = strsplit(position,'-');
            if strcmp(position_split{end},'centroid')
                % Plots square window around centroid; make sure
                % mask_slice is converted to logical first (to threshold).
                % Then it must be converted to double so that region props
                % calculates a single centroid.
                rp_slice = regionprops(double(mask_slice > 0),'Centroid','Area');
                i_centroid = round(rp_slice.Centroid(2));
                j_centroid = round(rp_slice.Centroid(1));

                % Possibly alter mask_slice based on "position" input
                switch position_split{1}
                    case 'centroid'
                        % Do nothing
                    case 'top'                    
                        [~,i_idx] = meshgrid(1:size(mask_slice,2),1:size(mask_slice,1));
                        mask_slice(i_idx > i_centroid) = 0;
                    case 'bottom'
                        [~,i_idx] = meshgrid(1:size(mask_slice,2),1:size(mask_slice,1));
                        mask_slice(i_idx < i_centroid) = 0;                    
                    case 'left'
                        [j_idx,~] = meshgrid(1:size(mask_slice,2),1:size(mask_slice,1));
                        mask_slice(j_idx > j_centroid) = 0;
                    case 'right'
                        [j_idx,~] = meshgrid(1:size(mask_slice,2),1:size(mask_slice,1));
                        mask_slice(j_idx < j_centroid) = 0;
                end                

                % Regrab centroid and also get bounding box
                rp_slice = regionprops(double(mask_slice > 0),'Centroid','BoundingBox','Area');
                i_centroid = round(rp_slice.Centroid(2));
                j_centroid = round(rp_slice.Centroid(1));

                % Determine window size based on bounding box (must account for
                % half pixel in bounding box)
                window = round(min([j_centroid-(rp_slice.BoundingBox(1)+0.5) ...
                                    (rp_slice.BoundingBox(1)+rp_slice.BoundingBox(3)-0.5)-j_centroid ...
                                    i_centroid-(rp_slice.BoundingBox(2)+0.5) ...
                                    (rp_slice.BoundingBox(2)+rp_slice.BoundingBox(4)-0.5)-i_centroid])*window_scale);
                                
                left = j_centroid-window;
                right = j_centroid+window;
                top = i_centroid-window;
                bottom = i_centroid+window;
            elseif strcmp(position,'slice')
                % Plot entire slice
                left = 1;
                right = size(mask_slice,2);
                top = 1;
                bottom = size(mask_slice,1);
            elseif strcmp(position,'boundingbox')                
                % Plot bounding box
                rp_slice = regionprops(double(mask_slice > 0),'BoundingBox','Area');

                % Determine window size based on bounding box (must account
                % for half pixel in bounding box)
                left = rp_slice.BoundingBox(1)+0.5;
                right = rp_slice.BoundingBox(1)+rp_slice.BoundingBox(3)-0.5;
                top = rp_slice.BoundingBox(2)+0.5;
                bottom = rp_slice.BoundingBox(2)+rp_slice.BoundingBox(4)-0.5;                                
            else
                error(['Unrecognized position input: ' position]);
            end
            
            % Plot background slice
            imagesc(bg_slice(top:bottom,left:right,:),'Parent',parent_axes,[0 1]);
            axis(parent_axes,'off');
            hold(parent_axes,'on');
            
            % Stuff before plotting
            if strcmp(type,'outlines')
                % Cycle over each outline
                for outline_num = 1:length(data_slice)                                    
                    % Get slice and plot the outline
                    outline_slice = data_slice{outline_num};
                    
                    % Get line color
                    line_color = info{1,outline_num};
                    
                    % Get line width
                    if size(info,1) >= 2
                        line_width = info{2,outline_num};
                    else
                        line_width = 0.25;
                    end

                    % Get boundaries
                    B = bwboundaries(outline_slice);
                    for k = 1:length(B)
                        boundary = B{k};
                        plot(parent_axes, ...
                             boundary(:,2) - (left-1), ...
                             boundary(:,1) - (top-1), ...
                             'Color', ...
                             line_color, ...
                             'LineWidth', ...
                             line_width)
                    end           
                end
            elseif strcmp(type,'sh_coefs')
                % Sample spherical harmonics using sphere() points
                [sphere_x,sphere_y,sphere_z] = sphere(info{2});
                [sh_sphere,~,~] = construct_SH_basis(info{1}, [sphere_x(:) sphere_y(:) sphere_z(:)], 2, 'real');
            elseif strcmp(type,'PAS')
                % Get sphere_coords for sampling the PAS
                [sphere_x,sphere_y,sphere_z] = sphere(info{4});
                sphere_coords_sample = cat(2,sphere_x(:),sphere_y(:),sphere_z(:));   
                
                % Get pointset from camino installation
                pointset_path = fullfile(info{3},'PointSets',['Elec' sprintf('%3.3d',info{2}) '.txt']);
                pointset = dlmread(pointset_path);
                pointset = reshape(pointset(2:end),3,[])'; % First value is number of elements, so skip it
            elseif strcmp(type,'DT')
                % Get sphere_coords for sampling the DT
                [sphere_x,sphere_y,sphere_z] = sphere(info{1});
                sphere_coords_sample = cat(2,sphere_x(:),sphere_y(:),sphere_z(:));   
            end
            
            % Plot data slice
            for i = top:bottom
                for j = left:right
                    if mask_slice(i,j)
                        switch type
                            case 'colorized_FA'
                                % Get line color
                                if length(info) >= 1
                                    line_color = info{1};
                                else
                                    line_color = 'w';
                                end
                                
                                % Get line_width
                                if length(info) >= 2
                                    line_width = info{2};
                                else
                                    line_width = 0.25;                                    
                                end
                                
                                % This plots v1 direction          
                                dwmri_visualizer.plot_direction(i-top+1, ...
                                                                j-left+1, ...
                                                                squeeze(data_slice(i,j,:)), ...
                                                                mask_slice(i,j), ...
                                                                line_color, ...
                                                                line_width, ...
                                                                parent_axes);
                            case 'directions'
                                % This plots multiple directions with an
                                % optional threshold and scale factor
                                
                                % Cycle over each direction
                                for outline_num = 1:size(data_slice,2)                                    
                                    % Get line color
                                    line_color = info{1,outline_num};

                                    % Get line_width
                                    if size(info,1) >= 2
                                        line_width = info{2,outline_num};
                                    else
                                        line_width = 0.25;                                    
                                    end
                                
                                    % Check threshold                                    
                                    if size(info,1) >= 3
                                        threshold_slice = data_slice{2,outline_num};
                                        if threshold_slice(i,j) < info{3,outline_num} % info{3} is the threshold
                                            % This is below threshold, so
                                            % just skip it
                                            continue
                                        end
                                    end
                                    
                                    % Get slice and plot the direction
                                    direction_slice = data_slice{1,outline_num};
                                    
                                    % plot direction
                                    if size(info,1) <= 3
                                        % Use mask to scale direction
                                        dwmri_visualizer.plot_direction(i-top+1, ...
                                                                        j-left+1, ...
                                                                        squeeze(direction_slice(i,j,:)), ...
                                                                        mask_slice(i,j), ...
                                                                        line_color, ...
                                                                        line_width, ...
                                                                        parent_axes);
                                    else
                                        % User scale factor in info{} and
                                        % the thresholds to scale direction
                                        dwmri_visualizer.plot_direction(i-top+1, ...
                                                                        j-left+1, ...
                                                                        squeeze(direction_slice(i,j,:)), ...
                                                                        info{4,outline_num}*threshold_slice(i,j), ...
                                                                        line_color, ...
                                                                        line_width, ...
                                                                        parent_axes);
                                    end
                                end
                            case 'sh_coefs'        
                                sh_coefs_slice = data_slice{1};
                                sphere_coords = data_slice{2};
                                
                                if any(~isfinite(sh_coefs_slice(i,j,:))) || any(~isreal(sh_coefs_slice(i,j,:)))
                                    % Some coefficients are either not
                                    % finite or not real, skip this
                                    continue
                                end
                                
                                % Multiply basis by sh coefficients to get values of OD on sphere
                                OD = reshape(sh_sphere * squeeze(sh_coefs_slice(i,j,:)),[info{2}+1, info{2}+1]);
                                                                
                                % Plot OD
                                dwmri_visualizer.plot_OD(i-top+1, ...
                                                         j-left+1, ...
                                                         OD, ...
                                                         sphere_coords, ...
                                                         orientation, ...
                                                         info{3}, ...
                                                         mask_slice(i,j), ...
                                                         [], ...
                                                         parent_axes);  
                            case 'OD'     
                                OD_slice = data_slice{1};
                                sphere_coords = data_slice{2};
                                
                                if any(~isfinite(OD_slice(i,j,:))) || any(~isreal(OD_slice(i,j,:)))
                                    % Some coefficients are either not
                                    % finite or not real, skip this
                                    continue
                                end
                                
                                % Reshape
                                OD = reshape(squeeze(OD_slice(i,j,:)),[sqrt(size(OD_slice,3)) sqrt(size(OD_slice,3))]);
                                
                                % Plot OD
                                dwmri_visualizer.plot_OD(i-top+1, ...
                                                         j-left+1, ...
                                                         OD, ...
                                                         sphere_coords, ...
                                                         orientation, ...
                                                         info, ...
                                                         mask_slice(i,j), ...
                                                         [], ...
                                                         parent_axes); 
                            case 'PAS'
                                pas_slice = data_slice{1};
                                sphere_coords = data_slice{2};
                                
                                % Here are the camino files used for
                                % sampling the PAS:
                                %
                                %   ~/camino/tools/CL_Initializer.java:
                                %       initMaxEnt() - this gets the
                                %       pointset used for reduced encoding.
                                %
                                %   ~/camino/misc/SphericalPoints.java:
                                %       getElecPointSet() - this gets
                                %       pointset file "Elec#.txt" and reads
                                %       it.
                                %
                                %   ~/camino/sphfunc/MaxEntProfile.java:
                                %       setReconDirs() - this will set the
                                %       pointset.
                                %
                                %   ~/camino/sphfunc/MaxEntProfile.java:
                                %       getRadius() - this gets the actual
                                %       radius used for the PAS. It loops
                                %       over the "lambda" directly output
                                %       when running pasmri from Camino and
                                %       also uses reconDirs which is the
                                %       pointset set from setReconDirs()
                                %       and equals the number of lambdas.
                                %       This calls 
                                %       SphDeconvKernels.kernel()
                                % 
                                %   ~/camino/mesd/SphDeconvKernels.java:
                                %       kernel() - this calls pasKernel()
                                % 
                                % If you open all these java files and look
                                % at the methods, you should be able to
                                % tell where the below code comes from.                                
                                
                                % Compute PAS from lambdas.
                                PAS = zeros(size(sphere_coords_sample,2),1);
                                for k = 1:size(sphere_coords_sample,1)
                                    % Get lambdas
                                    lambdas = squeeze(pas_slice(i,j,:));

                                    % Calculate PAS
                                    expo = lambdas(1);
                                    for m = 2:length(lambdas)
                                        pasKernel = cos(info{1} * sum(pointset(m-1,:).*sphere_coords_sample(k,1:3)));
                                        expo = expo + lambdas(m) * pasKernel;
                                    end
                                    PAS(k) = exp(expo);
                                end

                                % Reshape to fit sphere coords  
                                PAS = reshape(PAS,[info{4}+1, info{4}+1]);

                                % Plot PAS
                                dwmri_visualizer.plot_OD(i-top+1, ...
                                                         j-left+1, ...
                                                         PAS, ...
                                                         sphere_coords, ...
                                                         orientation, ...
                                                         info{5}, ...
                                                         mask_slice(i,j), ...
                                                         [], ...
                                                         parent_axes); 
                            case 'DT'
                                DT_slice = data_slice{1};
                                sphere_coords = data_slice{2};
                                
                                if any(~isfinite(DT_slice(i,j,:))) || any(~isreal(DT_slice(i,j,:)))
                                    % Some DT elements are either not
                                    % finite or not real, skip this
                                    continue
                                end
                                
                                % DT surface is r = (1/sqrt(2*g'*D^-1*g))^p
                                DT = squeeze(DT_slice(i,j,:));
                                DT_mat = zeros(3);
                                DT_mat(1,1:3) = DT(1:3);
                                DT_mat(2,2:3) = DT(4:5);
                                DT_mat(3,3) = DT(6);
                                % Fill other symmetric half
                                DT_mat(2,1) = DT_mat(1,2);
                                DT_mat(3,1:2) = DT_mat(1:2,3)';    
                                
                                % Diffusion tensor must be positive
                                % definite
                                [~,p] = chol(DT_mat);
                                if p > 0
                                    continue
                                end
                                          
                                DT_mat_inv = DT_mat^-1;
                                DT_inv = [DT_mat_inv(1,1:3) DT_mat_inv(2,2:3) DT_mat_inv(3,3)]';                            
                                
                                % Get exponent
                                if length(info) >= 2
                                    exponent = info{2};
                                else
                                    exponent = 1;
                                end
                                
                                % Sample (this is vectorized version of
                                % (1/sqrt(2*g'*D^-1*g))^p
                                OD = reshape(1./sqrt(2 * [sphere_coords_sample(:,1).^2 2*sphere_coords_sample(:,1).*sphere_coords_sample(:,2) 2*sphere_coords_sample(:,1).*sphere_coords_sample(:,3) sphere_coords_sample(:,2).^2 2*sphere_coords_sample(:,2).*sphere_coords_sample(:,3) sphere_coords_sample(:,3).^2] * DT_inv), ...
                                             [info{1}+1, info{1}+1]).^exponent;
                                
                                % Get DT color (if V1 was provided)
                                dt_color = [];                                
                                if length(data) >= 2
                                    v1_slice = data_slice{3};
                                    dt_color = squeeze(v1_slice(i,j,:))';
                                end
                                         
                                % Plot OD
                                dwmri_visualizer.plot_OD(i-top+1, ...
                                                         j-left+1, ...
                                                         OD, ...
                                                         sphere_coords, ...
                                                         orientation, ...
                                                         false, ...
                                                         mask_slice(i,j), ...
                                                         dt_color, ...
                                                         parent_axes);  
                        end
                    end
                end
                drawnow % For debugging
            end
                              
            % Stuff after plotting
            if strcmp(type,'sh_coefs') || strcmp(type,'OD') || strcmp(type,'PAS') || strcmp(type,'DT')
                % Remove grid around surf
                shading(parent_axes,'interp');                    
            end
        end
    end    
end