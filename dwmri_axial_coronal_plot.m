function f = dwmri_axial_coronal_plot(dv,info)
    % Inputs:
    %   info - see dwmri_info_plot for specification of "info"
               
    % Get dwmri info plot ------------------------------------------------%
    [f,pos_header,pos_info,pos_footer] = dwmri_info_plot(info); %#ok<ASGLU>
            
    % Set up axes --------------------------------------------------------%
    padding = 0.01;
    axes_header_height = 0.015;
    axes_area_height = pos_info(2)-(pos_footer(2)+pos_footer(4))-4*padding - axes_header_height;
    axes_area_width = 1-2*padding;    
    plot_height = 0.90*(axes_area_height/3);
    plot_width = 0.60*(axes_area_width/2);
    
    % Do Axial axes    
    for i = 1:3
        pos_axial(i,:) = [padding+(1/4)*axes_area_width-plot_width/2 pos_footer(2)+pos_footer(4)+padding+(3-i+.5)/3*axes_area_height-plot_height/2 plot_width plot_height]; %#ok<AGROW>
        axes_axial(i) = axes('Position',pos_axial(i,:),'Parent',f); %#ok<AGROW>
    end
    
    % Do Coronal axes
    for i = 1:3
        pos_coronal(i,:) = [padding+(3/4)*axes_area_width-plot_width/2 pos_footer(2)+pos_footer(4)+padding+(3-i+0.5)/3*axes_area_height-plot_height/2 plot_width plot_height]; %#ok<AGROW>
        axes_coronal(i) = axes('Position',pos_coronal(i,:),'Parent',f); %#ok<AGROW>
    end
    
    % End of axes --------------------------------------------------------%
        
    % Do header for axial and coronal
    pos_header_axial = [pos_axial(1,1) pos_axial(1,2)+pos_axial(1,4)+2*padding pos_axial(1,3) axes_header_height];
    uicontrol('style','text','units','normalized','String',{'Axial'},...
        'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
        'Position',pos_header_axial, ...
        'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
    
    pos_header_coronal = [pos_coronal(1,1) pos_coronal(1,2)+pos_coronal(1,4)+2*padding pos_coronal(1,3) axes_header_height];
    uicontrol('style','text','units','normalized','String',{'Coronal'},...
        'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
        'Position',pos_header_coronal, ...
        'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
        
    % Do plots -----------------------------------------------------------%
    
    % Get centroid of mask for location of axial and coronal slices; must 
    % do this in RAS orientation
    mask_vol_RAS = vol_utils.convert_to_xform(dv.get_mask_vol(),dv.get_xform_RAS());         
    rp_mask_vol = regionprops(double(mask_vol_RAS),'Centroid','Area');
    centroid = round(rp_mask_vol.Centroid);
    if length(centroid) == 2 % vol might be 1D or 2D by matlabs standards if trailing dimensions are 1
        centroid(3) = 1;
    end
    centroid = centroid([2 1 3]); % (x,y,z) => (j,i,k)
    
    % Plot axial
    dv.plot_slice(centroid(3),'axial','centroid',0.75,axes_axial(1));
    dv.plot_slice(centroid(3),'axial','top-centroid',0.85,axes_axial(2));
    dv.plot_slice(centroid(3),'axial','bottom-centroid',0.85,axes_axial(3));
    
    % Plot Coronal 
    dv.plot_slice(centroid(2),'coronal','centroid',1,axes_coronal(1));
    dv.plot_slice(centroid(2),'coronal','top-centroid',1,axes_coronal(2));
    dv.plot_slice(centroid(2),'coronal','bottom-centroid',1,axes_coronal(3));
        
    % set all to axis image
    axis(axes_axial(1),'image');
    axis(axes_axial(2),'image');
    axis(axes_axial(3),'image');
    axis(axes_coronal(1),'image');
    axis(axes_coronal(2),'image');
    axis(axes_coronal(3),'image');
    
    % Plot R-L, A-P, S-I labels ------------------------------------------%
    % Must update position of axes since they are scaled
    
    % Do Axial axes    
    for i = 1:3
        pos_axial(i,:) = matlab_utils.plotboxpos(axes_axial(i));
    end
    
    % Do Coronal axes
    for i = 1:3
        pos_coronal(i,:) = matlab_utils.plotboxpos(axes_coronal(i));
    end
    
    % Plot Left/Right for both axial and coronal slices
    text_height = 0.015;
    text_width = 0.015;
    for i = 1:3
        pos_L_axial = [pos_axial(i,1)-text_width/2 pos_axial(i,2)+pos_axial(i,4)/2-text_height/2 text_width text_height];
        uicontrol('style','text','units','normalized','String',{'L'},...
            'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
            'Position',pos_L_axial, ...
            'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
                
        pos_R_axial = [pos_axial(i,1)+pos_axial(i,3)-text_width/2 pos_axial(i,2)+pos_axial(i,4)/2-text_height/2 text_width text_height];
        uicontrol('style','text','units','normalized','String',{'R'},...
            'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
            'Position',pos_R_axial, ...
            'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
        
        pos_L_coronal = [pos_coronal(i,1)-text_width/2 pos_coronal(i,2)+pos_coronal(i,4)/2-text_height/2 text_width text_height];
        uicontrol('style','text','units','normalized','String',{'L'},...
            'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
            'Position',pos_L_coronal, ...
            'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
                
        pos_R_coronal = [pos_coronal(i,1)+pos_coronal(i,3)-text_width/2 pos_coronal(i,2)+pos_coronal(i,4)/2-text_height/2 text_width text_height];
        uicontrol('style','text','units','normalized','String',{'R'},...
            'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
            'Position',pos_R_coronal, ...
            'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
    end
    
    % Plot Anterior/Posterior for axial slices and Superior/Inferior for
    % coronal slices
    for i = 1:3
        pos_A_axial = [pos_axial(i,1)+pos_axial(i,3)/2-text_width/2 pos_axial(i,2)+pos_axial(i,4)-text_height/2 text_width text_height];
        uicontrol('style','text','units','normalized','String',{'A'},...
            'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
            'Position',pos_A_axial, ...
            'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
                
        pos_P_axial = [pos_axial(i,1)+pos_axial(i,3)/2-text_width/2 pos_axial(i,2)-text_height/2 text_width text_height];
        uicontrol('style','text','units','normalized','String',{'P'},...
            'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
            'Position',pos_P_axial, ...
            'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
        
        pos_S_coronal = [pos_coronal(i,1)+pos_coronal(i,3)/2-text_width/2 pos_coronal(i,2)+pos_coronal(i,4)-text_height/2 text_width text_height];
        uicontrol('style','text','units','normalized','String',{'S'},...
            'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
            'Position',pos_S_coronal, ...
            'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
                
        pos_I_coronal = [pos_coronal(i,1)+pos_coronal(i,3)/2-text_width/2 pos_coronal(i,2)-text_height/2 text_width text_height];
        uicontrol('style','text','units','normalized','String',{'I'},...
            'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
            'Position',pos_I_coronal, ...
            'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
    end
end

