function [f,pos_header,pos_footer] = dwmri_base_plot(info)  
    % This creates the header and footer sections of the plot.
    %
    % Inputs:
    %   info{1} - title
    %   info{2} - pipeline writer
    %   info{3} - pipeline writer email
    %   info{4} - paper
    %   info{5} - paper contact
    %   info{6} - version    
    
    % Using these settings seems to get the figure to save and print in a
    % correct way.
    f = figure('Units', 'Inches', 'Position', [0 0 8.5 11], ...
               'PaperPosition', [0 0 8.5 11], 'PaperUnits', ...
               'Inches', 'PaperSize', [8.5 11], 'Color', 'w', 'Resize', 'off');
    
    % Set main figure components
    padding = 0.01;
    header_height = 0.05;
    footer_height = 0.05;

    % title
    pos_header = [padding 1-padding-header_height 1-2*padding header_height];
    uicontrol('style','text','units','normalized','String',{'',info{1},''},...
        'FontUnits','Normalized','FontSize',2/7,'FontWeight','bold',...
        'Position',pos_header, ...
        'BackgroundColor',[0.85 0.85 0.85],'Parent',f);

    % footer
    pos_footer = [padding padding 1-2*padding footer_height];
    uicontrol('style','text','units','normalized','String',{'',['   ' info{1} ' pipeline written by: ' info{2} '; Contact: ' info{3}],['   Paper: ' info{4} '; Contact: ' info{5}], ['   Version: ' info{6} '; Date: ' char(datetime)]},...
        'FontUnits','Normalized','FontSize',1/6,...
        'Position',pos_footer,'HorizontalAlignment','left',...
        'BackgroundColor',[0.95 0.95 0.95],'Parent',f); 
end

