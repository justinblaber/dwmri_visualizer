function [f,pos_header,pos_info,pos_footer] = dwmri_info_plot(info)  
    % This creates the header, info, and  footer sections of the plot.
    %
    % Inputs:
    %   info{1} - title
    %   info{2} - header field 1
    %   info{3} - header field 2
    %   info{4} - header field 3
    %   info{5} - header field 4
    %   info{6} - header field 5
    %   info{7} - pipeline writer
    %   info{8} - pipeline writer email
    %   info{9} - paper
    %   info{10} - paper contact
    %   info{11} - version    
    
    % get base plot first
    [f,pos_header,pos_footer] = dwmri_base_plot([info(1) info(7:end)]);
    
    % Set main figure components
    padding = 0.01;
    info_height = 0.08;

    % info 
    pos_info = [padding pos_header(2)-padding-info_height 1-2*padding info_height];
    uicontrol('style','text','units','normalized','String',{'',info{2},info{3},info{4},info{5},info{6},''},...
        'FontUnits','Normalized','FontSize',13/110,...
        'Position',pos_info,'HorizontalAlignment','left',...
        'BackgroundColor',[0.95 0.95 0.85],'Parent',f);
end