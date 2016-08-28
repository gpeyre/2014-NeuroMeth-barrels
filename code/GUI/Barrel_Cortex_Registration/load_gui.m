function [ ] = load_gui( varargin )
   %Page 1: load data
    
    S.images = [];
    S.hImages = [];
    
   S.f = figure('visible','off',...
              'Units','normalized',...
              'outerposition',[0 0 1 1],...
              'menubar','none',...
              'resize','on',...
              'numbertitle','off',...
              'name','Barrel cortex histological registration');
    
    S.hPanel = uipanel('Position',[0.05 0.1 0.9 0.8], 'Parent',S.f,'tag','hPanel'); 
          
    S.load_btn = uicontrol('style','push',...
                 'units','normalized',...
                 'position',[0.05 0.9 0.9 0.07],...
                 'HorizontalAlign','left',...
                 'string','Load Images',...
                 'fontsize',11,'fontweight','bold',...
                 'callback',{@load_btn_call,S});
             
    S.next_btn = uicontrol('style','push',...
                 'units','normalized',...
                 'position',[0.75 0.02 0.2 0.05],...
                 'HorizontalAlign','left',...
                 'string','Next (segmentation)',...
                 'fontsize',11,...
                 'callback',{@next_btn_call});
    
    set(S.f,'Visible','on')
    
    % Global variables
    hMainImg = -1;
    hThumImg = [];
    numImgs = 0;
    hcbSelectSection = [];
    selected = [];
    
    %# callback functions
    function [] = load_btn_call(varargin)
        % Load images and associate a check box to each images for
        % selection
        %S = varargin{3};
            arrayfun(@delete,allchild(S.hPanel));

            [File, Folder] = uigetfile('/Users/Mellowmany/Documents/Stage_2013/Data/Histology/*.*', 'MultiSelect', 'on');
            S.images = [];%cell(1, length(File));
            ncol = ceil(length(File)/2);
            nrow = 2;
            if ~iscell(File)
                % There is only one file selected
                File = {File};
            end
            numImgs = length(File);
            for iFile = 1:numImgs
               filename = strcat(Folder, File{iFile});
               image = mat2gray((rgb2gray(imread(filename))));
    % %           S.hImages(iFile) = axes('Parent',S.hPanel, ...
    % %             'Position',[mod(iFile-1,ncol)/ncol 1-(floor((iFile-1)/ncol)+1)/nrow 1/ncol 1/nrow],'Tag','toto');
    % %           %axes(hImages(iFile));  % Better use a vector to store handles
    % %           imagesc(image);axis image;colormap gray
    % %           set(S.hImages(iFile),'Visible','on');
    % %           set(S.hImages(iFile),'xtick',[]);
    % %           set(S.hImages(iFile),'ytick',[]);
    % %           title(File{iFile});
               S.images(iFile).im = image;
               tmp = strsplit(File{iFile},'.');
               S.images(iFile).name = char(tmp(1));
             end

            %# Design GUI
            % Set the name of the window with the name of the data folder
            set(S.f,'name',Folder);
            hPanelViewer = uipanel('Position',[0 0 1 0.2], 'Parent',S.hPanel);
            %# main axis, and show first frame
            hAx = axes('Position',[0 0.2 1 0.8], 'Parent',S.hPanel);
            hMainImg = imagesc(S.images(1).im, 'Parent',hAx);axis image;colormap gray
            %axis(hAx, 'normal')
            set(hAx,'xtick',[]);
            set(hAx,'ytick',[]);
            %# thumbnail axes
            hThumImg = zeros(numImgs,1);
            selected = ones(numImgs,1);
            for i=1:numImgs
                %# create axis, show frame, hookup click callback
                hAx = axes('Parent',hPanelViewer, ...
                    'Position',[(i-1)/numImgs 0.2 1/numImgs 0.9]);
                hThumImg(i) = imagesc(S.images(i).im, 'Parent',hAx);axis image;colormap gray
                set(hThumImg(i), 'ButtonDownFcn',@click_thumb_call)
                set(hAx,'xtick',[]);
                set(hAx,'ytick',[]);
                title(File{i});
                % Check box for selection
                hcbSelectSection(i) = uicontrol(...
                'Parent',hPanelViewer,'Units','normalized','CData',[],...
                'Position',[(i-1)/numImgs+1/(2*numImgs) 0.05 0.1 0.1],...
                'String','','Style','checkbox','UserData',[],'Value',true,...
                'Callback',{@checkbox_call,i});
            end

            assignin('base','images',S.images); % Save data
            %set(S.hPanel,'Visible','on');
        
    end

    function [] = next_btn_call(varargin)
        set(S.f,'visible','off');
        segment_gui(S.images(find(selected==1)),S); % Call next page GUI
    end

    function click_thumb_call(src,~)
        %# update the main image
        set(hMainImg, 'CData',get(src,'CData'));
        drawnow
    end

    function checkbox_call(varargin) 
        selected(varargin{3}) = get(varargin{1},'Value');
    end
end

