function merge_barrels_gui(images, prevGui,imMerge,imPostProc)
    [M,N] = size(images(1).im);
    P = length(images);
    if isempty(prevGui)
          prevGui = -1;
    end
 
    nameWindow = '';
    S.exportDir = '';
    if isstruct(prevGui)
        nameWindow = get(prevGui.f,'name');
        S.exportDir = prevGui.exportDir;
    end
    
   %% Global variables
    fileCount = 1; % to allow saving multiple files (appends a number to the filename)
   hPlotPoly = -1;
   maskVtx = [];
   %imBarrels = ones(M,N); % image containing contour of the barrels
   selected = ones(P,1);   
   maskRect = []; % coordinates [xmin,ymin,width,height] of rectangle of 
                    % interest (selected on post-processed image)
   
    for i=1:P
        images(i).imReg = mat2gray(images(i).imReg);
    end
    

S.f = figure(...
'Units','normalized',...
'outerposition',[0 0 1 1],...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name',nameWindow,...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'PaperSize',get(0,'defaultfigurePaperSize'),...
'PaperType',get(0,'defaultfigurePaperType'),...
'UserData',[],...
'Tag','figure1',...
'Visible','on',...
'resize','on');

hPanelSelect = uipanel(...
'Parent',S.f,...
'Title',{  'Panel' },...
'Clipping','on',...
'Position',[0.0397967823878069 0.767384155455904 0.892464013547841 0.224215246636771],...
'Tag','uipanel1',...
'BackgroundColor',[1 1 1]);

hAxesMerge = axes(...
'Parent',S.f,...
'Position',[0.0414902624894157 0.0669058295964126 0.268416596104996 0.449925261584454],...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'LooseInset',[0.263344768439108 0.173152941176471 0.192444253859348 0.118058823529412],...
'XColor',get(0,'defaultaxesXColor'),...
'YColor',get(0,'defaultaxesYColor'),...
'ZColor',get(0,'defaultaxesZColor'),...
'Tag','axes1',...
'xtick',[],...
'xtick',[],...
'visible','off');

himMerge = -1;
if ~isempty(imMerge)
    himMerge = imagesc(imMerge); axis image;colormap gray
    title(['Merged']);
    set(hAxesMerge,'xtick',[]);
    set(hAxesMerge,'ytick',[]);
end

hAxesPostProc = axes(...
'Parent',S.f,...
'Position',[0.354784081287044 0.0669058295964126 0.268416596104996 0.449925261584454],...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'LooseInset',[0.263344768439108 0.173152941176471 0.192444253859348 0.118058823529412],...
'XColor',get(0,'defaultaxesXColor'),...
'YColor',get(0,'defaultaxesYColor'),...
'ZColor',get(0,'defaultaxesZColor'),...
'Tag','axes2',...
'xtick',[],...
'xtick',[],...
'visible','off');

himPostProc = -1;
if ~isempty(imPostProc)
    himPostProc = imagesc(imPostProc); axis image;colormap gray
    title(['Post-Processed']);
    set(hAxesPostProc,'xtick',[]);
    set(hAxesPostProc,'ytick',[]);
end

hAxesBarrels = axes(...
'Parent',S.f,...
'Position',[0.662997459779847 0.0654110612855007 0.268416596104996 0.449925261584454],...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'LooseInset',[0.263344768439108 0.173152941176471 0.192444253859348 0.118058823529412],...
'XColor',get(0,'defaultaxesXColor'),...
'YColor',get(0,'defaultaxesYColor'),...
'ZColor',get(0,'defaultaxesZColor'),...
'Tag','axes3',...
'xtick',[],...
'xtick',[],...
'visible','off');

h18 = uicontrol(...
'Parent',S.f,...
'Units','normalized',...
'Callback',@merge_call,...
'Position',[0.116850127011008 0.53 0.12108382726503 0.0762331838565022],...
'String','Merge',...
'Tag','pushbutton1');

hbgPostProc = uibuttongroup(...
'Parent',S.f,...
'Title','Post-process',...
'UserData',[],...
'Clipping','on',...
'Position',[0.410668924640135 0.523273542600897 0.147332768839966 0.227204783258595],...
'Tag','uipanel3',...
'SelectedObject',[],...
'SelectionChangeFcn',[],...
'OldSelectedObject',[]);

hcbLowerDrift = uicontrol(...
'Parent',hbgPostProc,...
'Units','normalized',...
'CData',[],...
'Position',[0.1 0.83136690647482 0.488235294117647 0.165467625899281],...
'String','Lower drift',...
'Style','checkbox',...
'UserData',[],...
'Tag','checkbox3');

heSizeGaussian = uicontrol(...
'Parent',hbgPostProc,...
'Units','normalized',...
'CData',[],...
'Position',[0.235294117647059 0.671151079136691 0.576470588235294 0.12661870503597],...
'String','10',...
'Style','edit',...
'UserData',[],...
'Tag','edit2');

h24 = uicontrol(...
'Parent',hbgPostProc,...
'Units','normalized',...
'CData',[],...
'Position',[0.223529411764706 0.779424460431655 0.594117647058824 0.0935251798561151],...
'String','Size of gaussian',...
'Style','text',...
'UserData',[],...
'Tag','text2');

hcbEnhance = uicontrol(...
'Parent',hbgPostProc,...
'Units','normalized',...
'CData',[],...
'Position',[0.1 0.5336690647482 0.788235294117647 0.165467625899281],...
'String','Enhance contrast',...
'Style','checkbox',...
'UserData',[],...
'Tag','checkbox3');

hcbMedianFilter = uicontrol(...
'Parent',hbgPostProc,...
'Units','normalized',...
'CData',[],...
'Position',[0.1 0.391294964028777 0.7 0.165467625899281],...
'String','Median filter',...
'Style','checkbox',...
'UserData',[],...
'Tag','checkbox4');

heBlockMedFilt = uicontrol(...
'Parent',hbgPostProc,...
'Units','normalized',...
'CData',[],...
'Position',[0.235294117647059 0.241151079136691 0.576470588235294 0.12661870503597],...
'String','10',...
'Style','edit',...
'UserData',[],...
'Tag','edit2');

hstatBlock = uicontrol(...
'Parent',hbgPostProc,...
'Units','normalized',...
'CData',[],...
'Position',[0.223529411764706 0.349424460431655 0.594117647058824 0.0935251798561151],...
'String','Block size',...
'Style','text',...
'UserData',[],...
'Tag','text2');

h25 = uicontrol(...
'Parent',hbgPostProc,...
'Units','normalized',...
'Callback',@postproc_call,...
'CData',[],...
'Position',[0.235294117647059 0.024748201438849 0.564705882352941 0.173381294964029],...
'String','Post-process',...
'UserData',[],...
'Tag','pushbutton6');

    iconSelectDir = imresize(imread('icons/folder_icon.png'),[16 16]);
    export_dir_btn = uicontrol('style','push',...
                 'units','normalized',...
                 'position',[0.75 0.62 0.03 0.03],...
                 'HorizontalAlign','left',...%'string','Select directory',...
                 'callback',{@select_dir_call},...
                 'TooltipString','Select export directory',...
                 'cdata',iconSelectDir);

    htExportDir = uicontrol(...
        'Parent',S.f,...
        'Units','normalized',...
        'CData',[],...
        'Position',[0.79 0.62 0.15 0.03],...
        'String',S.exportDir,...
        'Style','edit',...
        'UserData',[],...
        'Tag','text10',...
        'ButtonDownFcn',@select_dir_call);

    save_btn = uicontrol('style','push',...
                 'units','normalized',...
                 'position',[0.75 0.53 0.12108382726503 0.0762331838565022],...
                 'HorizontalAlign','left',...
                 'string','Save data',...
                 'callback',{@save_call});
             
    prev_btn = uicontrol('style','push',...
                 'units','normalized',...
                 'position',[0.05 0.02 0.2 0.05],...
                 'HorizontalAlign','left',...
                 'string','Previous (Registration)',...
                 'callback',{@prev_btn_call,S,prevGui});             
             
    % Create the toolbar
    th = uitoolbar(S.f);
    set(S.f,'Toolbar','figure');
       
   % Load images on the left side
   hImages = [];
   hcbSelectSections = [];
   for i = 1:P
          hImages(i) = axes('Parent',hPanelSelect,...
            'Position',[[(i-1)/P 0.2 1/P 0.9]]);
          %axes(hImages(iFile));  % Better use a vector to store handles
          imagesc(images(i).imReg);axis image;colormap gray
          set(hImages(i),'xtick',[]);
          set(hImages(i),'ytick',[]);
          hcbSelectSections(i) = uicontrol(...
            'Parent',hPanelSelect,'Units','normalized','CData',[],...
            'Position',[(i-1/2)/P 0.05 0.1 0.1],...
            'String','','Style','checkbox','UserData',[],'Value',true,...
            'Callback',{@checkbox_call,i});
   end
        
      
   %% Callback functions
   function [] = prev_btn_call(varargin)
        delete(varargin{3}.f); % Going back to previous GUI deletes the current GUI
        set(varargin{4}.f,'visible','on');
   end

    function [] = save_call(varargin)
        watchon;
        % Export images
        date = strsplit(images(1).name,'_');
        date = date(1);
        szGauss = str2double(get(heSizeGaussian,'String'));
        szBlock = str2double(get(heBlockMedFilt,'String'));
        imwrite(imMerge,char(strcat(S.exportDir,'/',date,'_merge.png')));
        filename = strcat(S.exportDir,'/',date,'_postproc[',...
                            num2str(szGauss),'_',num2str(szBlock),']');
        filename = char(strcat(filename,sprintf('(%d).png',fileCount)));
        imwrite(imPostProc,filename);
%         filename = strcat(S.exportDir,'/',date,'_barrels');
%         filename = char(strcat(filename,sprintf('(%d).png',fileCount)));
%         imwrite(imBarrels,filename,'Alpha',1-imBarrels);
        fileCount = fileCount+1;
        watchoff;
    end

    function [] = merge_call(varargin)
        watchon;drawnow;
        imMerge = reconstruct_barrels({images(find(selected==1)).imReg},{images(find(selected==1)).maskReg},false);
        axes(hAxesMerge)
        himMerge = imagesc(imMerge); axis image;colormap gray
        title(['Merged']);
        set(hAxesMerge,'xtick',[]);
        set(hAxesMerge,'ytick',[]);
        assignin('base','imMerge',imMerge); % Save data
       	watchoff;drawnow;
    end

    function [] = postproc_call(varargin)
        watchon;drawnow;
        tmp = imMerge;
        if get(hcbLowerDrift,'Value')
            tmp = lower_drift(tmp,str2double(get(heSizeGaussian,'String')),false);
        end
        if get(hcbEnhance,'Value')
            tmp = adapthisteq(tmp);
        end
        if get(hcbMedianFilter,'Value')
            sz = str2double(get(heBlockMedFilt,'String'))
            tmp = mat2gray(medfilt2(tmp,[sz sz]));
            %adapthisteq(tmp,'NumTiles',[20 20],'ClipLimit',0.05);
        end
        imPostProc = tmp;
        %imPostProcSmall = resize(imPostProc,scale);
        axes(hAxesPostProc)
        himPostProc = imagesc(imPostProc); axis image;colormap gray
        if (size(maskVtx,1) >0 )
                    hold on; plot([maskVtx(:,2);maskVtx(1,2)],[maskVtx(:,1);maskVtx(1,1)]);hold off;
                end
        %hImMask = imagesc(mask);axis image; colormap gray;
        %set(hImMask, 'AlphaData', 0.2); hold off
        title(['Post-processed']);
        set(hAxesPostProc,'xtick',[]);
        set(hAxesPostProc,'ytick',[]);
        assignin('base','imPostProc',imPostProc); % Save data
        watchoff;drawnow;
    end

    function checkbox_call(varargin) 
        selected(varargin{3}) = get(varargin{1},'Value');
    end

    function[] = select_dir_call(varargin)
        tmp = uigetdir();
        if tmp ~= 0
            S.exportDir = tmp;
            set(htExportDir, 'String', S.exportDir);
        end
    end
end