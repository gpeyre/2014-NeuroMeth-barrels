function patchAlphaSliderGUI(I1, I2)
    %# setup GUI
    scrsz = get(0,'ScreenSize');
    hFig = figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)],'toolbar','figure');
    hAx = axes('Parent',hFig);
    hSlider = uicontrol('Parent',hFig, 'Style','slider', ...
        'Value',0.5, 'Min',0, 'Max',1, 'SliderStep',[1 10]./100, ...
        'Units','pixels', 'Position',[150 5 300 20], ...
        'Callback',@slider_callback);
    hTxt = uicontrol('Parent',hFig, 'Style','text', 'String','0.5', ...
        'Units','pixels', 'Position',[290 28 30 15]);
    
    % overlay the two images
    imagesc(I1); axis image; colormap gray; hold on;
    hIm = imagesc(I2);axis image; colormap gray;
    set(hIm, 'AlphaData', 0.5);
    hold off;

    function slider_callback(hObj,eventdata)
        %# get new alpha value
        alpha = get(hObj,'Value');

        %# update patches transparency and label
        set(hTxt, 'String',num2str(alpha,'%.02f'))
        %set(hPatch, 'FaceAlpha',alpha)
        set(hIm, 'AlphaData', alpha);
    end
end