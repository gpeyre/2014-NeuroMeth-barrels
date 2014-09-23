function [I2r,c2r,mask2r,c1,TR,TT] = icp_register_vessels(I1,I2,sigma,...
                                        nvess,thresh,thetas,alignCentroids,...
                                        roi1,roi2,mask1,mask2,debug)
    [N,M,~] = size(I1);
    BW1 = roi1&mask1;
    BW2 = roi2&mask2;
    % 1) Detect vertical vessels with normalized cross correlation
    [c1,corr1,radius1] = detect_perpendicular_vessels(I1.*BW1,sigma,thresh,false);
    [c2,corr2,radius2]  = detect_perpendicular_vessels(I2.*BW2,sigma,thresh,false);
    % Sort by increasing correlations
    [~,idx1] = sort(corr1,'descend');
    c1 = c1(idx1,:);
    [~,idx2] = sort(corr2,'descend');
    c2 = c2(idx2,:);
    % Keep the best vessels
    nc = min([nvess,size(c1,1),size(c2,1)]); % number of vessels kept
    c1 = c1(1:nc,:);
    radius1 = radius1(idx1(1:nc));
    c2 = c2(1:nc,:);
    radius2 = radius2(idx2(1:nc));

    if debug
        patchAlphaSliderGUI(I1,I2);
        hold on;
        plot(c1(:,1), c1(:,2), 'g+');
        hold on;
        plot(c2(:,1), c2(:,2), 'r+');
    end

    [TR,TT] = pre_register(I1,c1,mask1,I2,c2,mask2,alignCentroids,thetas);

    c2r = [TR * c2' + repmat(TT, 1, size(c2,1))]';
    T2D = [TR(1,1) TR(2,1) 0; TR(1,2) TR(2,2) 0; TT(1) TT(2) 1];
    tform = maketform('affine',T2D);
    I2r = imtransform(I2, tform,...
                        'XData',[1 M],...
                        'YData',[1 N]);
    mask2r = imtransform(mask2, tform,...
                        'XData',[1 M],...
                        'YData',[1 N]);

    if debug
        % Display results
        patchAlphaSliderGUI(I1,I2pre);
        hold on;
        plot(c1(:,1), c1(:,2), 'r.');
        hold on;
        plot(c2r(:,1), c2r(:,2), 'g+');
    end
end