function [TR,TT] = pre_register(I1,centroids1,mask1,...
                                        I2,centroids2,mask2,center,thetas)
    [N,M,~] = size(I1);
    P = size(thetas,2);
    
    TR = eye(2,2); % rotation
    TT = zeros(2,1); % translation
    c1 = regionprops(mask1~=0,'Centroid');
    c1 = c1(1).Centroid;
    c = c1;%[floor(M/2) floor(N/2)];
    c2 = regionprops(mask2~=0,'Centroid');
    c2 = c2(1).Centroid;
    if center
        %% Center the centroids of the two slices
        % Second slice
        TMP = [1 0 0; 0 1 0; c1(1)-c2(1) c1(2)-c2(2) 1];
        TT = TT + [c1(1)-c2(1);c1(2)-c2(2)];
        centroids22 = move_points(centroids2,TMP,[0 0]);
    else
        c = c2;
        centroids22 = centroids2;
    end
        
    if P > 0
        %% Find the best initialization
        r0 = 50;
        nrjs = zeros(P,1);
        nrjstart = zeros(P,1);
        %TT0 = repmat(TT,[1,1,P]);
        %TR0 = repmat(eye(2,2), [1,1,P]);

        for j=1:P
            %I2pre = imrotate(I2,thetas(j),'crop'); % rotates around center of image
            th = thetas(j)*pi/180; % angle in radians
            %TR0 = TR*[cos(th) -sin(th); sin(th) cos(th)];
            TMP = [cos(th) sin(th) 0; -sin(th) cos(th) 0; 0 0 1]; % In the corresponding image the axis Y is directed to the bottom
            centroids222 = move_points(centroids22,TMP,c);
            [TR1, TT1, E] =  robust_icp(centroids1', centroids222',...
                @(d,r) r^2/6*(1-(1-(d/r).^2).^3).*(abs(d)<=r) + r^2/6.*(abs(d)>r),...
                @(d,r) ((1-(d/r).^2).^2).*(abs(d)<=r),...
                @(r) r, r0,5,5);
            %TR0(:,:,j) = TR1*TR0(:,:,j);
            %TT0(:,:,j) = TR1*TT0(:,:,j)+ TT1;
            nrjs(j) = E(end);
            nrjstart(j) = E(1);
        end
        [~,theta_idx] = min(nrjs); 
        %I2pre = imrotate(I2,thetas(theta_idx),'crop');
        %thetas(theta_idx)
        th = thetas(theta_idx)*pi/180; % angle in radians
        TT = TT - c'; % move centroid to the center so that the rotation is done about the center
        TR = TR*[cos(th) sin(th); -sin(th) cos(th)];
        TT = TR*TT + c'; % move back centroid to its place
        TMP = [cos(th) sin(th) 0; -sin(th) cos(th) 0;  0 0  1];
        centroids222 = move_points(centroids22,TMP,c);
    
        % rerun ICP
        [TR1, TT1, E] =  robust_icp(centroids1', centroids222',...
            @(d,r) r^2/6*(1-(1-(d/r).^2).^3).*(abs(d)<=r) + r^2/6.*(abs(d)>r),...
            @(d,r) ((1-(d/r).^2).^2).*(abs(d)<=r),...
            @(r) r,r0,5,5);
        TMP = [TR1(1,1) TR1(1,2) 0; TR1(2,1) TR1(2,2) 0; TT1(1) TT1(2) 1];
        centroids2222 = move_points(centroids222,TMP,[0 0]);%[TR1 * centroids222' + repmat(TT1, 1, size(centroids222,1))]';
        
        TR = TR1*TR;
        TT = TR1*TT + TT1;
    end
    centroids2r = [TR * centroids2' + repmat(TT, 1, size(centroids2,1))]';
    T2D = [TR(1,1) TR(2,1) 0; TR(1,2) TR(2,2) 0; TT(1) TT(2) 1];
    tform = maketform('affine',T2D);
    I2r = imtransform(I2, tform,...
                        'XData',[1 M],...
                        'YData',[1 N]);
end