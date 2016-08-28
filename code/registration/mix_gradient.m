function Ifus = mix_gradient(I1,mask1,I2,mask2,debug)
    [N,M] = size(I1);
    % Differential operators
    sy = [N 1:N-1];
    sx = [M 1:M-1];
    grad = @(f)cat(3, f-f(sy,:), f-f(:,sx));
    ty = [2:N 1];
    tx = [2:M 1];
    div = @(v)v(ty,:,1)-v(:,:,1) + v(:,tx,2)-v(:,:,2);
    delta = @(f)div(grad(f));

    sigma = [2.5 3 3.5]; % detect details of different sizes
    thresh = 0.75; % threshold on the correlation values
    niterSobInpaint = 150;
    erode = 60;
    I1clean = remove_vessels(I1,sigma,thresh,niterSobInpaint,false);
    I1cleand = medfilt2(I1clean,[5,5]); % smooth image
    % Set all contour pixel to 0 (to make sure object do not touch a border)
    I1cleand(1,:) = 0; I1cleand(N,:) = 0; I1cleand(:,1) = 0; I1cleand(:,N) = 0;
    % I1cleanEnh = double(int16(adapthisteq(uint8(I1clean))));
    % figure
    % imagesc(I1cleanEnh);axis image;colormap gray;
    % Erode the contour of the mask
    se = strel('disk',erode);
    maskShrink1 = imerode(mask1,se);
    % Compute gradient
    g1 = grad(I1cleand).*repmat(maskShrink1,[1,1,2]);
    ng1 = sqrt(sum(g1.^2,3));
    ng1l = medfilt2(ng1,[5,5]); % smooth gradient

    I2clean = remove_vessels(I2,sigma,thresh,niterSobInpaint,false);
    I2cleand = medfilt2(I2clean,[5,5]);
    % Set all contour pixel to 0 (to make sure object do not touch a border)
    I2cleand(1,:) = 0; I2cleand(N,:) = 0; I2cleand(:,1) = 0; I2cleand(:,N) = 0;
    % Create a mask to erode the contour of the slice
    se = strel('disk',erode);
    maskShrink2 = imerode(mask2,se);
    % Compute gradient
    g2 = grad(I2cleand).*repmat(maskShrink2,[1,1,2]);
    ng2 = sqrt(sum(g2.^2,3));
    ng2l = medfilt2(ng2,[5,5]);

    if debug
        figure
        subaxis(2,2,1); imagesc(I1cleand); colormap gray; axis image
        subaxis(2,2,2); imagesc(I2cleand); colormap gray; axis image
        subaxis(2,2,3); imagesc(ng1l); colormap gray; axis image
        subaxis(2,2,4); imagesc(ng2l); colormap gray; axis image
    end
    % Mix gradients
    tmpfun = @(x) atan(x)/pi + 1/2;
    %w = repmat(tmpfun(ng1l-ng2l),[1 1 2]); % weights
    w = repmat(ng1l>ng2l, [1 1 2]);
%     figure
%     imagesc(w(:,:,1));axis image;colormap gray
    mg = w.*g1 + (1-w).*g2; % We want a smooth filter
    Ifus = poisson_solver_function(mg(:,:,2),mg(:,:,1),zeros(N,M));
    % figure
    % subaxis(1,2,1); imagesc(mg(:,:,1)); colormap gray; axis image
    % subaxis(1,2,2); imagesc(mg(:,:,2)); colormap gray; axis image
    if debug
        figure;
        subaxis(2,2,1); imagesc(I1clean); colormap gray; axis image
        subaxis(2,2,2); imagesc(I2clean); colormap gray; axis image
        subaxis(2,2,3);imagesc(Ifus);axis image;colormap gray;colorbar;title('Fusion');
    end
end