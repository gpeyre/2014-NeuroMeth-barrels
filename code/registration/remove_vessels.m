function I1clean = remove_vessels(I1, sigma, thresh, niter, debug)

[N,M] = size(I1);
sy = [N 1:N-1];
sx = [M 1:M-1];
grad = @(f)cat(3, f-f(sy,:), f-f(:,sx));
ty = [2:N 1];
tx = [2:M 1];
div = @(v)v(ty,:,1)-v(:,:,1) + v(:,tx,2)-v(:,:,2);
delta = @(f)div(grad(f));

    I1 = double(I1);

    % 0) Detect vessels (black and white)
    [c1,corr1,radius1] = detect_perpendicular_vessels(I1,sigma,thresh,true);
    maskVessels1 = ones(N,M);
    c1 = round(c1);
    if debug
        figure
        title('Vessels detection');
        imagesc(I1); axis image; colormap gray;
        hold on;
        plot(c1(:,1), c1(:,2), 'r+');
    end
    for k = 1 : size(c1, 1),
        %r = ceil((ceil((radius1(k)*6 + 1)/2)*2 - 1)/2)
        r = ceil(radius1(k)*2);
        maskVessels1(max(c1(k,2)-r,1):min(c1(k,2)+r,N),...
            max(c1(k,1)-r,1):min(c1(k,1)+r,M)) = 0;
        if debug
         DrawCircle(c1(k,1), c1(k,2), r*2, 32, 'b-');
        end
    end
    hold off;
    
    % 1) Remove vessels by Sobolev inpainting
    Phi = @(f,mask)f.*mask + 0.5.*~mask;
    y = Phi(I1,maskVessels1);
    Pi = @(f,mask)f.*(1-mask) + y.*mask;
    tau = 0.99/4; % for convergence tau must be < 1/4
    E = zeros(niter,1);
    I1clean = y;
    % Projected gradient descent
    for i=1:niter
        gradSob = delta(I1clean);
        I1clean = Pi(I1clean + tau*gradSob,maskVessels1);
        E(i) = sqrt(sum(gradSob(:).^2));
    end
    if debug
    %     figure
    %     plot(E);
        figure
        subaxis(1,2,1)
        imagesc(I1);axis image; colormap gray
        title('Before removing vessels');
        subaxis(1,2,2)
        imagesc(I1clean);axis image; colormap gray
        title('After removing vessels');
    end
end