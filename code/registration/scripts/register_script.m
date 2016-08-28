cd '../Stage_2013/Code/Registration'
%%
addpath(genpath('../lib'))
clear all;
close all;

getd = @(p)path(p,path);
getd('../../../Matlab/toolbox_signal/');
getd('../../../Matlab/toolbox_general/');

%% 0) Load and display images
datafolder = 'VSD_20130306';
%pathname = '../../Data/Histology/HistologyEugenia/VSD_20130306/';
pathname = fullfile('../../Data/Histology/Histo',datafolder);
%pathname = fullfile('../../Data/Histology/Reconstruction~ok',datafolder);
fimages = dir(fullfile(pathname,'*.jpg'));
fimages = sort_nat({fimages.name});
subset = 1:8; % we do not consider the last layers that may not contain the barrels of interest 
fimages = fimages(subset);
L = size(fimages,2); % number of layers

tmp = fullfile(pathname,fimages{1});
tmp = imread(tmp);
[N,M,~] = size(tmp);
images = zeros(N,M,L,'double');
for i=1:L
   tmp = fullfile(pathname,fimages{i});
   tmp = imread(tmp);
   images(:,:,i) = rgb2gray(tmp);
end

figure
for i=1:L
    subaxis(2,4,i,'Spacing',0.01);
    imagesc(images(:,:,i)); axis image; colormap gray
end

if false
    barrelmap = imread(fullfile(pathname,'20130201_barrelmap.png'));
    figure
    imagesc(barrelmap);axis image;colormap gray
end

% figure
% for i=1:L
%     subaxis(2,4,i,'Spacing',0.01);
%     tmp = images(:,:,i);
%     hist(tmp(:),64);
% end

%% 1) Segment and extract cortex slices
imagesSeg = images;
imagesBin = images;
nReg = 2; % number of regions, typically 2 or 3
for i=1:L
    fprintf([num2str(i) '\n']);
    [imagesSeg(:,:,i),imagesBin(:,:,i)] = segment_cortex(images(:,:,i),nReg);
end

figure
for i=1:L
    subaxis(2,4,i,'Spacing',0.01);
    imagesc(imagesSeg(:,:,i)); axis image; colormap gray
end

%% 2) Go for registration !
close all
imagesReg = imagesSeg;
preciseReg = false;
for i = 4:L-1 % the first 3 layers are sometimes hard to register
    %%
    %close all
    i1 = i;
    i2 = i+1;
    I1 = imagesReg(:,:,i1);
    I2 = imagesSeg(:,:,i2);

    % 1) Detect vertical vessels with normalized cross correlation
    sigma = [2.5 3 3.5 4 4.5]; % detect details of different sizes
    thresh = 0.75; % threshold on the correlation values [#USR_INPUT]
    [c1,corr1,radius1] = detect_perpendicular_vessels(I1,sigma,thresh,false);
    [c2,corr2,radius2]  = detect_perpendicular_vessels(I2,sigma,thresh,false);
    % Sort by increasing correlations
    [~,idx1] = sort(corr1,'descend');
    c1 = c1(idx1,:);
    [~,idx2] = sort(corr2,'descend');
    c2 = c2(idx2,:);
    % Keep the best vessels
    nmax = 45;
    nc = min([nmax,size(c1,1),size(c2,1)]); % number of vessels kept
    c1 = c1(1:nc,:);
    radius1 = radius1(idx1(1:nc));
    c2 = c2(1:nc,:);
    radius2 = radius2(idx2(1:nc));

    % 1) Pre-registration
    % 1.1) Register the centroids of the two slices
    % 1.2) Coarse registration: try different initializations (light 
    % rotations) followed by fast ICP and pick the best result
    thetas = [-20:5:20]; % angles in gradians [#USR_INPUT]
    [I1pre,c1pre,I2pre,c2pre] = pre_register(I1,c1,I2,c2,thetas,false);
  
    % Display results
%     patchAlphaSliderGUI(I1,I2);
%     hold on;
%     plot(c1(:,1), c1(:,2), 'r.');
%     hold on;
%     plot(c2(:,1), c2(:,2), 'g+');
%     title(['L' num2str(i1) '-' num2str(i2) ' (original)']);
    
    patchAlphaSliderGUI(I1pre,I2pre);
    hold on;
    plot(c1pre(:,1), c1pre(:,2), 'r.');
    hold on;
    plot(c2pre(:,1), c2pre(:,2), 'g+');
    title(['L' num2str(i1) '-' num2str(i2) ' (registered)']);

   	% 2) Precise registration
    if preciseReg
        r0 = 10; % initial value for the parameter of the robust function
        [TR, TT, E] = robust_icp(c1pre', c2pre',...
                    @(d,r) r^2/6*(1-(1-(d/r).^2).^3).*(abs(d)<=r) + r^2/6.*(abs(d)>r),...
                    @(d,r) ((1-(d/r).^2).^2).*(abs(d)<=r),...
                    @(r) r,...
                    r0,5,5,I1pre,I2pre,false);
        c2r = [TR * c2pre' + repmat(TT, 1, size(c2pre,1))]';
        T2D = [TR(1,1) TR(2,1) 0; TR(1,2) TR(2,2) 0; TT(1) TT(2) 1];
        tform = maketform('affine',T2D);
        I2r = imtransform(I2pre, tform,...
                                        'XData',[1 M],...
                                        'YData',[1 N]);

        patchAlphaSliderGUI(I1pre,I2r);
        hold on;
        plot(c1pre(:,1), c1pre(:,2), 'r.');
        hold on;
        plot(c2r(:,1), c2r(:,2), 'g+');

        figure
        plot(E);
        title(['L' num2str(i1) '-' num2str(i2) ' : NRJ']);
    else
        I2r = I2pre;
    end
    
    imagesReg(:,:,i+1) = I2r;
end

%% Visualize registered images
figure
for i=1:L
    subaxis(2,4,i,'Spacing',0.01);
    imagesc(imagesReg(:,:,i)); axis image; colormap gray
end

%% Aggregate barrels (Gradient domain fusion)
istart = 4; % first slice where barrels appear
Ifusion = imagesReg(:,:,istart);
%figure
%count = 1;
for i=istart+1:8
   Ifusion = mix_gradient(Ifusion,imagesReg(:,:,i),false);
   %subaxis(1,3,count); imagesc(Ifusion);colormap gray; axis image
   %count = count+1;
end

barrelMap(end+1).histoRec = Ifusion;
barrelMap(end).name = datafolder;

%save('VSD_20130212.mat','imagesReg','Ifusion','barrelmap');
%load('VSD_20130219.mat');
%% Register with barrel map template
%patchAlphaSliderGUI(Ifusion,barrelmap);

figure
for i=1:size(barrelMap,2)
    subaxis(3,5,i); imagesc(barrelMap(i).histoRec); title(barrelMap(i).name);
    axis image; colormap gray
end

%% Lower drift by low-pass filtering
I = Ifusion19;
sigma = 10;
sz = ceil((sigma*6 + 1)/2)*2 - 1;
f = fspecial('gaussian',[sz sz], sigma);
figure
subaxis(2,2,1); imagesc(I); axis image; colormap gray
tmp = imfilter(I,f);
subaxis(2,2,2); imagesc(tmp); axis image; colormap gray
tmp1 = medfilt2(I - tmp,[10 10]);
subaxis(2,2,3); imagesc(tmp1); axis image; colormap gray
subaxis(2,2,4); imagesc(adapthisteq(tmp1)); axis image; colormap gray
%subaxis(2,2,4); imagesc(edge((tmp1),'canny',[])); axis image; colormap gray

%% Segmentation of the barrels
figure
subaxis(2,2,1); imagesc(tmp1); axis image; colormap gray
tmp2 = reshape(kmeans(double(tmp1(:)),2,'distance','sqEuclidean',...
                        'Replicates',3,'EmptyAction','drop'),N,M);
subaxis(2,2,2); imagesc(tmp2); axis image; colormap gray
tmp3 = reshape(kmeans(double(tmp1(:)),3,'distance','sqEuclidean',...
                        'Replicates',3,'EmptyAction','drop'),N,M);
subaxis(2,2,3); imagesc(tmp3); axis image; colormap gray