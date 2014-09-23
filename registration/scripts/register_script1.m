cd '../Stage_2013/Code/Registration'
%%
addpath(genpath('../lib'))
clear all;
close all;

getd = @(p)path(p,path);
getd('../../../Matlab/toolbox_signal/');
getd('../../../Matlab/toolbox_general/');

datas_recok = ['VSD_20130201';'VSD_20130212';...
                'VSD_20130214';'VSD_20130219'];
            
datas_3regions = ['VSD_20121107';'VSD_20121109';...
                'VSD_20130206';...
                'VSD_20130227';'VSD_20130228'];

datas_2regions = ['VSD_20121128';'VSD_20121218';'VSD_20121219';...
                'VSD_20121220';'VSD_20130129';'VSD_20130130';...
                'VSD_20130304';'VSD_20130305';'VSD_20130306'
                'VSD_20130304';'VSD_20130305'];            

datas = datas_3regions;

%%
for k=3:size(datas,1)
    close all
    %% 0) Load and display images
    datafolder = datas(k,:);
    disp(['Processing ',datafolder,'...']);
    
    %pathname = '../../Data/Histology/HistologyEugenia/VSD_20130306/';
    pathname = fullfile('../../Data/Histology/Histo',datafolder);
    %pathname = fullfile('../../Data/Histology/Reconstruction~ok',datafolder);
    fimages = dir(fullfile(pathname,'*.jpg'));
    fimages = sort_nat({fimages.name});
    subset = 1:min(size(fimages,2),8); % we do not consider the last layers that may not contain the barrels of interest 
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

    %% 1) Segment and extract cortex slices
    disp(['Segmenting slices ...']);
    imagesSeg = images;
    imagesBin = images;
    nReg = 3; % number of regions, typically 2 or 3
    tic
    for i=1:L
        fprintf([num2str(i) '\n']);
        [imagesSeg(:,:,i),imagesBin(:,:,i)] = segment_cortex(images(:,:,i),nReg);
    end
    
    figure
    for i=1:L
        subaxis(2,4,i,'Spacing',0.01);
        imagesc(imagesSeg(:,:,i)); axis image; colormap gray
    end
    
    disp(['Elapsed time : ',num2str(toc)]);

%     a = input('Ready to continue (y/n)? ','s');
%     if strcmpi(a,'y')
%         close all
%     else
%         return
%     end
    
    %% 2) Go for registration !
    imagesReg = imagesSeg;
    preciseReg = false;
    for i = 3:L-1 % the first 3 layers are often hard to register
        %%
        tic
        %close all
        i1 = i;
        i2 = i+1;
        I1 = imagesReg(:,:,i1);
        I2 = imagesSeg(:,:,i2);
        disp(['Registering slices ',num2str(i1),' and ',num2str(i2),' ...']);

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
        nmax = 30;
        nc = min([nmax,size(c1,1),size(c2,1)]); % number of vessels kept
        c1 = c1(1:nc,:);
        radius1 = radius1(idx1(1:nc));
        c2 = c2(1:nc,:);
        radius2 = radius2(idx2(1:nc));
        if false
            % Keep vessels inside the given bounding box
            xmin = 400; xmax = 800; ymin = 250; ymax = 650;
            c1 = c1((c1(:,1)>xmin & c1(:,1)<xmax) & (c1(:,2)>ymin & c1(:,2)<ymax),:) ;
            c2 = c2((c2(:,1)>xmin & c2(:,1)<xmax) & (c2(:,2)>ymin & c2(:,2)<ymax),:) ;
        end
        
        if false
            patchAlphaSliderGUI(I1,I2);
            hold on;
            plot(c1(:,1), c1(:,2), 'r.');
            hold on;
            plot(c2(:,1), c2(:,2), 'g+');
            title(['L' num2str(i1) '-' num2str(i2) ' (original)']);
        end
        
        % 1) Pre-registration
        % 1.1) Register the centroids of the two slices
        % 1.2) Coarse registration: try different initializations (light 
        % rotations) followed by fast ICP and pick the best result
        thetas = [-20:5:20]; % angles in gradians [#USR_INPUT]
        [I1pre,c1pre,I2pre,c2pre] = pre_register(I1,c1,I2,c2,true,thetas,false);

        % Display results
        patchAlphaSliderGUI(I1pre,I2pre);
        hold on;
        plot(c1pre(:,1), c1pre(:,2), 'r.');
        hold on;
        plot(c2pre(:,1), c2pre(:,2), 'g+');
        title(['L' num2str(i1) '-' num2str(i2) ' (registered)']);

        I2r = I2pre;

        imagesReg(:,:,i+1) = I2r;
        disp(['Elapsed time : ',num2str(toc)]);
        
    end

    a = input('Ready to continue (y/n)? ','s');
    if strcmpi(a,'y')
        
    else
        return
    end
    

    
    %% Aggregate barrels (Gradient domain fusion)
    disp(['Aggregating slices ...']);
    istart = 4; % first slice where barrels appear
    Ifusion = imagesReg(:,:,istart);
    tic
    for i=istart+1:L
       Ifusion = mix_gradient(Ifusion,imagesReg(:,:,i),false);
    end
    disp(['Elapsed time : ',num2str(toc)]);

    if exist('barrelMap')
        barrelMap(end+1).histoRec = Ifusion;
        barrelMap(end).imagesRec = imagesReg;
        barrelMap(end).name = datafolder;
    else
        barrelMap.histoRec = Ifusion;
        barrelMap(end).imagesRec = imagesReg;
        barrelMap.name = datafolder;
    end

    save('barrelMapRecok.mat','barrelMap');

    %% Lower drift by low-pass filtering
    sigma = 40;
    sz = ceil((sigma*6 + 1)/2)*2 - 1;
    f = fspecial('gaussian',[sz sz], sigma);
    figure
    subaxis(2,2,1); imagesc(Ifusion); axis image; colormap gray
    tmp = imfilter(Ifusion,f);
    subaxis(2,2,2); imagesc(tmp); axis image; colormap gray
    tmp1 = medfilt2(Ifusion - tmp,[10 10]);
    subaxis(2,2,3); imagesc(tmp1); axis image; colormap gray
    %subaxis(2,2,4); imagesc(adapthisteq(tmp1)); axis image; colormap gray
    subaxis(2,2,4); imagesc(edge((tmp1),'canny',[])); axis image; colormap gray

end


%% Register with barrel map template
clear all
load barrelMap3regions;
[N,M,~] = size(barrelMap(1).histoRec);

BMAP = zeros(N,M,size(barrelMap,2));
for i=1:size(barrelMap,2)
    BMAP(:,:,i) = mat2gray(barrelMap(i).histoRec);
end

figure
for i=1:size(BMAP,3)
    subaxis(3,4,i); imagesc(BMAP(:,:,i));
    axis image; colormap gray
end

%% Preprocessing
% lower drift by low-pass filtering
BMAPPOST = zeros(N,M,size(barrelMap,2));
sigma = 10;
sz = ceil((sigma*6 + 1)/2)*2 - 1;
f = fspecial('gaussian',[sz sz], sigma);
for i=1:size(barrelMap,2)
    tmp = imfilter(barrelMap(i).histoRec,f);
    BMAPPOST(:,:,i) = mat2gray(medfilt2(barrelMap(i).histoRec - tmp,[10 10]));
end

figure
for i=1:size(BMAPPOST,3)
    subaxis(3,4,i); imagesc(BMAPPOST(:,:,i));
    axis image; colormap gray
end

%
BMAPPOST2 = zeros(N,M,size(barrelMap,2));
for i=1:size(BMAPPOST2,3)
    BMAPPOST2(:,:,i) = imhmin(BMAPPOST(:,:,i),0.1);
end

figure
for i=1:size(BMAPPOST2,3)
    subaxis(3,4,i); imagesc(BMAPPOST2(:,:,i));
    axis image; colormap gray
end


%% 
I1 = BMAPPOST(:,:,6);
figure
BW1 = roipoly(I1);
I1(~BW1) = 0.5;
%figure;imagesc(I1);axis image; colormap gray
I2 = BMAPPOST(:,:,8);
BW2 = roipoly(I2);
I2(~BW2) = 0.6;
%figure;imagesc(I2);axis image; colormap gray

%% Segmentation with KMEANS
I11 = imhmin(I1,0.05);
I22 = imhmin(I2,0.05);
I11seg = reshape(kmeans(I11(:),2,'distance','sqEuclidean',...
                        'Replicates',3,'start','uniform'),N,M);
I22seg = reshape(kmeans(I22(:),2,'distance','sqEuclidean',...
    'Replicates',3,'start','uniform'),N,M);
figure;
subaxis(2,2,1);imagesc(I11);axis image; colormap gray
subaxis(2,2,2);imagesc(I11seg);axis image; colormap gray
subaxis(2,2,3);imagesc(I22);axis image; colormap gray
subaxis(2,2,4);imagesc(I22seg);axis image; colormap gray


%%
se = strel('disk', 6);
I11seg1 = imdilate(imerode(imfill(I11seg,'holes'),se),se);
figure;imshowpair(I11,I11seg1,'montage');
se = strel('disk', 6);
I22seg1 = imdilate(imerode(imfill(I22seg,'holes'),se),se);
figure;imshowpair(I22,I22seg1,'montage');
%%
I3 = imhmin(I1,0.2); %20 is the height threshold for suppressing shallow minima
figure;imagesc(I3);axis image;colormap gray
L = watershed(I3);
figure
imagesc(L);axis image;


%% Try different registration methods...
I1 = BMAPPOST(:,:,1);
I2 = BMAPPOST(:,:,6);
[optimizer,metric] = imregconfig('multimodal');
I2reg = imregister(I2, I1, 'rigid', optimizer, metric);
%figure, imshowpair(I1,I2);
figure, imshowpair(I1,I2reg);

figure, imshowpair(BMAPPOST(:,:,2), BMAPPOST(:,:,4), 'montage')
figure, imshowpair(BMAP(:,:,1),BMAP(:,:,6), 'montage')