load barrelMap3regions;
[N,M,~] = size(barrelMap(1).histoRec);

figure
for i=1:size(barrelMap,2)
    subaxis(3,4,i); imagesc(barrelMap(i).histoRec);
    axis image; colormap gray
end

%%
%clear HistoRec
for i=1:size(barrelMap,2)
    %%
    S = barrelMap(i).imagesRec;
    S(1,:,:) = 0; % in case, for the mix_gradient to erode the mask
    figure
    for j=1:size(S,3)
        S(:,:,j) = mat2gray(S(:,:,j));
        subaxis(3,3,j); imagesc(S(:,:,j)); axis image; colormap gray
    end
    istart = 3;
    HistoRec(:,:,i) = reconstruct_barrels1(S(:,:,istart:end),false);
    figure, imagesc(tmp(:,:,end)); axis image; colormap gray
end

%%
for i=1:size(barrelMap,2)
    %%
    tmp = mat2gray(barrelMap(i).histoRec);
    X = tmp(200:end,:);
    [N,M,~] = size(tmp);%size(X);
    X1 = mat2gray(lower_drift(X,10,false));
    X2 = adapthisteq(X1,'NumTiles',[20 20],'ClipLimit',0.05);
    X3 = reshape(kmeans(X2(:),2,'distance','sqEuclidean','Replicates',3,'start','uniform'),M,N);
    figure
    subaxis(2,2,1);imagesc(X);axis image; colormap gray
    subaxis(2,2,2);imagesc(X1);axis image; colormap gray
    subaxis(2,2,3);imagesc(X2);axis image; colormap gray
    subaxis(2,2,4);imagesc(X3);axis image; colormap gray
end

%%
idx = 2;
imagesReg = barrelMap(idx).imagesRec;
[N,M,~] = size(imagesReg(:,:,1));

for i=1:size(imagesReg,3)
    subaxis(3,4,i); imagesc(imagesReg(:,:,i)); axis image; colormap gray
end

figure
for i=1:size(imagesReg,3)
    tmp = mat2gray(lower_drift(adapthisteq(imagesReg(:,:,i)),20,false));
    tmp(imagesReg(:,:,i) == 0) = 0;
    tmp(1,:) = 0;
    Itmp(:,:,i) = tmp;
    subaxis(3,4,i); imagesc(Itmp(:,:,i));
    axis image; colormap gray
end

figure
for i=1:size(imagesReg,3)
    subaxis(3,4,i); imagesc((Itmp(:,:,i)));
    axis image; colormap gray
end

figure
for i=1:size(imagesReg,3)
    subaxis(3,4,i); imagesc(imagesReg(:,:,i));
    axis image; colormap gray
end
%% Aggregate barrels (Gradient domain fusion)
istart = 5; % first slice where barrels appear
tic
clear Ifusiontmp
Ifus(:,:,1) = I3(:,:,1);
figure
for i=2:size(I3,3)
   Ifus(:,:,i) = mat2gray(mix_gradient(Ifus(:,:,i-1),I3(:,:,i),false));
   subaxis(2,4,i-1);imagesc(Ifus(:,:,i));axis image; colormap gray
end
disp(['Elapsed time : ',num2str(toc)]);

cpt = 1;
figure
for i=istart+1:size(Itmp,3)
   subaxis(2,3,cpt);imagesc(Ifusiontmp(:,:,cpt+1));axis image; colormap gray
   cpt = cpt+1;
end
%%
istart = 5; % first slice where barrels appear
tic
cpt = 1;
clear Ifusion
Ifusion(:,:,1) = imagesReg(:,:,istart);
figure
for i=istart+1:size(Itmp,3)
   Ifusion(:,:,cpt+1) = mat2gray(mix_gradient(Ifusion(:,:,cpt),imagesReg(:,:,i),false));
   subaxis(2,3,cpt);imagesc(Ifusion(:,:,cpt+1));axis image; colormap gray
   cpt = cpt+1;
end
disp(['Elapsed time : ',num2str(toc)]);

cpt = 1;
figure
for i=istart+1:size(Itmp,3)
   subaxis(2,3,cpt);imagesc(Ifusion(:,:,cpt+1));axis image; colormap gray
   cpt = cpt+1;
end
%%
tmp1 = imagesReg(:,:,6);

tmp2 = adapthisteq(tmp1);

tmp3 = mat2gray(lower_drift(tmp2,10,false));
tmp4 = mat2gray(lower_drift(tmp2,20,false));
tmp5 = mat2gray(lower_drift(tmp2,30,false));
tmp6 = mat2gray(lower_drift(tmp2,40,false));

figure
subaxis(2,3,1);imagesc(tmp1);axis image; colormap gray
subaxis(2,3,2);imagesc(tmp2);axis image; colormap gray
subaxis(2,3,3);imagesc(tmp3);axis image; colormap gray
subaxis(2,3,4);imagesc(tmp4);axis image; colormap gray
subaxis(2,3,5);imagesc(tmp5);axis image; colormap gray
subaxis(2,3,6);imagesc(tmp6);axis image; colormap gray

%%

figure
for i=5:size(imagesReg,3)
    tmp = mat2gray(adapthisteq(imagesReg(:,:,i)));
    tmp(imagesReg(:,:,i) == 0) = 0;
    tmp(1,:) = 0;
   	tmp1(:,:,i) = tmp;
    wSize= 5;sigma = [3, 0.1];
    tmp2(:,:,i) = medfilt2(imagesReg(:,:,i),[5,5]);
    subaxis(3,4,i-4); imagesc(imagesReg(:,:,i));axis image; colormap gray
    subaxis(3,4,i); imagesc(tmp1(:,:,i));axis image; colormap gray
    subaxis(3,4,i+4); imagesc(tmp2(:,:,i));axis image; colormap gray
end

%%
tmp1 = imagesReg(:,:,8);
tmp2 = mat2gray(lower_drift(tmp1,20,false));tmp2(tmp1 == 0) = 0;
tmp3 = adapthisteq(tmp1,'NumTiles',[20 20],'ClipLimit',0.005);
tmp4 = mat2gray(lower_drift(tmp3,20,false));tmp4(tmp1 == 0) = 0;
figure
subaxis(2,3,1);imagesc(tmp1);axis image; colormap gray
subaxis(2,3,2);imagesc(tmp2);axis image; colormap gray
subaxis(2,3,3);imagesc(tmp3);axis image; colormap gray
subaxis(2,3,4);imhist(tmp1(:));
subaxis(2,3,5);imhist(tmp2(:));
subaxis(2,3,6);imhist(tmp3(:));

%%
istart = 5;
clear I1 I2 I3 I4
cpt = 1;
for i=istart:size(imagesReg,3)
    I1(:,:,cpt) = imagesReg(:,:,i);
    tmp2 = mat2gray(lower_drift(I1(:,:,cpt),20,false));tmp2(I1(:,:,cpt) == 0) = 0; I2(:,:,cpt) = tmp2;
    tmp3 = adapthisteq(I1(:,:,cpt),'NumTiles',[20 20],'ClipLimit',0.005); tmp3(I1(:,:,cpt) == 0) = 0; I3(:,:,cpt) = tmp3;
    tmp4 = mat2gray(lower_drift(I3(:,:,cpt),20,false));tmp4(I1(:,:,cpt) == 0) = 0; I4(:,:,cpt) = tmp4;
    cpt = cpt+1;
end

%%
clear Ifus1 Ifus2 Ifus3 Ifus4;
Ifus1(:,:,1) = I1(:,:,1);
Ifus2(:,:,1) = I2(:,:,1);
Ifus3(:,:,1) = I3(:,:,1);
Ifus4(:,:,1) = I4(:,:,1);
for i=2:size(I1,3)
    disp(num2str(i-1));
    Ifus1(:,:,i) = mat2gray(mix_gradient(Ifus1(:,:,i-1),I1(:,:,i),false));
	Ifus2(:,:,i) = mat2gray(mix_gradient(Ifus2(:,:,i-1),I2(:,:,i),false));
    Ifus3(:,:,i) = mat2gray(mix_gradient(Ifus3(:,:,i-1),I3(:,:,i),false));
    Ifus4(:,:,i) = mat2gray(mix_gradient(Ifus4(:,:,i-1),I4(:,:,i),false));
end

figure
for i=1:size(Ifus1,3)
   subaxis(4,4,i);imagesc(Ifus1(:,:,i));axis image; colormap gray
   subaxis(4,4,i+4);imagesc(Ifus2(:,:,i));axis image; colormap gray
   subaxis(4,4,i+8);imagesc(Ifus3(:,:,i));axis image; colormap gray
   subaxis(4,4,i+12);imagesc(Ifus4(:,:,i));axis image; colormap gray
end

%% Non uniform background removal
X1 = Ifus1(:,:,end);
X2 = Ifus2(:,:,end);
X3 = Ifus3(:,:,end);
X4 = Ifus4(:,:,end);
X5 = reshape(kmeans(X1(:),2,'distance','sqEuclidean','Replicates',3,'start','uniform'),N,M);
X6 = reshape(kmeans(X2(:),2,'distance','sqEuclidean','Replicates',3,'start','uniform'),N,M);
X7 = reshape(kmeans(X3(:),2,'distance','sqEuclidean','Replicates',3,'start','uniform'),N,M);
X8 = reshape(kmeans(X4(:),2,'distance','sqEuclidean','Replicates',3,'start','uniform'),N,M);
figure
subaxis(2,2,1);imagesc(X1);axis image; colormap gray
subaxis(2,2,2);imagesc(X2);axis image; colormap gray
subaxis(2,2,3);imagesc(X3);axis image; colormap gray
subaxis(2,2,4);imagesc(X4);axis image; colormap gray
figure
subaxis(2,2,1);imagesc(X5);axis image; colormap gray
subaxis(2,2,2);imagesc(X6);axis image; colormap gray
subaxis(2,2,3);imagesc(X7);axis image; colormap gray
subaxis(2,2,4);imagesc(X8);axis image; colormap gray
%%
X = mat2gray(lower_drift(Ifus1(:,:,end),10,false));
X1 = adapthisteq(X,'NumTiles',[20 20],'ClipLimit',0.05);
X2 = adapthisteq(X,'NumTiles',[20 20],'ClipLimit',0.01);
X3 = adapthisteq(X,'NumTiles',[20 20],'ClipLimit',0.02);
X4 = adapthisteq(X,'NumTiles',[20 20],'ClipLimit',0.03);
X5 = adapthisteq(X,'NumTiles',[10 10],'ClipLimit',0.03);
X6 = reshape(kmeans(X1(:),2,'distance','sqEuclidean','Replicates',3,'start','uniform'),N,M);
X7 = reshape(kmeans(X2(:),2,'distance','sqEuclidean','Replicates',3,'start','uniform'),N,M);
X8 = reshape(kmeans(X3(:),2,'distance','sqEuclidean','Replicates',3,'start','uniform'),N,M);
X9 = reshape(kmeans(X4(:),2,'distance','sqEuclidean','Replicates',3,'start','uniform'),N,M);
X10 = reshape(kmeans(X5(:),2,'distance','sqEuclidean','Replicates',3,'start','uniform'),N,M);
figure
subaxis(2,3,1);imagesc(X1);axis image; colormap gray
subaxis(2,3,2);imagesc(X2);axis image; colormap gray
subaxis(2,3,3);imagesc(X3);axis image; colormap gray
subaxis(2,3,4);imagesc(X4);axis image; colormap gray
subaxis(2,3,5);imagesc(X5);axis image; colormap gray
subaxis(2,3,6);imagesc(X);axis image; colormap gray
figure
subaxis(2,3,1);imagesc(X6);axis image; colormap gray
subaxis(2,3,2);imagesc(X7);axis image; colormap gray
subaxis(2,3,3);imagesc(X8);axis image; colormap gray
subaxis(2,3,4);imagesc(X9);axis image; colormap gray
subaxis(2,3,5);imagesc(X10);axis image; colormap gray
subaxis(2,3,6);imagesc(X);axis image; colormap gray
%%
% X1 = adapthisteq(X,'NumTiles',[20 20],'ClipLimit',0.03);
% X2 = lower_drift(X1,10,false);
% X3 = imadjust(imtophat(X1,strel('disk',15)));
% X4 = imadjust(imtophat(X1,strel('disk',50)));
% X5 = imadjust(imtophat(X1,strel('disk',200)));
% X6 = mat2gray(X1-X3);
% X7 = mat2gray(X1-X5);
% X8 = mat2gray(X5+X3);

% X7 = imhmin(X6,0.05);
X9 = reshape(kmeans(X1(:),2,'distance','sqEuclidean','Replicates',3,'start','uniform'),N,M);
figure
subaxis(2,3,1);imagesc(X1);axis image; colormap gray
subaxis(2,3,2);imagesc(X3);axis image; colormap gray
subaxis(2,3,3);imagesc(X5);axis image; colormap gray
subaxis(2,3,4);imagesc(X6);axis image; colormap gray
subaxis(2,3,5);imagesc(X8);axis image; colormap gray
subaxis(2,3,6);imagesc(X9);axis image; colormap gray

%%
X1 = imopen(X,strel('disk',5));
X2 = imadjust(X-X1);


Itmp = imadjust(imtophat(I,strel('disk',25)));
figure,imshowpair(I,Itmp,'montage');

level = graythresh(I3);
bw = im2bw(I3,level);
bw = bwareaopen(bw, 50);
imshow(bw)

% Display the Background Approximation as a Surface
figure, surf(double(Iback(1:8:end,1:8:end))),zlim([0 1]);
set(gca,'ydir','reverse');

figure, imagesc(mat2gray(mix_gradient(S(:,:,6),S(:,:,7),false)));axis image; colormap gray

%% Empiler les sections recalées pour mieux identifier les vaisseaux
tmp = imagesReg;
for j=1:length(imagesReg)
    tmp{j} = mat2gray(imagesReg{j});
end
%%
stack = tmp{1};
[M,N] = size(stack);

se = strel('disk', 5);
for j=2:length(tmp)
    mask = zeros(M,N);
    mask(stack ~= 0) = 1;
    mask = imerode(mask,se);
    stack = stack.*mask + tmp{j}.*~mask;
end

figure,imshow(stack);