function [tvmask, hypomask, nmask, img, dat] = cellSegCFos(filename)
%CELLSEG_RNASCOPE Summary of this function goes here
% given a particular czi file, read it and its meta data

%% Import data & get x,y,z volume matrices by signal 
% filename = '/dcl01/lieber/ajaffe/Keri/Imaging/cFos/Stitched_for_Badoi/e1_KI_8150/KI_8150_Slide1_8_1_Stitch.czi';
data = bfOpen3DVolume(filename);
omeMeta = data{1, 4};
metadata = data{1, 2};

%% OME data 
L = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
M = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
N = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
C = omeMeta.getPixelsSizeC(0).getValue(); % number of channels

voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0)...
    .value(ome.units.UNITS.MICROM).doubleValue(); 
voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0)...
    .value(ome.units.UNITS.MICROM).doubleValue(); 
voxelSizeZ = omeMeta.getPixelsPhysicalSizeZ(0)...
    .value(ome.units.UNITS.MICROM).doubleValue();

%% Get channel information
nameC1 = metadata.get(['Global Experiment|AcquisitionBlock|',...
    'MultiTrackSetup|TrackSetup|Detector|Folder #1']);
nameC2 = metadata.get(['Global Experiment|AcquisitionBlock|',...
    'MultiTrackSetup|TrackSetup|Detector|Folder #2']);
nameC3 = metadata.get(['Global Experiment|AcquisitionBlock|',...
    'MultiTrackSetup|TrackSetup|Detector|Folder #3']);
if(C==4) % if there's a tDtomato channel
    nameC4 = metadata.get(['Global Experiment|AcquisitionBlock|',...
        'MultiTrackSetup|TrackSetup|Detector|Folder #4']);
    cName = {nameC1,nameC2,nameC3,nameC4};
else
    cName = {nameC1,nameC2,nameC3};
end

%% Load in images 
out = ReadImage6D(filename);
image6d = out{1}; %dims = series,time, z, c, x, y 
img.c1 = permute(squeeze(image6d(:,:,:,find(strcmp(cName,'Cy5')),:,:)),[2,3,1]); % oxytocin
img.c2 = permute(squeeze(image6d(:,:,:,find(strcmp(cName,'EGFP')),:,:)),[2,3,1]); % cFos 
img.nuc = permute(squeeze(image6d(:,:,:,find(strcmp(cName,'DAPI')),:,:)),[2,3,1]); % DAPI
clear out image6d data;

%% Find third ventrical & epithelial cell mask
maxPlane = squeeze(max(max(img.nuc,[],2),[],1));
nuc2 = sqrt(img.nuc./reshape(repelem(maxPlane,L*M),[L M N]));
tvmask = nuc2 > graythresh(nuc2(:))-std(nuc2(:));
tvmask = imclose(tvmask,strel('disk',10));
tvmask = imerode(tvmask,strel('disk',25));
plane = 10;
% imshow(mat2gray(tvmask(:,:,plane)));
 imshow(mat2gray(img.nuc(:,:,plane)));

%% Find PVN region of OXT(+) interneurons
maxPlane = squeeze(max(max(img.nuc,[],2),[],1));
oxt2 = sqrt(img.c1./reshape(repelem(maxPlane,L*M),[L M N]));

oxtMaxZ = mean(oxt2,3);
hypomask = mat2gray(imgaussfilt(oxtMaxZ,25));
hypomask = (hypomask > graythresh(hypomask)) .* (mean(tvmask,3)>0);
hypomask = repmat(imclose(hypomask,strel('disk',25)),[1 1 N]);
% imshow(mat2gray(oxtMaxZ));
% imshow(mat2gray(hypomask));

%% 3D nuclei segmentation by adaptive thresholding
prm.method = 'adth';	% adaptive thresholding
prm.adth.filtrad = 11;
prm.adth.th = 0.00001;
prm.splitth = 1;
prm.h=[voxelSizeX,voxelSizeY,voxelSizeZ]; %x,y,z pixel to um 
[nucBW,nmask,~,~] = cellsegm.segmct(img.nuc.*hypomask,0.002,2000,'prm',prm);
% cellbw = cellsegm.splitcells(cellbw1,prm.splitth,prmout.splitvolvox,prm.h);

% imshow(mat2gray(img.nuc(:,:,plane)));
% imshow(mat2gray(nuc2(:,:,plane)));
% imshow(label2rgb(nmask(:,:,plane)));

%% split cells using watershed algorithm
% see paper 'Segmentation and detection of fluorescent 3D spots'
conn = 6;
nucBlur = imgaussfilt3(nuc2,2.5);
dist = bwdist(~nucBW).*nucBlur;
counts = histcounts(nmask(:),max(nmask(:))+1);
[~,I] = sort(counts);
tmpmask = eq(nmask,I(1)-1);
h = max(max(max(dist.*tmpmask)));
mdist = -imhmax(dist.*nucBlur,h,conn);
mdist(~nucBW) = -Inf;
nmask2 = double(watershed(mdist,conn));

% imshow(mat2gray(nucBlur(:,:,plane)));
% imshow(label2rgb(nmask(:,:,plane)));
% imshow(label2rgb(nmask2(:,:,plane)));
% imshow(nmask2(:,:,plane)==100);

%% 3D oxytocin segmentation
c1 = img.c1;
for i = 1:N c1(:,:,i) = sqrt(c1(:,:,i)/max(max(max(c1(:,:,i))))); end

c1mask = c1 > mean2(c1) + std(c1(:));
c1mask = bwlabeln(imopen(c1mask,strel('disk',5)));

dist = bwdist(~c1mask);
counts = histcounts(c1mask(:),max(c1mask(:))+1);
[~,I] = sort(counts);
tmpmask = eq(c1mask,I(3)-1);
h = max(max(max(dist.*tmpmask)));
mdist = -imhmax(dist,h,conn);
mdist(~c1mask) = -Inf;
c1mask = double(watershed(mdist,conn));

% imshow(label2rgb(c1mask(:,:,plane)));

%% extract nuclei information
dat.cName = cName;
dat.h = prm.h;
dat.nuc = [];

for i = 3:max(unique(nmask2(:)))
    tmp = find(nmask2==i);
    [x, y, z] = ind2sub([L M N],tmp); % coordinates
    vol = length(tmp)*prod(prm.h); % volume in um^3
    rad = (vol*3/4*pi()).^(1/3); % radius in um
    c1 = mean(img.c1(tmp)); c2 = mean(img.c2(tmp));
    dat.nuc=[dat.nuc;double([mean(x),mean(y),mean(z),i,vol,rad,c1,c2])];
end

%% extract oxytocin information
dat.c1 = [];

for i = 3:max(unique(c1mask(:)))
    tmp = find(c1mask==i);
    [x, y, z] = ind2sub([L M N],tmp); % coordinates
    vol = length(tmp)*prod(prm.h); % volume in um^3
    diameter = 2*(vol*3/4*pi()).^(1/3); % radius in um
    intensity = mean(img.c1(tmp));
    
    dat.c1 = [dat.c1;double([mean(x),mean(y),mean(z),i,vol,diameter,intensity])];
end


dat.PVNvol = prod(prm.h)*sum(hypomask(:));

