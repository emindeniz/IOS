%% This program is used to display and analyze single IOS maps. Requires pre and post stimulus
% averages to be obtained first. These could be done using Intrinsic
% Signal Processor VSD1 macro for ImageJ. 
% http://neuroscience.ubc.ca/faculty/murphy_software.html


%Read the images
base = imread('cone_base.tif');%read prestimulus average
base = im2double(base);
resp = imread('cone_resp.tif');%read poststimulus average
resp = im2double(resp);
grayResponse=imsubtract(resp,base);
grayResponse=imdivide(grayResponse,base);
%Platform specific, The area represented by each pixel in square microns
areaconversionfactor=9.766; 
ROIDidNotExist=0;

%Display parameters scaling, drawing boundaries etc.
scalingLow=0.00040;
scalingHigh=0.00020;
displayRstart=100;
displayRend=840;
displayCstart=200;
displayCend=988;
drawBoundaries=1; %Make this 0 if you do not want boundaries to be drawn

%determine genotype of files in this folder
%Genotypes are represented as text files in the same folder
if exist('wt.txt')
    genotype='wt';
    genotypecode=0;
elseif exist('het.txt')
    genotype='het';
    genotypecode=1;
else
    genotype='NaN';
    genotypecode=2;
end

%determine whisker type of files in this folder
if size(ls('ctwo*'))
    whisker='ctwo';
    whiskercode=0;
elseif size(ls('beta*'))
    whisker='beta';
    whiskercode=1;
else
    whisker='NaN';
    whiskercode='NaN';
end

%determine the animal id files in this folder
%Animal ID is a simple number in a txt file
if size(ls('animal.txt*'))
    animalid=csvread('animal.txt');
else
    animalid='NaN';
end

%ROI mat file represents all the parameters saved from a previous analysis
%If you want to start over just delete ROIs.mat file and restart
%You will be asked to crop the images in order to determine baseline
%and response regions
if exist('ROIs.mat','file')
    load('ROIs.mat','grayResponseBaseline');
    load('ROIs.mat','grayResponseCropped');
    load('ROIs.mat','rect');
    load('ROIs.mat','rectBaseline');
else
    ROIDidNotExist=1;
    %CROP TO DETERMINE POST STIMULUS BASELINE
    figure, imshow(grayResponse, []);
    [grayResponseBaseline,rectBaseline]=imcrop;

    %CROP TO DETERMINE THE AREA OF ACTIVATION(to get rid of artifact)
    [grayResponseCropped,rect]=imcrop;
    save('ROIs.mat','grayResponseBaseline','grayResponseCropped','rect','rectBaseline');
    close all;
end

%FILTER THE IMAGE
h=fspecial('average',[30 30]); %Images need to be filtered before analysis
grayResponse = imfilter(grayResponse,h);
grayResponseBaseline = imfilter(grayResponseBaseline,h);
grayResponseCropped = imfilter(grayResponseCropped,h);

%determine the response size
medianratio = median(grayResponseBaseline(:));
minratio = min(grayResponseCropped(:));
responsesize = medianratio - minratio ;
fprintf('Response Amplitude is = %e' , ...
		responsesize);

%the method of absolute thresholding
figure, imshow(grayResponse(100:840,200:988), [medianratio-scalingLow medianratio+scalingHigh]);
thresholdinglevels = [0.00020 0.00030 0.00040];
IOSareasAbsolute=zeros(size(thresholdinglevels));
IOScentroidsAbsolute=zeros(length(thresholdinglevels),2);

for l = 1:length(thresholdinglevels) %threshold seperately for each threshold level
binaryImage = (grayResponseCropped > -1000) & (grayResponseCropped < (medianratio-thresholdinglevels(l)));
if bwarea(binaryImage)~=0
    binaryImage = bwareafilt(binaryImage,1);
    %FIND BOUNDARIES AND DRAW AREAS
    [B,L,N,A] = bwboundaries(binaryImage,4,'noholes');

           for k = 1:length(B)
            boundary = B{k};
            BW2 = bwselect(binaryImage,boundary(:,2), boundary(:,1));
            disp(bwarea(BW2)*areaconversionfactor);%This conversion factor should have been explained!!!
            IOSareasAbsolute(l)=bwarea(BW2)*areaconversionfactor;
            s = regionprops(BW2, 'centroid');%find the center of the area
            centroid = cat(1, s.Centroid);
            if drawBoundaries
                axis off;
                hold on;
                plot(boundary(:,2)+rect(1)-200, boundary(:,1)+rect(2)-100, 'w', 'LineWidth', 3);
            else
            end
            IOScentroidsAbsolute(l,1)=centroid(:,1)+rect(1);
            IOScentroidsAbsolute(l,2)=centroid(:,2)+rect(2);
        end
end
end

if drawBoundaries
plot(centroid(:,1)+rect(1)-200, centroid(:,2)+rect(2)-100, 'w*');
else
end
set(gcf,'InvertHardCopy','off');
print -dbmp absolutethresholdfig;

%the method of relative thresholding
figure, imshow(grayResponse(100:840,200:988), [medianratio-scalingLow medianratio+scalingHigh]);
thresholdinglevels = [responsesize*0.5 responsesize*0.6 responsesize*0.7 responsesize*0.8 ];
IOSareasRelative=zeros(size(thresholdinglevels));
IOScentroidsRelative=zeros(length(thresholdinglevels),2);


for m = 1:length(thresholdinglevels) %threshold seperately for each threshold level
binaryImage = (grayResponseCropped > -1000) & (grayResponseCropped < (medianratio-thresholdinglevels(m)));
if bwarea(binaryImage)~=0
    binaryImage = bwareafilt(binaryImage,1);
    %FIND BOUNDARIES AND DRAW AREAS
    [B,L,N,A] = bwboundaries(binaryImage,4,'noholes');

          for n = 1:length(B)
            boundary = B{n};
            BW2 = bwselect(binaryImage,boundary(:,2), boundary(:,1));
            disp(bwarea(BW2)*areaconversionfactor);
            IOSareasRelative(m)=bwarea(BW2)*areaconversionfactor;
            s = regionprops(BW2, 'centroid');%find the center of the area
            centroid = cat(1, s.Centroid);
            axis off;
            hold on;
            plot(boundary(:,2)+rect(1)-200, boundary(:,1)+rect(2)-100, 'w', 'LineWidth', 3);
            IOScentroidsRelative(m,1)=centroid(:,1)+rect(1);
            IOScentroidsRelative(m,2)=centroid(:,2)+rect(2);
        end
end
end
plot(centroid(:,1)+rect(1)-200, centroid(:,2)-100+rect(2), 'w*');
set(gcf,'InvertHardCopy','off');
print -dbmp relativethresholdfig;

%{
Save the analysis results in a txt file. Results are:
Response size: Instrinsic imaging response size
IOSareasAbsolute: Area of the respnse as determined by absolute
thresholding methods
IOSareasAbsolute: Area of the respnse as determined by relative
thresholding methods
%}
fileID = fopen('IOSresults.txt','w');
fprintf(fileID,'%0.8f\t',responsesize, IOSareasAbsolute,IOSareasRelative,IOScentroidsAbsolute,IOScentroidsRelative);
fclose(fileID);
 


