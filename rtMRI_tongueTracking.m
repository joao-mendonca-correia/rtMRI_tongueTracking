dataPath = '..\DATA\DICOMS\';
outputPath = '..\DATA\OUTPUT\';

if ~exist(outputPath)
    mkdir(outputPath);
end

mydir = pwd;

d=dir([dataPath '\' '*.dcm']);
if ~length(d)
    cd(dataPath);
    system('rename *.ima *.dcm');
    cd(mydir);
    d=dir([dataPath '\' '*.dcm']);
end


%% The method is:
%1-dilate
%2-clear small clusters
%3-erode


%% options to adjust
firstFrame = 40;

BWfactor = 2;
SE_dilMask_strel = 5; %Radius used to dilate a disk shape, 5
clearSmallClusters = 50; %Threshold of pixels to clear small clusters , 50

conn = 4    ; %Number of connected pixels used to binarize image = 4 or 8., 4
erodeFilterKernel = [3 3]; %kernel used to filter the data for erosion, [3 3]

SE_dilCOI_strel = 1; %how much to dilate the final COI (contour of interest), 1
erodeNrPixelsTh = 20;
erodeConn = 26;

LipThickness = 7; %in pixels, 4
maxHorizontalDiscontinuity = 1;

frames = [firstFrame:length(d)];


%% Loop through frames
t=1;

for frame = frames
    
    templateFrame = dicomread([dataPath '\' d(frame).name]);
    I = imadjust(templateFrame,stretchlim(templateFrame));
    
    if (t==1)
        %define region of interest ROI
        fprintf('*** Use the mouse the demark the Region Of Interest (ROI) for tongue segmentation\nOnly min and max will be used to make a rectangular ROI\n');
        SelectedMargin = roipoly(I);
        sumRows = sum(SelectedMargin,2);
        sumCols = sum(SelectedMargin,1);
        
        UpperMarginOut = find(sumRows,1,'first');
        LowerMarginOut = find(sumRows,1,'last');
        
        LeftMarginOut = find(sumCols,1,'first');
        RightMarginOut = find(sumCols,1,'last');
        
        ROImask = ones(size(templateFrame));
        ROImask(LowerMarginOut:end,:) = 0;
        ROImask(1:UpperMarginOut,:) = 0;
        ROImask(:,RightMarginOut:end) = 0;
        ROImask(:,1:LeftMarginOut) = 0;
        
    end
    I(LowerMarginOut:end,:) = [];
    I(1:UpperMarginOut,:) = [];
    I(:,RightMarginOut:end) = [];
    I(:,1:LeftMarginOut) = [];
    
    %calculates histogram
    [counts, x] = imhist(I);
    %cumulative sum of the histogram
    cumcounts = cumsum(counts);
    %find threshold on histogram to binarize
    [minCounts, minCountsIdx] = min(abs(cumcounts-length(I(:))/BWfactor));
    %binarize image
    BW = im2bw(I, graythresh(I)*.5);
    %BW = im2bw(I,minCountsIdx/length(cumcounts));
    
    BW = bwareaopen(BW,clearSmallClusters,conn);
    
    %get the anterior contour of the lips
    clear antPoints lowAntPoint LLantPoints ULantPoints HorizontalDisp
    lowAntPoint = find(BW(end,:),1,'first');
    for ii=1:size(BW,1)
        antPoints(ii,:) = [ii find(BW(ii,:),1,'first')];
    end
    
    %   antPoints(antPoints(:,2)>lowAntPoint,:) = [];
    HorizontalDisp = diff(antPoints(:,2));
    
    LLantPoints = antPoints;
    ULantPoints = antPoints;
    
    LLantPoints(1:find(-1*HorizontalDisp>maxHorizontalDiscontinuity,1,'last'),:) = [];
    ULantPoints(find(HorizontalDisp>maxHorizontalDiscontinuity,1,'first'):end,:) = [];
    
    tmp1 = min(LLantPoints,[],1);
    tmp2 = max(ULantPoints,[],1);
    if ~isempty(tmp1) & ~isempty(tmp2)
        lipApperture(t) = tmp1(1)-tmp2(1);
    else
        lipApperture(t) = nan;
    end
    
    LLmask = zeros(size(I));
    ULmask = zeros(size(I));
    for ii=1:size(LLantPoints,1)
        LLmask(LLantPoints(ii,1),LLantPoints(ii,2):LLantPoints(ii,2)+LipThickness) = 1;
    end
    for ii=1:size(ULantPoints,1)
        ULmask(ULantPoints(ii,1),ULantPoints(ii,2):ULantPoints(ii,2)+LipThickness) = 1;
    end
    
    
    
    if(t==1)
        %define tongue mask
        fprintf('*** Use the mouse the demark the contour of the tongue\n');
        templateMask = roipoly(I);
        templateMask0 = templateMask;
        
        %define other structures that you don't ever want the tongue to invade
        fprintf('*** Use the mouse the demark other structures you don''t want the tongue to invade, like palate or velum\n');
        ExternalMask = roipoly(I);
        
    end
    
    
    
    
    %1-dilate
    SE_dilMask = strel('disk', SE_dilMask_strel);
    dilatedTemplateMask = imdilate(templateMask,SE_dilMask);
    maskedBW = BW;
    maskedBW(~dilatedTemplateMask) = 0;
    
    %2-clear small clusters
    maskedBW = bwareaopen(maskedBW,clearSmallClusters,conn);
    %spur the image twice
    maskedBW = bwmorph(maskedBW,'spur');
    maskedBW = bwmorph(maskedBW,'spur');
    
    %3-erode
    h = fspecial('gaussian', erodeFilterKernel);
    maskedBW = filter2(h, maskedBW);
    maskedBW(maskedBW<1) = 0;
    maskedBW = bwareaopen(maskedBW,erodeNrPixelsTh,conn);
    
    %remove overlap with Lips
    maskedBW(LLmask==1) = 0;
    maskedBW(ULmask==1) = 0;
    
    %remove external mask
    maskedBW(ExternalMask==1) = 0;
    
    
    %erodes and then dilates the BW image
    maskedBW = bwmorph(maskedBW,'open');
    
    
    CC = bwconncomp(maskedBW);
    COI = zeros(size(maskedBW));
    for CCIdx=1:length(CC.PixelIdxList)
        tmpCC = zeros(size(BW));
        tmpCC(CC.PixelIdxList{1,CCIdx}) = 1;
        tmpDiff(CCIdx) = sum(abs(templateMask0(:)-tmpCC(:)));
    end
    
    [minDiff minDiffIdx] = min(tmpDiff);
    tmpDiff = [];
    COI(CC.PixelIdxList{1,minDiffIdx}) = 1;
    
    %erodes and then dilates the BW image
    COI = bwmorph(COI,'open');
    
    
    templateMask = COI;
    templateMask = bwmorph(templateMask,'fill');
    
    
    %get contour of mask
    [~, threshold] = edge(templateMask, 'sobel');
    fudgeFactor = .3;
    EdgeDilatedCOI = edge(templateMask,'sobel', threshold * fudgeFactor);
    
    
    LLmask = bwmorph(LLmask,'open');
    ULmask = bwmorph(ULmask,'open');
    
    [L N] = bwlabel(templateMask,4);
    for ii=1:N
        nrPixels(ii) = sum(L(:)==ii);
    end
    [tmp biggestROI] = max(nrPixels);
    templateMask = L==biggestROI;
    
    
    %plot
    Output = imadjust(templateFrame,stretchlim(templateFrame));
    ROImaskIndices = find(ROImask);
    Output(ROImaskIndices(find(EdgeDilatedCOI))) = 2^16;
    subplot(1,2,1);
    imshow(Output);
    
    Output = imadjust(templateFrame,stretchlim(templateFrame));
    ROImaskIndices = find(ROImask);
    Output(ROImaskIndices(find(maskedBW))) = 2^16;
    
    Output(ROImaskIndices(find(LLmask))) = (2^16)/4;
    Output(ROImaskIndices(find(ULmask))) = (2^16)/2;
    subplot(1,2,2);
    imshow(Output);
    
    tmpMask = zeros(size(templateFrame));
    tmpMask(ROImaskIndices(find(maskedBW))) = 1;
    TongueMask(t,:,:) = tmpMask;
    LowerLipMask(t,:,:) = LLmask;
    UpperLipMask(t,:,:) = ULmask;
    
    Movie(t) = getframe(gcf);
    t=t+1;
end

lipApperture(lipApperture<0) = 0;
implay(Movie);
save([outputPath 'outputMasks.mat'], 'TongueMask', 'LowerLipMask', 'UpperLipMask', 'lipApperture');
save([outputPath 'outputMovie.mat'], 'Movie', '-v7.3');

figure('color', [1 1 1]); plot(smooth(lipApperture));
xlabel('Frames'); ylabel('pixels');
title('Time-course of lip apperture');

