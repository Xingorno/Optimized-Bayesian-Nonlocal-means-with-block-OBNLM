clc;
clear all;

%% 
% Input Image Data
% path ='D:\Files\Projects\HISKY\SpeckleDenoising\MyContribution\HiskyData\';

% % path = 'D:\Files\Projects\HISKY\SpeckleDenoising\MyContribution\BayesianNLM\'
% % 
% % file = dir(fullfile(path,'*.png')); % (*.dat)
% % 
% % fileNames = {file.name}';
% % 
% % numFiles = size(fileNames,1);
% % 
% % for i = 1:1
% %     singleImgName = strcat(path, fileNames(i));
% % %     fid = fopen(singleImgName{1});
% % %     tline = fread(fid, [512, 512], 'int16');
% % %     R = tline';
% % %     R(R<0) = 0; 
% % %     mask = logical(R);
% % %     R = uint8(R);
% % %     f = @(X) imadjust(X,[],[],1);
% % %     img = roifilt2(R, mask, f);
% % %     figure
% % %     imshow(img,[0, 150])
% %     img = imread(singleImgName{1});
% %     figure
% %     imshow(img)
% % end

img = imread('noisyImage.png')
blockSize = 5; % size of the block
windowSize = 21; % size of the search window
gapBwnBlock = 2; % gap between the search block (in order to solve computational burden)
h = 5; % filtering parameter controlling the decay of the exponential function

img = ImgNormalize(img);
processedImg = BayesianNLM(img, blockSize, windowSize, gapBwnBlock, h)


figure
subplot 131
imshow(img)
title('Origin Image')
subplot 132
imshow(processedImg)
title('Despecked Image')
subplot 133
delta = ~logical(img - processedImg);
imshow(double(delta))
title('Subtraction Image')

imwrite(processedImg, 'despeckledImage.png')