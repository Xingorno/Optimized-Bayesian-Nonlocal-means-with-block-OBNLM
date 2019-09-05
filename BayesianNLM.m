
function processedImg = BayesianNLM(img, blockSize, windowSize, gapBwnBlock, h)

% Function: Bayesian-based nonlocal means algorithm
% Input: 
%   img: the origin image [data type: double or integer]
%   blockSize: size of the block
%   windowSize: size of the search window
%   gapBwnBlock: gap between the search blocks (in order to solve computational burden)
%   h: filtering parameter controlling the decay of the exponential function
%   
%
%   _ _ _ _ _ _ _ imgSize_ _ _ _ _ _ _ __
%   |* * * * * * * * * *                 |
%   |*                 *                 |
%   |*    + + +        *                 |
%   |*    +   +        *                 |
%   |*    + + +        *                 |
%   |*  blockSize      *                 |  
%   |*                 *                 |
%   |*                 *                 |
%   |* * windowSize  * *                 |
%   |                                    | 
%   |                                    |
%   |_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ |
%
%
% Output: the despeckled image

% Author: Shuwei Xing
% Date: 2019-09-03
% Reference: Coup¨¦, Pierrick, et al. "Nonlocal means-based speckle filtering for ultrasound images." IEEE transactions on image processing 18.10 (2009): 2221-2229.

%%
% Define variables for algorithm

originImg = double(img); % the raw image with noises
imgSize = size(originImg,1); % Assumption: the shape of image is suqare
processingImg = zeros(size(originImg)); % the value of processing
indexImg = zeros(size(originImg)); % show the number of each pixel
denoisedImg = zeros(size(originImg)); % the resotred(denoised) image
epsilon = 10^(-13); % handle 0/0 case 
sideBlock = fix(blockSize/2);

a = ceil(blockSize/2);
b = imgSize - a + 1;

A = ceil(windowSize/2);
B = imgSize - A + 1;

numCompletedPixel = 0;

%%

% Optimizaing the search space

rowTotal = ceil(a / gapBwnBlock)* gapBwnBlock : gapBwnBlock : b ;
colTotal = ceil(a / gapBwnBlock) * gapBwnBlock : gapBwnBlock : b;

indexImgTotal = zeros(imgSize);
indexImgTotal(rowTotal, colTotal) = 1;

mask = logical(img);
searchIndexMatrix =double(mask).* indexImgTotal;

[rows, cols] = find(searchIndexMatrix);

% Bayesian-based NonLocal Means Algorithm

%
tic
% for row = ceil(a / gapBwnBlock)* gapBwnBlock : gapBwnBlock : b 
%     for col = ceil(a / gapBwnBlock) * gapBwnBlock : gapBwnBlock : b

for i = 1:size(cols, 1)
    col = cols(i);
    row = rows(i);
%        row = 12;
%        col = 200;
        
       % Get the block matrix of origin image
       % To do list: if the alogrithm can run fastly, write the specific
       % function to get the block.       
        originBlock = originImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]);
        
        %
        % Discuss this solution within different conditions
        %
        
        % CONDITION 1        
        if row >= A && row <= B            
            if col >= A && col <= B                
                weightList1 = [];                
                NLBlockWithoutNormalization1 = zeros(blockSize, blockSize);                
                for row_NB = ceil((row - A + a )/ gapBwnBlock)* gapBwnBlock : gapBwnBlock : (row + A -a)                    
                    for col_NB = ceil((col - A + a)/ gapBwnBlock) * gapBwnBlock : gapBwnBlock : (col + A - a)                        
                        neighborBlock = originImg([(row_NB - sideBlock) : (row_NB + sideBlock) ], [(col_NB - sideBlock) : (col_NB + sideBlock)]);
                        
                        % Compute Pearson Distance                        
                        distance = getPearsonDistance(originBlock, neighborBlock);              
                        
                        % Compute the weight without normalization , just for single
                        % neighborBlock
                        weightWithoutNormalization = exp(-distance/(h^2));                        
                        weightList1 = [weightList1, weightWithoutNormalization];
                        
                        % Compute the restored neighbor block NL_u_B without normalization
                        subNLBlock1 = neighborBlock*weightWithoutNormalization;                        
                        NLBlockWithoutNormalization1 = NLBlockWithoutNormalization1 + subNLBlock1;             
                    end                    
                end
                
                % Sum these weightes and get the final restored neighbor
                % block NL_u_B                
                NLBlock = NLBlockWithoutNormalization1/sum(weightList1);
                
                % pixels'value of the NL_u_B is put in the processingImg                
                processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + NLBlock;             
                
                % pixels' logic value of the NL_u_B is put in the indexImg                
                indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + ones(blockSize);                
%                 disp('run to this position')
            end
        
            
            % CONDITION 2            
            if col < A               
                weightList2 = [];                
                NLBlockWithoutNormalization2 = zeros(blockSize, blockSize);                
                for row_NB = ceil((row - A + a )/ gapBwnBlock)* gapBwnBlock : gapBwnBlock : (row + A -a)                    
                    for col_NB = ceil(a / gapBwnBlock)* gapBwnBlock : gapBwnBlock : (windowSize + 1 -a)                        
                        neighborBlock = originImg([(row_NB - sideBlock) : (row_NB + sideBlock) ], [(col_NB - sideBlock) : (col_NB + sideBlock)]);
                        
                        % Compute Pearson Distance                        
                        distance = getPearsonDistance(originBlock, neighborBlock);              
                        
                        % Compute the weight without normalization , just for single
                        % neighborBlock
                        weightWithoutNormalization = exp(-distance/(h^2));                        
                        weightList2 = [weightList2, weightWithoutNormalization];
                        
                        % Compute the restored neighbor block NL_u_B without normalization
                        subNLBlock2 = neighborBlock*weightWithoutNormalization;                        
                        NLBlockWithoutNormalization2 = NLBlockWithoutNormalization2 + subNLBlock2;                        
                    end                    
                end
                
                % Sum these weights and get the final restored neighbor
                % block NL_u_B                
                NLBlock = NLBlockWithoutNormalization2/sum(weightList2);
                
                % pixels'value of the NL_u_B is put in the processingImg                
                processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + NLBlock;             
                
                % pixels' logic value of the NL_u_B is put in the indexImg                
                indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + ones(blockSize);            
            end
            
            
            % CONDITION 3            
            if col > B                
                weightList3 = [];                
                NLBlockWithoutNormalization3 = zeros(blockSize, blockSize);                
                for row_NB = ceil((row - A + a )/ gapBwnBlock)* gapBwnBlock : gapBwnBlock : (row + A -a)                    
                    for col_NB = ceil((imgSize - windowSize + a) / gapBwnBlock)* gapBwnBlock : gapBwnBlock : b                         
                        neighborBlock = originImg([(row_NB - sideBlock) : (row_NB + sideBlock) ], [(col_NB - sideBlock) : (col_NB + sideBlock)]);
                        
                        % Compute Pearson Distance                       
                        distance = getPearsonDistance(originBlock, neighborBlock);              
                        
                        % Compute the weight without normalization , just for single
                        % neighborBlock
                        weightWithoutNormalization = exp(-distance/(h^2));                        
                        weightList3 = [weightList3, weightWithoutNormalization];
                        
                        % Compute the restored neighbor block NL_u_B without normalization
                        subNLBlock = neighborBlock*weightWithoutNormalization;                        
                        NLBlockWithoutNormalization3 = NLBlockWithoutNormalization3 + subNLBlock;                                                
                    end
                    
                end
                
                % Sum these weights and get the final restored neighbor
                % block NL_u_B                
                NLBlock = NLBlockWithoutNormalization3/sum(weightList3);
                
                % pixels'value of the NL_u_B is put in the processingImg                
                processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + NLBlock;             
                
                % pixels' logic value of the NL_u_B is put in the indexImg                
                indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + ones(blockSize);                
            end                        
        end
                              
        if row < A            
            % CONDITION 4            
            if col >= A && col <= B                                
                weightList4 = [];                
                NLBlockWithoutNormalization4 = zeros(blockSize, blockSize);                
                 for row_NB = ceil(a / gapBwnBlock)* gapBwnBlock : gapBwnBlock : (windowSize + 1 -a)                    
                    for col_NB = ceil((col - A + a) / gapBwnBlock)* gapBwnBlock : gapBwnBlock : (col + A - a)                         
                        neighborBlock = originImg([(row_NB - sideBlock) : (row_NB + sideBlock) ], [(col_NB - sideBlock) : (col_NB + sideBlock)]);
                        
                        % Compute Pearson Distance                        
                        distance = getPearsonDistance(originBlock, neighborBlock);              
                        
                        % Compute the weight without normalization , just for single
                        % neighborBlock
                        weightWithoutNormalization = exp(-distance/(h^2));                        
                        weightList4 = [weightList4, weightWithoutNormalization];
                        
                        % Compute the restored neighbor block NL_u_B without normalization
                        subNLBlock = neighborBlock*weightWithoutNormalization;                        
                        NLBlockWithoutNormalization4 = NLBlockWithoutNormalization4 + subNLBlock;                                                
                    end                    
                 end
                 
                % Sum these weights and get the final restored neighbor
                % block NL_u_B                
                NLBlock = NLBlockWithoutNormalization4/sum(weightList4);
                
                % pixels'value of the NL_u_B is put in the processingImg                
                processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + NLBlock;             
                
                % pixels' logic value of the NL_u_B is put in the indexImg                
                indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + ones(blockSize);     
                 
            end
            
            
            % CONDITION 5
            if col < A                
                weightList5 = [];                
                NLBlockWithoutNormalization5 = zeros(blockSize, blockSize);                
                for row_NB = ceil(a / gapBwnBlock)* gapBwnBlock : gapBwnBlock : (windowSize + 1 - a)                    
                    for col_NB = ceil(a / gapBwnBlock)* gapBwnBlock : gapBwnBlock : (windowSize + 1 - a)                         
                        neighborBlock = originImg([(row_NB - sideBlock) : (row_NB + sideBlock) ], [(col_NB - sideBlock) : (col_NB + sideBlock)]);
                        
                        % Compute Pearson Distance                        
                        distance = getPearsonDistance(originBlock, neighborBlock);              
                        
                        % Compute the weight without normalization , just for single
                        % neighborBlock
                        weightWithoutNormalization = exp(-distance/(h^2));                        
                        weightList5 = [weightList5, weightWithoutNormalization];
                        
                        % Compute the restored neighbor block NL_u_B without normalization
                        subNLBlock = neighborBlock*weightWithoutNormalization;                        
                        NLBlockWithoutNormalization5 = NLBlockWithoutNormalization5 + subNLBlock;                                                
                    end                    
                end
                
                % Sum these weights and get the final restored neighbor
                % block NL_u_B                
                NLBlock = NLBlockWithoutNormalization5/sum(weightList5);
                
                % pixels'value of the NL_u_B is put in the processingImg                
                processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + NLBlock;             
                
                % pixels' logic value of the NL_u_B is put in the indexImg                
                indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + ones(blockSize);     
                 
            end
            
            
            
            % CONDITION 6            
            if col > B                
                weightList6 = [];                
                NLBlockWithoutNormalization6 = zeros(blockSize, blockSize);                
                for row_NB = ceil(a / gapBwnBlock)* gapBwnBlock : gapBwnBlock : (windowSize + 1 - a)                    
                    for col_NB = ceil((imgSize - windowSize + a) / gapBwnBlock)* gapBwnBlock : gapBwnBlock : b                         
                        neighborBlock = originImg([(row_NB - sideBlock) : (row_NB + sideBlock) ], [(col_NB - sideBlock) : (col_NB + sideBlock)]);
                        
                        % Compute Pearson Distance                        
                        distance = getPearsonDistance(originBlock, neighborBlock);              
                        
                        % Compute the weight without normalization , just for single
                        % neighborBlock
                        weightWithoutNormalization = exp(-distance/(h^2));                        
                        weightList6 = [weightList6, weightWithoutNormalization];
                        
                        % Compute the restored neighbor block NL_u_B without normalization
                        subNLBlock = neighborBlock*weightWithoutNormalization;                        
                        NLBlockWithoutNormalization6 = NLBlockWithoutNormalization6 + subNLBlock;                                                
                    end                    
                end
                
                % Sum these weights and get the final restored neighbor
                % block NL_u_B                
                NLBlock = NLBlockWithoutNormalization6/sum(weightList6);
                
                % pixels'value of the NL_u_B is put in the processingImg                
                processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + NLBlock;             
                
                % pixels' logic value of the NL_u_B is put in the indexImg                
                indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + ones(blockSize);                                                    
            end
            
        end
        
                
        if row > B
            
            % CONDITION 7            
            if col >= A && col <= B                
                weightList7 = [];                
                NLBlockWithoutNormalization7 = zeros(blockSize, blockSize);                
                for row_NB = ceil((imgSize - windowSize + a) / gapBwnBlock)* gapBwnBlock : gapBwnBlock : b                    
                    for col_NB = ceil((col - A + a) / gapBwnBlock)* gapBwnBlock : gapBwnBlock : (col + A - a)                         
                        neighborBlock = originImg([(row_NB - sideBlock) : (row_NB + sideBlock) ], [(col_NB - sideBlock) : (col_NB + sideBlock)]);
                        
                        % Compute Pearson Distance                       
                        distance = getPearsonDistance(originBlock, neighborBlock);              
                        
                        % Compute the weight without normalization , just for single
                        % neighborBlock
                        weightWithoutNormalization = exp(-distance/(h^2));                        
                        weightList7 = [weightList7, weightWithoutNormalization];
                        
                        % Compute the restored neighbor block NL_u_B without normalization
                        subNLBlock = neighborBlock*weightWithoutNormalization;                        
                        NLBlockWithoutNormalization7 = NLBlockWithoutNormalization7 + subNLBlock;                                                
                    end                    
                end
                
                % Sum these weights and get the final restored neighbor
                % block NL_u_B                
                NLBlock = NLBlockWithoutNormalization7/sum(weightList7);
                
                % pixels'value of the NL_u_B is put in the processingImg                
                processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + NLBlock;             
                
                % pixels' logic value of the NL_u_B is put in the indexImg                
                indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + ones(blockSize);     
                                 
            end
            
            
            % CONDITION 8            
            if col < A                
                weightList8 = [];                
                NLBlockWithoutNormalization8 = zeros(blockSize, blockSize);                
                for row_NB = ceil((imgSize - windowSize + a) / gapBwnBlock)* gapBwnBlock : gapBwnBlock : b                    
                    for col_NB = ceil(a / gapBwnBlock)* gapBwnBlock : gapBwnBlock : (windowSize - a + 1)
                         
                        neighborBlock = originImg([(row_NB - sideBlock) : (row_NB + sideBlock) ], [(col_NB - sideBlock) : (col_NB + sideBlock)]);
                       
                        % Compute Pearson Distance                        
                        distance = getPearsonDistance(originBlock, neighborBlock);              
                        
                        % Compute the weight without normalization , just for single
                        % neighborBlock
                        weightWithoutNormalization = exp(-distance/(h^2));                        
                        weightList8 = [weightList8, weightWithoutNormalization];
                        
                        % Compute the restored neighbor block NL_u_B without normalization
                        subNLBlock = neighborBlock*weightWithoutNormalization;                        
                        NLBlockWithoutNormalization8 = NLBlockWithoutNormalization8 + subNLBlock;                                                
                    end                    
                end
                
                % Sum these weights and get the final restored neighbor
                % block NL_u_B                
                NLBlock = NLBlockWithoutNormalization8/sum(weightList8);
                
                % pixels'value of the NL_u_B is put in the processingImg                
                processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + NLBlock;             
                
                % pixels' logic value of the NL_u_B is put in the indexImg                
                indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + ones(blockSize);     
                           
            end
            
            
            % CONDITION 9            
            if col > B                
                weightList9 = [];                
                NLBlockWithoutNormalization9 = zeros(blockSize, blockSize);                
                for row_NB = ceil((imgSize - windowSize + a) / gapBwnBlock)* gapBwnBlock : gapBwnBlock : b                    
                    for col_NB = ceil((imgSize - windowSize + a) / gapBwnBlock)* gapBwnBlock : gapBwnBlock : b                        
                        neighborBlock = originImg([(row_NB - sideBlock) : (row_NB + sideBlock) ], [(col_NB - sideBlock) : (col_NB + sideBlock)]);
                        
                        % Compute Pearson Distance                        
                        distance = getPearsonDistance(originBlock, neighborBlock);              
                        
                        % Compute the weight without normalization , just for single
                        % neighborBlock
                        weightWithoutNormalization = exp(-distance/(h^2));                        
                        weightList9 = [weightList9, weightWithoutNormalization];
                        
                        % Compute the restored neighbor block NL_u_B without normalization
                        subNLBlock = neighborBlock*weightWithoutNormalization;                        
                        NLBlockWithoutNormalization9 = NLBlockWithoutNormalization9 + subNLBlock;                                                
                    end                    
                end
                
                % Sum these weights and get the final restored neighbor
                % block NL_u_B                
                NLBlock = NLBlockWithoutNormalization9/sum(weightList9);
                
                % pixels'value of the NL_u_B is put in the processingImg                
                processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = processingImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + NLBlock;             
                
                % pixels' logic value of the NL_u_B is put in the indexImg                
                indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) = indexImg([(row - sideBlock) : (row + sideBlock) ], [(col - sideBlock) : (col + sideBlock)]) + ones(blockSize);                     
            end                
        end
        
        numCompletedPixel = numCompletedPixel + 1;        
        disp(strcat('The completed number of pixel: ', num2str(numCompletedPixel)))
    end    
% end

%%

% Restore the pixel value of image

maskImg = logical(indexImg);
NLMUndonePixel = double(~maskImg).*originImg;
fullIndexImg = double(~maskImg) + indexImg;
processedImg = (processingImg + NLMUndonePixel)./fullIndexImg;
processedImg = uint8(processedImg);
toc

end