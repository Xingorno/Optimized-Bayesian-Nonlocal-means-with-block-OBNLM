function normalizedImg = ImgNormalize(originalImg, method)

% Input: 
%   originalImg: origin image [data type: double or integer]
%   method: if method == 1, histogram linear stretching; if method ==2, histogram nonlinear stertching

% Output: normalized image [data type: integer(uint8)]

% Author: Shuwei Xing
% Date: 2019-09-04

if nargin == 1
    method = 1; % histogram linear stretch
end

switch method
    
    case 1
        originalImg = double(originalImg);
        mappedMax = 255;
        mappedMin = 0;
        originalMax = max(max(originalImg));
        %originalMin = min(min(originalImg)) % reduce the disturbtion from noise
        originalMin = 0;
        deltaOrigin = originalMax - originalMin;
        scalar = (mappedMax - mappedMin)/deltaOrigin;
        normalizedImg = (originalImg - originalMin) * scalar + mappedMin;
        normalizedImg = uint8(normalizedImg);
        
    case 2
        disp('Undone!!')
end
end
