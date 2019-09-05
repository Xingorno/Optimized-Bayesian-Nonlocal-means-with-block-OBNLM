function [ distance ] = getPearsonDistance( originBlock, neighborBlock, gama)

% originBlock: 2D, n*n
% neighborBlock: 2D, n*n
% gama: double type

% gama = 0.5
% 
% originBlock = [1 2; 3 4];
% 
% neighborBlock = [5 6; 8 7];
epsilon = 10^(-13); % handle 0/0 case 
if nargin <=1
    
    disp('The input does not meet demands, please recheck it!!!')
    
end

if nargin ==2 
    
    gama = 0.5; 
    
end

    temp1 = (originBlock - neighborBlock).*(originBlock - neighborBlock);

    temp2 = neighborBlock.^(2*gama) + epsilon;

    distance = sum(sum(temp1./temp2));


end

