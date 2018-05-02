% THIS FUNCTION NEEDS BLIST TO WORK!!!
%
% INPUT:
% img = the image to binarize.
%
% OUTPUT:
% Ims = The segmented image.
% 
function [Ims] = ADBLIST(img)
 
img = im2uint8(img);

[nrows, ncols] = size(img);
Ims = ones(nrows, ncols);

% use 20% of the total image width as window
w = ncols*0.02;

maxind = (ncols/w);
lvls = zeros(1, maxind);
rs = zeros(1, maxind);
r=1;
h=w/2;

% apply global BLIST to each window to determine local thresholds
for i=1:w:ncols-w
    window = img(:, i:i+w);
    [~, topt] = BLIST(window, 'pairs');
    lvl = topt/256;
    lvls(r) = lvl;
    rs(r) = r*w-h+1;
    r=r+1;
end

% call a 3rd level polinomial regression on the local thresholds
tam = [1:ncols]';
finalTs = poly_regression(rs', lvls', tam, 3);

for i=1:ncols
    
	if finalTs(i) >= 0 & finalTs(i) <= 1
		Ims(:, i) = im2bw(img(:,i), finalTs(i));
    end
    
    if finalTs(i) < 0
		Ims(:, i) = im2bw(img(:,i), 0);
    end
    
    if finalTs(i) > 1
		Ims(:, i) = im2bw(img(:,i), 1);
    end
end
