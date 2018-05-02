% INPUT:
% img = grey scale image
% window = window size
% k = standard deviation coeficient
% 
% OUTPUT
% Ims = Segmented image
%
function [Ims] = NIBLACK(img, window, k)


[nrows ncols] = size(img);
Ims = ones([nrows ncols]);


% determine halfwindow to be used in iterations
halfwindow = round((window-1)/2);

finalrow = nrows - halfwindow;
finalcol = ncols - halfwindow;

% run through limits of image, findind a threhsold for each window
for j=halfwindow+1:finalcol-1
    for i=halfwindow+1:finalrow-1
        curwindow = img(i-halfwindow:i+halfwindow, j-halfwindow:j+halfwindow);
        cur = double(curwindow);
        curwindow(:);
        soma = sum(curwindow(:));
        average = soma/(window*window);
        variance = var(cur(:));
        t = average + k * sqrt(variance);
        if img(i,j)<t
            Ims(i,j) = 0;
        end
    end
end
toc
