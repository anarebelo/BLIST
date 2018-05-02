% INPUT:
% img = grey scale image
% w = window size
% c = minimum contrast between pixels in a window to apply threshold
% 
% OUTPUT
% Ims = Segmented image
%
function [Ims] = BERNSEN(img, w, c)

img = im2uint8(img);
[nrows ncols] = size(img);
Ims=ones(nrows, ncols);

% determine half window h to be used in iterations
h = w/2;
h = int16(h);

for i=h:ncols-h
    for j=h:nrows-h
        window = img(j-h+1:j+h, i-h+1:i+h);
        iwindow = double(window);
        maxG = max(iwindow(:));
        minG = min(iwindow(:));
        medG = (maxG+minG)/2;
        if maxG-minG > c
            if img(j,i) < medG
                Ims(j,i) = 0;
            end
         else
             if medG < 128
                 Ims(j,i) = 0;
             end
        end
    end
end
