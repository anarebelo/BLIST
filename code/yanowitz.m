% INPUT:
% img = grey scale image
% 
% OUTPUT
% Ims = Segmented image
%
function [Ims] = YANOWITZ(img)

[nrows ncols] = size(img);

% edges using sobel
e_img = edge(img, 'sobel');
ind=find(e_img == 1);

% smoothing with a 3x3 average filter
h=fspecial('average',3);

p = im2double(imfilter(img,h));

% iter will be only 20% of N 
iter = round(ncols*0.05);

% iterative method to determine threshold surface p
for c=1:iter

% calculate r
    r = del2(p);

% calculate next p
    p1 = p + r;
    p1(ind) = img(ind);
    p1([1 end],:) = p1([2 end-1],:);
    p1(:,[1 end]) = p1(:,[2 end-1]);
    p = p1;
end


% use threshold surface p to threshold image img into Ims
Ims = zeros(nrows, ncols);
for i=1:nrows
    for j=1:ncols
        if img(i, j) > p(i, j)
            Ims(i, j) = 1;
        end
    end
end