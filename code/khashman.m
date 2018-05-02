% INPUT:
% img = grey scale image
% 
% OUTPUT
% Ims = Segmented image
% topt = optimum threshold
%
function [Ims, topt] = KHASHMAN(img)



Gmax = max(img(:));


% INCORRECT VERSION used in the project :(
%
M = median(img(:));
D = Gmax - M;
topt = M - D;


% CORRECT VERSION
%
% [y x] = size(img);
% M = sum(img(:))/(x*y);
% D = Gmax - M;
% topt = M - D;

topt = cast(topt, 'double');
lvl = topt/256;
Ims = im2bw(img, lvl);

close all