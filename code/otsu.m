% INPUT:
% img = grey scale image
% 
% OUTPUT
% Ims = Segmented image
% topt = optimum threshold
%
function [Ims, topt] = OTSU(img)

lvl = graythresh(img);

topt = lvl * 256;

Ims = im2bw(img, lvl);

close all