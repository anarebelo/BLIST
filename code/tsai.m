% INPUT:
% img = grey scale image
% 
% OUTPUT
% Ims = Segmented image
% topt = optimum threshold
%
function [Ims, topt] = TSAI(img)

I = img;
img = double(img);
% Calculate the histogram.
H = hist(img(:),0:255);

% Find Prob
Prob = zeros(1,256);
for t = 0:255
  Prob(t+1) = Aux(H,t)/Aux(H,255);
end

% Find x0
x2 = (m1(H)*m2(H)-Aux(H, 255)*m3(H)) / (Aux(H, 255)*m2(H)-m1(H)^2);
x1 = (m1(H)*m3(H)-m2(H)^2) / (Aux(H,255)*m2(H)-m1(H)^2);
x0 = .5 - (m1(H)/Aux(H,255)+x2/2) / sqrt(x2^2-4*x1);

% The threshold is chosen such that Prob is closest to x0
[~,ind] = min(abs(Prob-x0));
topt = ind-1;
lvl = topt/256;
Ims = im2bw(I, lvl);

%%% AUXILIAR FUNCTIONS
%
%
function x = Aux(H, j)

x = sum(H(1:j+1));
%
%
function x = m1(H)

ind = 0:255;
x = ind*H(1:256)';
%
%
function x = m2(H)

ind = 0:255;
x = ind.^2*H(1:256)';
%
%
function x = m3(H)

ind = 0:255;
x = ind.^3*H(1:256)';
