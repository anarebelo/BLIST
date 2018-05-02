%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
% img = grey scale image
%
% OUTPUT:
% Ims = binarization
% topt = optimum threshold
%
function [Ims, topt] = KAPUR(img)

% Use greyHist to determine the grey level histogram of the image
H = hist(img(:),0:255);

% Uses kapurThresh to determine optimum threshold topt
topt = kapurThresh(H);

% Determine lvl to be used in im2bw
lvl = topt/256;
Ims = im2bw(img, lvl);

close all

%
% INPUT:
% H = grey level histogram
%
% OUTPUT:
% topt = optimum threshold
%
function [topt] = kapurThresh(H)

% Call entropy function
v = zeros(1,256);
for t=1:256
   v(t) = Ent(H,t);
end 

% Calculate minimum entropy
[~,topt] = min(v);


% Entropy function
%
function v = Ent(H,t)

n = 256;

v = E(H,t)/Aux(H,t) - log10(Aux(H,t)) + (E(H,n)-E(H,t))/(Aux(H,n)-Aux(H,t)) - log10(Aux(H,n)-Aux(H,t));
%
%
%  
% Other functions
function x = E(H,j)

H = H(1:j);
H = H(H~=0);
x = sum(H.*log10(H));
%
function x = Aux(H,j)

x = sum(H(1:j));