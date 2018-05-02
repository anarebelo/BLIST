
%
%
% INPUT:
% img = grey scale image
% q = real value such as q > 0; if q = 1, the algorithm uses Shannon's
% Entropy
% 
% OUTPUT
% Ims = Segmented image
% topt = optimum threshold
%
function [Ims, topt] = TSALLIS(img,q)

% Use greyHist to determine the grey level histogram of the image
H = hist(img(:),0:255);

% Calculate optimum threhsold topt
topt = tsallisThresh(H,q);

% Determine lvl to be used in im2bw
lvl = topt/256;
Ims = im2bw(img, lvl);


close all

%
%
% INPUT:
% H = grey level histogram
% q = Tsallis entropic parameter
%
% OUTPUT:
% topt = optimum threshold
%
function topt = tsallisThresh(H,q)


% calculate tsallis for each class
r = 1:256;
for t=0:255
   [sA,sB] = Tsallis(H,q,t);
   r(t+1) = sA + sB + ((1-q) * sA * sB);
end 


% calculate maximum entropy
[~,topt] = max(r);




%
%
% Function to calculate Tsallis entropy
%
% INPUT: 
% H = grey scale image
% q = entropic parameter
% t = threshold that separates classes
%
% INPUT:
% eA -> left class entropy
% eB -> right class entropy
%
%
function [eA,eB] = Tsallis(H,q,t)

eA = 0;
eB = 0;

% calculate left class values
tot1 = 0;
for i=0:t
   tot1 = tot1 + H(1,i+1);
end

sA = 0;
if (tot1 ~= 0)
 for i=0:t
    sA = sA + (H(1,i+1)/tot1)^q;   
 end
 eA = (1-sA)/(q-1);
end

% calculate right class values
tot2 = 0;
for i=t+1:255
   tot2 = tot2 + H(1,i);
end

sB = 0;
if (tot2 ~= 0)
 for i=t+1:255
    sB = sB + (H(1,i)/tot2)^q;   
 end
 eB = (1-sB)/(q-1);
end

if ((sA == 0) | (sB == 0))
   eA = 0;
   eB = 0;
end


