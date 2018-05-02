%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% INPUT:
% img = grey scale image
% a1, a3 -> real values such as: 0<a1<1 e a3>1
%
% OUTPUT:
% Ims = binarized image
% topt = optimum threshold
%
function [Ims, topt] = SAHOO(img,a1,a3)

% Use greyHist to determine the grey level histogram of the image
H = hist(img(:),0:255);


% Calculate optimum threhsold topt
topt = sahooThresh(H,a1,a3);

% Determine lvl to be used in im2bw

lvl = topt/256;
Ims = im2bw(img, lvl);

close all
%
% This function calculates Renyi's entropy for 3 diferent alpha values (a1,
% a2, a3) and then generates their weighted average as topt
%
% INPUT:
% H = grey level histogram
% a1, a3 -> real values such as: 0<a1<1 e a3>1
%
% OUTPUT:
% topt = optimum threshold
%
function [topt] = sahooThresh(H,a1,a3)

% Calculate Renyi's entropy for 0<a1<1
r1 = 1:256;
for t1=0:255
   [sA,sB] = Renyi(H,a1,t1);
   r1(t1+1) = sA + sB;
end 
% calculate maximum entropy
[~,t1opt] = max(r1);

% Calculate Renyi's entropy for a=1
r2 = 1:256;
for t2=0:255
   [sA,sB] = Renyi(H,1,t2);
   r2(t2+1) = sA + sB;
end 
% calculate maximum entropy
[~,t2opt] = max(r2);

% Calculate Renyi's entropy for a3>1
r3 = 1:256;
for t3=0:255
   [sA,sB] = Renyi(H,a3,t3);
   r3(t3+1) = sA + sB;
end 
% calculate maximum entropy
[~,t3opt] = max(r3);


% Order topt by lowest to highest
if t1opt < t2opt
    if  t1opt < t3opt
        t1 = t1opt;
        if t2opt < t3opt
            t2 = t2opt;
            t3 = t3opt;
        else
            t2 = t3opt;
            t3 = t2opt;
        end
    else
        t1 = t3opt;
        t2 = t1opt;
        t3 = t2opt;
    end
else
    if t2opt < t3opt
        t1 = t2opt;
        if t1opt < t3opt
            t2 = t1opt;
            t3 = t3opt;
        else
            t2 = t3opt;
            t3 = t1opt;
        end
    else
        t1 = t3opt;
        t2 = t2opt;
        t3 = t1opt;
    end
end

% Determine b1, b2, b3 based on t1, t2, t3
if abs(t1 - t2) <= 5 && abs(t2 - t3) <= 5
    b1 = 1;
    b2 = 2;
    b3 = 1;
end
if abs(t1 - t2) > 5 && abs(t2 - t3) > 5
    b1 = 1;
    b2 = 2;
    b3 = 1;
end
if abs(t1 - t2) <= 5 && abs(t2 - t3) > 5
    b1 = 0;
    b2 = 1;
    b3 = 3;
end
if abs(t1 - t2) > 5 && abs(t2 - t3) <= 5
    b1 = 3;
    b2 = 1;
    b3 = 0;
end

tot1 = 0;
for i=0:t1-1
   tot1 = tot1 + H(1,i+1);
end
tot3 = 0;
for i=0:t3-1
   tot3 = tot3 + H(1,i+1);
end

pt1 = H(1, t1)/tot1;
pt3 = H(1, t3)/tot3;

w = pt3 - pt1;

% calculate topt as a weighted average of t1opt, t2opt and t3opt
topt = t1 * (pt1 + (w * b1)/4) + (t2 * w * b2)/4 + t3 * (1 - pt3 + (w * b3)/4);


% Renyi's entropy function
%
function [eA,eB] = Renyi(H,a,t)

eA = 0;
eB = 0;

% first class (before t)
tot1 = 0;
for i=0:t
   tot1 = tot1 + H(1,i+1);
end

sA = 0;
if (tot1 ~= 0)
 for i=0:t
    sA = sA + (H(1,i+1)/tot1)^a;   
 end
 eA = (log(sA))/(1-a);
end

% second class (after t)
tot2 = 0;
for i=t+1:255
   tot2 = tot2 + H(1,i);
end

sB = 0;
if (tot2 ~= 0)
 for i=t+1:255
    sB = sB + (H(1,i)/tot2)^a;   
 end
 eB = (log(sB))/(1-a);
end

if ((sA == 0) || (sB == 0))
   eA = 0;
   eB = 0;
end