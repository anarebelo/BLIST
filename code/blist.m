% INPUT:
% img = the image to binarize.
% seq = sequence to be used ('pairs' or 'triplets')
%
% OUTPUT:
% Ims = The segmented image.
% topt = optimum threshold
%
function [Ims, topt] = BLIST(img, seq)

    img = im2uint8(img);

    topt = globalMusicThreshold(img, seq);
    
    lvl = topt/256;
    Ims = im2bw(img, lvl);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function threshold = globalMusicThreshold(img, method)
%method: 'pairs' or 'triplets'
    switch lower(method)
       case 'pairs'
          methodID=1;
       case 'triplets'
          methodID=2;
       otherwise
          error('Unknown method.');
    end
    img = im2uint8(img);

    [nrows ncols] = size(img);
    
    medvalue = median(img(:));
    minValue = min(img(:));
    
    acumHist = zeros(nrows,1);      %the maximum value of the sum of two or three consecutive runs is nrows
    indHists = zeros(nrows,256);    %histogram for each threshold
    
    for th=minValue+1:medvalue

        bw = (img>th);
        for col=1:ncols
            data = bw(:,col); 
            data = rle(data);
            values1=data{2};            
            values2=values1(2:end);
            values3=values1(3:end);
            if methodID==1  % if using pairs
                sumConsecutiveRuns=values1(1:end-1)+values2;
            else            % if using triplets
                sumConsecutiveRuns=values1(1:end-2)+values2(1:end-1)+values3;
                % keep only triplets (black,white,black)
                sumConsecutiveRuns=sumConsecutiveRuns((bw(1,col)+1):2:end);
            end         
            for i=1:length(sumConsecutiveRuns)
                indHists(sumConsecutiveRuns(i),th)=indHists(sumConsecutiveRuns(i), th)+1;
                acumHist(sumConsecutiveRuns(i))=acumHist(sumConsecutiveRuns(i))+1;
            end
        end
    end
    [~, referenceLength]=max(acumHist);
    if nargin==3
        referenceLength=extReferenceLength;
    end
    
    %work with absolute indHists
    value = max(indHists, [], 1);
    i=0;
    idx=[];
    while isempty(idx)
        ref = min(nrows, referenceLength+i);        
        idx = find(value==indHists(ref,:)); 
        if ~isempty(idx)
            break;
        end
        ref = max(1, referenceLength-i);
        idx = find(value==indHists(ref,:));       
        i=i+1;
    end    
    
    [~, bestThr] = max(value(idx));
    threshold = idx(bestThr);
    return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = rle(x)
% data = rle(x) (de)compresses the data with the RLE-Algorithm
%   Compression:
%      if x is a numbervector data{1} contains the values
%      and data{2} contains the run lengths
%   Decompression:
%      if x is a cell array, data contains the uncompressed values
%      Version 1.0 by Stefan Eireiner (<a href="mailto:stefan-e@web.de?subject=rle">stefan-e@web.de</a>)
%      based on Code by Peter J. Acklam
%      last change 14.05.2004
if iscell(x) % decoding
	i = cumsum([ 1 x{2} ]);
	j = zeros(1, i(end)-1);
	j(i(1:end-1)) = 1;
	data = x{1}(cumsum(j));
else % encoding
	if size(x,1) > size(x,2), x = x'; end % if x is a column vector, transpose
    i = [ find(x(1:end-1) ~= x(2:end)) length(x) ];
	data{2} = diff([ 0 i ]);
	data{1} = x(i);
end


