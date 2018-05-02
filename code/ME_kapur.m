function ME_kapur(directory)

% directory of grey scale images
dt = strcat(directory, '/images/binarizations/kapur');
d=dir(dt);
d=struct2cell(d);
names=d(1,3:end,:);

% initiate counter
counter =1;

%initiate file id to write ME 
file = strcat(directory, '/results/ME/kapur.txt');
fid = fopen(file,'wt');

MEmat = zeros(1,65);

for i=1:size(names,2)
    
    % add image name to str
    str=strcat(dt, '/');
    image=strcat(str, char(names(i)));
    % check if file is .DS_Store
    ds_store = strcat(str, '.DS_Store');
    if(strcmpi(image, ds_store) == 0)

        % can only calculate ME for images with ground-truths
        if counter==1 || counter == 2 || counter == 10 || counter==12 ...
                || counter==27 || counter ==36 || counter==46 || ...
                counter==50 || counter==53 || counter==63    
          
            % print image number
            sprintf('Calculating ME for image number %d', counter)
            
            % read grey scale image to variable img
            img=imread(image);
            
            % fetch ground-truth image name
            if counter < 10
                gt_image=sprintf('img0%d.png', counter);
            else
                gt_image=sprintf('img%d.png', counter);
            end
            
            gt_image=strcat('/images/ground_truth/',gt_image);
            gt_image=strcat(directory, gt_image);
            gt=imread(gt_image);

%%%%%%%%%%% RUN ME %%%%%%%%%%%
            ME = sum(xor(img(:), gt(:)));
            % in percentage
            dim = size(img(:));
            ME = (ME/dim(1))*100;

            MEmat(counter) = ME;
            
            % print ME
            sprintf('ME: %d', ME)

            % save ME in file
            if counter < 10
                fprintf(fid, '0%d: %d\n', counter, ME);
            else
                fprintf(fid, '%d: %d\n', counter, ME);
            end

        end
        % next image
        counter = counter + 1;
    end
end

% calculate average ratio
len = length(find(MEmat~=0));
ratio = sum(MEmat)/len;

% print ratio
sprintf('ME Ratio: %d', ratio)

% write in file
fprintf(fid, 'ME Ratio: %d\n', ratio);


fclose(fid);