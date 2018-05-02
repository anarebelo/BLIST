function ME_niblack(directory)

% directory of grey scale images
dt = strcat(directory, '/images/binarizations/niblack');
d=dir(dt);
d=struct2cell(d);
names=d(1,3:end,:);

% initiate counter
counter =1;

%initiate file ids
fileME = strcat(directory, '/results/ME/niblack.txt');
fileMOPx = strcat(directory, '/results/MOPx/niblack.txt');
fileFOPx = strcat(directory, '/results/FOPx/niblack.txt');

fid = fopen(fileME,'wt');
fidm = fopen(fileMOPx,'wt');
fidf = fopen(fileFOPx,'wt');

MEmat = zeros(1,65);
MOPxmat = zeros(1,65);
FOPxmat = zeros(1,65);


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
            sprintf('Calculating ME, MOPX and FOPx for image number %d',...
                counter)
            
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

            
            Fbin = find(img(:) == 0)';
            Fgt = find(gt(:) == 0)';
            

%%%%%%%%%%% RUN ME %%%%%%%%%%%
            ME = sum(xor(img(:), gt(:)));
            % in percentage
            dim = size(img(:));
            ME = (ME/dim(1))*100;

%%%%%%%%%%% RUN MOPx %%%%%%%%%%%
            m1 = length(Fgt)-length(intersect(Fbin, Fgt));
            m2 = length(Fgt);
            % in percentage
            MOPx = (m1/m2)*100;

%%%%%%%%%%% RUN FOPx %%%%%%%%%%%
            f1 = length(Fbin)-length(intersect(Fbin, Fgt));
            f2 = length(Fbin);
            % in percentage
            FOPx = (f1/f2)*100;

            MEmat(counter) = ME;
            MOPxmat(counter) = MOPx;
            FOPxmat(counter) = FOPx;
            
            % print ME, MOPx and FOPx
            sprintf('ME: %d', ME)
            sprintf('MOPx: %d', MOPx)
            sprintf('FOPx: %d', FOPx)
            

            % save values in file
            if counter < 10
                fprintf(fid, '0%d: %d\n', counter, ME);
                fprintf(fidm, '0%d: %d\n', counter, MOPx);
                fprintf(fidf, '0%d: %d\n', counter, FOPx);
            else
                fprintf(fid, '%d: %d\n', counter, ME);
                fprintf(fidm, '%d: %d\n', counter, MOPx);
                fprintf(fidf, '%d: %d\n', counter, FOPx);
            end

        end
        % next image
        counter = counter + 1;
    end
end

% calculate average ratio
len = length(find(MEmat~=0));
ratio = sum(MEmat)/len;

len = length(find(MOPxmat~=0));
ratiom = sum(MOPxmat)/len;

len = length(find(FOPxmat~=0));
ratiof = sum(FOPxmat)/len;

% print ratio
sprintf('ME Ratio: %d', ratio)
sprintf('MOPx Ratio: %d', ratiom)
sprintf('FOPx Ratio: %d', ratiof)

% write in file
fprintf(fid, 'ME Ratio: %d\n', ratio);
fprintf(fidm, 'MOPx Ratio: %d\n', ratiom);
fprintf(fidf, 'FOPx Ratio: %d\n', ratiof);


fclose(fid);
fclose(fidm);
fclose(fidf);