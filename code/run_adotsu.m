function run_adotsu(directory)

% create output dir
filedir=strcat(directory, '/images/binarizations/adotsu/');
% if dir exists replace it with a new one
if isdir(filedir) == 1
    rmdir(filedir, 's');
    mkdir(filedir);
else
	mkdir(filedir);
end

% directory of grey scale images
dt = strcat(directory, '/images/grey');
d=dir(dt);
d=struct2cell(d);
names=d(1,3:end,:);

% initiate counter
counter =1;

for i=1:size(names,2)
    
    % add image name to str
    str=strcat(dt, '/');
    image=strcat(str, char(names(i)));
    % check if file is .DS_Store
    ds_store = strcat(str, '.DS_Store');
    if(strcmpi(image, ds_store) == 0)
        
        % print image number
        sprintf('Running for image number %d', counter)
        
        % read grey scale image to variable img
        img=imread(image);
        
%%%%%%% RUN METHOD %%%%%%%
        Ims = ADOTSU(img);
                
        % determine name of image to save
        if counter < 10
            filename=sprintf('img0%d.png', counter);
        else
            filename=sprintf('img%d.png', counter);
        end
        filename=strcat(filedir, filename);
        % save image
        imwrite(Ims, filename,'png');

        % next image
        counter = counter + 1;
    end
end