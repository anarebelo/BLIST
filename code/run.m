function RUN()

clc
clear all


% path to RR root
directory = ('..');

%%%%% Create output dirs %%%%%


dir=strcat(directory, '/images/binarizations/');
% if dir does not exist create a new one
if isdir(dir) == 0
	mkdir(dir);
end

dir=strcat(directory, '/results/');
% if dir does not exist create a new one
if isdir(dir) == 0
	mkdir(dir);
end

dir=strcat(directory, '/results/global_thresholds/');
% if dir does not exist create a new one
if isdir(dir) == 0
	mkdir(dir);
end

%%%%% Run each method for all 65 images %%%%%
% run_adblist(directory);
% run_adotsu(directory);
% run_bernsen(directory);
% run_blist_pairs(directory);
% run_blist_triplets(directory);
% run_chen(directory);
% run_kapur(directory);
% run_khashman(directory);
% run_niblack(directory);
 run_otsu(directory);
% run_sahoo(directory);
% run_tsai(directory);
% run_tsallis(directory);
% run_yanowitz(directory);


%%%%% Create output dirs %%%%%

dir=strcat(directory, '/results/ME');
% if dir does not exist create a new one
if isdir(dir) == 0
	mkdir(dir);
end

dir=strcat(directory, '/results/MOPx');
% if dir does not exist create a new one
if isdir(dir) == 0
	mkdir(dir);
end

dir=strcat(directory, '/results/FOPx');
% if dir does not exist create a new one
if isdir(dir) == 0
	mkdir(dir);
end

%%%%% Run Misclassification Error for all ground-truth images %%%%%
% ME_adblist(directory);
% ME_adotsu(directory);
% ME_bernsen(directory);
% ME_blist_pairs(directory);
% ME_blist_triplets(directory);
% ME_chen(directory);
% ME_kapur(directory);
% ME_khashman(directory);
% ME_niblack(directory);
 ME_otsu(directory);
% ME_sahoo(directory);
% ME_tsai(directory);
% ME_tsallis(directory);
% ME_yanowitz(directory);