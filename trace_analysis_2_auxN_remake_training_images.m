% chenzhe, 2018-03-29
% get (almost) the same training images, when you have new cluster result,
% so made new images.
% do manually:
% (a) dir_ref = folder of twin/notwin, training images used previously. The
% purpose is to look at the image, and find the image name.
% (b) dir_source = image folder, for [all/categorized] your new images.
% (c) dir_target = new [twin/notwin] destination folder
%
% chenzhe, 2018-08-31, change script name, make it an auxiliary function
% The purpose of this code is to make sure that the training images are
% kept the same, if we want to compare the effect of image processing.
% 

%% copy twin images for training
dir_ref = 'D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\train_img_color\twin';
names = dir(dir_ref);

%% copy twin images
dir_source = 'D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab_after_realign\all_img_colored';
dir_target = 'D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab_after_realign\train_img_colored\twin';
ii = 0;
%% copy one-by-one. Click, observe (and correct), then continue
ii = ii + 1
copyfile(fullfile(dir_source,names(ii+2).name),fullfile(dir_target,names(ii+2).name));
    






%% non-twin
dir_ref = 'D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\train_img_color\notwin';
names = dir(dir_ref);

%% copy non-twin images
dir_source = 'D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab_after_realign\all_img_colored';
dir_target = 'D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab_after_realign\train_img_colored\notwin';
ii = 0;
%% copy one-by-one. Click, observe (and correct), then continue
ii = ii + 1
copyfile(fullfile(dir_source,names(ii+2).name),fullfile(dir_target,names(ii+2).name));