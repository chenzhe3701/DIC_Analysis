% chenzhe, 2018-03-29
% get (almost) the same training images, when you have new cluster result,
% so made new images.
% do manually:
% (a) dira = folder of twin/notwin, training images used previously. The
% purpose is to look at the image, and find the image name.
% (b) dirb = image folder, for [all/categorized] your new images.
% (c) dirc = new [twin/notwin] destination folder

dira = 'D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\train_img_color\twin';
names = dir(dira);
dirb = 'D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\A_all_img_standard_silhouette\twin';
dirc = 'D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\B_train_img_standard_silhouette\twin';
ii = 0;
%%
ii = ii + 1
copyfile(fullfile(dirb,names(ii+2).name),fullfile(dirc,names(ii+2).name));
    
