p = uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\train_img_color\twin','training image Dir for twins');
p2 = uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\train_img_color\notwin','training image Dir for twins');

list = dir(p);
ii = 1;
filenames = [];
while ii<= length(list)
   if 1==list(ii).isdir
      list(ii) = [];
   else
       filenames{ii} = fullfile(list(ii).folder, list(ii).name);
       ii = ii+1;
   end
end
figure; montage(filenames,'size',[10,5]); %title('Traning Images for Twins');


list = dir(p2);
ii = 1;
filenames = [];
while ii<= length(list)
   if 1==list(ii).isdir
      list(ii) = [];
   else
       filenames{ii} = fullfile(list(ii).folder, list(ii).name);
       ii = ii+1;
   end
end
figure; montage(filenames,'size',[10,10]); %title('Traning Images for Non-Twins');