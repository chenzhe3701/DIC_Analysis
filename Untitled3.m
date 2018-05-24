% for debug, plot the labeled segment


img = zeros(256);
img(ind_list_pre(p1:p2))=p1:p2;
figure;
imagesc(img);

img = zeros(256);
img(ind_list(q1:q2))=q1:q2;
figure;
imagesc(img);

img = zeros(256);
img(ind_list(skl_list==1)) = 10;
img(ind_list_pre(:)) = img(ind_list_pre(:)) + 1;
% img(ind_list_pre(642))=2;
img(ind_list(:)) = img(ind_list(:)) + 3;
% img(ind_list(641)) = img(ind_list(641)) + 1 ;
figure;
imagesc(img);


img = zeros(256);
img(ind_list_pre)=1;
% img(ind_list_pre([646,648,649,651:658]))=5;
img(ind_list_pre([1,10]))=5;
img(ind_list)=img(ind_list)+2;
figure;
imagesc(img);


ind_list_pre(p1-1),
ind_list_pre(p2+1),
ind_list(q1-1),
ind_list(q2+1),