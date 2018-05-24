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
img(ind_list_pre(p1:p2))=1;
img(ind_list(q1:q2))=img(ind_list(q1:q2))+2;
figure;
imagesc(img);


ind_list_pre(p1-1),
ind_list_pre(p2+1),
ind_list(q1-1),
ind_list(q2+1),