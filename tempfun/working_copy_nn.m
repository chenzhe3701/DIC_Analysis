
clc;
sx = 3;
sy = 2;
sz = 1;
sp = 1;

layers = [imageInputLayer([sy,sx,sz],'name','inputa','Normalization','none');
    fullyConnectedLayer(1,'name','fca');
    regressionLayer('name','rga')];
lgraph = layerGraph(layers);

x = ones(sy,sx,sz,sp);
y = 2*ones(sp,1);
option = trainingOptions('sgdm');

img = [1 0 0; 0 0 0];


% nett = trainNetwork(x,y,layers,option);
% nett.predict(img);



% layers(2) = nett.Layers(2);   % this is required if want to modify to use, because parameteres were not properly inferred. 
% however, using the following workaround, it's OK:
weights = [1 2 3; 4 5 6];
% layers(2).Weights = weights(:)';
% layers(2).Bias = [0];
net = trainNetwork_init(layers);
a = net.predict(img)


net.activations(img,2)
net.activations(img,3)