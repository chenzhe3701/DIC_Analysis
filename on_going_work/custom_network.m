%% 2-layer
close all;
clc;

% build training data
x = -2:0.1:2;
t = tansig(x);
figure;
plot(x,t,'xr'); hold on;

% define network
net = network;
net.numInputs = 1;
net.numLayers = 2;
net.biasConnect = [1;0];
net.inputConnect = [1;0];
net.layerConnect = [0 0; 1 0];
net.outputConnect = [0 1];

net.inputs{1}.size = 1; % I think if input is a col vector, size = vector length  

net.layers{1}.size = 1;
net.layers{1}.transferFcn = 'tansig';
net.layers{1}.initFcn = 'initnw';
net.layers{2}.size = 1;

net.trainFcn = 'trainlm';
net.divideFcn = 'dividerand';
net.plotFcns = {'plotperform','plottrainstate'};

net = init(net);
sim(net,x)

net = train(net,x,t)

plot(x,net(x),'ok');
legend('truth','predicted','location','best');

%% 1 layer
close all; clc;

% build training data
x = -2:0.1:2;
t = tansig(x);
figure;
plot(x,t,'xr'); hold on;

% define network
net = network;
net.numInputs = 1;
net.numLayers = 1;
net.biasConnect = [1;];
net.inputConnect = [1;];
net.layerConnect = [0];
net.outputConnect = [1];

net.inputs{1}.size = 1; % I think if input is a col vector, size = vector length  

net.layers{1}.size = 1;
net.layers{1}.transferFcn = 'logsig';
net.layers{1}.initFcn = 'initnw';


net.trainFcn = 'trainlm';
net.divideFcn = 'dividerand';
net.plotFcns = {'plotperform','plottrainstate'};

net = init(net);
sim(net,x)

net = train(net,x,t)
plot(x,net(x),'ok');
legend('truth','predicted','location','best');

%% input is 2x1 vector, dimension=2
close all; clc;

% build training data
x = -2:0.1:2;
t = tansig(x) + (rand(1,41)-0.5)/10;

xx = [x;x];
xt = tansig(x) + (rand(1,41)-0.5)/10;
yt = -tansig(x) + (rand(1,41)-0.5)/10;
tt = [xt;yt];

figure;
plot(x,tt,'xr'); hold on;

% define network
net = network;
net.numInputs = 1;
net.numLayers = 1;
net.biasConnect = [1;];
net.inputConnect = [1;];
net.layerConnect = [0];
net.outputConnect = [1];

net.inputs{1}.size = 2; % I think if input is a col vector, size = vector length  

net.layers{1}.size = 2;
net.layers{1}.transferFcn = 'tansig';
net.layers{1}.initFcn = 'initnw';


net.trainFcn = 'trainlm';
net.divideFcn = 'dividerand';
net.plotFcns = {'plotperform','plottrainstate'};

net = init(net);

net = train(net,xx,tt);
plot(x,net(xx),'ok');
legend('truth','predicted','location','best');

net.iw{1}   % row 1 sum = 1, row 2 sum = -1
net.b{1}
%% try 2 inputs. For a layer, if have 2 neurons, the output vector is of size 2x1
close all; clc;

% build training data
x = -2:0.1:2;
t = tansig(x) + (rand(1,41)-0.5)/10;

xt = tansig(x) + (rand(1,41)-0.5)/10;
yt = -tansig(x) + (rand(1,41)-0.5)/10;


figure;
plot(x,[xt;yt],'xr'); hold on;

% define network
net = network;
net.numInputs = 2;
net.numLayers = 1;
net.biasConnect = [1;];
net.inputConnect = [1 1;];
net.layerConnect = [0];
net.outputConnect = [1];

net.inputs{1}.size = 1; % I think if input is a col vector, size = vector length  
net.inputs{2}.size = 1;

net.layers{1}.size = 2;
net.layers{1}.transferFcn = 'tansig';
net.layers{1}.initFcn = 'initnw';

% --> one layer can only have one output (although it can be multi dimension ?)  
net.trainFcn = 'trainlm';
net.divideFcn = 'dividerand';
net.plotFcns = {'plotperform','plottrainstate'};

net.iw{1} = [3.82; -2.8];
net.iw{2} = [-1.82; 0.81];

XX = [x;x];
XX = mat2cell(XX, [1,1], ones(1,length(x)));
TT = [xt;yt];
TT = mat2cell(TT, 2, ones(1,length(xt)));
net = train(net,XX,TT);
pp = cell2mat(net(TT));
plot(x,pp,'ok');
legend('truth','predicted','location','best');




