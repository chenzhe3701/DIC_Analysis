% chenzhe, 2018-02-12
% based on alexnet tutorial, try cnn twin area

%% Get Started with Transfer Learning
% This example shows how to use transfer learning to retrain AlexNet, a
% pretrained convolutional neural network, to classify a new set of images.
% Try this example to see how simple it is to get started with deep
% learning in MATLAB(R).

%%
% Unzip and load the new images as an image datastore. Divide the data into
% training and validation data sets. Use 70% of the images for training and
% 30% for validation.
saveDataPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','choose a path [to save the]/[of the saved] processed data, or WS, or etc.'),'\'];

img_size = 227; % 227 for alexnet, 224 for vgg, googlenet
trainFolder = ['D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\train_img_',num2str(img_size),'_labeled\'];
trainFolder = uigetdir(trainFolder,'select parent folder of training images');

trainImages = imageDatastore(trainFolder,'IncludeSubfolders',true,'LabelSource','foldernames');
[trainingImages,validationImages] = splitEachLabel(trainImages,0.7,'randomized');

allImgPath = ['D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\twin_img_',num2str(img_size),'_labeled\'];
allImgPath = uigetdir(allImgPath,'choose all image path');
imgs = dir([allImgPath,'\**\*.tif']);

%%
% Load the pretrained AlexNet network. If Neural Network Toolbox(TM) Model
% _for AlexNet Network_ is not installed, then the software provides a
% download link. AlexNet is trained on more than one million images and can
% classify images into 1000 object categories.

net = alexnet;
% net = vgg16;
%%
% To retrain AlexNet to classify new images, replace the last three layers
% of the network. Set the final fully connected layer to have the same size
% as the number of classes in the new data set (2, in this example). To
% learn faster in the new layers than in the transferred layers, increase
% the learning rate factors of the fully connected layer.
layersTransfer = net.Layers(1:end-3);
numClasses = numel(categories(trainingImages.Labels));
layers = [
    layersTransfer
    fullyConnectedLayer(numClasses,'WeightLearnRateFactor',5,'BiasLearnRateFactor',5)
    softmaxLayer
    classificationLayer];

%%
% Specify the training options, including learning rate, mini-batch size,
% and validation data.
options = trainingOptions('sgdm',...
    'MiniBatchSize',32,...
    'MaxEpochs',100,...
    'InitialLearnRate',1e-4,...
    'Verbose',true,...
    'Plots','training-progress',...
    'ValidationData',validationImages,...
    'ValidationFrequency',5,...
    'ValidationPatience',10,...
    'Shuffle','every-epoch',...
    'ExecutionEnvironment','gpu');      % can try cpu

%%
% Train the network using the training data.
netTransfer = trainNetwork(trainingImages,layers,options);

%%
% Classify the validation images using the fine-tuned network, and
% calculate the classification accuracy.
predictedLabels = classify(netTransfer,validationImages);
accuracy = mean(predictedLabels == validationImages.Labels)
%%
timeStr = datestr(now,'yyyymmdd_HHMM');
save([saveDataPath,'\trainedAlexNet.mat'],'netTransfer');
save([saveDataPath,'\trainedAlexNet_',timeStr,'.mat'],'netTransfer');

%%
close all;
ind = randsample(length(imgs),1);
imshow(fullfile(imgs(ind).folder,imgs(ind).name));
I = imread(fullfile(imgs(ind).folder,imgs(ind).name));
predictedLabels = classify(netTransfer,I)
predict(netTransfer,I)






