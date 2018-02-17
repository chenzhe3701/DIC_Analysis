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
images = imageDatastore('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\train_img\','IncludeSubfolders',true,'LabelSource','foldernames');
[trainingImages,validationImages] = splitEachLabel(images,0.7,'randomized');

%%
% Load the pretrained AlexNet network. If Neural Network Toolbox(TM) Model
% _for AlexNet Network_ is not installed, then the software provides a
% download link. AlexNet is trained on more than one million images and can
% classify images into 1000 object categories.
net = alexnet;

%%
% To retrain AlexNet to classify new images, replace the last three layers
% of the network. Set the final fully connected layer to have the same size
% as the number of classes in the new data set (5, in this example). To
% learn faster in the new layers than in the transferred layers, increase
% the learning rate factors of the fully connected layer.
layersTransfer = net.Layers(1:end-3);
numClasses = numel(categories(trainingImages.Labels));
layers = [
    layersTransfer
    fullyConnectedLayer(numClasses,'WeightLearnRateFactor',20,'BiasLearnRateFactor',20)
    softmaxLayer
    classificationLayer];

%%
% Specify the training options, including learning rate, mini-batch size,
% and validation data.
options = trainingOptions('sgdm',...
    'MiniBatchSize',10,...
    'MaxEpochs',100,...
    'InitialLearnRate',1e-4,...
    'Verbose',false,...
    'Plots','training-progress',...
    'ValidationData',validationImages,...
    'ValidationFrequency',5);

%%
% Train the network using the training data.
netTransfer = trainNetwork(trainingImages,layers,options);

%%
% Classify the validation images using the fine-tuned network, and
% calculate the classification accuracy.
predictedLabels = classify(netTransfer,validationImages);
accuracy = mean(predictedLabels == validationImages.Labels)

%%
% For a more detailed transfer learning example, see
% <docid:nnet_examples.mw_10e960e5-a91b-48cb-8ea4-f3cd716c7e2d Transfer
% Learning Using AlexNet>.

imgPath = 'D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab\twin_img\'
imgs = dir(imgPath);
imgs = imgs(3:end);
%%
close all;
ind = randsample(length(imgs),1);
imshow(fullfile(imgs(ind).folder,imgs(ind).name));
I = imread(fullfile(imgs(ind).folder,imgs(ind).name));
predictedLabels = classify(netTransfer,I)





