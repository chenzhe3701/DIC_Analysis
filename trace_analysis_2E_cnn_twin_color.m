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
    fullyConnectedLayer(numClasses,'WeightLearnRateFactor',10,'BiasLearnRateFactor',10)
    softmaxLayer
    classificationLayer];

%%
% Specify the training options, including learning rate, mini-batch size,
% and validation data.
options = trainingOptions('sgdm',...
    'MiniBatchSize',25,...
    'MaxEpochs',100,...
    'InitialLearnRate',1e-4,...
    'Verbose',true,...
    'Plots','training-progress',...
    'ValidationData',validationImages,...
    'ValidationFrequency',5,...
    'ValidationPatience',5,...
    'Shuffle','every-epoch',...
    'ExecutionEnvironment','gpu');      % can try cpu

%%
% Train the network using the training data.
rng('default');
netTransfer = trainNetwork(trainingImages,layers,options);

%%
% Classify the validation images using the fine-tuned network, and
% calculate the classification accuracy.
predictedLabels = classify(netTransfer,validationImages);
accuracy = mean(predictedLabels == validationImages.Labels)

%% can also test on the groud truth data
truthFolder_twin = uigetdir(allImgPath,'select parent folder 1 of twin ground truth images');
truthImages_twin = imageDatastore({truthFolder_twin},'IncludeSubfolders',true,'LabelSource','foldernames');

truthFolder_notwin = uigetdir(allImgPath,'select parent folder 1 of nontwin ground truth images');
truthImages_notwin = imageDatastore({truthFolder_notwin},'IncludeSubfolders',true,'LabelSource','foldernames');

predictedLabels = classify(netTransfer,truthImages_twin);
disp('accuracy of twin:');
accuracy = mean(predictedLabels == truthImages_twin.Labels)

predictedLabels = classify(netTransfer,truthImages_notwin);
disp('accuracy of notwin:');
accuracy = mean(predictedLabels == truthImages_notwin.Labels)
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

%% This part copies all images into different folders according to classification result, to evaluate
% modify / or keep an eye on these settings for the specific sample to analyze  ------------------------------------------------------------------------------------
addChenFunction;
[fileSetting,pathSetting] = uigetfile('','select setting file which contains sampleName, stopNames, FOVs, translations, etc');
load_settings([pathSetting,fileSetting],'sampleName','cpEBSD','cpSEM','sampleMaterial','stressTensor');
STOP = {'0','1','2','3','4','5','6','7'};
B=1;    % 0-based B=1.  1-based B=0.
iE_start = 2;   % elongation levels to analyze. 0-based.
iE_stop = 5;

classifyImgPath = [uigetdir('D:\WE43_T6_C1_insitu_compression\Analysis_by_Matlab','Choose/make a parent path for output classified images.'),'\'];
mkdir(classifyImgPath,'Classified_twin');
mkdir(classifyImgPath,'Classified_notwin');


[fileNet,pathNet] = uigetfile([saveDataPath,'\trainedAlexNet.mat'],'select the trained cnn net');
net = load([pathNet,fileNet],'netTransfer');
net = net.netTransfer;
kk = find(cellfun(@(x) strcmp(x,'twin'),net.Layers(end).ClassNames));

for iE = iE_start:iE_stop
    fName_c2t_result = [sampleName,'_s',num2str(STOP{iE+B}),'_cluster_to_twin_result.mat'];
    load([saveDataPath,fName_c2t_result],'stru');
    
    hWaitbar = waitbar(0,'Matching cluster with twin system ...');
    for iS =1:length(stru)
        
        ID_current = stru(iS).gID;
        % ==== change cluster number into twin system number, or 0
        nCluster = length(stru(iS).cLabel);
        
        for iCluster = 1:nCluster
            cNum = stru(iS).cLabel(iCluster);
            
            % cnn identify twin
            imgName = ([num2str(iE*100000 + ID_current*10 + cNum),'.tif']);
            try
                I = imread(fullfile(allImgPath,imgName));
                predictedLabels = classify(net,I);
                predictedProb = predict(net,I);
                
                prob = predictedProb(kk);
                
                
                % move classified twin images
                if prob>0.5
                    copyfile(fullfile(allImgPath,imgName),fullfile(classifyImgPath,'Classified_twin',imgName));
                else
                    copyfile(fullfile(allImgPath,imgName),fullfile(classifyImgPath,'Classified_notwin',imgName));
                end
            end
            
        end

        
        waitbar(iS/length(stru), hWaitbar);
    end
    
    
    try
        close(hWaitbar);
    catch
    end
    
end






