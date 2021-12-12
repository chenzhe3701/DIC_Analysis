% chenzhe, 2018--5-13
% Use the following as a reference for algorithm. 
% But I think there are a few errors there.
% Also note that his definition of the matrices is different from that in
% matlab built-in entropy. So, maybe it's better to make them the same.
% https://stackoverflow.com/questions/23691398/mutual-information-and-joint-entropy-of-two-images-matlab


function [mutualInformation,jointEntropy,jointHistogram,entropy1,entropy2] = mutual_info(img1,img2)

% need to preprocess img1 and img2 into unit8 with values within [0,255]

% Find unique (pixel grayscale) values of input image --> binC.
% And the id(indC) of the pixels in the uniqueValueVector(binC)
[binC,~,indC] = unique(img1(:));    % [C,IA,IC] = unique(A), so that A(IA)=C, and C(IC)=A
[binR,~,indR] = unique(img2(:));

nElements = numel(indR);    % number of elements in input image

% [A] count jointHistogram
jointHistogram = accumarray([indR,indC], 1);

% normalize to get probability. Eliminate zeros to to log2().
jointProb = jointHistogram/nElements;
indNonZero = (jointHistogram~=0);
jointProbNonZero = jointProb(indNonZero);

% [B] calculate joint entropy according to definition
jointEntropy = -sum(jointProbNonZero.*log2(jointProbNonZero));

% histogram and entropy of image 1
histogram1 = sum(jointHistogram, 1);
prob1 = histogram1/nElements;
indNonZero = (histogram1~=0);
prob1NonZero = prob1(indNonZero);
entropy1 = -sum(prob1NonZero.*log2(prob1NonZero));

% histogram and entropy of image 2
histogram2 = sum(jointHistogram, 2);
prob2 = histogram2/nElements;
indNonZero = (histogram2~=0);
prob2NonZero = prob2(indNonZero);
entropy2 = -sum(prob2NonZero.*log2(prob2NonZero));

% calculate mutual information
mutualInformation = entropy1 + entropy2 - jointEntropy;
end