clear;
dir=('C:\Users\Windows 10\Desktop\data_set\flower_data\train');
testdir=('C:\Users\Windows 10\flower_data\val');
trainingSet = imageSet(dir,'recursive');
testSet = imageSet(testdir,'recursive');


[trainingFeatures,trainingLabels,testFeatures,testLabels]=extractFeature(trainingSet,testSet);


classifier = fitcecoc(trainingFeatures, trainingLabels);
save classifier.mat classifier;


predictedLabels = predict(classifier, testFeatures);



confMat=confusionmat(testLabels, predictedLabels)
accuracy=(confMat(1,1)/sum(confMat(1,:))+confMat(2,2)/sum(confMat(2,:))+...
    confMat(3,3)/sum(confMat(3,:)))/3


