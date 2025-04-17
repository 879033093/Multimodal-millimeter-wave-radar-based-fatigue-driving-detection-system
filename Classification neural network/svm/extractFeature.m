function [trainingFeatures,trainingLabels,testFeatures,testLabels]=extractFeature(trainingSet,testSet)

img = read(trainingSet(1), 1);

img=rgb2gray(img);

lvl = graythresh(img);
img = im2bw(img, lvl);
img=imresize(img,[256 256]);
cellSize = [4 4];
[hog_feature, vis_hog] = extractHOGFeatures(img,'CellSize',cellSize);
glcm_feature = getGLCMFeatures(img);
SizeOfFeature = length(hog_feature)+ length(glcm_feature);


trainingFeatures = [];
trainingLabels   = [];
for digit = 1:numel(trainingSet)       
    numImages = trainingSet(digit).Count;
    features  = zeros(numImages, SizeOfFeature, 'single');
    for i = 1:numImages
        img = read(trainingSet(digit), i);
        
        img=rgb2gray(img);                
        glcm_feature = getGLCMFeatures(img);  
       
        lvl = graythresh(img);            
        img = im2bw(img, lvl);            
        img=imresize(img,[256 256]);
        
        [hog_feature, vis_hog] = extractHOGFeatures(img,'CellSize',cellSize);
        
        features(i, :) = [hog_feature glcm_feature];
    end
    
    labels = repmat(trainingSet(digit).Description, numImages, 1);  
    
    trainingFeatures = [trainingFeatures; features];
    trainingLabels   = [trainingLabels; labels];       
end



testFeatures = [];
testLabels   = [];
for digit = 1:numel(testSet)
           
    numImages = testSet(digit).Count;
    
    features  = zeros(numImages, SizeOfFeature, 'single');
    
    for i = 1:numImages
        
        img = read(testSet(digit), i);
        
        img=rgb2gray(img);
        glcm_feature = getGLCMFeatures(img);
        
        lvl = graythresh(img);
        img = im2bw(img, lvl);
        img=imresize(img,[256 256]);
        [hog_4x4, vis4x4] = extractHOGFeatures(img,'CellSize',cellSize);
        features(i, :) = [hog_4x4 glcm_feature];
    end
    
    
    labels = repmat(testSet(digit).Description, numImages, 1);
        
    testFeatures = [testFeatures; features];
    testLabels=[testLabels; labels];
        
end
end
