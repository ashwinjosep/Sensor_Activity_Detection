% get list of all users in groundtruth folder
groundTruthUsers = dir('groundTruth/');
indices = ismember({groundTruthUsers.name},{'.', '..', '.DS_Store'});
groundTruthUsers = groundTruthUsers(~indices);

% get list of all users in myodata folder
myoDataUsers = dir('MyoData/');
indices = ismember({myoDataUsers.name},{'.', '..', '.DS_Store'});
myoDataUsers = myoDataUsers(~indices);

%Global Data tables
forkEMG = [];
forkIMU = [];
spoonEMG = [];
spoonIMU = [];
activityRow = zeros(20,1);
featureEatingUser = [];
featureNonEatingUser = [];
featureArray = [];
for i = 1 : length(groundTruthUsers)
    
    % Reading users fork data from GroundTruth
    folderName = "groundTruth/"+groundTruthUsers(i).name+"/fork/*.txt";
    files = dir(folderName);
    fileName = files(1).folder + "/" + files(1).name;
    
    % Fork Ground Data
    forkDataTruth = csvread(fileName);

    % Reading user's fork data from MyoData
    folderName = "MyoData/"+groundTruthUsers(i).name+"/fork/*.txt";
    files = dir(folderName);
    
    %Reading fork EMG Sensor Data
    fileName = files(1).folder + "/" + files(1).name;
    forkDataEMG = csvread(fileName);
    
    %Reading fork IMU Sensor Data
    fileName = files(2).folder + "/" + files(2).name;
    forkDataIMU = csvread(fileName);    
    
    % Reading users spoon data from GroundTruth
    folderName = "groundTruth/"+groundTruthUsers(i).name+"/spoon/*.txt";
    files = dir(folderName);
    fileName = files(1).folder + "/" + files(1).name;
    
    % Spoon Ground Data
    spoonDataTruth = csvread(fileName);
    
    % Reading user's spoon data from MyoData
    folderName = "MyoData/"+groundTruthUsers(i).name+"/spoon/*.txt";
    files = dir(folderName);
    
    % Spoon EMG Data
    fileName = files(1).folder + "/" + files(1).name;
    spoonDataEMG = csvread(fileName);
    
    % Spoon IMU Data
    fileName = files(2).folder + "/" + files(2).name;
    spoonDataIMU = csvread(fileName);
    
    %Get Durations of food eating
    forkFoodDuration = (forkDataTruth(:,2)-forkDataTruth(:,1))*1000/30;
    spoonFoodDuration = (spoonDataTruth(:,2)-spoonDataTruth(:,1))*1000/30;
    
    %Get Durations of Non Food Eating Activity
    forkNoFoodDuration = ([forkDataTruth(:,1);0]-[0;forkDataTruth(:,2)])*1000/30;
    forkNoFoodDuration = forkNoFoodDuration(2:end-1);
    spoonNoFoodDuration = ([spoonDataTruth(:,1);0]-[0;spoonDataTruth(:,2)])*1000/30;
    spoonNoFoodDuration = spoonNoFoodDuration(2:end-1);
    
    %Comparing groundtruth and sensor data for forks
    forkDataEMG(:,10) = 0;
    forkDataIMU(:,12)=0;
    [forkFoodDurationRows,forkFoodDurationColumns] = size(forkFoodDuration);
    [forkNoFoodDurationRows,forkNoFoodDurationColumns] = size(forkNoFoodDuration);
    [forkDataEMGRows, forkDataEMGColumns] = size(forkDataEMG);
    [forkDataIMURows, forkDataIMUColumns] = size(forkDataIMU);
    j = forkFoodDurationRows;
    beginTimeStamp = forkDataEMG(forkDataEMGRows,1);
    
    % Get all food activity time ranges for fork data
    while(j >= 1)
        endTimeStamp = beginTimeStamp - forkFoodDuration(j,1);
        forkFoodDuration(j,2) = beginTimeStamp;
        forkFoodDuration(j,3) = endTimeStamp;
        if(j>1)
            beginTimeStamp = endTimeStamp - forkNoFoodDuration(j-1,1);
        end
        j = j-1;
    end
    
    % Get all non food activity time ranges for fork data
    beginTimeStamp = forkFoodDuration(forkFoodDurationRows,3);
    j = forkNoFoodDurationRows;
    while(j >= 1)
        endTimeStamp = beginTimeStamp - forkNoFoodDuration(j,1);
        forkNoFoodDuration(j,2) = beginTimeStamp;
        forkNoFoodDuration(j,3) = endTimeStamp;
        if(j>1)
            beginTimeStamp = endTimeStamp - forkFoodDuration(j+1,1);
        end
        j = j-1;
    end
    
    % Now that fork food durations have been recorded classify each
    % eating activity and assign their class labels
    for k=1:forkFoodDurationRows
        
        % Classify each instance of the activity data for EMG sensors
        idx = (and((forkDataEMG(:,1)>=forkFoodDuration(k,3)),(forkDataEMG(:,1)<=forkFoodDuration(k,2))));
        forkEMGActivity = forkDataEMG(idx,:);
        
        % Classify each instance of the activity data for EMG sensors    
        idx = (and((forkDataIMU(:,1)>=forkFoodDuration(k,3)),(forkDataIMU(:,1)<=forkFoodDuration(k,2))));
        forkIMUActivity = forkDataIMU(idx,:);
        
        %Apply feature transformation function
        if(and((~isempty(forkEMGActivity)),(~isempty(forkIMUActivity))))
            activityRow = featureExtract(forkEMGActivity,forkIMUActivity);
        end

        %Assign class label and user for the activity
        activityRow(:,19) = str2double(extractBetween(groundTruthUsers(i).name,5,6));
        activityRow(:,20) = 1;
        
        %Append to feature array
        featureArray = [featureArray;activityRow];
    end
    
    % Now that fork food durations have been recorded classify each
    % activity and assign their class labels
    for k=1:forkNoFoodDurationRows
        
        % Classify each instance of the activity data for EMG sensors
        idx = (and((forkDataEMG(:,1)>=forkNoFoodDuration(k,3)),(forkDataEMG(:,1)<=forkNoFoodDuration(k,2))));
        forkEMGActivity = forkDataEMG(idx,:);
        
        % Classify each instance of the activity data for EMG sensors    
        idx = (and((forkDataIMU(:,1)>=forkNoFoodDuration(k,3)),(forkDataIMU(:,1)<=forkNoFoodDuration(k,2))));
        forkIMUActivity = forkDataIMU(idx,:);
        
        %Apply feature transformation function
        if(and((~isempty(forkEMGActivity)),(~isempty(forkIMUActivity))))
            activityRow = featureExtract(forkEMGActivity,forkIMUActivity);
        end

        %Assign class label and user for the activity
        activityRow(:,19) = str2double(extractBetween(groundTruthUsers(i).name,5,6));
        activityRow(:,20) = 0;
        
        %Append to feature array
        featureArray = [featureArray;activityRow];
    end
    
% Performing the same for each activity. So no longer required 
%     % Assign Class label in data
%     for j = 1:forkDataEMGRows
%         for k = 1:forkFoodDurationRows
%             if(and((forkDataEMG(j,1)>=forkFoodDuration(k,3)),(forkDataEMG(j,1)<=forkFoodDuration(k,2))))
%                 forkDataEMG(j,10) = 1;
%                 break;
%             end
%         end
%     end
%     
%     % Assign Class label in data
%     for j = 1:forkDataIMURows
%         for k = 1:forkFoodDurationRows
%             if(and((forkDataIMU(j,1)>=forkFoodDuration(k,3)),(forkDataIMU(j,1)<=forkFoodDuration(k,2))))
%                 forkDataIMU(j,12) = 1;
%                 break;
%             end
%         end
%     end
    
%     %Comparing groundtruth and sensor data for spoons
%     spoonDataEMG(:,10) = 0;
%     spoonDataIMU(:,12) = 0;
%     [spoonFoodDurationRows,spoonFoodDurationColumns] = size(spoonFoodDuration);
%     [spoonDataEMGRows, spoonDataEMGColumns] = size(spoonDataEMG);
%     [spoonDataIMURows, spoonDataIMUColumns] = size(spoonDataIMU);
%     
%     j = spoonFoodDurationRows;
%     beginTimeStamp = spoonDataEMG(spoonDataEMGRows,1);
%     
%     % Get all food activity time ranges for spoon data
%     while(j >= 1)
%         endTimeStamp = beginTimeStamp - spoonFoodDuration(j,1);
%         spoonFoodDuration(j,2) = beginTimeStamp;
%         spoonFoodDuration(j,3) = endTimeStamp;
%         if(j>1)
%             beginTimeStamp = endTimeStamp - spoonNoFoodDuration(j-1,1);
%         end
%         j = j-1;
%     end
%     
%     % Assign Class label in data
%     for j = 1:spoonDataEMGRows
%         for k = 1:spoonFoodDurationRows
%             if(and((spoonDataEMG(j,1)>=spoonFoodDuration(k,3)),(spoonDataEMG(j,1)<=spoonFoodDuration(k,2))))
%                 spoonDataEMG(j,10) = 1;
%                 break;
%             end
%         end
%     end
%     
%     % Assign Class label in data
%     for j = 1:spoonDataIMURows
%         for k = 1:spoonFoodDurationRows
%             if(and((spoonDataIMU(j,1)>=spoonFoodDuration(k,3)),(spoonDataIMU(j,1)<=spoonFoodDuration(k,2))))
%                 spoonDataIMU(j,12) = 1;
%                 break;
%             end
%         end
%     end

end

% PCA
[coeff, score, latent] = pca(featureArray(:,1:18));
% Score contains the new feature matrix

%User dependent 
metrics = [];
users = unique(featureArray(:,19));
for i=1:size(users,1)
    %Prepare training and test data
    data = [score(:,1:3) featureArray(:,19) featureArray(:,20)];
    data = data(data(:,4)==users(i),:);
    % Randomize and partition
    data = data(randperm(size(data,1)),:);
    train_data = data(1:ceil(0.6*size(data,1)),:); 
    test_data = data(ceil(0.6*size(data,1))+1:end,:);
   
    %SVM
    SVMModel = fitcsvm(train_data(:,1:3),train_data(:,5));
    label = predict(SVMModel,test_data(:,1:3));
    %Get metrics
    [fscoreSVM, PrecisionSVM, RecallSVM] = compareResults(test_data(:,5), label);

    %Decision Tree
    tree = fitctree(train_data(:,1:3),train_data(:,5));
    label = predict(tree,test_data(:,1:3));
    %Get metrics
    [fscoreDT, PrecisionDT, RecallDT] = compareResults(test_data(:,5), label);    
    
    metrics = [metrics; [fscoreSVM, PrecisionSVM, RecallSVM, users(i), "SVM"]];
    metrics = [metrics; [fscoreDT, PrecisionDT, RecallDT, users(i), "DT"]];
    
    nnX = data(:,1:3);
    nnY = data(:,5);
  
end

%User Independent testing
%Prepare training and test data
data = [score(:,1:3) featureArray(:,20)];
% Randomize and partition
data = data(randperm(size(data,1)),:);
train_data = data(1:ceil(0.6*size(data,1)),:); 
test_data = data(ceil(0.6*size(data,1))+1:end,:);

%SVM
SVMModel = fitcsvm(train_data(:,1:3),train_data(:,4));
label = predict(SVMModel,test_data(:,1:3));
%Get metrics
[fscoreSVM, PrecisionSVM, RecallSVM] = compareResults(test_data(:,4), label);

%Decision Tree
tree = fitctree(train_data(:,1:3),train_data(:,4));
label = predict(tree,test_data(:,1:3));
%Get metrics
[fscoreDT, PrecisionDT, RecallDT] = compareResults(test_data(:,4), label);

%Neural Network
nnX = data(:,1:3);
nnY = data(:,4);

% gscatter(train_data(:,1),train_data(:,2),train_data(:,4));
% sv = SVMModel.SupportVectors;
% hold on
% plot(sv(:,1 ),sv(:,2),'ko','MarkerSize',10)
% legend('eating','non-eating','Support Vector')
% hold off

function [fscore, Precision, Recall] = compareResults(target, output)
    C = confusionmat(target, output);
    true_positive = C(1,1);
    false_positive = C(2,1);
    false_negative = C(1,2);
    Precision = true_positive / (true_positive + false_positive);
    Recall = true_positive / (true_positive + false_negative);
    fscore = 2 * Precision * Recall / (Precision + Recall);
    
end

% Function for Feature Extraction
function transformedRow = featureExtract(ActivityDataEMG, ActivityDataIMU)

    % Scaling data
    for j=2:9
        ActivityDataEMG(:,j) = ActivityDataEMG(:,j)/(abs(max(ActivityDataEMG(:,j))-min(ActivityDataEMG(:,j))));
    end
    for j=2:11
        ActivityDataIMU(:,j) = ActivityDataIMU(:,j)/(abs(max(ActivityDataIMU(:,j))-min(ActivityDataIMU(:,j))));
    end
    
    %Feature Extraction  
    
    %Transforming Features
    fftEMG1 = fft(ActivityDataEMG(:,2),50);
    fftEMG1(1)=[];
    powerfftEMG1=fftEMG1.*conj(fftEMG1)/50;
    EMG1 = powerfftEMG1(3); %FFT 8 3  

    fftEMG2 = fft(ActivityDataEMG(:,3),50);
    fftEMG2(1)=[];
    powerfftEMG2=fftEMG2.*conj(fftEMG2)/50;    
    EMG2 = powerfftEMG2(6); %FFT 20 6
    
    EMG3 = rms(ActivityDataEMG(:,4));
    EMG4 =rms(ActivityDataEMG(:,5));
    EMG5 = rms(ActivityDataEMG(:,6));
    
    fftEMG6 = fft(ActivityDataEMG(:,7),50);
    fftEMG6(1)=[];
    powerfftEMG6=fftEMG6.*conj(fftEMG6)/50;    
    EMG6 = powerfftEMG6(7);%FFT 24 7
    
    fftEMG7 = fft(ActivityDataEMG(:,8),50);
    fftEMG7(1)=[];
    powerfftEMG7=fftEMG7.*conj(fftEMG7)/50;     
    EMG7 = powerfftEMG7(17);%FFT 64 17
    
    
    fftEMG8 = fft(ActivityDataEMG(:,9),50);
    fftEMG8(1)=[];
    powerfftEMG8=fftEMG8.*conj(fftEMG8)/50;
    EMG8 = powerfftEMG8(13);%FFT 48 13
    
    IMU1 = var(ActivityDataIMU(:,2));
    IMU2 = var(ActivityDataIMU(:,3));
    IMU3 = var(ActivityDataIMU(:,4));
    IMU4 = var(ActivityDataIMU(:,5));
    IMU5 = entropy(ActivityDataIMU(:,6));
    IMU6 = mean(ActivityDataIMU(:,7));
    IMU7 = entropy(ActivityDataIMU(:,8));
    IMU8 = entropy(ActivityDataIMU(:,9));
    IMU9 = rms(ActivityDataIMU(:,10));
    IMU10 = rms(ActivityDataIMU(:,11));
    IMU = [IMU1 IMU2 IMU3 IMU4 IMU5 IMU6 IMU7 IMU8 IMU9 IMU10];
    EMG = [EMG1 EMG2 EMG3 EMG4 EMG5 EMG6 EMG7 EMG8];
    transformedRow = [EMG IMU];
end

