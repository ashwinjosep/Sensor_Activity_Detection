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
featureEatingUser = [];
featureNonEatingUser = [];
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
    
    % Assign Class label in data
    for j = 1:forkDataEMGRows
        for k = 1:forkFoodDurationRows
            if(and((forkDataEMG(j,1)>=forkFoodDuration(k,3)),(forkDataEMG(j,1)<=forkFoodDuration(k,2))))
                forkDataEMG(j,10) = 1;
                break;
            end
        end
    end
    
    % Assign Class label in data
    for j = 1:forkDataIMURows
        for k = 1:forkFoodDurationRows
            if(and((forkDataIMU(j,1)>=forkFoodDuration(k,3)),(forkDataIMU(j,1)<=forkFoodDuration(k,2))))
                forkDataIMU(j,12) = 1;
                break;
            end
        end
    end
    
    %Comparing groundtruth and sensor data for spoons
    spoonDataEMG(:,10) = 0;
    spoonDataIMU(:,12) = 0;
    [spoonFoodDurationRows,spoonFoodDurationColumns] = size(spoonFoodDuration);
    [spoonDataEMGRows, spoonDataEMGColumns] = size(spoonDataEMG);
    [spoonDataIMURows, spoonDataIMUColumns] = size(spoonDataIMU);
    
    j = spoonFoodDurationRows;
    beginTimeStamp = spoonDataEMG(spoonDataEMGRows,1);
    
    % Get all food activity time ranges for spoon data
    while(j >= 1)
        endTimeStamp = beginTimeStamp - spoonFoodDuration(j,1);
        spoonFoodDuration(j,2) = beginTimeStamp;
        spoonFoodDuration(j,3) = endTimeStamp;
        if(j>1)
            beginTimeStamp = endTimeStamp - spoonNoFoodDuration(j-1,1);
        end
        j = j-1;
    end
    
    % Assign Class label in data
    for j = 1:spoonDataEMGRows
        for k = 1:spoonFoodDurationRows
            if(and((spoonDataEMG(j,1)>=spoonFoodDuration(k,3)),(spoonDataEMG(j,1)<=spoonFoodDuration(k,2))))
                spoonDataEMG(j,10) = 1;
                break;
            end
        end
    end
    
    % Assign Class label in data
    for j = 1:spoonDataIMURows
        for k = 1:spoonFoodDurationRows
            if(and((spoonDataIMU(j,1)>=spoonFoodDuration(k,3)),(spoonDataIMU(j,1)<=spoonFoodDuration(k,2))))
                spoonDataIMU(j,12) = 1;
                break;
            end
        end
    end
    
    %Append data to global data matrix
    forkEMG=[forkEMG;forkDataEMG];
    forkIMU=[forkIMU;forkDataIMU];
    spoonEMG=[spoonEMG;spoonDataEMG];
    spoonIMU=[spoonIMU;spoonDataIMU];

    %Scaling Data
    for j=2:9
        forkEMG(:,j) = forkEMG(:,j)/(abs(max(forkEMG(:,j))-min(forkEMG(:,j))));
        spoonEMG(:,j) = spoonEMG(:,j)/(abs(max(spoonEMG(:,j))-min(spoonEMG(:,j))));
    end
    for j=2:11
        forkIMU(:,j) = forkIMU(:,j)/(abs(max(forkIMU(:,j))-min(forkIMU(:,j))));
        spoonIMU(:,j) = spoonIMU(:,j)/(abs(max(spoonIMU(:,j))-min(spoonIMU(:,j))));
    end
    
    %Dividing into Eating and Non Eating Data
    idx=(forkEMG(:,10)==1);
    EatingDataEMG = forkEMG(idx,:);
    NonEatingDataEMG = forkEMG(~idx,:);

    idx=(forkIMU(:,12)==1);
    EatingDataIMU = forkIMU(idx,:);
    NonEatingDataIMU = forkIMU(~idx,:);

    %Feature Extraction  
    %Transforming Features
    fftEatingEMG1 = fft(EatingDataEMG(:,2),50);
    fftEatingEMG1(1)=[];
    powerfftEatingEMG1=fftEatingEMG1.*conj(fftEatingEMG1)/50;
    EMG1Eating = powerfftEatingEMG1(3); %FFT 8 3  
    
    fftNonEatingEMG1 = fft(NonEatingDataEMG(:,2),50);
    fftNonEatingEMG1(1)=[];
    powerfftNonEatingEMG1=fftNonEatingEMG1.*conj(fftNonEatingEMG1)/50;    
    EMG1NonEating = powerfftNonEatingEMG1(3); %FFT 8 3
    
    fftEatingEMG2 = fft(EatingDataEMG(:,3),50);
    fftEatingEMG2(1)=[];
    powerfftEatingEMG2=fftEatingEMG2.*conj(fftEatingEMG2)/50;    
    EMG2Eating = powerfftEatingEMG2(6); %FFT 20 6
    
    fftNonEatingEMG2 = fft(NonEatingDataEMG(:,3),50);
    fftNonEatingEMG2(1)=[];
    powerfftNonEatingEMG2=fftNonEatingEMG2.*conj(fftNonEatingEMG2)/50;        
    EMG2NonEating = powerfftNonEatingEMG2(6); %FFT 20 6
    
    EMG3Eating = rms(EatingDataEMG(:,4));
    EMG3NonEating = rms(NonEatingDataEMG(:,4));
    EMG4Eating =rms(EatingDataEMG(:,5));
    EMG4NonEating = rms(NonEatingDataEMG(:,5));
    EMG5Eating = rms(EatingDataEMG(:,6));
    EMG5NonEating = rms(NonEatingDataEMG(:,6));
    
    fftEatingEMG6 = fft(EatingDataEMG(:,7),50);
    fftEatingEMG6(1)=[];
    powerfftEatingEMG6=fftEatingEMG6.*conj(fftEatingEMG6)/50;    
    EMG6Eating = powerfftEatingEMG6(7);%FFT 24 7
    
    fftNonEatingEMG6 = fft(NonEatingDataEMG(:,7),50);
    fftNonEatingEMG6(1)=[];
    powerfftNonEatingEMG6=fftNonEatingEMG6.*conj(fftNonEatingEMG6)/50;
    EMG6NonEating = powerfftNonEatingEMG6(7);%FFT 24 7
    
    fftEatingEMG7 = fft(EatingDataEMG(:,8),50);
    fftEatingEMG7(1)=[];
    powerfftEatingEMG7=fftEatingEMG7.*conj(fftEatingEMG7)/50;     
    EMG7Eating = powerfftEatingEMG7(17);%FFT 64 17
    
    fftNonEatingEMG7 = fft(NonEatingDataEMG(:,8),50);
    fftNonEatingEMG7(1)=[];
    powerfftNonEatingEMG7=fftNonEatingEMG7.*conj(fftNonEatingEMG7)/50;
    EMG7NonEating = powerfftNonEatingEMG7(17);   %FFT 64 17
    
    fftEatingEMG8 = fft(EatingDataEMG(:,9),50);
    fftEatingEMG8(1)=[];
    powerfftEatingEMG8=fftEatingEMG8.*conj(fftEatingEMG8)/50;
    EMG8Eating = powerfftEatingEMG8(13);%FFT 48 13
    
    fftNonEatingEMG8 = fft(NonEatingDataEMG(:,9),50);
    fftNonEatingEMG8(1)=[];
    powerfftNonEatingEMG8=fftNonEatingEMG8.*conj(fftNonEatingEMG8)/50;
    EMG8NonEating = powerfftNonEatingEMG8(13);%FFT 48 13
    
    IMU1Eating = var(EatingDataIMU(:,2));
    IMU1NonEating = var(NonEatingDataIMU(:,2));
    IMU2Eating = var(EatingDataIMU(:,3));
    IMU2NonEating = var(NonEatingDataIMU(:,3));
    IMU3Eating = var(EatingDataIMU(:,4));
    IMU3NonEating = var(NonEatingDataIMU(:,4));
    IMU4Eating = var(EatingDataIMU(:,5));
    IMU4NonEating = var(NonEatingDataIMU(:,5));
    IMU5Eating = entropy(EatingDataIMU(:,6));
    IMU5NonEating = entropy(NonEatingDataIMU(:,6));
    IMU6Eating = mean(EatingDataIMU(:,7));
    IMU6NonEating = mean(NonEatingDataIMU(:,7));
    IMU7Eating = entropy(EatingDataIMU(:,8));
    IMU7NonEating = entropy(NonEatingDataIMU(:,8));
    IMU8Eating = entropy(EatingDataIMU(:,9));
    IMU8NonEating = entropy(NonEatingDataIMU(:,9));
    IMU9Eating = rms(EatingDataIMU(:,10));
    IMU9NonEating = rms(NonEatingDataIMU(:,10));
    IMU10Eating = rms(EatingDataIMU(:,11));
    IMU10NonEating = rms(NonEatingDataIMU(:,11));
    IMUEating = [IMU1Eating IMU2Eating IMU3Eating IMU4Eating IMU5Eating IMU6Eating IMU7Eating IMU8Eating IMU9Eating IMU10Eating];
    IMUNonEating = [IMU1NonEating IMU2NonEating IMU3NonEating IMU4NonEating IMU5NonEating IMU6NonEating IMU7NonEating IMU8NonEating IMU9NonEating IMU10NonEating];
    EMGEating = [EMG1Eating EMG2Eating EMG3Eating EMG4Eating EMG5Eating EMG6Eating EMG7Eating EMG8Eating];
    EMGNonEating = [EMG1NonEating EMG2NonEating EMG3NonEating EMG4NonEating EMG5NonEating EMG6NonEating EMG7NonEating EMG8NonEating];
    Eating = [EMGEating IMUEating];
    featureEatingUser = [featureEatingUser;Eating];
    NonEating = [EMGNonEating IMUNonEating];
    featureNonEatingUser = [featureNonEatingUser;NonEating];
end

featureUser = featureEatingUser;

% Eigen Vector Analysis
% eat_length = length(featureUser);
% featureUser = [featureUser;featureNonEatingUser];
% C = cov(featureUser);
% [V,D] = eig(C);
% [d,ind] = sort(diag(D));
% Ds = D(ind,ind);
% Vs = V(:,ind);

% PCA
[coeff, score, latent] = pca(featureUser);
% Score contains the new feature matrix

% Figure Plot
% vbls = {'','','','','','','','','','','','','IMU5','','IMU7','IMU8','',''};
% figure()
% plot(score((1:eat_length),1), score((1:eat_length),2),'b+')
% hold on
% plot(score(eat_length+1:length(featureUser),1), score(eat_length+1:length(featureUser),2),'ro')
% hold off
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')
% zlabel('3rd Principal Component')
% biplot(coeff(:,1:3),'Scores',score(:,1:3),'Varlabels',vbls);

%User dependent testing

%User Independent testing

