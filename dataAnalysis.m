% get list of all users in groundtruth folder
groundTruthUsers = dir('groundTruth/');
indices = ismember({groundTruthUsers.name},{'.', '..'});
groundTruthUsers = groundTruthUsers(~indices);

% get list of all users in myodata folder
myoDataUsers = dir('MyoData/');
indices = ismember({myoDataUsers.name},{'.', '..'});
myoDataUsers = myoDataUsers(~indices);

for i = 1 : 4%length(groundTruthUsers)
    
    forkDataTruth=[];
    spoonDataTruth = [];
    forkDataEMG = [];
    forkDataIMU = [];
    spoonDataEMG = [];
    
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
        if(j>1)spoon
            beginTimeStamp = endTimeStamp - forkNoFoodDuration(j-1,1);
        end
        j = j-1;
    end
    
    % Assign Class label in data
    activityFlag=0;
    for j = 1:forkDataEMGRows
        for k = 1:forkFoodDurationRows
            if(and((forkDataEMG(j,1)>=forkFoodDuration(k,3)),(forkDataEMG(j,1)<=forkFoodDuration(k,2))))
                forkDataEMG(j,10) = 1;
                flag = 1;
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
    
    %Limiting Data to one activity of each type
    prevClass = forkDataEMG(1,10);
    count=0;
    for j=2:forkDataEMGRows
        if prevClass ~= forkDataEMG(j,10)
            count=count+1;
            prevClass = forkDataEMG(j,10);
        end    
        if count==2
           forkDataEMG(j+1:forkDataEMGRows,:)=[]; 
           break;
        end
    end
    
    prevClass = forkDataIMU(1,12);
    count=0;
    for j=2:forkDataIMURows
        if prevClass ~= forkDataIMU(j,12)
            count=count+1;
            prevClass = forkDataIMU(j,12);
        end    
        if count==2
           forkDataIMU(j+1:forkDataIMURows,:)=[]; 
           break;
        end
    end

    %Feature Extraction
    %Scaling Data
    for j=2:9
        forkDataEMG(:,j) = forkDataEMG(:,j)/(abs(max(forkDataEMG(:,j))-min(forkDataEMG(:,j))));
        spoonDataEMG(:,j) = spoonDataEMG(:,j)/(abs(max(spoonDataEMG(:,j))-min(spoonDataEMG(:,j))));
    end
    for j=2:11
        forkDataIMU(:,j) = forkDataIMU(:,j)/(abs(max(forkDataIMU(:,j))-min(forkDataIMU(:,j))));
        spoonDataIMU(:,j) = spoonDataIMU(:,j)/(abs(max(spoonDataIMU(:,j))-min(spoonDataIMU(:,j))));
    end

    %Dividing into Eating and Non Eating Data
    idx=(forkDataEMG(:,10)==1);
    EatingDataEMG = forkDataEMG(idx,:);
    NonEatingDataEMG = forkDataEMG(~idx,:);

    idx=(forkDataIMU(:,12)==1);
    EatingDataIMU = forkDataIMU(idx,:);
    NonEatingDataIMU = forkDataIMU(~idx,:);

    %EMG Sensor Data Transformations
    % Calculating x axis values for plotting
    EatingDataEMGLength = size(EatingDataEMG);
    xValuesEatingEMG = zeros(EatingDataEMGLength(1),1);
    for j=1:EatingDataEMGLength(1)
           xValuesEatingEMG(j) = j/EatingDataEMGLength(1);  
    end

    NonEatingDataEMGLength = size(NonEatingDataEMG);
    xValuesNonEatingEMG = zeros(NonEatingDataEMGLength(1),1);
    for j=1:NonEatingDataEMGLength(1)
           xValuesNonEatingEMG(j) = j/NonEatingDataEMGLength(1);  
    end

    %Calculating RMS values 
    rmsEatingEMG = zeros(EatingDataEMGLength(1),8);
    rmsNonEatingEMG = zeros(NonEatingDataEMGLength(1),8);

    %Calculating Means
    meanEatingEMG = zeros(EatingDataEMGLength(1),8);
    meanNonEatingEMG = zeros(NonEatingDataEMGLength(1),8);

    %Caculating entropy
    entropyEatingEMG = zeros(EatingDataEMGLength(1),8);
    entropyNonEatingEMG = zeros(NonEatingDataEMGLength(1),8);

    %Calculating Variance
    varEatingEMG = zeros(EatingDataEMGLength(1),8);
    varNonEatingEMG = zeros(NonEatingDataEMGLength(1),8);

    for sensor = 2:9
        %Mean
        meanValue = mean(EatingDataEMG(:,sensor));
        meanEatingEMG(:,sensor) = meanValue;

        meanValue = mean(NonEatingDataEMG(:,sensor));
        meanNonEatingEMG(:,sensor) = meanValue;

        %RMS
        rmsValue = rms(EatingDataEMG(:,sensor));
        rmsEatingEMG(:,sensor) = rmsValue;

        rmsValue = rms(NonEatingDataEMG(:,sensor));
        rmsNonEatingEMG(:,sensor) = rmsValue;

        %Entropy
        entropyValue = entropy(EatingDataEMG(:,sensor));
        entropyEatingEMG(:,sensor) = entropyValue;

        entropyValue = entropy(NonEatingDataEMG(:,sensor));
        entropyNonEatingEMG(:,sensor) = entropyValue;

        %Variance
        varValue = var(EatingDataEMG(:,sensor));
        varEatingEMG(:,sensor) = varValue;

        varValue = var(NonEatingDataEMG(:,sensor));
        varNonEatingEMG(:,sensor) = varValue;
    end

    for sensor = 2:9
%         plotFigure(xValuesEatingEMG,xValuesNonEatingEMG,EatingDataEMG(:,sensor),NonEatingDataEMG(:,sensor),meanEatingEMG(:,sensor),meanNonEatingEMG(:,sensor),"mean","EMG",sensor-1,groundTruthUsers(i).name);
%         plotFigure(xValuesEatingEMG,xValuesNonEatingEMG,EatingDataEMG(:,sensor),NonEatingDataEMG(:,sensor),rmsEatingEMG(:,sensor),rmsNonEatingEMG(:,sensor),"RMS","EMG",sensor-1,groundTruthUsers(i).name);
%         plotFigure(xValuesEatingEMG,xValuesNonEatingEMG,EatingDataEMG(:,sensor),NonEatingDataEMG(:,sensor),entropyEatingEMG(:,sensor),entropyNonEatingEMG(:,sensor),"Entropy","EMG",sensor-1,groundTruthUsers(i).name);
%         plotFigure(xValuesEatingEMG,xValuesNonEatingEMG,EatingDataEMG(:,sensor),NonEatingDataEMG(:,sensor),varEatingEMG(:,sensor),varNonEatingEMG(:,sensor),"Variance","EMG",sensor-1,groundTruthUsers(i).name);
    end

    %IMU Sensor Data Transformations
    % Calculating x axis values for plotting
    EatingDataIMULength = size(EatingDataIMU);
    xValuesEatingIMU = zeros(EatingDataIMULength(1),1);
    for j=1:EatingDataIMULength(1)
           xValuesEatingIMU(j) = j/EatingDataIMULength(1);  
    end

    NonEatingDataIMULength = size(NonEatingDataIMU);
    xValuesNonEatingIMU = zeros(NonEatingDataIMULength(1),1);
    for j=1:NonEatingDataIMULength(1)
           xValuesNonEatingIMU(j) = j/NonEatingDataIMULength(1);  
    end

    %Calculating RMS values 
    rmsEatingIMU = zeros(EatingDataIMULength(1),8);
    rmsNonEatingIMU = zeros(NonEatingDataIMULength(1),8);

    %Calculating Means
    meanEatingIMU = zeros(EatingDataIMULength(1),8);
    meanNonEatingIMU = zeros(NonEatingDataIMULength(1),8);

    %Caculating entropy
    entropyEatingIMU = zeros(EatingDataIMULength(1),8);
    entropyNonEatingIMU = zeros(NonEatingDataIMULength(1),8);

    %Calculating Variance
    varEatingIMU = zeros(EatingDataIMULength(1),8);
    varNonEatingIMU = zeros(NonEatingDataIMULength(1),8);

    for sensor = 2:11
        %Mean
        meanValue = mean(EatingDataIMU(:,sensor));
        meanEatingIMU(:,sensor) = meanValue;

        meanValue = mean(NonEatingDataIMU(:,sensor));
        meanNonEatingIMU(:,sensor) = meanValue;

        %RMS
        rmsValue = rms(EatingDataIMU(:,sensor));
        rmsEatingIMU(:,sensor) = rmsValue;

        rmsValue = rms(NonEatingDataIMU(:,sensor));
        rmsNonEatingIMU(:,sensor) = rmsValue;

        %Entropy
        entropyValue = entropy(EatingDataIMU(:,sensor));
        entropyEatingIMU(:,sensor) = entropyValue;

        entropyValue = entropy(NonEatingDataIMU(:,sensor));
        entropyNonEatingIMU(:,sensor) = entropyValue;

        %Variance
        varValue = var(EatingDataIMU(:,sensor));
        varEatingIMU(:,sensor) = varValue;

        varValue = var(NonEatingDataIMU(:,sensor));
        varNonEatingIMU(:,sensor) = varValue;
    end

    for sensor = 2:11
%         plotFigure(xValuesEatingIMU,xValuesNonEatingIMU,EatingDataIMU(:,sensor),NonEatingDataIMU(:,sensor),meanEatingIMU(:,sensor),meanNonEatingIMU(:,sensor),"mean","IMU",sensor-1,groundTruthUsers(i).name);
%         plotFigure(xValuesEatingIMU,xValuesNonEatingIMU,EatingDataIMU(:,sensor),NonEatingDataIMU(:,sensor),rmsEatingIMU(:,sensor),rmsNonEatingIMU(:,sensor),"RMS","IMU",sensor-1,groundTruthUsers(i).name);
%         plotFigure(xValuesEatingIMU,xValuesNonEatingIMU,EatingDataIMU(:,sensor),NonEatingDataIMU(:,sensor),entropyEatingIMU(:,sensor),entropyNonEatingIMU(:,sensor),"Entropy","IMU",sensor-1,groundTruthUsers(i).name);
%         plotFigure(xValuesEatingIMU,xValuesNonEatingIMU,EatingDataIMU(:,sensor),NonEatingDataIMU(:,sensor),varEatingIMU(:,sensor),varNonEatingIMU(:,sensor),"Variance","IMU",sensor-1,groundTruthUsers(i).name);
    end


    %Calculating FTT for EMG Sensors
    for sensor=2:9
        %Calculating Fourier Transforms
        fftEatingEMG = fft(EatingDataEMG(:,sensor),50);
        fftEatingEMG(1)=[];
        powerfftEatingEMG=fftEatingEMG.*conj(fftEatingEMG)/50;

        fftNonEatingEMG = fft(NonEatingDataEMG(:,sensor),50);
        fftNonEatingEMG(1)=[];
        powerfftNonEatingEMG=fftNonEatingEMG.*conj(fftNonEatingEMG)/50;

        frequency=1000/250*(0:2001);

        figure('visible','off');
        plot(frequency(1:24),powerfftEatingEMG(1:24),'b',frequency(1:24),powerfftNonEatingEMG(1:24),'r');
        legend("Eating", "Non-Eating")
        sensorNumber = sensor-1;
        chartTitle = "Power Spectral Density - FFT - EMG "+sensorNumber+" :  "+groundTruthUsers(i).name;
        title(chartTitle);
        xlabel('Frequency (Hz)');
        ylabel('Power Spectral Density');
        filename=chartTitle+".png";
        saveas(gcf,filename);
    end
    
    %Calculating FTT for IMU Sensors    
    for sensor=2:11
        %Calculating Fourier Transforms
        fftEatingIMU = fft(EatingDataIMU(:,sensor),50);
        fftEatingIMU(1)=[];
        powerfftEatingIMU=fftEatingIMU.*conj(fftEatingIMU)/50;

        fftNonEatingIMU = fft(NonEatingDataIMU(:,sensor),50);
        fftNonEatingIMU(1)=[];
        powerfftNonEatingIMU=fftNonEatingIMU.*conj(fftNonEatingIMU)/50;

        frequency=1000/250*(0:2000);

        figure('visible','off');
        plot(frequency(1:24),powerfftEatingIMU(1:24),'b',frequency(1:24),powerfftNonEatingIMU(1:24),'r');
        legend("Eating", "Non-Eating")
        sensorNumber = sensor-1;
        chartTitle = "Power Spectral Density - FFT - IMU "+sensorNumber+" :  "+groundTruthUsers(i).name;
        title(chartTitle);
        xlabel('Frequency (Hz)');
        ylabel('Power Spectral Density')
        filename=chartTitle+".png";
        saveas(gcf,filename);
    end
    
end
    
function plotFigure(x1,x2,a1,a2,b1,b2,property,sensorType,sensorNo,user)
    legendb1 = "Eating "+property;
    legendb2 = "Non-Eating "+property;
    chartTitle = sensorType+sensorNo+" Sensor "+property+" Values : "+user;
    filename= chartTitle+".png";
    
    figure('visible','off');
    plot(x2,a2,'b',x1,a1,'r',x1,b1,'k',x2,b2,'g');
    legend('Non-Eating','Eating', legendb1, legendb2);
    title(chartTitle);
    saveas(gcf,filename);    
end


