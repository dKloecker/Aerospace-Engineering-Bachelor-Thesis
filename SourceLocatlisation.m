clear;
clc;

youngsModulus = 73e9;
rho = 2785;
poisonsRatio = 0.33;
thickness = 0.002;

%distancesInMM = [75, 70, 80, 60, 65, 55, 50, 40, 45, 35, 30, 20, 25, 19, 18, 17, 16, 15, 14, 13, 12, 10, 11, 8, 9, 7, 6, 5, 4, 3, 2, 1, 0];
distancesInMM = [85,80,70,75,65,60,50,55,45,40,38,39,37,36,35,34,33,32,31,30,29,28,26,27,25,24,23,22,21,20,18,19,17,16,15,14,13,11,12,10,8,9,7,6,5,4,3,2,1,0];
[distancesInMM, rightIndexes] = sort(distancesInMM,'ascend');

load("A0Times.mat");
load("A0InputOrig.mat")

load("CombinedInputOrig.mat");
load("CombinedTimes.mat")

load("S0InputOrig.mat");
load("S0Times.mat")

signalAmpsS0 = S0InputOrig;
signalTimesS0 = S0Times; 


signalAmpsS0 = signalAmpsS0(:,rightIndexes);

nrOfSignals = length(distancesInMM);

signalDataSize = size(signalAmpsS0);
normalisedSignalsS0 = zeros(signalDataSize);
FFTfrequenciesS0 = zeros(ceil(signalDataSize(1)/2), signalDataSize(2));
P1AmpsS0 = zeros(ceil(signalDataSize(1)/2), signalDataSize(2));

timeDifferencesS0 = diff(signalTimesS0);
meanSampleTimeS0 = mean(timeDifferencesS0);
FsS0 = 1/meanSampleTimeS0;

indexOfMaxAmpS0 = [];
singalStartTimesS0 = [];
signalStartIndexS0 = [];
maxAmpsS0 = [];
maxNormedAmpsS0 = [];
riseTimeS0 = [];
riseTimeRatioS0 = [];
dominantFrequenciesS0 = [];

%Normalise all signals
for i = 1:nrOfSignals
    
    currentSignalS0 = signalAmpsS0(:,i);
    currentMaxS0 = max(abs(currentSignalS0));
    
    %Find first none zero value in the signal and its corresponding index
    signalStartIndexS0(i) = find(currentSignalS0~=0, 1, 'first');
    singalStartTimesS0(i) = signalTimesS0(signalStartIndexS0(i));
    
    indexOfMaxAmpS0(i) = find(abs(currentSignalS0)== currentMaxS0, 1, 'first');
    maxAmpsS0(i) = currentMaxS0;
        
    maxAmpTimesS0(i) = signalTimesS0(indexOfMaxAmpS0(i));
    
    riseTimeS0(i) = signalTimesS0(indexOfMaxAmpS0(i)) - signalTimesS0(signalStartIndexS0(i));
    
    [normalisedSignalsS0(:,i), normalisedWithReferenceS0(:,i)] = normSignal(currentSignalS0, signalAmpsS0(:, 2));
    
    maxNormedAmpsS0(i) =  normalisedWithReferenceS0(indexOfMaxAmpS0(i),i);
    
    absRiseTimeRatioS0(i) = riseTimeS0(i)/abs(maxNormedAmpsS0(i));
    riseTimeRatioS0(i) = riseTimeS0(i)/maxNormedAmpsS0(i);
    
    [FFTfrequenciesS0(:,i),P1AmpsS0(:,i)] = signalFFT(FsS0, normalisedSignalsS0(:,i));
    
    [FFTfrequenciesReferenceNorm(:,i), P1AmpsReferenceNormS0(:,i)] = signalFFT(FsS0,normalisedWithReferenceS0(:,i));
    
    [powerS0(:,i),freqSpecS0(:,i)] = powerSpectrum(normalisedWithReferenceS0(:,i),FsS0);
    
    indexOfDominantFreqS0 = find(P1AmpsS0(:,i) == max(P1AmpsS0(:,i)), 1, 'first');
    
    dominantFrequenciesS0(i) = freqSpecS0(indexOfDominantFreqS0,i);
    
    meanFrequenciesS0(i) = meanfreq(powerS0(:,i),FsS0);
    
    energyOfSignalS0(i) = trapz(signalTimesS0, abs(normalisedWithReferenceS0(:,i))); 
    
end
    
%CALCULATE ARRIVAL TIMES OF DIFFERENT WAVE MODES
cLongitudinalS0 = sqrt(youngsModulus/rho);
cTransverseS0 = sqrt(youngsModulus/(2*rho*(1+poisonsRatio)));
cExtensionalS0 = sqrt(youngsModulus/(rho*(1-poisonsRatio^2)));
cFlexuralS0 = ((youngsModulus*thickness^2)/(12*rho*(1-poisonsRatio^2)))^0.25 ...
    * sqrt(dominantFrequenciesS0*2*pi);

arrivalTimesLongitudinalS0 = maxAmpTimesS0(1) + (distancesInMM/1000)/cLongitudinalS0;
arrivalTimesTransverseS0 =  maxAmpTimesS0(1) + (distancesInMM/1000)/cTransverseS0;
arrivalTimesExtensionalS0 =  maxAmpTimesS0(1) + (distancesInMM/1000)/cExtensionalS0;
arrivalTimesFlexuralS0 =  maxAmpTimesS0(1) + (distancesInMM/1000)./cFlexuralS0;





%--------------------------------------------------------------------------%
%                       A0 MODE INPUT PROCESSING
%                       
%--------------------------------------------------------------------------%


signalAmpsA0 = A0InputOrig;
signalTimesA0 = A0Times; 


signalAmpsA0 = signalAmpsA0(:,rightIndexes);

nrOfSignals = length(distancesInMM);

signalDataSize = size(signalAmpsA0);
normalisedSignalsA0 = zeros(signalDataSize);
FFTfrequenciesA0 = zeros(ceil(signalDataSize(1)/2), signalDataSize(2));
P1AmpsA0 = zeros(ceil(signalDataSize(1)/2), signalDataSize(2));

timeDifferencesA0 = diff(signalTimesA0);
meanSampleTimeA0 = mean(timeDifferencesA0);
FsA0 = 1/meanSampleTimeA0;

indexOfMaxAmpA0 = [];
singalStartTimesA0 = [];
signalStartIndexA0 = [];
maxAmpsA0 = [];
maxNormedAmpsA0 = [];
riseTimeA0 = [];
riseTimeRatioA0 = [];
dominantFrequenciesA0 = [];

%Normalise all signals
for i = 1:nrOfSignals
    
    currentSignalA0 = signalAmpsA0(:,i);
    currentMaxA0 = max(abs(currentSignalA0));
    
    %Find first none zero value in the signal and its corresponding index
    signalStartIndexA0(i) = find(currentSignalA0 ~=0, 1, 'first');
    singalStartTimesA0(i) = signalTimesA0(signalStartIndexA0(i));
    
    indexOfMaxAmpA0(i) = find(abs(currentSignalA0)== currentMaxA0, 1, 'first');
    maxAmpsA0(i) = currentMaxA0;
        
    maxAmpTimesA0(i) = signalTimesA0(indexOfMaxAmpA0(i));
    
    riseTimeA0(i) = signalTimesA0(indexOfMaxAmpA0(i)) - signalTimesA0(signalStartIndexA0(i));
    
    [normalisedSignalsA0(:,i), normalisedWithReferenceA0(:,i)] = normSignal(currentSignalA0, signalAmpsA0(:, 2));
    
    maxNormedAmpsA0(i) =  normalisedWithReferenceA0(indexOfMaxAmpA0(i),i);
    
    absRiseTimeRatioA0(i) = riseTimeA0(i)/abs(maxNormedAmpsA0(i));
    riseTimeRatioA0(i) = riseTimeA0(i)/maxNormedAmpsA0(i);
    
    [FFTfrequenciesA0(:,i),P1AmpsA0(:,i)] = signalFFT(FsA0, normalisedSignalsA0(:,i));
    
    [FFTfrequenciesReferenceNorm(:,i), P1AmpsReferenceNormA0(:,i)] = signalFFT(FsA0,normalisedWithReferenceA0(:,i));
    
    [powerA0(:,i),freqSpecA0(:,i)] = powerSpectrum(normalisedWithReferenceA0(:,i),FsA0);
    
    indexOfDominantFreqA0 = find(P1AmpsA0(:,i) == max(P1AmpsA0(:,i)), 1, 'first');
    
    dominantFrequenciesA0(i) = freqSpecA0(indexOfDominantFreqA0,i);
    
    meanFrequenciesA0(i) = meanfreq(powerA0(:,i),FsA0);

    energyOfSignalA0(i) = trapz(signalTimesA0, abs(normalisedWithReferenceA0(:,i))); 
end
    
%CALCULATE ARRIVAL TIMES OF DIFFERENT WAVE MODES
cLongitudinalA0 = sqrt(youngsModulus/rho);
cTransverseA0 = sqrt(youngsModulus/(2*rho*(1+poisonsRatio)));
cExtensionalA0 = sqrt(youngsModulus/(rho*(1-poisonsRatio^2)));
cFlexuralA0 = ((youngsModulus*thickness^2)/(12*rho*(1-poisonsRatio^2)))^0.25 ...
    * sqrt(dominantFrequenciesA0*2*pi);

arrivalTimesLongitudinalA0 = maxAmpTimesA0(1) + (distancesInMM/1000)/cLongitudinalA0;
arrivalTimesTransverseA0 =  maxAmpTimesA0(1) + (distancesInMM/1000)/cTransverseA0;
arrivalTimesExtensionalA0 =  maxAmpTimesA0(1) + (distancesInMM/1000)/cExtensionalA0;
arrivalTimesFlexuralA0 =  maxAmpTimesA0(1) + (distancesInMM/1000)./cFlexuralA0;


%--------------------------------------------------------------------------%
%                       Combined MODE INPUT PROCESSING
%                       
%--------------------------------------------------------------------------%


signalAmpsComb = CombinedInputOrig;
signalTimesComb = CombinedTimes; 


signalAmpsComb = signalAmpsComb(:,rightIndexes);

nrOfSignals = length(distancesInMM);

signalDataSize = size(signalAmpsComb);
normalisedSignalsComb = zeros(signalDataSize);
FFTfrequenciesComb = zeros(ceil(signalDataSize(1)/2), signalDataSize(2));
P1AmpsComb = zeros(ceil(signalDataSize(1)/2), signalDataSize(2));

timeDifferencesComb = diff(signalTimesComb);
meanSampleTimeComb = mean(timeDifferencesComb);
FsComb = 1/meanSampleTimeComb;

indexOfMaxAmpComb = [];
singalStartTimesComb = [];
signalStartIndexComb = [];
maxAmpsComb = [];
maxNormedAmpsComb = [];
riseTimeComb = [];
riseTimeRatioComb = [];
dominantFrequenciesComb = [];

%Normalise all signals
for i = 1:nrOfSignals
    
    currentSignalComb = signalAmpsComb(:,i);
    currentMaxComb = max(abs(currentSignalComb));
    
    %Find first none zero value in the signal and its corresponding index
    signalStartIndexComb(i) = find(currentSignalComb ~=0, 1, 'first');
    singalStartTimesComb(i) = signalTimesComb(signalStartIndexComb(i));
    
    indexOfMaxAmpComb(i) = find(abs(currentSignalComb)== currentMaxComb, 1, 'first');
    maxAmpsComb(i) = currentMaxComb;
        
    maxAmpTimesComb(i) = signalTimesComb(indexOfMaxAmpComb(i));
    
    riseTimeComb(i) = signalTimesComb(indexOfMaxAmpComb(i)) - signalTimesComb(signalStartIndexComb(i));
    
    [normalisedSignalsComb(:,i), normalisedWithReferenceComb(:,i)] = normSignal(currentSignalComb, signalAmpsComb(:, 2));
    
    maxNormedAmpsComb(i) =  normalisedWithReferenceComb(indexOfMaxAmpComb(i),i);
    
    absRiseTimeRatioComb(i) = riseTimeComb(i)/abs(maxNormedAmpsComb(i));
    riseTimeRatioComb(i) = riseTimeComb(i)/maxNormedAmpsComb(i);
    
    [FFTfrequenciesComb(:,i),P1AmpsComb(:,i)] = signalFFT(FsComb, normalisedSignalsComb(:,i));
    
    [FFTfrequenciesReferenceNorm(:,i), P1AmpsReferenceNormComb(:,i)] = signalFFT(FsComb,normalisedWithReferenceComb(:,i));
    
    [powerComb(:,i),freqSpecComb(:,i)] = powerSpectrum(normalisedWithReferenceComb(:,i),FsComb);
    
    indexOfDominantFreqComb = find(P1AmpsComb(:,i) == max(P1AmpsComb(:,i)), 1, 'first');
    
    dominantFrequenciesComb(i) = freqSpecComb(indexOfDominantFreqComb,i);
    
    meanFrequenciesComb(i) = meanfreq(powerComb(:,i),FsComb);
    
    energyOfSignalComb(i) = trapz(signalTimesComb, abs(normalisedWithReferenceComb(:,i))); 
end


distancesOfInterest = [ 1, 5, 10, 15, 20, 25, 40, 60, 80];
indexOfClosestDist = [];
for dist = 1:length(distancesOfInterest)
    [ ~, indexOfClosestDist(dist) ] = min(abs(distancesInMM - distancesOfInterest(dist))); 
end
    
distancesUsed = distancesInMM(indexOfClosestDist);

str = {};
for i = 1:length(distancesUsed)
    
    str = [str , strcat('x= ' , num2str(distancesUsed(i)),'mm' )]
    
end


%CALCULATE ARRIVAL TIMES OF DIFFERENT WAVE MODES
cLongitudinalComb = sqrt(youngsModulus/rho);
cTransverseComb = sqrt(youngsModulus/(2*rho*(1+poisonsRatio)));
cExtensionalComb = sqrt(youngsModulus/(rho*(1-poisonsRatio^2)));
cFlexuralComb = ((youngsModulus*thickness^2)/(12*rho*(1-poisonsRatio^2)))^0.25 ...
    * sqrt(dominantFrequenciesComb*2*pi);

arrivalTimesLongitudinalComb = maxAmpTimesComb(1) + (distancesInMM/1000)/cLongitudinalComb;
arrivalTimesTransverseComb =  maxAmpTimesComb(1) + (distancesInMM/1000)/cTransverseComb;
arrivalTimesExtensionalComb =  maxAmpTimesComb(1) + (distancesInMM/1000)/cExtensionalComb;
arrivalTimesFlexuralComb =  maxAmpTimesComb(1) + (distancesInMM/1000)./cFlexuralComb;


sensor_x = zeros(50, 1);

sensor_y = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,...
    21, 22, 23, 24, 25,26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, ...
    45, 50, 55, 60, 65, 70, 75, 80, 85]' * 1e-3;

sensors  = [sensor_x, sensor_y]; 

sensor1 = sensors(10,:);
sensor2 = sensors(20,:);
sensor3 = sensors(40,:);

convergence_criteria = 0.001;

t1 = singalStartTimesS0(10);
t2 = singalStartTimesS0(20);
t3 = singalStartTimesS0(40);

t1 = maxAmpTimesS0(10);
t2 = maxAmpTimesS0(5);
t3 = maxAmpTimesS0(40);

velocity  = cExtensionalS0;
[xs_exten, ys_exten, totalError_exten] = calculateDistance( sensor1, sensor2, sensor3, t1, t2, t3, velocity, convergence_criteria );

velocity  = cLongitudinalS0;
[xs_long, ys_long, totalError_long] = calculateDistance( sensor1, sensor2, sensor3, t1, t2, t3, velocity, convergence_criteria );

velocity  = cTransverseS0;
[xs_trans, ys_trans, totalError_trans] = calculateDistance( sensor1, sensor2, sensor3, t1, t2, t3, velocity, convergence_criteria );

results = [xs_exten, ys_exten, totalError_exten;...
    xs_long, ys_long, totalError_long;...
    xs_trans, ys_trans, totalError_trans]';

