%Code used for signal plots:
%Aslak Grinsted (2021). Subaxis - Subplot (https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot),
%MATLAB Central File Exchange. Retrieved May 6, 2021

clear;
clc;



youngsModulus = 73e9;
rho = 2785;
poisonsRatio = 0.33;
thickness = 0.002;

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


distancesOfInterest = [ 35];
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




calibrationForce = 1;

for i = 2:50
    
    summedNormedStrainsS0(i) =  sum(normalisedWithReferenceS0(:,i));
    
end

TransferFunctions = summedNormedStrainsS0/calibrationForce;



for i = 2:50
    
    multiplier = TransferFunctions(i);
    currentSignal = normalisedWithReferenceS0(:,i);
    d = multiplier * currentSignal;
    reconstructed(:, i) = d;
    
end





figure(1)
plot(distancesInMM, riseTimeRatioComb, 'mo--')
hold on
plot(distancesInMM, riseTimeRatioA0, 'r+-.');
plot(distancesInMM, riseTimeRatioS0, 'kx-');
grid on 
xlabel("Distance from source (mm)")
ylabel("RA (Max Amp/Rise Time")
legend("Complete Input", "A0 Mode Input", "S0 Mode Input")

% 
figure(2)
plot(distancesInMM, absRiseTimeRatioComb, 'mo--')
hold on 
plot(distancesInMM, absRiseTimeRatioA0, 'r+-.');
plot(distancesInMM, absRiseTimeRatioS0, 'kx-');
grid on 
xlabel("Distance from source (mm)")
ylabel("RA Max Amp/Rise Time")
legend("Complete Input", "A0 Mode Input", "S0 Mode Input")



figure(3)
plot(distancesInMM, riseTimeComb, 'mo--')
hold on
plot(distancesInMM, riseTimeA0, 'r+-.');
plot(distancesInMM, riseTimeS0, 'kx-');
grid on 
xlabel("Distance from source (mm)")
ylabel("Rise Time (s)")
legend("Combined Input", "A0 Mode Input", "S0 Mode Input")


figure(4)
plot(distancesInMM(2:end), abs(maxNormedAmpsComb(2:end)), 'mo--')
hold on
plot(distancesInMM(2:end), abs(maxNormedAmpsA0(2:end)), 'r+-.');
plot(distancesInMM(2:end), abs(maxNormedAmpsS0(2:end)), 'kx-');
grid on 
xlabel("Distance from source (mm)")
ylabel("Maximum normalised displacement")
legend("Complete Input", "A0 Mode Input", "S0 Mode Input")

figure(5)
plot(distancesInMM,dominantFrequenciesComb, 'mo--')
hold on 
plot(distancesInMM,dominantFrequenciesA0, 'r+-.');
plot(distancesInMM,dominantFrequenciesS0, 'kx-')
xlabel("Distance from source (mm)")
ylabel("Dominant Frequency (Hz)")
grid on
legend("Complete Input", "A0 Mode Input", "S0 Mode Input")


figure(6)
plot(distancesInMM,meanFrequenciesComb, 'mo--')
hold on 
plot(distancesInMM,meanFrequenciesA0, 'r+-.');
plot(distancesInMM,meanFrequenciesS0, 'kx-')
grid on
xlabel("Distance from source (mm)")
ylabel("Mean Frequency (Hz)")
legend("Complete Input", "A0 Mode Input", "S0 Mode Input")

figure(7)
subplot(2,2,1)
plot(signalTimesS0, normalisedWithReferenceS0(:,6))
xline(arrivalTimesLongitudinalS0(6), 'r--')
xline(arrivalTimesTransverseS0(6), 'b--') 
xline(arrivalTimesExtensionalS0(6), 'k--')
xline(arrivalTimesFlexuralS0(6), 'g--')

subplot(2,2,2)
plot(signalTimesS0, normalisedWithReferenceS0(:,11))
xline(arrivalTimesLongitudinalS0(11), 'r--')
xline(arrivalTimesTransverseS0(11), 'b--')
xline(arrivalTimesExtensionalS0(11), 'k--')
xline(arrivalTimesFlexuralS0(11), 'g--')

subplot(2,2,3)
plot(signalTimesS0, normalisedWithReferenceS0(:,21))
xline(arrivalTimesLongitudinalS0(21), 'r--')
xline(arrivalTimesTransverseS0(21), 'b--')
xline(arrivalTimesExtensionalS0(21), 'k--')
xline(arrivalTimesFlexuralS0(21), 'g--')

subplot(2,2,4)
plot(signalTimesS0, normalisedWithReferenceS0(:,30))
xline(arrivalTimesLongitudinalS0(30), 'r--')
xline(arrivalTimesTransverseS0(30), 'b--')
xline(arrivalTimesExtensionalS0(30), 'k--')
xline(arrivalTimesFlexuralS0(30), 'g--')

figure(8)
subplot(2,2,1)
plot(signalTimesS0, normalisedWithReferenceA0(:,1))
xline(arrivalTimesLongitudinalA0(1), 'r--')
xline(arrivalTimesTransverseA0(1), 'b--') 
xline(arrivalTimesExtensionalA0(1), 'k--')
xline(arrivalTimesFlexuralA0(1), 'g--')

subplot(2,2,2)
plot(signalTimesS0, normalisedWithReferenceA0(:,2))
xline(arrivalTimesLongitudinalA0(2), 'r--')
xline(arrivalTimesTransverseA0(2), 'b--')
xline(arrivalTimesExtensionalA0(2), 'k--')
xline(arrivalTimesFlexuralA0(2), 'g--')

subplot(2,2,3)
plot(signalTimesS0, normalisedWithReferenceA0(:,3))
xline(arrivalTimesLongitudinalA0(3), 'r--')
xline(arrivalTimesTransverseA0(3), 'b--')
xline(arrivalTimesExtensionalA0(3), 'k--')
xline(arrivalTimesFlexuralA0(3), 'g--')

subplot(2,2,4)
plot(signalTimesS0, normalisedWithReferenceA0(:,5))
xline(arrivalTimesLongitudinalA0(5), 'r--')
xline(arrivalTimesTransverseA0(5), 'b--')
xline(arrivalTimesExtensionalA0(5), 'k--')
xline(arrivalTimesFlexuralA0(5), 'g--')


figure(9)
plot(signalTimesA0, normalisedWithReferenceA0(:,indexOfClosestDist))
grid on
xlabel("Time (s)")
ylabel("Normalised Amplitude")
legend(str{:})

figure(10)
plot(signalTimesS0, normalisedWithReferenceS0(:,indexOfClosestDist))
grid on
xlabel("Time (s)")
ylabel("Normalised Amplitude")
legend(str{:})

figure(11)
plot(signalTimesComb, normalisedWithReferenceComb(:,indexOfClosestDist))
grid on
xlabel("Time (s)")
ylabel("Normalised Amplitude")
legend(str{:})

energyOfSignalComb = energyOfSignalComb/energyOfSignalComb(2);
energyOfSignalS0 = energyOfSignalS0 /energyOfSignalS0(2);
energyOfSignalA0 = energyOfSignalA0 /energyOfSignalA0(2);

figure(12)
plot(distancesInMM(1:end), energyOfSignalA0(1:end), 'mo--')
hold on
plot(distancesInMM(1:end), energyOfSignalComb(1:end), 'r+-.');
plot(distancesInMM(1:end), energyOfSignalS0(1:end), 'kx-')
xlabel("Distance from Source (mm)")
ylabel("Normalised Energy") 
grid on
legend("Complete Input", "A0 Mode Input", "S0 Mode Input")


maxNormedAmpsComb = abs(maxNormedAmpsComb);
maxNormedAmpsS0 = abs(maxNormedAmpsS0);
maxNormedAmpsA0 = abs(maxNormedAmpsA0);

[xDataComb , yDataComb ] = prepareCurveData( distancesInMM, maxNormedAmpsComb );
% Set up fittype and options.
ftComb  = fittype( 'exp2' );
optsComb  = fitoptions( 'Method', 'NonlinearLeastSquares' );
optsComb .Display = 'Off';
opts.StartPoint = [2.47222495251658 -0.152084633769866 0.0564221595072211 -0.0144753927412974];


% Fit model to data.
[fitresultComb , gofComb ] = fit( xDataComb , yDataComb , ftComb , optsComb  );


[xDataS0 , yDataS0 ] = prepareCurveData( distancesInMM, maxNormedAmpsS0 );
% Set up fittype and options.
ftS0  = fittype( 'exp2' );
optsS0  = fitoptions( 'Method', 'NonlinearLeastSquares' );
optsS0 .Display = 'Off';
opts.StartPoint = [2.47222495251658 -0.152084633769866 0.0564221595072211 -0.0144753927412974];


% Fit model to data.
[fitresultS0 , gofS0 ] = fit( xDataS0 , yDataS0 , ftS0 , optsS0  );


[xDataA0 , yDataA0 ] = prepareCurveData( distancesInMM, maxNormedAmpsA0 );
% Set up fittype and options.
ftA0  = fittype( 'exp2' );
optsA0  = fitoptions( 'Method', 'NonlinearLeastSquares' );
optsA0 .Display = 'Off';
opts.StartPoint = [2.47222495251658 -0.152084633769866 0.0564221595072211 -0.0144753927412974];

% Fit model to data.
[fitresultA0 , gofA0 ] = fit( xDataA0 , yDataA0 , ftA0 , optsA0  );

figure(13)
subplot(3,1,1)
plot( fitresultComb, xDataComb, yDataComb );
grid on
xlabel("Distance from source (mm)")
ylabel("Normalised Displacement")
title("Combined Input")

subplot(3,1,2)
plot( fitresultS0, xDataS0, yDataS0 );
grid on
xlabel("Distance from source (mm)")
ylabel("Normalised Displacement")
title("S0 Mode Input")

waveLengthsTransLongS0 =  cFlexuralS0./dominantFrequenciesS0;
%Wave prediction 

[fitresultS0 , gofS0 ] = fit( xDataS0 , yDataS0 , ftS0 , optsS0  );

decayAmountExpected = fitresultS0(distancesInMM)

sensor1Signal = abs(normalisedWithReferenceS0(:,2));

expectedSignals = sensor1Signal'.*decayAmountExpected;
expectedSignals = expectedSignals';


timeIncreases = (distancesInMM/1000)/cExtensionalS0

signalSize = size(expectedSignals);
expectedTimes = zeros(signalSize(1), signalSize(2))';

subplot(3,1,3)
plot( fitresultA0, xDataA0, yDataA0 );
grid on
xlabel("Distance from source (mm)")
ylabel("Normalised Displacement")
title("A0 Mode Input")

figure(14)
plot(FFTfrequenciesA0(1:50,1), P1AmpsReferenceNormA0(1:50,indexOfClosestDist))
grid on
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend(str{:})


figure(15)
plot(FFTfrequenciesS0(1:50,1), P1AmpsReferenceNormS0(1:50,indexOfClosestDist))
grid on
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend(str{:})


figure(16)
plot(FFTfrequenciesComb(1:50,1), P1AmpsReferenceNormComb(1:50,indexOfClosestDist))
grid on
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend(str{:})




distancesOfInterest = [0, 1, 2, 3, 5, 8, 10, 15, 20, 25, 30, 40, 60, 80];
indexOfClosestDist = [];
for dist = 1:length(distancesOfInterest)
    [ ~, indexOfClosestDist(dist) ] = min(abs(distancesInMM - distancesOfInterest(dist))); 
end

str = {};
for i = 1:length(distancesUsed)
    
    str = [str , strcat('x= ' , num2str(distancesUsed(i)),'mm' )]
    
end

distancesUsed = distancesInMM(indexOfClosestDist);

figure(17)

for sig = 1:length(indexOfClosestDist)
    str =  strcat('x= ' , num2str(distancesUsed(sig)),'mm' );
        
    subaxis(ceil(length(indexOfClosestDist)/2), 2,sig, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
    plot(signalTimesComb, normalisedSignalsComb(:,indexOfClosestDist(sig)))
    hold all
    if indexOfClosestDist(sig) ~= 1
        hold on 
        plot(signalTimesComb, normalisedWithReferenceComb(:,indexOfClosestDist(sig)))
    end
    xline(arrivalTimesLongitudinalComb(indexOfClosestDist(sig)), 'r--')
    xline(arrivalTimesTransverseComb(indexOfClosestDist(sig)), 'b--')
    xline(arrivalTimesExtensionalComb(indexOfClosestDist(sig)), 'k--')
    xline(arrivalTimesFlexuralComb(indexOfClosestDist(sig)), 'm--')
    axis tight
    axis off
    ylim([-1, 1])
    xlim([0, signalTimesComb(end)])
end
    
subaxis(8,2,1, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
plot(signalTimesS0, normalisedWithReferenceS0(:,1))
axis tight
axis off

subaxis(8,2,2, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
plot(signalTimesS0, normalisedWithReferenceS0(:,2))
axis tight
axis off

subaxis(8,2,3, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
plot(signalTimesS0, normalisedWithReferenceS0(:,10))
axis tight
axis off

subaxis(8,2,4, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
plot(signalTimesS0, normalisedWithReferenceS0(:,15))
axis tight
axis off

subaxis(8,2,5, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
plot(signalTimesS0, normalisedWithReferenceS0(:,17))
axis tight
axis off

subaxis(8,2,6, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
plot(signalTimesS0, normalisedWithReferenceS0(:,18))
axis tight
axis off

subaxis(8,2,7, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
plot(signalTimesS0, normalisedWithReferenceS0(:,20))
axis tight
axis off

subaxis(8,2,8, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
plot(signalTimesS0, normalisedWithReferenceS0(:,22))
axis tight
axis off

subaxis(8,2,9, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
plot(signalTimesS0, normalisedWithReferenceS0(:,24))
axis tight
axis off

subaxis(8,2,10, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
plot(signalTimesS0, normalisedWithReferenceS0(:,26))
axis tight
axis off

subaxis(8,2,11, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
plot(signalTimesS0, normalisedWithReferenceS0(:,28))
axis tight
axis off

subaxis(8,2,12, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
plot(signalTimesS0, normalisedWithReferenceS0(:,29))
axis tight
axis off

subaxis(8,2,13, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
plot(signalTimesS0, normalisedWithReferenceS0(:,30))
axis tight
axis off

subaxis(8,2,14, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
plot(signalTimesS0, normalisedWithReferenceS0(:,31))
axis tight
axis off

subaxis(8,2,15, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
plot(signalTimesS0, normalisedWithReferenceS0(:,32))
axis tight
axis off

subaxis(8,2,16, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
plot(signalTimesS0, normalisedWithReferenceS0(:,33))
axis tight
axis off


%STFT for the various modes

distancesOfInterest = [1, 5, 10, 15, 20, 30, 40, 60, 80];
indexOfClosestDist = [];
for dist = 1:length(distancesOfInterest)
    [ ~, indexOfClosestDist(dist) ] = min(abs(distancesInMM - distancesOfInterest(dist))); 
end

str = {};
for i = 1:length(distancesUsed)
    
    str = [str , strcat('x = ' , num2str(distancesUsed(i)),'mm' )]
    
end

distancesUsed = distancesInMM(indexOfClosestDist);


figure(18)

for i = 1:length(distancesUsed)
    
    subplot(3,3,i)
    stft(normalisedWithReferenceS0(:,indexOfClosestDist(i)),FsS0,'Window',kaiser(256,5), 'OverlapLength',220,'FFTLength',512);
    ylim([0 2])
    caxis([-40 40])
    title(str{i})
    
end

%A0 mode results
for i = 1:length(distancesUsed)
    
    subplot(3,3,i)
    stft(normalisedWithReferenceA0(:,indexOfClosestDist(i)),FsS0,'Window',kaiser(256,5), 'OverlapLength',220,'FFTLength',512);
    ylim([0 2])
    caxis([-40 40])
    title(str{i})
    
end

%Combined STFT

for i = 1:length(distancesUsed)
    
    subplot(3,3,i)
    stft(normalisedWithReferenceComb(:,indexOfClosestDist(i)),FsS0,'Window',kaiser(256,5), 'OverlapLength',220,'FFTLength',512);
    ylim([0 2])
    caxis([-40 40])
    title(str{i})
    
end

figure(19)
subaxis(1,3,1, 'Spacing', 0.03, 'Padding', 0.01, 'Margin', 0.01);
plot(signalTimesComb,normalisedSignalsComb(:,1))
axis tight
ylim([-1 1])
% axis off
subaxis(1,3,2, 'Spacing', 0.03, 'Padding', 0.01, 'Margin', 0.01);
plot(signalTimesS0,normalisedSignalsS0(:,1))
axis tight
ylim([-1 1])
% axis off
subaxis(1,3,3, 'Spacing', 0.03, 'Padding', 0.01, 'Margin', 0.01);
plot(signalTimesA0,normalisedSignalsA0(:,1))
axis tight
ylim([-1 1])
%axis off
