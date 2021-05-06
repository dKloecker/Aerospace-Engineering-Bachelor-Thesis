function [normalisedSignal,normalisedWithReference] = normSignal(signal,referenceSignal)

maxValue = max(abs(signal));
normalisedSignal = signal/maxValue;

if isnan(referenceSignal) 
    maxReferenceValue = maxValue;
else 
    maxReferenceValue = max(abs(referenceSignal));
end

normalisedWithReference = signal/maxReferenceValue;

end

