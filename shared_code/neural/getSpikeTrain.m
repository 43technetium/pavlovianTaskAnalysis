% This function gets spike timestamps and converts it into binned spike
% trains
% This funciton is called by spikeCrossCorrelation and getShuffleCorrector
%
% taquino sep/18

function newTrain = getSpikeTrain(train,binSize,spikesRange,binStep)
% Defining time bins for spikes
binStarts = spikesRange(1):binStep:spikesRange(2)-binSize;
binEnds = binStarts+binSize;
newTrain = zeros(length(binStarts),1);
for bI = 1:length(binStarts)
    newTrain(bI) = sum(train>=binStarts(bI)&train<=binEnds(bI));
end
end