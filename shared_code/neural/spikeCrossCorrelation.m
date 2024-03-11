% This function computes the cross-correlation between two spike trains,
% given a bin size (in s) and a time range [timeMin, timeMax]
% see Dayan & Abbott p. 28
%
% If convertTrain == 1, take spike timestamps and convert it into spike
% train with 1s and 0s
%
% taquino aug/18

function [crossCorrelation, timeLags, delay] = spikeCrossCorrelation(train1, train2, params)
binSize = params.binSize; spikesRange = params.spikesRange; binStep = params.binStep;
convertTrain = params.convertTrain;
% Getting spikes within time range
rangeTrain1 = train1(train1>=spikesRange(1)&train1<=spikesRange(2));
rangeTrain2 = train2(train2>=spikesRange(1)&train2<=spikesRange(2));
if convertTrain
    rangeTrain1 = getSpikeTrain(rangeTrain1,binSize,spikesRange,binStep);
    rangeTrain2 = getSpikeTrain(rangeTrain2,binSize,spikesRange,binStep);
end

[crossCorrelation,delayTrial,~,lags] = compute_xcorr(rangeTrain1,rangeTrain2,'none');
rangeSize = spikesRange(2) - spikesRange(1);
timeLags = linspace(-rangeSize,rangeSize,length(lags));
delay = timeLags(lags == delayTrial);
% timeLags = lags*binSize;

end

function [X,delayTrial,delayAverage,lags] = compute_xcorr(x1,x2,type)
[X,lags] = xcorr(x1,x2,type);
delayTrial = finddelay(x1,x2);
delayAverage = finddelay(nanmean(x1),nanmean(x2));
end




