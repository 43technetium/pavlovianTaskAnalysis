% This function computes the shuffle corrector between two spike trains,
% given a bin size (in s) and a time range [timeMin, timeMax]
% see Brody, 1999
%
% taquino aug/18

function [crossCorrelation] = getShuffleCorrector(spikes1, spikes2, params)
binSize = params.binSize; spikesRange = params.spikesRange; binStep = params.binStep;
[allTrains1,~] = getPSTH_corr(spikes1,binSize,spikesRange,binStep);
[allTrains2,~] = getPSTH_corr(spikes2,binSize,spikesRange,binStep);
[crossCorrelation,~,~,~] = compute_xcorr(allTrains1,allTrains2,'none');
end


