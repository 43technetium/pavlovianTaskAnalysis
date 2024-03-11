function [psth,variance] = getPSTH_corr(spikes,binSize,spikesRange,binStep)
nTrials = length(spikes);
psthMatrix = [];
for tI = 1:nTrials
    psthMatrix = [psthMatrix getSpikeTrain(spikes{tI},binSize,spikesRange,binStep)];
end
psth = mean(psthMatrix,2);
variance = var(psthMatrix,[],2);
end
