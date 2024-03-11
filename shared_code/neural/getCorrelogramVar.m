function V = getCorrelogramVar(spikes1,spikes2,N,params)
binSize = params.binSize; spikesRange = params.spikesRange; binStep = params.binStep;
[allTrains1,V1] = getPSTH_corr(spikes1,binSize,spikesRange,binStep);
[allTrains2,V2] = getPSTH_corr(spikes2,binSize,spikesRange,binStep);
V = (compute_xcorr(V1,V2,'none')+compute_xcorr(allTrains1.^2,V2,'none')+...
    compute_xcorr(V1,allTrains2.^2,'none'))./N;
end