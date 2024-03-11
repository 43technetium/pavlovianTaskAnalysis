function pavlovCrossCorrelation()
dbstop if error;
warning('off','all')
% All sessions
params.sessions = {'P49CS','P51CS_1','P51CS_2','P53CS','P54CS','P58CS',...
    'P60CS','P61CS','P62CS','P63CS','P70CS','P71CS','P76CS'};

blockTrials = [1:24;25:48;49:72;73:96];

nSessions = length(params.sessions);
sessionPath = '..\patientData\';
sessionCell = cell(nSessions, 1);

% Loading sessionParams files and creating cell of all units
allUnitCell = {};
sessionNames = {};
latents = {};
sessionIds = {};
for sI = 1:nSessions
    %% Loading data
    session = params.sessions{sI};
    behavior_folder = [sessionPath '\allBehavior\'];
    behavior_data = load([behavior_folder 'sessionBehavior_' session '.mat']);
    sessionFolder = [sessionPath session '\'];
    load([sessionFolder 'sessionData.mat']) 
    unitCell = sessionData.neuralData.unitCell;  
    nTrials = length(unitCell{1,1}.trialReferencedSpikes);
    nUnits = length(unitCell);   
    sessionCell{sI} = sessionData.neuralData;
    nameCell = cell(length(sessionData.neuralData.unitCell), 1);
    nameCell(:) = {session};
    sIdxCell = cell(length(sessionData.neuralData.unitCell), 1);
    sIdxCell(:) = {sI};    
    sessionNames = vertcat(sessionNames, nameCell);
    allUnitCell = vertcat(allUnitCell, sessionData.neuralData.unitCell); 
    sessionIds = vertcat(sessionIds, sIdxCell);
    % Getting behavior data
    latents = [latents behavior_data];
end
allUnitCell = horzcat(allUnitCell, sessionNames,sessionIds);

%% Getting brain areas
nUnits = length(allUnitCell);
areaVec = zeros(nUnits,1);
for uI = 1:nUnits
    brainArea = allUnitCell{uI,1}.unitInfo(4);
    areaVec(uI) = brainArea;
end
% Merging left/right side for all brain areas
mergedAreas = areaVec;
mergedAreas(mod(mergedAreas,2)==0)=mergedAreas(mod(mergedAreas,2)==0)-1;
% Hippocampus/Amygdala/Anterior cingulate/Supplementary motor areas
hipUnits = find(mergedAreas==1); amyUnits = find(mergedAreas==3);
accUnits = find(mergedAreas==5); smaUnits = find(mergedAreas==7);
ofcUnits = find(mergedAreas==9);
unitsByArea = {ofcUnits; amyUnits; hipUnits; accUnits; smaUnits};

% Selecting all valid area 1 units
corrUnits1 = allUnitCell(accUnits,:);
% Selecting all valid area 2 units
corrUnits2 = allUnitCell(smaUnits,:);

%% ========================================================================
% Defining cross-correlation parameters
params.binSize = 0.05;
% params.binSize = 0.015;
params.binStep = 0.005;
params.binRange = [-2,2];
% params.spikesRange = [4 7];
params.spikesRange = [7 10];
params.convertTrain = 1;
params.nTimes = 2*length(params.spikesRange(1):params.binSize:params.spikesRange(2)-params.binSize)-1;
params.nCorrTimes = 2*length(params.spikesRange(1):params.binStep:params.spikesRange(2)-params.binSize)-1;

%% ========================================================================
% Getting unit pairs between selected unit groups, within sessions
pairs = [];
learningSessions = params.sessions;

for sI = 1:length(learningSessions)
    ids1 = find(cell2mat(cellfun(@(x) strcmp(x,learningSessions{sI}),corrUnits1(:,2),'UniformOutput',false)));
    ids2 = find(cell2mat(cellfun(@(x) strcmp(x,learningSessions{sI}),corrUnits2(:,2),'UniformOutput',false)));
    [X,Y] = meshgrid(ids1,ids2);
    pairs = [pairs; [X(:) Y(:)]];
end

%% ========================================================================
% Performing cross-correlation in each pair
nBins = 3; 
unitAvgXCorr = zeros(params.nCorrTimes,nBins,length(pairs(:,1)));
V = zeros(params.nCorrTimes,nBins,length(pairs(:,1)));
shuffleCorrectorMat = zeros(params.nCorrTimes,nBins,length(pairs(:,1)));
pairDelay = zeros(nBins,length(pairs(:,1)));
nBinTrials = zeros(nBins,length(pairs(:,1)));
for bI = 1:nBins
    for pI = 1:length(pairs(:,1))
        display(['Pair ' num2str(pI) ' of ' num2str(length(pairs(:,1))) ...
            '/ Regressor bin ' num2str(bI) ' of ' num2str(nBins)])
        unitSpikes1 = corrUnits1{pairs(pI,1),1}.trialReferencedSpikes;
        unitSpikes2 = corrUnits2{pairs(pI,2),1}.trialReferencedSpikes;
        % Selecting trials within the desired bin 
%         binTrials = corrUnits1{pairs(pI,1),1}.groupedLabels{trialType}.binnedAbsPE==bI;
        sIdx = corrUnits1{pairs(pI,1),3};
        
        %% Getting regressors
        behavior_data = latents{sIdx};
        CSp_options = unique(behavior_data.stim2Vec);
        T_CSp = zeros(length(behavior_data.T_cell),2);
        for tI = 1:length(behavior_data.T_cell)
            CSp = behavior_data.stim2Vec(tI);        
            idx = find(CSp == CSp_options);
            T = behavior_data.T_cell{tI};
            [block, ~, ~] = find(blockTrials == tI);
            % Getting valence of CSd
            if behavior_data.trialTypes_flat(tI)==1||behavior_data.trialTypes_flat(tI)==3
                CSd_idx = 1;
            else
                CSd_idx = 2;
            end

            if CSp == behavior_data.CSpPlus(block)
                T_CSp(tI,idx) = T(CSd_idx,3);
                T_CSp(tI,3-idx) = T(CSd_idx,4);
            else
                T_CSp(tI,idx) = T(CSd_idx,4);
                T_CSp(tI,3-idx) = T(CSd_idx,3);
            end
        end
        
        x = behavior_data.SPE2_vec; regressorName = 'uPE 2';
%         x = behavior_data.SPE1_vec; regressorName = 'uPE 1'; 
        quantiles = getQuantileLabels(x,3);        
        binTrials = quantiles == bI;
        % Getting averages over trials
        xcorrMat = [];
        delayMat = [];
        trialIds = find(binTrials);
        for i = 1:length(trialIds)
            uIdx = trialIds(i);
            train1 = unitSpikes1{uIdx};
            train2 = unitSpikes2{uIdx};
            [crossCorrelation, timeLags, delay] = spikeCrossCorrelation(train1, train2, params);
            xcorrMat = [xcorrMat crossCorrelation];
            delayMat = [delayMat delay];
        end
        % Averaging across trials
        N = size(xcorrMat,2);
        nBinTrials(bI,pI) = N;
        shuffleCorrector = getShuffleCorrector(unitSpikes1(trialIds),unitSpikes2(trialIds),params);
        V(:,bI,pI) = (getCorrelogramVar(unitSpikes1(trialIds),unitSpikes2(trialIds),N,params));
        shuffleCorrectorMat(:,bI,pI) = shuffleCorrector;
        trialAvgXCorr = nanmean(xcorrMat,2);
        unitAvgXCorr(:,bI,pI) = trialAvgXCorr-shuffleCorrector;            
        pairDelay(bI,pI) = nanmean(delayMat);
    end
end

colorCell = {'r','b','k'};
semCorr = std(unitAvgXCorr,[],3)./sqrt(length(unitAvgXCorr(1,1,:)));
meanCorr = nanmean(unitAvgXCorr,3);
figure; hold on; 
% Plotting xcorr and error bars
handleCell = cell(nBins,1);
for bI = [1,3]
    [handleCell{bI}, ~] = boundedline(timeLags, meanCorr(:,bI), semCorr(:,bI), colorCell{bI}, 'alpha');
end
xlabel('Time (s)')
ylabel('Correlation')
%% Plotting significant times
nTimes = length(trialAvgXCorr);
group = [];
for bI = 1:nBins
    group = [group; bI.*ones(length(pairs),1)];
end
p=[];
for i = 1:nTimes
    y = squeeze(unitAvgXCorr(i,:,:)).';
    y = y(:);    
    [p(i),~,~]=anova1(y,group,'off');
end
legend([handleCell{1} handleCell{3}],{'low','high'});


%% Getting stats
nPairs = size(pairs,1);
summary_sigma = sqrt(sum(V,3))./nPairs;
timeLagLimit = 0.2;
timeLagIdx = [timeLags>=-1*timeLagLimit&timeLags<=timeLagLimit];
alpha = 0.05; 
figure; hold on; 
titleCell = {'low','med','high'};
bN = 0;
plottedBins = [1,3];
for bI = plottedBins
    bN = bN + 1;
    % Get significance
    binMean = meanCorr((timeLagIdx),bI);
    binSigma = summary_sigma(timeLagIdx,bI);
    p = zeros(length(binMean),1);
    for i = 1:length(find(timeLagIdx))
        p(i) = 1-normcdf(binMean(i),0,binSigma(i));
    end
    
    [sorted_p,sort_i] = sort(p);
    k_comparison = zeros(length(p),1);
    for k = 1:length(p)
        k_comparison(k) = sorted_p(k) <= k*alpha/length(p);
    end
    selectedTimes = timeLags(timeLagIdx);
    sigIdxs = sort_i(k_comparison==1);
    sigTimes = selectedTimes(sigIdxs);    
    subplot(1,length(plottedBins),bN); hold on;    
    bar(selectedTimes,meanCorr((timeLagIdx),bI),'hist')
    plot(selectedTimes,2.*summary_sigma(timeLagIdx,bI),['--r'])
    plot(selectedTimes,-2.*summary_sigma(timeLagIdx,bI),['--r'])
    title(titleCell{bI})
    xlabel('Time(s)')
    ylabel('Mean counts')
    % Plot sig stars
    selectedCorr = meanCorr((timeLagIdx),bI);
    for sigI = 1:length(sigTimes)
        sigIdx = sigIdxs(sigI);
        text(sigTimes(sigI),1.1*selectedCorr(sigIdx),'*','HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle')
    end
    xlim([selectedTimes(1)-0.016 selectedTimes(end)+0.016])
end



%% Plot integrals (divide positive and negative)
positive_times = [timeLags>=-1*timeLagLimit&timeLags<=timeLagLimit&timeLags>0];
negative_times = [timeLags>=-1*timeLagLimit&timeLags<=timeLagLimit&timeLags<0];
time_sums_p = squeeze(nansum(unitAvgXCorr(positive_times,:,:),1)).*params.binStep;
time_sums_n = squeeze(nansum(unitAvgXCorr(negative_times,:,:),1)).*params.binStep;
sum_mean_p = mean(time_sums_p,2);
sum_mean_n = mean(time_sums_n,2);
sum_sem_p = std(time_sums_p,[],2)./sqrt(length(time_sums_p));
sum_sem_n = std(time_sums_n,[],2)./sqrt(length(time_sums_n));

numel = size(time_sums_p,1)*size(time_sums_p,2);
y = [reshape(time_sums_p.',[numel,1]); reshape(time_sums_n.',[numel,1])];
% Grouping for ANOVA
binVector = [];
for k = 1:2
    for i = 1:nBins
        binVector = [binVector; i.*ones(size(time_sums_p,2),1)];
    end
end
pn_vector = [ones(numel,1);-1.*ones(numel,1)];

%% Integral figure

figure; hold on;
barY = [sum_mean_n(plottedBins) sum_mean_p(plottedBins)];
errY = [sum_sem_n(plottedBins) sum_sem_p(plottedBins)];
bar(1:length((plottedBins)),barY)                

ngroups = size(barY, 1);
nbars = size(barY, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, barY(:,i), errY(:,i), 'k', 'linestyle', 'none');
end
ylabel('Crosscorrelogram integral')
set(gca,'XTick',1:nBins)
set(gca,'XTickLabels',{'low','high'})

title(['Correlogram integrals by ' regressorName], 'Interpreter', 'none')
% legend({'vmPFC leads','Amygdala leads'})
legend({'dACC leads','preSMA leads'})

%% ANOVA stats
middleTrials = binVector==2;
highAndLowBin = binVector; highAndLowBin(middleTrials) = [];
highAndLowArea = pn_vector; highAndLowArea(middleTrials) = [];
highAndLow_y = y; highAndLow_y(middleTrials) = [];

[p,tbl,stats] = anovan(highAndLow_y,{highAndLowBin.' highAndLowArea.'},'model','interaction','varnames',{'bin','P/N'},'display','off');


end
