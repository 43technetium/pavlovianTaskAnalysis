% This function aggrgates rating and RT analyses for the Pavlovian
% Conditioning Task
%
% Trial type labeling:
% 1: CSd+/CSp+
% 2: CSd-/CSp-
% 3: CSd+/CSp-
% 4: CSd-/CSp+
%
% taquino jul/18

function behaviorAnalysis
dbstop if error

sessions = {'P49CS','P51CS_1','P51CS_2','P53CS','P54CS','P58CS',...
    'P60CS','P61CS','P62CS','P63CS','P70CS','P71CS','P76CS'};
subIDs = {'49','51','51','53','54','58','60','61','62','63','70','71','76'};
% Which blocks are valid for each session (out of 1,2,3,4)
blocksAvailable = {[1 2 3 4],[1 2 3 4],[1],[1 2 3 4],...
    [1 2 3 4],[1 2 3 4],[1,2,3,4],[1,2,3,4],[1,2,3,4],...
    [1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3,4]};
% Excluded TTL trials (to account for correct RT labeling)
excludedTTLTrials = {[],[1],[1],[],[],[1],[],[],[],[],[],[]}; %#ok<NBRAK,NASGU>

dataBasepath = '..\patientData\';
nSessions = length(sessions);
nStimuli = 22;
nTrialsPerBlock = 24;

% Creating matrix with all rating changes across sessions
allMeanChanges = [];
allAbsoluteChanges = [];
% Creating matrix with the mean of all CS ratings across sessions
allCSRatings = [];
allRTs = [];

% Looping over selected sessions
changeCell = cell(nSessions,1);
allCSdPlusRatings = [];
allCSdMinusRatings = [];
for sI = 1:nSessions    
    behaviorPath = [dataBasepath sessions{sI} '\behavior\'];        
    % Defining data structures
    nAvailableBlocks = length(blocksAvailable{sI});
    allRatings = NaN(nAvailableBlocks+1, nStimuli);
    Z = NaN(nAvailableBlocks+1, nStimuli);
    sideChangeMatrix = zeros(nAvailableBlocks,nTrialsPerBlock);
    trialTypes =  NaN(nAvailableBlocks, nTrialsPerBlock);
    CSMatrix = zeros(nAvailableBlocks,4);
    changeMatrix  = zeros(nAvailableBlocks,4);
    absoluteChangeMatrix  = zeros(nAvailableBlocks,4);
    CSRatingMatrix = zeros(nAvailableBlocks,4);
    blockStructCell = cell(1,4);
    % Analyzing each block separately
    sessionParams.blocksAvailable = blocksAvailable{sI}; %#ok<STRNU>
    for bI = blocksAvailable{sI}
        behaviorFile = [behaviorPath subIDs{sI} '_Sub_' num2str(bI) '_Block.mat']; 
        taskStruct = load(behaviorFile);
        blockStructCell{bI} = taskStruct.blockStructCell{bI};
        % Getting which stimuli were rated
        if bI == 1
            initialRatingStimuli = taskStruct.initialRatings.ratingVector;
            initialRatings = taskStruct.initialRatings.initialRatingVector;            
            % Filling in rating matrix with initial ratings
            for rI = 1:length(initialRatingStimuli)
                allRatings(1,initialRatingStimuli(rI)) = initialRatings(rI);
            end
            blockMeans = nanmean(allRatings,2);
            blockStds = nanstd(allRatings,[],2);
            Z(1,:) = (allRatings(1,:) - blockMeans(1))./blockStds(1);
        end
        if isfield(taskStruct.blockStructCell{1,bI},'finalRatings')
            finalRatingStimuli = taskStruct.blockStructCell{1,bI}.finalRatings.ratingVector;
            % Getting stimulus ratings        
            finalRatings = taskStruct.blockStructCell{1,bI}.finalRatings.finalRatingVector;
            % Filling in rating matrix with final ratings
            for rI = 1:length(finalRatingStimuli)
                allRatings(bI+1,finalRatingStimuli(rI)) = finalRatings(finalRatingStimuli(rI));
            end
        end
        % Getting which stimuli were used in each condition for this block
        CSdPlus = taskStruct.stimuliStructCell{1,bI}.CSdPlus;
        CSdMinus = taskStruct.stimuliStructCell{1,bI}.CSdMinus;
        CSpPlus = taskStruct.stimuliStructCell{1,bI}.CSpPlus;
        CSpMinus = taskStruct.stimuliStructCell{1,bI}.CSpMinus;   
        
        CSMatrix(bI,1) = CSdPlus; CSMatrix(bI,2) = CSdMinus;
        CSMatrix(bI,3) = CSpPlus; CSMatrix(bI,4) = CSpMinus;
        
        allCSdPlusRatings = [allCSdPlusRatings; allRatings(1,CSdPlus)];
        allCSdMinusRatings = [allCSdMinusRatings; allRatings(1,CSdMinus)];
        
        % Getting rating changes given task structure
        % Z-scoring all ratings (zscore can't deal with nans)
        blockMeans = nanmean(allRatings,2);
        blockStds = nanstd(allRatings,[],2);
        Z(bI+1,:) = (allRatings(bI+1,:) - blockMeans(bI+1))./blockStds(bI+1);
        % Getting rating changes for each stimulus after every block
        % Getting the last block in which the stimulus was rated        
        changeMatrix(bI,1) = Z(bI+1,CSdPlus)-findLastRating(Z,bI,CSdPlus);
        changeMatrix(bI,2) = Z(bI+1,CSdMinus)-findLastRating(Z,bI,CSdMinus);
        changeMatrix(bI,3) = Z(bI+1,CSpPlus)-findLastRating(Z,bI,CSpPlus);
        changeMatrix(bI,4) = Z(bI+1,CSpMinus)-findLastRating(Z,bI,CSpMinus);               
        % Getting absolute rating changes for each stimulus after every block
        absoluteChangeMatrix(bI,1) = allRatings(bI+1,CSdPlus)-findLastRating(allRatings,bI,CSdPlus);
        absoluteChangeMatrix(bI,2) = allRatings(bI+1,CSdMinus)-findLastRating(allRatings,bI,CSdMinus);
        absoluteChangeMatrix(bI,3) = allRatings(bI+1,CSpPlus)-findLastRating(allRatings,bI,CSpPlus);
        absoluteChangeMatrix(bI,4) = allRatings(bI+1,CSpMinus)-findLastRating(allRatings,bI,CSpMinus); 
        % Getting z-scored ratings for each stimulus type
        CSRatingMatrix(bI,1) = Z(bI+1,CSdPlus);
        CSRatingMatrix(bI,2) = Z(bI+1,CSdMinus);
        CSRatingMatrix(bI,3) = Z(bI+1,CSpPlus);
        CSRatingMatrix(bI,4) = Z(bI+1,CSpMinus);        
        % Initializing the first last CSp side to the first trial
        lastSide = taskStruct.blockStructCell{1,bI}.CSpPosition(1); 
        % Trial loop
        for tI = 1:nTrialsPerBlock
            % Getting RTs for every trial
            % This is aligned to the end of button press
            RT = taskStruct.blockStructCell{1,bI}.trialStructCell{1,tI}.CSpResponseTime;
            allRTs = [allRTs RT];
            % Determining whether there was a side change on this trial            
            currentSide = taskStruct.blockStructCell{1,bI}.CSpPosition(tI);
            if lastSide == currentSide
                sideChangeMatrix(bI,tI) = 0;
            else
                sideChangeMatrix(bI,tI) = 1;
            end
            lastSide = currentSide;
        end        
        % Getting trial types
        trialTypes(bI,:) = taskStruct.blockStructCell{1,bI}.trialTypes;
    end
    % Getting mean change matrix for each session
    changeCell{sI} = changeMatrix;
    allMeanChanges = [allMeanChanges; changeMatrix];
    allAbsoluteChanges = [allAbsoluteChanges; absoluteChangeMatrix];
    allCSRatings = [allCSRatings; CSRatingMatrix];
    
    meanChange = mean(absoluteChangeMatrix,1);
    s = size(absoluteChangeMatrix); nRatings = s(1);
    semChange = std(absoluteChangeMatrix,[],1)./sqrt(nRatings);
    hVector = zeros(1,4);
       
end

% Testing whether initial ratings was different for CSdPlus vs. CSdMinus
[h,p] = ttest2(allCSdPlusRatings,allCSdMinusRatings);

% Plotting summary of rating changes
% CSdPlus, CSdMinus, CSpPlus, CSpMinus
% Change matrix: CSd+/CSd-/CSp+/CSp-
nRatings = size(allMeanChanges,1);
[h,p] = ttest2(allMeanChanges(:,1),allMeanChanges(:,2),'Tail','right');
[h,p] = ttest2(allMeanChanges(:,3),allMeanChanges(:,4),'Tail','right');
mdl = fitglm([[ones(nRatings,1);zeros(nRatings,1);ones(nRatings,1);zeros(nRatings,1)] ...
    [ones(nRatings,1);ones(nRatings,1);zeros(nRatings,1);zeros(nRatings,1)]],allMeanChanges(:));   

meanChange = nanmean(allMeanChanges,1);
s = size(allMeanChanges); nRatings = s(1);
semChange = nanstd(allMeanChanges,[],1)./sqrt(nRatings);
hVector = zeros(1,4);
figure; hold on;
ylabel('Mean rating change')
title('Change in Z-scored ratings (all sessions)')
% Getting bar centers
hBar = bar(1:4,meanChange);
for i = 1:4
    if i == 1 || i == 3 % Hypothesizing increment for +
        [hVector(i),p(i)] = ttest(allMeanChanges(:,i),0,'Tail','right');
    else % Hypothesizing shrinking for -
        [hVector(i),p(i)] = ttest(allMeanChanges(:,i),0,'Tail','left');
    end
    if hVector(i)
        ctr = bsxfun(@plus, hBar.XData, [hBar.XOffset]');
        plot(ctr(i), meanChange(1,i)*1.5, '*k')
    end
end
errorbar(1:4,meanChange,semChange,'.')
labelNames = {'CSd+';'CSd-';'CSp+';'CSp-'};
set(gca, 'XTickLabel',labelNames, 'XTick',1:numel(labelNames))

% Plotting Z-score rating aggregating positive vs. negative CS
aggZChange = [[allMeanChanges(:,1); allMeanChanges(:,3)] ...
    [allMeanChanges(:,2); allMeanChanges(:,4)]];
meanChange = nanmean(aggZChange,1);
semChange = nanstd(aggZChange,[],1)./sqrt(size(aggZChange,1));
figure; hold on;
hBar = bar(1:2,meanChange);
errorbar(1:2,meanChange,semChange,'.')
[~,p_plus] = ttest(aggZChange(:,1),0,'Tail','right');
[~,p_minus] = ttest(aggZChange(:,1),0,'Tail','left');
[~,p] = ttest2(aggZChange(:,1),aggZChange(:,2),'Tail','right')
labelNames = {'CS+';'CS-'};
ylabel('Mean rating change')
title('Change in Z-scored ratings (all sessions)')
set(gca, 'XTickLabel',labelNames, 'XTick',1:numel(labelNames))

% Plotting summary of absolute rating changes
% CSdPlus, CSdMinus, CSpPlus, CSpMinus
meanChange = nanmean(allAbsoluteChanges,1);
[h,p] = ttest(allMeanChanges);
s = size(allAbsoluteChanges); nRatings = s(1);
semChange = nanstd(allAbsoluteChanges,[],1)./sqrt(nRatings);
hVector = zeros(1,4);
figure; hold on;
ylabel('Mean rating change')
title('Change in raw ratings (all sessions)')
labelNames = {'CSd+';'CSd-';'CSp+';'CSp-'};
% Getting bar centers
hBar = bar(1:4,meanChange);
for i = 1:4
    hVector(i) = ttest(allAbsoluteChanges(:,i));
    if hVector(i)
        ctr = bsxfun(@plus, hBar.XData, [hBar.XOffset]');
        plot(ctr(i), meanChange(1,i)*1.5, '*k')
    end
end
errorbar(1:4,meanChange,semChange,'.')
set(gca, 'XTickLabel',labelNames, 'XTick',1:numel(labelNames))

% Plotting summary of ratings by stimulus type
meanRatings = nanmean(allCSRatings,1);
semRatings = nanstd(allCSRatings,[],1)./sqrt(nRatings);
figure; hold on;
bar(1:4,meanRatings)
errorbar(1:4,meanRatings,semRatings,'.')
set(gca, 'XTickLabel',labelNames, 'XTick',1:numel(labelNames))
ylabel('Mean rating')
title('Mean CS Z-scored ratings across subjects')

% Getting absolute (normed) rating change
normZChange = abs([[allMeanChanges(:,1); allMeanChanges(:,2)] ...
    [allMeanChanges(:,3); allMeanChanges(:,4)]]);
[~,p_abs] = ttest(normZChange(:,1),normZChange(:,2));

end

% Get the last valid rating for comparisons
function r = findLastRating(Z,bI,stimID)
stimRatings = Z(:,stimID);
possibleRatings = stimRatings(1:bI);
idx = find(~isnan(possibleRatings),1,'last');
if isempty(idx)
    r = nan;
else
    r = possibleRatings(idx);
end
end