% This function implements SVM decoding for Pavlovian conditioning data,
% such as in Pauli et al. (2019)

% This is the main script
% taquino oct/21

% v2: lasso regression
function svmDecoding()
% Trial type labeling:
% 1: CSd+/CSp+
% 2: CSd-/CSp-
% 3: CSd+/CSp-
% 4: CSd-/CSp+

dbstop if error
% sessions = {'P47CS'}; % Different version

% Without P51CS_2 (only 1 block)
sessions = {'P49CS','P51CS_1','P53CS','P54CS','P58CS',...
    'P60CS','P61CS','P62CS','P63CS','P70CS','P71CS','P76CS'};

blockTrials = [1:24;25:48;49:72;73:96];

% Leftout blocks 
CV_block_folds = {1,2,3,4};
nFolds = length(CV_block_folds);
basefolder = 'C:\Users\Tomas\Documents\PhD\OLab\pavlovianLearning\pavlovianConditioningTask\patientData\';

%% Which design matrix to use
% design = 'outcome'; nClasses = 2;
% design = 'csd';
% design = 'csp';
% design = 'selected_stimEVs';
% design = 'selectedStimEVs_block';
% design = 'EV_CSd'; nClasses = 3;
% design = 'EV_CSp'; nClasses = 3;
% design = 'csp_identity';
% design = 'csp_identity_saccade';
design = 'csp_presumed_identity'; nClasses = 2;
% design = 'csp_identity_saccade_distal';
% design = 'csp_presumed_EV'; nClasses = 3;
% design = 'SPE1'; nClasses = 3;
% design = 'SPE2'; nClasses = 3;
% design = 'RPE'; nClasses = 3;
% design = 'EV_CSd_MB'; nClasses = 3;
% design = 'EV_CSp_MB'; nClasses = 3;

% % Which time alignment to use for spikes
reference = 'trial';
% reference = 'decision';
% reference = 'outcome';
% reference = 'CSd_saccades';
% reference = 'CSp_saccades';

%% Time window arrangement
% windowing = 'standard_decision'; 
% windowing = 'standard_outcome';
% windowing = 'pre_outcome';
windowing = 'standard_trial'; 
% windowing = 'csp_window';
% windowing = 'post-saccade';
% windowing = 'saccade_windowed';
% windowing = 'pre_saccade';
% windowing = 'whole_trial';

%% Tested brain areas (6 = all)
brainAreas = [1,2,3,4,5,6];
nAreas = length(brainAreas);
nPermutations = 0;
nSessions = length(sessions);
train_accuracy = nan(nSessions,1,nFolds,nAreas,nPermutations+1);
test_accuracy = nan(nSessions,1,nFolds,nAreas,nPermutations+1);

for sI = 1:nSessions
     
    %% Loading data
    session = sessions{sI};
    behavior_folder = [basefolder '\allBehavior\'];
    behavior_data = load([behavior_folder 'sessionBehavior_' session '.mat']);
    sessionFolder = [basefolder session '\'];
    load([sessionFolder 'sessionData.mat']) 
    unitCell = sessionData.neuralData.unitCell;  
    nTrials = length(unitCell{1,1}.trialReferencedSpikes);
    nUnits = length(unitCell);

    %% Setting up windowing
    switch windowing
        case 'standard_outcome'
            % Bandit-related time windows
            binSize = 2;
            windowSpacing = 0.01;
            % For post-reward windows
            period = [0.25 2.25];
            windowStarts = period(1):windowSpacing:period(2)-binSize;
            windowEnds = windowStarts + binSize;
            nBins = length(windowStarts);
        case 'pre_outcome'
            % Bandit-related time windows
            binSize = 3;
            windowSpacing = 0.01;
            % For post-reward windows
            period = [-3 0];
            windowStarts = period(1):windowSpacing:period(2)-binSize;
            windowEnds = windowStarts + binSize;
            nBins = length(windowStarts);
        case 'standard_trial'
            % Bandit-related time windows
            binSize = 2.75;
            windowSpacing = 0.01;
            % For post-reward windows
            period = [0.25 3];
            windowStarts = period(1):windowSpacing:period(2)-binSize;
            windowEnds = windowStarts + binSize;
            nBins = length(windowStarts);
        case 'csp_window'
            % Bandit-related time windows
            binSize = 2.75;
            windowSpacing = 0.01;
            % For post-reward windows
            period = [4.25 7];
            windowStarts = period(1):windowSpacing:period(2)-binSize;
            windowEnds = windowStarts + binSize;
            nBins = length(windowStarts);
        case 'post_saccade'
            % Bandit-related time windows
            binSize = 1;
            windowSpacing = 0.01;
            % For post-reward windows
            period = [0 1];
            windowStarts = period(1):windowSpacing:period(2)-binSize;
            windowEnds = windowStarts + binSize;
            nBins = length(windowStarts);
        case 'saccade_windowed'
            % Bandit-related time windows
            binSize = 0.5;
            windowSpacing = 0.01;
            % For post-reward windows
            period = [-1 1];
            windowStarts = period(1):windowSpacing:period(2)-binSize;
            windowEnds = windowStarts + binSize;
            nBins = length(windowStarts);        
        case 'pre_saccade'
            % Bandit-related time windows
            binSize = 1;
            windowSpacing = 0.01;
            % For post-reward windows
            period = [-1 0];
            windowStarts = period(1):windowSpacing:period(2)-binSize;
            windowEnds = windowStarts + binSize;
            nBins = length(windowStarts);
        case 'whole_trial'
            % Bandit-related time windows
            binSize = 0.02;
            windowSpacing = 0.02;
            % For post-reward windows
            period = [0 11];
            windowStarts = period(1):windowSpacing:period(2)-binSize;
            windowEnds = windowStarts + binSize;
            nBins = length(windowStarts);
    end
    %% Setting up regressors
    trialTypes = behavior_data.trialTypes_flat;
    outcome = ones(nTrials,1);
    outcome(trialTypes==2) = 0;    
    outcome(trialTypes==3) = 0;
    outcome(trialTypes==4) = 1;
    % Getting CSp identity
    CSp_options = unique(behavior_data.stim2Vec);
    csp_identity = (behavior_data.stim2Vec == CSp_options(1)).';    
    % The most likely proximal identity
    csp_presumed_identity = csp_identity;
    csp_presumed_identity(trialTypes==3|trialTypes==4) = 1-csp_presumed_identity(trialTypes==3|trialTypes==4);
    % Getting EVs for CSp options (stimulus locked)
    CSp_stims = unique(behavior_data.stim2Vec);
    stimLocked_stimEVs = behavior_data.stimEVs(:,CSp_stims);
    id = csp_presumed_identity;
    csp_presumed_EV = zeros(length(id),1);
    for tI = 1:length(id)
        csp_presumed_EV(tI) = stimLocked_stimEVs(tI,id(tI)+1);
    end
    SPE1 = behavior_data.SPE1_vec;
    SPE2 = behavior_data.SPE2_vec;
    RPE = behavior_data.RPE;
    EV_CSd_MB = behavior_data.stimEVs_MB(:,1);
    EV_CSp_MB = behavior_data.stimEVs_MB(:,2);
    EV_CSd = behavior_data.EV_CSd;
    EV_CSp = behavior_data.EV_CSp;
    
    %% Getting brain areas
    areaVec = zeros(length(unitCell),1);
    for uI = 1:length(unitCell)
        brainArea = unitCell{uI,1}.unitInfo;
        areaVec(uI) = brainArea(4);
    end
    % Merging left/right side for all brain areas
    mergedAreas = areaVec;
    mergedAreas(mod(mergedAreas,2)==0)=mergedAreas(mod(mergedAreas,2)==0)-1;
    % Hippocampus/Amygdala/Anterior cingulate/Supplementary motor areas
    hipUnits = find(mergedAreas==1); amyUnits = find(mergedAreas==3);
    accUnits = find(mergedAreas==5); smaUnits = find(mergedAreas==7);
    ofcUnits = find(mergedAreas==9); 
    allUnits = find(mergedAreas==1|mergedAreas==3|mergedAreas==5|mergedAreas==7|mergedAreas==9);
    unitsByArea = {ofcUnits; amyUnits; hipUnits; accUnits; smaUnits; allUnits};
    
    switch design
        case 'outcome'
            y = outcome;
        case 'csp_identity'
            y = csp_identity;
        case 'csp_presumed_identity'
            y = csp_presumed_identity;
        case 'csp_identity_saccade'
            y = csp_identity_saccade;
        case 'csp_identity_saccade_distal'
            y = csp_identity_saccade_distal;
        case 'csp_presumed_EV'
            y = getQuantileLabels(csp_presumed_EV,nClasses);
        case 'SPE1'
            y = getQuantileLabels(SPE1,nClasses);
        case 'SPE2'
            y = getQuantileLabels(SPE2,nClasses);
        case 'RPE'
            y = getQuantileLabels(RPE,nClasses);
        case 'EV_CSd'
            y = getQuantileLabels(EV_CSd,nClasses);
        case 'EV_CSp'
            y = getQuantileLabels(EV_CSp,nClasses);
        case 'EV_CSd_MB'
            y = getQuantileLabels(EV_CSd_MB,nClasses);
        case 'EV_CSp_MB'
            y = getQuantileLabels(EV_CSp_MB,nClasses);       
    end
    for pI = 1:nPermutations+1   
        if pI > 1
            y = y(randperm(length(y)));
        end
        for aI = 1:nAreas
            display(['Session number: ' num2str(sI) ...
            '/ Permutation number: ' num2str(pI) ...
            '/ Area number: ' num2str(aI)])
            % Decoding independently for time bins
            for bI = 1:nBins
                % Iterating over time bins and testing GLM
                x = zeros(length(y),nUnits,nBins);
                for uI = 1:nUnits                         
                    if strcmp(reference,'outcome')
                        spikes = sessionData.neuralData.unitCell{uI}.outcomeReferencedSpikes;
                    elseif strcmp(reference,'trial')
                        spikes = sessionData.neuralData.unitCell{uI}.trialReferencedSpikes;    
                    elseif strcmp(reference,'decision')
                        spikes = sessionData.neuralData.unitCell{uI}.decisionReferencedSpikes; 
                    elseif strcmp(reference,'CSd_saccades')
                        spikes = sessionData.neuralData.unitCell{uI}.saccadeReferencedSpikes_CSd;
                    elseif strcmp(reference,'CSp_saccades')
                        spikes = sessionData.neuralData.unitCell{uI}.saccadeReferencedSpikes_CSp;
                    end
                    % Neural data
                    x(:,uI,bI) = cellfun(@(x) nnz(x>windowStarts(bI)&...
                        x<=windowEnds(bI)), spikes);     
                end
                x = x(:,unitsByArea{aI},bI);
                if ~isempty(x)
                    % Setting up cross-validation
                    nTrials = length(y);
                    for fI = 1:nFolds
                        fold = CV_block_folds{fI};
                        testBlocks = fold;
                        trainBlocks = setdiff(1:4,testBlocks);
                        % Getting train/test trials
                        train = []; test = [];
                        for blockI = 1:length(trainBlocks)
                            blockID = trainBlocks(blockI);
                            train = [train blockTrials(blockID,:)];
                        end
                        for blockI = 1:length(testBlocks)
                            blockID = testBlocks(blockI);
                            test = [test blockTrials(blockID,:)];
                        end

                        % Train SVM
                        if nClasses > 2
                            SVMModel = fitcecoc(x(train,:),y(train));
                        else 
                            SVMModel = fitcsvm(x(train,:),y(train));
%                             
                        end
                        [train_label] = predict(SVMModel,x(train,:));
                        [test_label] = predict(SVMModel,x(test,:));            
                        train_accuracy(sI,bI,fI,aI,pI) = sum(y(train) == train_label)./length(y(train));
                        test_accuracy(sI,bI,fI,aI,pI) = sum(y(test) == test_label)./length(y(test));
                        
                    end
                end
            end
        end
    end
end 

savefolder = [basefolder 'decodingResults/' design '/' reference '/' windowing '/'];
% Make folder if it doesn't exist
if exist(savefolder,'dir')~=7
    mkdir(savefolder);
end

save([savefolder 'decodingResults'],'test_accuracy')

end
