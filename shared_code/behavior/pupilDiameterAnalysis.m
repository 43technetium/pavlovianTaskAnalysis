% pupilDiameterAnalysis is the main file for analyzing pupil diameter data
%
% taquino/oct21

% Trial type labeling:
% 1: CSd+/CSp+
% 2: CSd-/CSp-
% 3: CSd+/CSp-
% 4: CSd-/CSp+

function pupilDiameterAnalysis()
sessions = {'P49CS','P51CS_1','P51CS_2','P54CS','P58CS',...
    'P60CS','P61CS','P62CS','P63CS','P70CS','P71CS'};
% Which blocks are valid for each session (out of 1,2,3,4)
blocksAvailable = {[1 2 3 4],[1 2 3 4],[1],...
    [1 2 3 4],[1 2 3 4],[1,2,3,4],[1,2,3,4],[1,2,3,4],...
    [1,2,3,4],[1,2,3,4],[2,3,4]};


blockTrials = [1:24;25:48;49:72;73:96];

basefolder = '..\pavlovianConditioningTask\';
load([basefolder 'patientData\allSessions_eyetracking\eyeTrackingData_allSessions.mat']);

% Ratio of datapoints removed at which we discard trial
nan_criterion = 0.2;

trialTypeVector = eyeTrackingData.trialTypeVector;
timesVector = eyeTrackingData.timesVector;
pupilVector = eyeTrackingData.pupilVector;
nanRatio = eyeTrackingData.allPupilNanRatio;

mdlCell_D = {};
mdlCell_P = {};
diameterDCell = {};
diameterPCell = {};
trialTypes_flat_cell = {};
all_trialTypes_flat = [];
all_RT = [];
all_PEd = [];
all_EV_CSp = [];
all_sessionID = [];
all_diameterP = [];
all_EV_CSd = [];
all_diameterD = [];
all_outcome = [];
removed_trial_ratio = [];
all_diameterO = [];
diameterOCell = {};
% Tracking correlation between pupil diameter and EVs
for sI = 1:length(sessions)
    trialsAvailable = reshape(blockTrials(blocksAvailable{sI},:).', ...
        [24*length(blocksAvailable{sI}),1]);
    session = sessions{sI};
    behavior_folder = [basefolder 'patientData\allBehavior\'];
    behavior_data = load([behavior_folder 'sessionBehavior_' session '.mat']);
    sessionFolder = [basefolder 'patientData\' session '\'];
    load([sessionFolder 'sessionData.mat']) 
    unitCell = sessionData.neuralData.unitCell;  
    nTrials = length(unitCell{1,1}.trialReferencedSpikes);
    nUnits = length(unitCell);
    IndexC = strfind(eyeTrackingData.sessionVector,session);
    sessionTrials = find(not(cellfun('isempty',IndexC)));
    sessionPupil = pupilVector(sessionTrials);
    sessionTimes = timesVector(sessionTrials);
    sessionNan = nanRatio(sessionTrials);
    removed_trial_ratio(sI) = sum(sessionNan > nan_criterion)/length(sessionNan);
    nTrials = length(sessionPupil);
    diameterD = zeros(nTrials,1);
    diameterP = zeros(nTrials,1);
    diameterO = zeros(nTrials,1);
    for tI = 1:nTrials
        selectedTimes_D = sessionTimes{tI}>500&sessionTimes{tI}<1300;
        selectedTimes_P = sessionTimes{tI}>4500&sessionTimes{tI}<5800;
        selectedTimes_O = sessionTimes{tI}>8000&sessionTimes{tI}<9300;
        if sessionNan(tI) < nan_criterion
            diameterD(tI) = nanmean(sessionPupil{tI}(selectedTimes_D));
            diameterP(tI) = nanmean(sessionPupil{tI}(selectedTimes_P));       
            diameterO(tI) = nanmean(sessionPupil{tI}(selectedTimes_O));  
        else
            diameterD(tI) = nan;
            diameterP(tI) = nan;
            diameterO(tI) = nan;
        end
    end
    mdlCell_D{sI} = fitglm(behavior_data.EV_CSd(trialsAvailable),diameterD);   
    mdlCell_P{sI} = fitglm(behavior_data.EV_CSp(trialsAvailable),diameterP);
    diameterDCell{sI} = diameterD;
    diameterPCell{sI} = diameterP;
    diameterOCell{sI} = diameterO;
    trialTypes_flat_cell{sI} = behavior_data.trialTypes_flat(trialsAvailable);    
    all_trialTypes_flat = [all_trialTypes_flat; behavior_data.trialTypes_flat(trialsAvailable)];
    all_RT = [all_RT; behavior_data.rtVec(trialsAvailable).'];
    
    all_PEd = [all_PEd; behavior_data.SPE1_vec(trialsAvailable)];
    all_EV_CSd = [all_EV_CSd; behavior_data.stimEVs_MB(trialsAvailable,1)];
    all_EV_CSp = [all_EV_CSp; behavior_data.stimEVs_MB(trialsAvailable,2)];   
    
    outcome = behavior_data.RPE;
    outcome(outcome>0) = 1;
    outcome(outcome<0) = 0;
    all_outcome = [all_outcome; outcome(trialsAvailable)];
    
    all_diameterD = [all_diameterD; diameterD];
    all_diameterO = [all_diameterO; diameterO];
    all_sessionID = [all_sessionID; sI.*ones(length(diameterD),1)];
    all_diameterP = [all_diameterP; diameterP];
end

% Saving pupil data
pupilData = struct();
pupilData.all_RT = all_RT;
pupilData.all_EV_CSd = all_EV_CSd;
pupilData.all_EV_CSp = all_EV_CSp;
pupilData.all_sessionID = all_sessionID;
pupilData.all_diameterD = all_diameterD;
pupilData.all_diameterP = all_diameterP;
pupilData.all_PEd = all_PEd;
pupilData.sessions = sessions;
save([basefolder 'patientData\pupilData\allPupilSessions.mat'],'pupilData')

%% Fitting linear mixed model for RT
tblRT = table(all_RT,all_EV_CSd, all_EV_CSp, all_sessionID, ...
    'VariableNames',{'RT','EV_CSd','EV_CSp','sessionID'});
formulaRT = 'RT~EV_CSd+EV_CSp+(1|sessionID)';
lmeRT = fitlme(tblRT,formulaRT);
%% Fitting linear mixed model for pupil diameter, PEd, EV_CSp
% 'y ~ treatment + (1|block)'

tblP = table(all_diameterP,all_EV_CSp,all_PEd,all_EV_CSp.*all_PEd, all_sessionID, ...
    'VariableNames',{'diameterP','EV_CSp','PEd','EV_CSp_PEd_interaction','sessionID'});
formulaP = 'diameterP~EV_CSp+PEd+EV_CSp_PEd_interaction+(1|sessionID)';
lmeP = fitlme(tblP,formulaP);

p_trialType_distal = [];
p_trialType_proximal = [];
for sI = 1:length(sessions)
    CSdPlus = [diameterDCell{sI}(trialTypes_flat_cell{sI}==1); diameterDCell{sI}(trialTypes_flat_cell{sI}==3)];
    CSdMinus = [diameterDCell{sI}(trialTypes_flat_cell{sI}==2); diameterDCell{sI}(trialTypes_flat_cell{sI}==4)];
    CSpPlus = [diameterPCell{sI}(trialTypes_flat_cell{sI}==1); diameterPCell{sI}(trialTypes_flat_cell{sI}==4)];
    CSpMinus = [diameterPCell{sI}(trialTypes_flat_cell{sI}==2); diameterPCell{sI}(trialTypes_flat_cell{sI}==3)];
    [~,p_trialType_distal(sI)] = ttest2(CSdPlus,CSdMinus,'Tail','right');
    [~,p_trialType_proximal(sI)] = ttest2(CSpPlus,CSpMinus,'Tail','right');
end


% Plot diameter change by trial type
d1 = all_diameterP(all_trialTypes_flat==1);
d2 = all_diameterP(all_trialTypes_flat==2);
d3 = all_diameterP(all_trialTypes_flat==3);
d4 = all_diameterP(all_trialTypes_flat==4);

barY = [nanmean(d1),nanmean(d2),nanmean(d3),nanmean(d4)];
errY = [nanstd(d1)/sqrt(sum(~isnan(d1))),nanstd(d2)/sqrt(sum(~isnan(d2))),...
    nanstd(d3)/sqrt(sum(~isnan(d3))),nanstd(d4)/sqrt(sum(~isnan(d4)))];

figure; hold on; 
bar(1:4,barY)                
ngroups = size(barY, 1);
nbars = size(barY, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(i, barY(:,i), errY(:,i), 'r', 'linestyle', 'none');
end
ylabel('% Pupil diameter change')
title('Pupil at CSp')
set(gca,'XTick',1:4)

avgTimeVector = 1:2:7000; % all

figure; hold on;
hCell = {};
colorCell = {'g','r','k','b'};
trialTypeGroups = {[1,4],[2,3]};

for i = 1:length(trialTypeGroups)
    idx = [];
    for j = 1:length(trialTypeGroups{i})
        jID = trialTypeGroups{i}(j);
        idx = union(idx,find(trialTypeVector==jID));
    end
    idx = idx.';
    [avgSizes, stdSizes] = getMeanPupilSize(timesVector, pupilVector, avgTimeVector, idx);
    
    % Plotting pupil diameter and error bars
    [hCell{i}, ~] = boundedline(avgTimeVector-avgTimeVector(1), avgSizes, ...
        stdSizes/sqrt(length(idx)), 'alpha', 'cmap', colorCell{i});
    % Making thicker lines
    set(hCell{i}, 'LineWidth', 2)

end
legend([hCell{1},hCell{2}],'CSd+','CSd-')

xlabel('Time (ms)')
ylabel('Normalized pupil response')
title('Pupil change from baseline')
end