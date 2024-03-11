% eyeTrackingAnalysis is the main file for setting up eye tracking data for
% all patients (looking at blink rate, saccade latency, pupil diameter)
%
% taquino/jan20

function eyeTrackingAnalysis_setup()
dbstop if error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Folder that contains analysis structcure
basefolder = 'pavlovianConditioningTask\';
cd([basefolder 'dataAnalysis\behavior']);
% The individual sessions to analyse
% P39 was halted for deplantation - old version
% P47 had ET problems - old version
% P53 data is not on this PC
sessions = {'P49CS','P51CS_1','P51CS_2','P54CS','P58CS','P60CS', ...
    'P61CS','P62CS','P63CS','P70CS','P71CS'};
subIDs = {'49','51','51','54','58','60','61','62','63','70','71'};
% Which blocks are valid for each session (out of 1,2,3,4)
blocksAvailable = {[1 2 3 4],[1 2 3 4],[1],[1 2 3 4],[1 2 3 4], ...
    [1 2 3 4],[1 2 3 4],[1 2 3 4],[1 2 3 4],[1 2 3 4],[2 3 4]};


trialTypeVector = [];
trialBlock = [];
CSdBlinks_allSessions = [];
CSpBlinks_allSessions = [];
allPupilNanRatio = [];
timesVector = {};
pupilVector = {};
sessionVector = {};
allRTs = [];
dMin_vector = [];
pMin_vector = [];
nSessions = length(sessions);
allSaccEndTimesVector = {};
allSaccTrialXendVector = {};
allSaccTrialYendVector = {};
allCalAcc = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looping over sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for sI = 1:nSessions    
    eyeTrackingData = struct();
    display(sessions{sI})
    behaviorPath = [basefolder 'patientData\' sessions{sI} '\behavior\'];                
    edfFiles = {};
    nTrialsVector = [];    
    blockStructCell_data = {};
    for bI = 1:length(blocksAvailable{sI})
        bID = blocksAvailable{sI}(bI);        
        % Loading behavior for this block
        behaviorFile = [behaviorPath subIDs{sI} '_Sub_' num2str(bID) '_Block.mat']; 
        load(behaviorFile)
        blockStructCell_data{bID} = blockStructCell{bID};
        trialTypeVector = [trialTypeVector trialTypes];
        nTrialsVector = [nTrialsVector length(trialTypes)];                
        trialBlock = [trialBlock bID*ones(1,nTrialsVector(bI))];         
        % Loading eyetracking for this block        
        edfFile = dir([behaviorPath '*b' num2str(bID) '.edf']);
        edfFiles{bID} = edfFile;        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setting up eyetracking variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timesVector_session = {};
    pupilVector_session = {};
    blinkVector = {};
    saccTrialTimesVector = {};
    saccEndTimesVector = {};
    saccTrialXVector = {};
    saccTrialXendVector = {};
    saccTrialYVector = {};
    saccTrialYendVector = {};
    fixationStartTrialVector = {};
    fixationEndTrialVector = {};
    fixationPosXTrialVector = {};
    fixationPosYTrialVector = {};
    fixationPupilTrialVector = {};
    firstSaccadesUSVector = [];
    trialNameCell = {};
    
    % Timestamp vectors
    trialTimeCSdONVec = [];
    trialTimeCSdOFFVec = [];
    trialTimeCSpONVec = [];
    trialTimeCSpOFFVec = [];
    trialTimeUSVec = [];
    trialTimeEndVec = [];    
    timeReferenceVec = [];
    timeReferenceDVec = [];
    timeReferencePVec = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Looping over blocks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for slice = 1:length(blocksAvailable{sI})
        sliceID = blocksAvailable{sI}(slice);
        display(['block: ' num2str(sliceID)])
        nTrials = nTrialsVector(slice);        
        filename = edfFiles{sliceID};
        file_path = [behaviorPath filename.name];
        % Setting useOldProcedure to false crashes MATLAB!
        % This function doesn't like filenames with spaces!
        edf = Edf2Mat(file_path, true);
        pupilSize = edf.Samples.pupilSize;
        
        % Getting calibration accuracy
        calIndex = find(contains(edf.Events.Messages.info,'!CAL VALIDATION'));        
        if ~isempty(calIndex)
            calIndex = calIndex(end);
            k = strfind(edf.Events.Messages.info{calIndex},'avg');
            calAcc = str2num(edf.Events.Messages.info{calIndex}(k-5:k-2));
        else
            calAcc = nan;
        end
        allCalAcc = [allCalAcc; calAcc];
        blinkMargin = 50; % Removing data points within a certain time margin from each blink
        blinkTimestamps = []; % Getting all the blink onsets
        % Discarding pupil sizes during blink events (as tracked by EyeLink)
        blinkTimesVector = zeros(length(edf.Events.Eblink.end),2);
        for blink = 1:length(edf.Events.Eblink.end)
            blinkStart = edf.Events.Eblink.start(blink);
            blinkTimestamps = [blinkTimestamps blinkStart];
            blinkEnd = edf.Events.Eblink.end(blink);
            blinkingIdx = edf.Samples.time > blinkStart-blinkMargin & edf.Samples.time < blinkEnd+blinkMargin;
            pupilSize(blinkingIdx) = nan;
            blinkTimesVector(blink,1) = blinkStart;
            blinkTimesVector(blink,2) = blinkEnd;
        end
        % Removing pupil sizes further than 5 s.d. from the average
        m = nanmean(pupilSize);
        s = nanstd(pupilSize);
%         pupilSize(pupilSize>m+3*s|pupilSize<m-3*s) = nan;
        pupilSize(pupilSize>m+5*s|pupilSize<m-5*s) = nan;

        % Interpolating and filtering with moving average
        ipPupilSize = pupilSize;
        %substitution of nan with closest value
        toBeRemoved = false(size(ipPupilSize,1),size(ipPupilSize,2));
        for lineIndex = 1:size(ipPupilSize,2)
            for rowIndex = size(ipPupilSize,1):-1:2
                CurrentValue = ipPupilSize(rowIndex, lineIndex);
                if ~isnan(CurrentValue) && isnan(ipPupilSize(rowIndex-1, lineIndex))
                    ipPupilSize(rowIndex-1, lineIndex) = CurrentValue;
                    toBeRemoved (rowIndex,lineIndex) = 1;
                end
            end
        end
        window = 50;
        b = ones(1, window)/window; % Rational transfer function
        filteredPupilSize = filter(b, 1, ipPupilSize);
        % Looping over trials
        for i = 1:nTrials
            % This is aligned to the end of button press
            RT = blockStructCell_data{1,sliceID}.trialStructCell{1,i}.CSpResponseTime; 
            allRTs = [allRTs RT];
            sessionVector = [sessionVector sessions{sI}];
            trialName = ['TRIALID ' num2str(i)];
            % Finding current trial on edf file
            ind = find(ismember(edf.Events.Messages.info,trialName));
            % Flagging trials with a different structure
            numTrialMessages = length(ind);
            ind = ind(1);
            % A trial has 8 events on the edf file:
            trialIdx = ind-1:ind+11;
            trialTimestamps = (edf.Events.Messages.time(trialIdx)).';

            % Organizing timestamps for each trial
            if numTrialMessages ~= 4
                trialTimeCSdON = trialTimestamps(1);
                trialTimeCSdOFF = trialTimestamps(3);
                trialTimeCSpON = trialTimestamps(5);
                trialTimeCSpOFF = trialTimestamps(11);
                trialTimeUS = trialTimestamps(13);
                trialTimeEnd = trialTimestamps(13)+2000; % Videos take 2s
            else
                trialTimeCSdON = trialTimestamps(1);
                trialTimeCSdOFF = trialTimestamps(3);
                trialTimeCSpON = trialTimestamps(5);
                trialTimeCSpOFF = trialTimestamps(7);
                trialTimeUS = trialTimestamps(9);
                trialTimeEnd = trialTimestamps(9)+2000; % Videos take 2s
            end
            trialTimeCSdONVec = [trialTimeCSdONVec trialTimeCSdON];
            trialTimeCSdOFFVec = [trialTimeCSdOFFVec trialTimeCSdOFF];
            trialTimeCSpONVec = [trialTimeCSpONVec trialTimeCSpON];
            trialTimeCSpOFFVec = [trialTimeCSpOFFVec trialTimeCSpOFF];
            trialTimeUSVec = [trialTimeUSVec trialTimeUS];
            trialTimeEndVec = [trialTimeEndVec trialTimeEnd];            
            timeReferenceD = trialTimeCSdON; % Distal stim baseline ref
            timeReferenceP = trialTimeCSpON; % Proximal stim baseline ref
            % Counting blinks in this trial
            nCSdBlinks = sum(blinkTimesVector(:,1)>trialTimeCSdON&blinkTimesVector(:,1)<trialTimeCSdOFF);
            nCSpBlinks = sum(blinkTimesVector(:,1)>trialTimeCSpON&blinkTimesVector(:,1)<trialTimeCSpOFF);
            CSdBlinks_allSessions = [CSdBlinks_allSessions nCSdBlinks];
            CSpBlinks_allSessions = [CSpBlinks_allSessions nCSpBlinks];
            % Taking one of the baseline references to align trial times
            timeReference = timeReferenceD;
            timeReferenceVec = [timeReferenceVec timeReference];
            timeReferenceDVec = [timeReferenceDVec timeReferenceD];
            timeReferencePVec = [timeReferencePVec timeReferenceP];
            pupilTimeIdx = edf.Samples.time > trialTimeCSdON & edf.Samples.time < trialTimeEnd;
            blinkTimeIdx = blinkTimestamps > trialTimeCSdON & blinkTimestamps < trialTimeEnd;
            saccTimeIdx = edf.Events.Esacc.start > trialTimeCSdON & edf.Events.Esacc.start < trialTimeEnd;
            saccUSIdx = edf.Events.Esacc.start > trialTimeUS & edf.Events.Esacc.start < trialTimeEnd;
            fixationIdx = edf.Events.Efix.start > trialTimeCSdON & edf.Events.Efix.start < trialTimeEnd;
            % Creating a window of 1s before the alignment time and
            % averaging pupil sizes to take a reference.
            referenceWindow = 1000;
            referenceDTimeIdx = edf.Samples.time > timeReferenceD-referenceWindow ...
                & edf.Samples.time < timeReferenceD & edf.Samples.time > 0;            
            trialRefDSize = nanmean(filteredPupilSize(referenceDTimeIdx));            
            referencePTimeIdx = edf.Samples.time > timeReferenceP-referenceWindow ...
                & edf.Samples.time < timeReferenceP & edf.Samples.time > 0;            
            trialRefPSize = nanmean(filteredPupilSize(referencePTimeIdx));                      
            % Aligning trials by reference time
            pupilTrialTimes = edf.Samples.time(pupilTimeIdx) - timeReference;
            blinkTrialTimes = blinkTimestamps(blinkTimeIdx) - timeReference;
            saccTrialTimes = edf.Events.Esacc.start(saccTimeIdx) - timeReference;
            saccEndTimes = edf.Events.Esacc.end(saccTimeIdx) - timeReference;
            saccTrialX = edf.Events.Esacc.posX(saccTimeIdx);
            saccTrialXend = edf.Events.Esacc.posXend(saccTimeIdx);
            saccTrialY = edf.Events.Esacc.posY(saccTimeIdx);
            saccTrialYend = edf.Events.Esacc.posYend(saccTimeIdx);
            rawPupilTrialSizes = pupilSize(pupilTimeIdx);
            pupilNanRatio = sum(isnan(rawPupilTrialSizes))/length(rawPupilTrialSizes);
            allPupilNanRatio = [allPupilNanRatio pupilNanRatio];
            
            pupilTrialSizes = filteredPupilSize(pupilTimeIdx);
            fixationStartTrial = edf.Events.Efix.start(fixationIdx);
            fixationEndTrial = edf.Events.Efix.end(fixationIdx);
            fixationPosXTrial = edf.Events.Efix.posX(fixationIdx);
            fixationPosYTrial = edf.Events.Efix.posY(fixationIdx);
            fixationPupilTrial = edf.Events.Efix.pupilSize(fixationIdx);
            % Selecting first saccade after onset of US
            saccUS = edf.Events.Esacc.start(saccUSIdx) - trialTimeUS;                        
            % Normalizing pupil size according to reference
            timesVector_session = [timesVector_session pupilTrialTimes];
            % Pupil normalization is done in two stages:
            % -> Prior to CSd fixation cross (D reference)
            % -> After CSp fixation cross (P reference)
            referenceChangeTime = trialTimeCSpON-timeReference-referenceWindow;
            refChangeIdx = find(pupilTrialTimes < referenceChangeTime, 1, 'last');
            % Section normalized by distal reference
            % Proportional normalizing
            normalizedPupilD = (pupilTrialSizes(1:refChangeIdx)-trialRefDSize)/trialRefDSize;
            % Absolute normalizing
%             normalizedPupilD = (pupilTrialSizes(1:refChangeIdx)-trialRefDSize);
            % Section normalized by proximal or distal reference
            normalizedPupilP = (pupilTrialSizes(refChangeIdx:end)-trialRefDSize)/trialRefDSize;
%             normalizedPupilP = (pupilTrialSizes(refChangeIdx:end)-trialRefPSize)/trialRefPSize;
%             normalizedPupilP = (pupilTrialSizes(refChangeIdx:end)-trialRefDSize);
            newTrialPupil = [normalizedPupilD; normalizedPupilP];
            [dMin,~] = min(normalizedPupilD);
            [pMin,~] = min(normalizedPupilP);
            dMin_vector = [dMin_vector dMin];
            pMin_vector = [pMin_vector pMin];
            pupilVector_session = [pupilVector_session newTrialPupil];
            % Saccade times at each trial (inserting nan if no saccades (highly unlikely))
            if ~isempty(saccTrialTimes)
                saccTrialTimesVector = [saccTrialTimesVector saccTrialTimes];
                saccEndTimesVector = [saccEndTimesVector saccEndTimes];
                saccTrialXVector = [saccTrialXVector saccTrialX];
                saccTrialXendVector = [saccTrialXendVector saccTrialXend];
                saccTrialYVector = [saccTrialYVector saccTrialY];
                saccTrialYendVector = [saccTrialYendVector saccTrialYend];
            else
                saccTrialTimesVector = [saccTrialTimesVector nan];
                saccEndTimesVector = [saccEndTimesVector nan];
                saccTrialXVector = [saccTrialXVector nan];
                saccTrialXendVector = [saccTrialXendVector nan];
                saccTrialYVector = [saccTrialYVector nan];
                saccTrialYendVector = [saccTrialYendVector nan];
            end
            % Blink times at each trial (inserting nan if no blinks)
            if ~isempty(blinkTrialTimes)
                blinkVector = [blinkVector blinkTrialTimes];
            else
                blinkVector = [blinkVector nan];
            end

            % Saccades after US (inserting nan if no saccades)
            if ~isempty(saccUS)
                saccUS = saccUS(1);
                firstSaccadesUSVector = [firstSaccadesUSVector saccUS];
            else
                firstSaccadesUSVector = [firstSaccadesUSVector nan];
            end

            % Fixation vectors for this trial (inserting nan if no fixation)
            if ~isempty(fixationStartTrial)
                fixationStartTrialVector = [fixationStartTrialVector fixationStartTrial];
                fixationEndTrialVector = [fixationEndTrialVector fixationEndTrial];
                fixationPosXTrialVector = [fixationPosXTrialVector fixationPosXTrial];
                fixationPosYTrialVector = [fixationPosYTrialVector fixationPosYTrial];
                fixationPupilTrialVector = [fixationPupilTrialVector fixationPupilTrial];
            else
                fixationStartTrialVector = [fixationStartTrialVector nan];
                fixationEndTrialVector = [fixationEndTrialVector nan];
                fixationPosXTrialVector = [fixationPosXTrialVector nan];
                fixationPosYTrialVector = [fixationPosYTrialVector nan];
                fixationPupilTrialVector = [fixationPupilTrialVector nan];
            end
            % Invalid trials generate empty vectors, so keep track of the label
            % for each valid trial here
            if ~isempty(pupilTrialSizes)
                trialNameCell = [trialNameCell trialName];
            end
        end % Trial loop          
      
    end % Block loop       
    timesVector = horzcat(timesVector,timesVector_session);
    pupilVector = horzcat(pupilVector,pupilVector_session);
    allSaccEndTimesVector{sI} = saccEndTimesVector;
    allSaccTrialXendVector{sI} = saccTrialXendVector;
    allSaccTrialYendVector{sI} = saccTrialYendVector;
end % Session loop

% Saving eyetracking data for this session
eyeTrackingData.allRTs = allRTs;
eyeTrackingData.allCalAcc = allCalAcc;
eyeTrackingData.dMin_vector = dMin_vector;
eyeTrackingData.pMin_vector = pMin_vector;
eyeTrackingData.CSdBlinks_allSessions = CSdBlinks_allSessions;
eyeTrackingData.CSpBlinks_allSessions = CSpBlinks_allSessions;
eyeTrackingData.allPupilNanRatio = allPupilNanRatio;
eyeTrackingData.sessionVector = sessionVector;
eyeTrackingData.trialTypeVector = trialTypeVector;
eyeTrackingData.trialBlock = trialBlock;
eyeTrackingData.timesVector = timesVector;
eyeTrackingData.pupilVector = pupilVector;
eyeTrackingData.allSaccEndTimesVector = allSaccEndTimesVector;
eyeTrackingData.allSaccTrialXendVector = allSaccTrialXendVector;
eyeTrackingData.allSaccTrialYendVector = allSaccTrialYendVector;

save([basefolder 'patientData\allSessions_eyetracking\eyeTrackingData_allSessions.mat'], 'eyeTrackingData')

end