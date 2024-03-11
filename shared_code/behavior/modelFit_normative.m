% This function fits some specified model to behavioral data from one 
% Pavlovian task session
% 
% This function is called by modelComparison.m
% taquino jul/18
%% Main function
function fitData = modelFit_normative(behaviorData,model)
% Getting EV matrix from fitted model
% Normative parameters from Pauli et al., 2016
fitParams = [0.22,0.89];
% [fitData.EV,fitData.PEd,fitData.PEp,fitData.stimEVs] = getEV(fitData.fitParams,behaviorData,model);
fitData = getEV(fitParams,behaviorData,model);
fitData.fitParams = fitParams;
end % End of main function

%% Getting expected value from each trial from outcome history
function fitData = getEV(params,b,model)
transformedParams = parameterTransform(params,model);
if strcmp(model,'transitionMatrix')
    eta = transformedParams(:,1);
    nAvailableBlocks = length(b.sessionParams.blocksAvailable);
    trialEVs = zeros(nAvailableBlocks*b.taskStruct.nTrials,2); 
    SPE1_vec = zeros(nAvailableBlocks*b.taskStruct.nTrials,1);
    SPE2_vec = zeros(nAvailableBlocks*b.taskStruct.nTrials,1);
    T_cell = cell(nAvailableBlocks*b.taskStruct.nTrials,1);
    oldT = zeros(6,6);
    trialNumber = 1;
    for bI = b.sessionParams.blocksAvailable               
        % Transition matrix between stimuli and rewards (A,B,X,Y,1,0)
        T = zeros(6,6);
        % Initializing T with uniform transition priors
        T(1,3) = 0.5; T(2,3) = 0.5; T(1,4) = 0.5; T(2,4) = 0.5;
        if bI == 1
            T(3,5) = 0.5; T(4,5) = 0.5; T(3,6) = 0.5; T(4,6) = 0.5;
        else
            % Taking into account that CSp may be reversed
            if b.taskStruct.stimuliStructCell{1,bI}.CSpPlus == b.taskStruct.stimuliStructCell{1,bI-1}.CSpPlus
                T(3,5) = oldT(3,5); T(4,5) = oldT(4,5); T(3,6) = oldT(3,6); T(4,6) = oldT(4,6);
            else
                T(3,5) = oldT(4,5); T(4,5) = oldT(3,5); T(3,6) = oldT(4,6); T(4,6) = oldT(3,6);
            end
        end
        nTrials = length(b.blockStructCell{1,bI}.trialStructCell);

        % Looping over trials        
        for tI = 1:nTrials        
            trialType = b.blockStructCell{1,bI}.trialTypes(tI);
            T_cell{trialNumber} = T;
            switch trialType
                case 1 % CSd+/CSp+
                    trialEVs(trialNumber,1) = T(1,3)*T(3,5)+T(1,4)*T(4,5);
                    SPE1 = 1-T(1,3); T(1,3) = T(1,3)+eta*SPE1;
                    T(1,[1,2,4,5,6])=T(1,[1,2,4,5,6])*(1-eta);                    
                    trialEVs(trialNumber,2) = T(3,6)*0+T(3,5)*1;
                    SPE2 = 1-T(3,5); T(3,5) = T(3,5)+eta*SPE2;
                    T(3,[1,2,3,4,6])=T(3,[1,2,3,4,6])*(1-eta);                    
                case 2 % CSd-/CSp-
                    trialEVs(trialNumber,1) = T(2,3)*T(3,5)+T(2,4)*T(4,5);
                    SPE1 = 1-T(2,4); T(2,4) = T(2,4)+eta*SPE1;
                    T(2,[1,2,3,5,6])=T(2,[1,2,3,5,6])*(1-eta);
                    trialEVs(trialNumber,2) = T(4,6)*0+T(4,5)*1;
                    SPE2 = 1-T(4,6); T(4,6) = T(4,6)+eta*SPE2;
                    T(4,[1,2,3,4,5])=T(4,[1,2,3,4,5])*(1-eta);
                case 3 % CSd+/CSp-
                    trialEVs(trialNumber,1) = T(1,3)*T(3,5)+T(1,4)*T(4,5);
                    SPE1 = 1-T(1,4); T(1,4) = T(1,4)+eta*SPE1;
                    T(1,[1,2,3,5,6])=T(1,[1,2,3,5,6])*(1-eta);
                    trialEVs(trialNumber,2) = T(4,6)*0+T(4,5)*1;
                    SPE2 = 1-T(4,6); T(4,6) = T(4,6)+eta*SPE2;
                    T(4,[1,2,3,4,5])=T(4,[1,2,3,4,5])*(1-eta);
                case 4 % CSd-/CSp+
                    trialEVs(trialNumber,1) = T(2,3)*T(3,5)+T(2,4)*T(4,5);
                    SPE1 = 1-T(2,3); T(2,3) = T(2,3)+eta*SPE1;
                    T(2,[1,2,4,5,6])=T(2,[1,2,4,5,6])*(1-eta);
                    trialEVs(trialNumber,2) = T(3,6)*0+T(3,5)*1;
                    SPE2 = 1-T(3,5); T(3,5) = T(3,5)+eta*SPE2;
                    T(3,[1,2,3,4,6])=T(3,[1,2,3,4,6])*(1-eta);
            end            
            SPE1_vec(trialNumber) = SPE1;
            SPE2_vec(trialNumber) = SPE2;
                        
            trialNumber = trialNumber+1;
        end        
        oldT = T;
    end
    fitData.trialEVs = trialEVs;
    fitData.SPE1_vec = SPE1_vec;
    fitData.SPE2_vec = SPE2_vec;
    fitData.T_cell = T_cell;
    if length(trialEVs) ~= trialNumber-1
        warning('Number of trials incompatible, recheck EV assignment')
    end
elseif strcmp(model,'normative')
    % Getting free parameters
    alpha = params(1); lambda = params(2);    
    nAvailableBlocks = length(b.sessionParams.blocksAvailable);
    trialEVs = 0.5*ones(nAvailableBlocks*b.taskStruct.nTrials,2);
    trialPEd = 0.5*ones(nAvailableBlocks*b.taskStruct.nTrials,1);
    trialPEp = 0.5*ones(nAvailableBlocks*b.taskStruct.nTrials,1);
    % Initializing EVs at zero for all stimuli (A,B,X,Y)
    EVs = 0.5*ones(nAvailableBlocks*b.taskStruct.nTrials,22);
    trialNumber = 1;
    for bI = b.sessionParams.blocksAvailable
        nTrials = length(b.blockStructCell{1,bI}.trialStructCell);
        if bI == max(b.sessionParams.blocksAvailable) % last block
            nTrials = nTrials - 1;
        end
        for tI = 1:nTrials
            trialStruct = b.blockStructCell{1,bI}.trialStructCell{1,tI};
            stim1 = trialStruct.stim1; stim2 = trialStruct.stim2;
            trialType = b.blockStructCell{1,bI}.trialTypes(tI);
            % Parsing previous EV values
            EVs(trialNumber+1,:) = EVs(trialNumber,:);
            % Updating value of CSd
            EVs(trialNumber+1,stim1) = EVs(trialNumber,stim1) + alpha*(lambda*EVs(trialNumber,stim2)-EVs(trialNumber,stim1));
            trialPEd(trialNumber) = (lambda*EVs(trialNumber,stim2)-EVs(trialNumber,stim1));
            if trialType == 1 || trialType == 4 % Reward
                % Updating value of CSp
                EVs(trialNumber+1,stim2) = EVs(trialNumber,stim2) + alpha*(1-EVs(trialNumber,stim2));
                trialPEp(trialNumber) = (1-EVs(trialNumber,stim2));
            elseif trialType == 2 || trialType == 3 % No reward
                % Updating value of CSp
                EVs(trialNumber+1,stim2) = EVs(trialNumber,stim2) + alpha*(-1*EVs(trialNumber,stim2));                
                trialPEd(trialNumber) = (-1*EVs(trialNumber,stim2));
            end              
            trialEVs(trialNumber+1,:) = [EVs(trialNumber+1,stim1) EVs(trialNumber+1,stim2)]; 
            trialNumber = trialNumber + 1;
        end
    end
    fitData.trialEVs = trialEVs;
    fitData.trialPEd = trialPEd;
    fitData.trialPEp = trialPEp;
    fitData.EVs = EVs;
elseif strcmp(model,'modelFree')
    % Getting free parameters
    alpha = params(1); lambda = params(2);    
    nAvailableBlocks = length(b.sessionParams.blocksAvailable);
    trialEVs = 0.5*ones(nAvailableBlocks*b.taskStruct.nTrials,2);
    trialPEd = 0.5*ones(nAvailableBlocks*b.taskStruct.nTrials,1);
    trialPEp = 0.5*ones(nAvailableBlocks*b.taskStruct.nTrials,1);
    % Initializing EVs at zero for all stimuli (A,B,X,Y)
    EVs = 0.5*ones(nAvailableBlocks*b.taskStruct.nTrials,22);
    trialNumber = 1;
    for bI = b.sessionParams.blocksAvailable
        nTrials = length(b.blockStructCell{1,bI}.trialStructCell);
        if bI == max(b.sessionParams.blocksAvailable) % last block
            nTrials = nTrials - 1;
        end
        for tI = 1:nTrials
            trialStruct = b.blockStructCell{1,bI}.trialStructCell{1,tI};
            stim1 = trialStruct.stim1; stim2 = trialStruct.stim2;
            trialType = b.blockStructCell{1,bI}.trialTypes(tI);
            % Parsing previous EV values
            EVs(trialNumber+1,:) = EVs(trialNumber,:);
            if trialType == 1 || trialType == 4 % Reward
                % Updating value of CSp
                EVs(trialNumber+1,stim2) = EVs(trialNumber,stim2) + alpha*(1-EVs(trialNumber,stim2));
                trialPEp(trialNumber) = (1-EVs(trialNumber,stim2));
            elseif trialType == 2 || trialType == 3 % No reward
                % Updating value of CSp
                EVs(trialNumber+1,stim2) = EVs(trialNumber,stim2) + alpha*(0-EVs(trialNumber,stim2));                
                trialPEp(trialNumber) = (0-EVs(trialNumber,stim2));
            end              
            % Updating value of CSd
            EVs(trialNumber+1,stim1) = EVs(trialNumber,stim1) + alpha*(lambda*EVs(trialNumber+1,stim2)-EVs(trialNumber,stim1));
            trialPEd(trialNumber) = (lambda*EVs(trialNumber+1,stim2)-EVs(trialNumber,stim1));
            
            trialEVs(trialNumber+1,:) = [EVs(trialNumber+1,stim1) EVs(trialNumber+1,stim2)]; 
            trialNumber = trialNumber + 1;
        end
    end
    fitData.trialEVs = trialEVs;
    fitData.trialPEd = trialPEd;
    fitData.trialPEp = trialPEp;
    fitData.EVs = EVs;
end

end

% Transforming surrogate parameters to real parameters if needed
function newParams = parameterTransform(params,model)
    newParams = params;
end

