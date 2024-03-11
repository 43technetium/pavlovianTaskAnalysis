import numpy as np
import pandas as pd
import statsmodels.api as sm
import os
from sklearn.preprocessing import OneHotEncoder
import statsmodels.formula.api as smf
import scipy.stats as sst
import pickle
from sklearn.linear_model import PoissonRegressor
import matplotlib.pyplot as plt
from custom_loadmat import custom_loadmat
from scipy.stats import f_oneway
import json

def fit_model(x,y):    
    try:            
        df = pd.DataFrame()
        poisson_results = sm.GLM(y, x, family=sm.families.Poisson(), missing='drop').fit()
        df['BB_COUNT'] = np.squeeze(y)
        df['BB_LAMBDA'] = poisson_results.mu
        df['AUX_OLS_DEP'] = df.apply(lambda x: ((x['BB_COUNT'] - x['BB_LAMBDA'])**2 - x['BB_LAMBDA']) / x['BB_LAMBDA'], axis=1)
        ols_expr = """AUX_OLS_DEP ~ BB_LAMBDA - 1"""
        aux_olsr_results = smf.ols(ols_expr, df).fit()
        print(aux_olsr_results.params.iloc[0])
        res = sm.GLM(y, x, family=sm.families.NegativeBinomial(alpha=aux_olsr_results.params.iloc[0])).fit()
        return res
    except: # Units that did not converge
        return np.nan

def fit_model_ols(x,y):        
    res = sm.OLS(y,x).fit()
    return res

def get_counts_from_spikes(spikes,windowStart,windowEnd):
    nTrials = len(spikes)
    y = np.zeros((nTrials,1))
    for tI in np.arange(nTrials):
        trial_spikes = np.array(spikes[tI])
        y[tI] = np.sum(np.logical_and(trial_spikes>windowStart,trial_spikes<=windowEnd))
    return y

sessions = ['P49CS','P51CS_1','P53CS','P54CS','P58CS','P60CS','P61CS','P62CS','P63CS','P70CS','P71CS','P76CS']

# Which blocks are valid for each session (out of 1,2,3,4)
blocksAvailable = [[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3]]

blockTrials = np.array([np.arange(0,24),np.arange(24,48),np.arange(48,72),np.arange(72,96)])
basefolder = 'C:/Users/Tomas/Documents/PhD/OLab/pavlovianLearning/pavlovianConditioningTask/patientData/'

# Which regression model to test
#design = 'csp_presumed_identity_EV_CSd_MB_csd_identity'
#design = 'csp_presumed_identity_continuous_EV_CSd_MB_csd_identity'
design = 'csp_presumed_identity_EV_CSd_MB'
#design = 'EV_CSd_MB_csd_identity'
#design = 'csp_identity_csp_EV_uPE1_button'
#design = 'outcome_spe2_rpe'

# Which time alignment to use for spikes
reference = 'trial'
#reference = 'decision'
#reference = 'outcome'

## Time window arrangement
# windowing = 'standard_decision'; 
#windowing = 'standard_outcome'
# windowing = 'pre_outcome';
windowing = 'standard_trial'
#windowing = 'csp_window'
# windowing = 'csp_windowed';

for sI in np.arange(len(sessions)):
    # Loading data
    session = sessions[sI]
    behavior_folder = os.path.join(basefolder, 'allBehavior')
    behavior_data = custom_loadmat(os.path.join(behavior_folder,'sessionBehavior_'+session+'.mat'))
    sessionFolder = os.path.join(basefolder, session)
    data = custom_loadmat(os.path.join(sessionFolder, 'sessionData.mat'))
    unitCell = data['sessionData']['neuralData']['unitCell']
    nTrials = len(unitCell[0]['trialReferencedSpikes'])
    nUnits = len(unitCell)
    # Getting trials available in this session
    trialsAvailable = np.squeeze(blockTrials[blocksAvailable[sI], :].reshape(-1, 1))
    # Setting up windowing 
    dt = 0.00001
    if windowing == 'standard_outcome':
        # Bandit-related time windows
        binSize = 2
        windowSpacing = 0.01
        # For post-reward windows
        period = [0.25, 2.25]
        windowStarts = np.arange(period[0],period[1]-binSize+dt,windowSpacing)
        windowEnds = windowStarts + binSize
        nBins = len(windowStarts)
    elif windowing == 'pre_outcome':
        # Bandit-related time windows
        binSize = 3
        windowSpacing = 0.01
        # For post-reward windows
        period = [-3, 0]
        windowStarts = np.arange(period[0],period[1]-binSize+dt,windowSpacing)
        windowEnds = windowStarts + binSize
        nBins = len(windowStarts)
    elif windowing == 'standard_trial':
        # Bandit-related time windows
        binSize = 2.75
        windowSpacing = 0.01
        # For post-reward windows
        period = [0.25, 3]
        windowStarts = np.arange(period[0],period[1]-binSize+dt,windowSpacing)
        windowEnds = windowStarts + binSize
        nBins = len(windowStarts)
    elif windowing == 'csp_windowed':
        # Bandit-related time windows
        binSize = 0.5
        windowSpacing = 0.05
        # For post-reward windows
        period = [3.75, 7.5]
        windowStarts = np.arange(period[0],period[1]-binSize+dt,windowSpacing)
        windowEnds = windowStarts + binSize
        nBins = len(windowStarts)
    elif windowing == 'csp_window':                                    
        # Bandit-related time windows
        binSize = 2.75
        windowSpacing = 0.01
        # For post-reward windows
        period = [4.25, 7]
        windowStarts = np.arange(period[0],period[1]-binSize+dt,windowSpacing)
        windowEnds = windowStarts + binSize
        nBins = len(windowStarts)

    ## Creating regressor vectors
    EV_CSd_MB = behavior_data['stimEVs_MB'][:,0]
    EV_CSp_MB = behavior_data['stimEVs_MB'][:,1]
    
    # Getting CSp identity
    CSp_options = np.unique(behavior_data['stim2Vec'])
    CSd_options = np.unique(behavior_data['stim1Vec'])    
    enc = OneHotEncoder(handle_unknown='ignore')
    CSd_onehot = enc.fit_transform(behavior_data['stim1Vec'].reshape(-1,1)).toarray()   
    CSd_onehot = np.delete(CSd_onehot, 0, 1)  # delete first column of onehot (redundant for GLM)
    # Unsigned prediction error
    uPE1 = np.squeeze(behavior_data['SPE1_vec'])
    uPE2 = np.squeeze(behavior_data['SPE2_vec'])
    RPE = np.squeeze(behavior_data['RPE'])
    outcome = np.copy(behavior_data['trialTypes'])
    outcome[outcome==2] = 0
    outcome[outcome==3] = 0
    outcome[outcome==4] = 1
    button_press = np.squeeze(behavior_data['responseVec'])-1

    # Constructing CSd matrix
    sessionBlocks = blocksAvailable[sI]
    csd_identity = np.zeros((nTrials,len(sessionBlocks)))
    for bI in sessionBlocks:
        block_unique = np.unique(behavior_data['stim1Vec'][blockTrials[bI]])
        block_ID = np.array((block_unique[0]==behavior_data['stim1Vec'])).astype(int)-np.array((block_unique[1]==behavior_data['stim1Vec'])).astype(int)
        csd_identity[:,bI] = block_ID        
        
    csp_identity = (behavior_data['stim2Vec'] == CSp_options[0]).astype(int)
    # The most likely proximal identity
    trialTypes = behavior_data['trialTypes_flat']
    nTrials = len(trialTypes)
    csp_presumed_identity = np.copy(csp_identity)
    csp_presumed_identity[np.logical_or(trialTypes==3,trialTypes==4)] = 1-csp_presumed_identity[np.logical_or(trialTypes==3,trialTypes==4)]
    # Building matrix of continuous presumed identities
    csp_presumed_identity_continuous = np.zeros((nTrials,1))
    for tI in np.arange(nTrials):
        # Getting transition probability matrix for this trial
        T = behavior_data['T_cell'][tI]
        trial_type = behavior_data['trialTypes_flat'][tI]
        if trial_type == 1 or trial_type == 3: # CSd+ / leaving from A->X or A->Y
            T_row = 0
        elif trial_type == 2 or trial_type == 4: # CSd- / leaving from B->Y or B->X
            T_row = 1
        # Probability of moving onto first X identity will be tracked (X's identity reverses between blocks)
        if tI in blockTrials[0] or tI in blockTrials[2]:
            T_column = 2
        elif tI in blockTrials[1] or tI in blockTrials[3]:
            T_column = 3
        csp_presumed_identity_continuous[tI] =T[T_row,T_column]
    # Getting design matrix (with constant column)
    if design == 'csp_presumed_identity_EV_CSd_MB_csd_identity':        
        x = np.concatenate((np.ones((nTrials,1)), csp_presumed_identity.reshape(-1,1), EV_CSd_MB.reshape(-1,1), CSd_onehot),axis=1)        
    elif design == 'csp_presumed_identity_continuous_EV_CSd_MB_csd_identity':        
        x = np.concatenate((np.ones((nTrials,1)), csp_presumed_identity_continuous, EV_CSd_MB.reshape(-1,1), CSd_onehot),axis=1)        
    elif design == 'csp_presumed_identity_EV_CSd_MB':
        x = np.concatenate((np.ones((nTrials,1)), csp_presumed_identity.reshape(-1,1), EV_CSd_MB.reshape(-1,1)),axis=1)
    elif design == 'csp_presumed_identity_continuous_EV_CSd_MB':
        x = np.concatenate((np.ones((nTrials,1)), csp_presumed_identity_continuous.reshape(-1,1), EV_CSd_MB.reshape(-1,1)),axis=1)
    elif design == 'EV_CSd_MB_csd_identity':
        x = np.concatenate((np.ones((nTrials,1)), EV_CSd_MB.reshape(-1,1), CSd_onehot),axis=1)    
    elif design ==  'csp_identity_csp_EV_uPE1_button':
        x = np.concatenate((np.ones((nTrials,1)), csp_identity.reshape(-1,1), EV_CSp_MB.reshape(-1,1), uPE1.reshape(-1,1), button_press.reshape(-1,1)),axis=1)    
    elif design == 'outcome_spe2_rpe':
        x = np.concatenate((np.ones((nTrials,1)), outcome.reshape(-1,1), uPE2.reshape(-1,1), RPE.reshape(-1,1)),axis=1)

    CSd_trials_all = []
    
    for stimI in np.arange(len(CSd_options)):
        CSd_trials_all.append(np.where(np.array(behavior_data['stim1Vec']==CSd_options[stimI]))[0])    
    
    # Enforce only trials from available blocks are selected
    x = x[trialsAvailable,:]
    n_regressors = x.shape[1]
    mdlCell=[]
    all_unit_info = []
    t_values = np.empty((0,n_regressors))
    p_values = np.empty((0,n_regressors))
    anova_p = np.zeros((nUnits,nBins))
    for uI in np.arange(nUnits):
        print('Session number: ', str(sI), ' / Unit number: ', str(uI))
        unit_info = data['sessionData']['neuralData']['unitCell'][uI]['unitInfo']
        all_unit_info.append(unit_info)
        # Looping over time bins        
        if reference == 'outcome':
            spikes = data['sessionData']['neuralData']['unitCell'][uI]['outcomeReferencedSpikes']            
        elif reference == 'trial':
            spikes = data['sessionData']['neuralData']['unitCell'][uI]['trialReferencedSpikes']            
        elif reference == 'decision':
            spikes = data['sessionData']['neuralData']['unitCell'][uI]['decisionReferencedSpikes']                    
        bin_results = []        
        for bI in np.arange(nBins):
            y = get_counts_from_spikes(spikes,windowStarts[bI],windowEnds[bI])       
            final_y = y[trialsAvailable]
            y_z = sst.zscore(y[trialsAvailable])
            mdl = fit_model(x,final_y)
            anova_res = f_oneway(y_z[CSd_trials_all[0]], y_z[CSd_trials_all[1]], y_z[CSd_trials_all[2]], y_z[CSd_trials_all[3]], y_z[CSd_trials_all[4]], y_z[CSd_trials_all[5]], y_z[CSd_trials_all[6]], y_z[CSd_trials_all[7]],axis=0)
            anova_p[uI,bI] = anova_res.pvalue
                  
            if type(mdl) is float or type(mdl) is np.float64: # Skipping nans
                filler = np.full([1,n_regressors], np.nan)
                p_values = np.vstack((p_values,filler))    
                t_values = np.vstack((t_values,filler))                     
            else:
                p_values = np.vstack((p_values,mdl.pvalues.reshape(1,-1)))    
                t_values = np.concatenate((t_values,mdl.tvalues.reshape(1,-1))) 
            bin_results.append(mdl)

        mdlCell.append(bin_results)
    
    # Creating a dictionary to store the values
    forwardData = {
        'mdlCell': np.squeeze(mdlCell),
        'all_unit_info': all_unit_info,
        'anova_p': anova_p,
        'reference': reference,
        'binSize': binSize,
        'windowSpacing': windowSpacing,
        'period': period,
        'windowStarts': windowStarts,
        'windowEnds': windowEnds,
        'nBins': nBins
    }

    # Generating the save folder path
    savefolder = os.path.join(sessionFolder, 'forward', design, reference, windowing)

    # Make folder if it doesn't exist
    if not os.path.exists(savefolder):
        os.makedirs(savefolder)
    
    # Saving the data using pickle in Python (similar to MATLAB's save function)
    with open(os.path.join(savefolder, 'GLMResults_linear.pkl'), 'wb') as file:
        pickle.dump(forwardData, file)
    
    # Saving data for plotting firing rates
    with open(os.path.join(savefolder, 'p_values.npy') , 'wb') as f:
        np.save(f,p_values)    
    with open(os.path.join(savefolder, 't_values.npy') , 'wb') as f:
        np.save(f,t_values)