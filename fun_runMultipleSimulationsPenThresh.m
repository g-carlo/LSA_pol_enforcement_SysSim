% run multiple simulations

function fun_runMultipleSimulationsPenThresh(n_sim,Mno11,Mno21,threspen21,threspen22,threspen23)
% clear all;
% delete(gcp);
% myPool = parpool();
%n_sim = 100;        % number of repetitions

% paramVector = [6:12];
% 
% n_sim = length(paramVector);

% saveFileLabel = 'GAIN_100_incumbents_power_9dB_121_sites';
% saveFileLabel = 'MISBEH_16_Dec_p_0-1_0-01_L120_T25_4s-50';
saveFileLabel = '180410_KV4_';
% PenaltyFunction = 'linear';           

%mbeh_probs =   [0.3 0.2 0];
%penaltyThreshold = [0 50 100 150 200];

mb_v1 = Mno11';
mb_v2 = Mno21';
mb_v3 = zeros(1,1);

mbeh_probs = [mb_v1, mb_v2, mb_v3];

threspen_v1 = threspen21;
threspen_v2 = threspen22;
threspen_v3 = threspen23;
% threspen_v4 = [threspen24];
% threspen_v5 = [threspen25];

penaltyThreshold = threspen_v1: threspen_v2: threspen_v3;

n_mbeh_probs = size(mbeh_probs,1);                
n_penaltyThreshold = numel(penaltyThreshold);
% tic
penalties = [0 0 0];
parfor idx_sim = 1:n_sim
% for idx_sim = 1:n_sim
    k1 = min(n_mbeh_probs,ceil(idx_sim/(n_sim/n_mbeh_probs)));
    k2 = min(n_penaltyThreshold,ceil(idx_sim/(n_sim/n_penaltyThreshold)));
    %%% specify parameters
    paramStructToChange = cell(2,2);
    paramStructToChange{1,1} = 'conf.operator.misBehProbability';
    paramStructToChange{1,2} = mbeh_probs(k1,:);
    paramStructToChange{2,1} = 'conf.operator.initialPenalty';
    paramStructToChange{2,2} = penalties;
    paramStructToChange{3,1} = 'conf.L1_scheduling.penaltyThreshold';
    paramStructToChange{3,2} = penaltyThreshold(k2);
%     paramStructToChange{4,1} = 'conf.L1_scheduling.penaltyCalculation';
%     paramStructToChange{4,2} = PenaltyFunction;
    
    %%% 
    tic
%     penalties = ...
    runSimulation_v3_fnc_MBEHPenThresh(idx_sim,saveFileLabel,paramStructToChange);  % specify file label
    disp([num2str(idx_sim) ' out of ' num2str(n_sim)]);
%     disp(penaltyThreshold(k2));
    %disp(penalties);
    toc       
end
% toc
%% 
processResults_fncPenThresh(n_sim,saveFileLabel,penaltyThreshold);
%

