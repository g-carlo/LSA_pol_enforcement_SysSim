function result = mainSim_v3_fnc_MBEH(conf,carrierArray,cellTypeArray,carrierTypeArray,type2TypeArray,sim_idx,sim_label)

%%% SON SIMULATOR: main simulation file
%
%   Author: Carlo Galiotto
%
%   Date:   July 2010, last modification: March, 2015
%
%   Description: this file is the main simulation file, where the main
%   variables and the simulation loop are defined.

%

global macroInd;
global picoInd;

varSaveVars = 1;

if varSaveVars  == 1
    varfnamevars = 'lsaSimulatorVars';
end
    
numOfCarriers = conf.general.numOfCarriers;
isEiRP = 0; 
%%% Simulator profiling
if conf.results.profiling
    profile on;
end

%random generators
% seedFileName = strcat(conf.storage.inputFolderPath,conf.storage.seedFileName);
% seedNumber = getSeedNumber(conf.randomness.randGeneratorState,isEiRP,...
%     seedFileName);
% RandStream.setGlobalStream(RandStream('mt19937ar','seed',seedNumber));
% RandStream.setGlobalStream(RandStream('mt19937ar','seed',0));


%% initialization of cells

% tic;

% % figure
% % hold on

if conf.general.SLS_activated
    
    [cls,conf] = initializeCells_v2(conf,cellTypeArray,type2TypeArray,...
        carrierTypeArray, [carrierArray.numOfAvailablePRB]);
    
else
    
    [cls,conf] = initializeCells_noSLS(conf,cellTypeArray,type2TypeArray,...
        carrierTypeArray, [carrierArray.numOfAvailablePRB],conf.deployment.N_sites_extedndedScenario);
    
end

%     if varSaveVars  == 1
%         varcellpos = [cls.pos];
%         varcellpos_x = [varcellpos.x];
%         varcellpos_y = [varcellpos.y];
%         save(varfnamevars,'varcellpos_x','-append');
%         save(varfnamevars,'varcellpos_y','-append');
%     end
%   
[carrierArray, cellsPerCarrier] = setCellsPerCarrier_v1(...
    carrierArray, [cls.complexPosition], ...
    conf.general.randomAssigCellsPerCarrier, conf.general.cellRatio,[cls.networkOP_ID]);

%load cells' power from file

varLoadTxFromFile = sum([carrierArray(:).loadPtxData]);
if varLoadTxFromFile>0
    filename = strcat(conf.storage.inputFolderPath,'cellsPower.dat');
    if exist(filename, 'file')
        fid = fopen(filename);
        A=fscanf(fid, '%d %f\n', [2 inf]);
        fclose(fid);

        tmp_err=0;

        counter=0;
        for ii=1:numOfCarriers
            if carrierArray(ii).loadPtxData==1
                for j=1:length(carrierArray(ii).cellsInCarrier)
                    counter=counter+1;
                    if A(1,counter)~=ii
                        tmp_err=1;
                        error('Structure of inputData/cellsPtx.dat is not OK!');
                    end
                    cellId = carrierArray(ii).cellsInCarrier(j);
                    cls(cellId).raVars(ii).maxTxPw = A(2,counter);                    
                end
            end
        end
    end
end

% tfin = toc;
% disp(strcat({'Initialize cells: '},num2str(tfin),{' s'}));


%% initialize event handler

% generate database and incuments 
paramVar.deployment  	= conf.deployment;
paramVar.database       = conf.database;
paramVar.incumbent      = conf.incumbent;
database_var            = LSA_dataBase_Cls (conf.database.numOfLSAChannels,conf.deployment.numOfIncumbents,paramVar);   
            
areaLimVector = [conf.deployment.incumbentEastBound conf.deployment.incumbentWestBound conf.deployment.incumbentSouthBound conf.deployment.incumbentNorthBound];

%% initialize LSA controller

paramLSAcontroller.general          = conf.general;
paramLSAcontroller.deployment       = conf.deployment;
paramLSAcontroller.incumbent        = conf.incumbent;
paramLSAcontroller.PHY              = conf.PHY;
paramLSAcontroller.channel          = conf.channel;
paramLSAcontroller.L1_scheduling    = conf.L1_scheduling;
paramLSAcontroller.database         = conf.database;

% LSA_controller = LSA_controller_Cls(conf.database.numOfLSAChannels,conf.deployment.numOfIncumbents,conf.deployment.numOfNetworkOperators,paramLSAcontroller,conf.carrier.frequency);         
LSA_controller = LSA_controller_Cls(conf.database.numOfLSAChannels,conf.deployment.numOfIncumbents,conf.deployment.numOfNetworkOperators,paramLSAcontroller,conf.carrier.frequency);         

% generate incumbents
paramIncumbent.incumbent = conf.incumbent;
paramIncumbent.PHY       = conf.PHY;

for nn = 1:conf.deployment.numOfIncumbents
    incumb(nn) = incumbent_Cls(nn,paramIncumbent,areaLimVector);                 
end

netOperatorEventGen = networkOperator_eventGenerator_Cls(conf.deployment.numOfNetworkOperators,conf.traffic.lam_newDem)  ;

% initialize database
database_var.initializeDataBase(incumb,netOperatorEventGen);
LSA_controller.initializeLSAcontroller(incumb,netOperatorEventGen,conf.operator.initialPenalty)

for nn = 1:conf.deployment.numOfIncumbents

    eventQueue(nn) =  eventSLS_Cls('incumbent',incumb(nn).activityCnt,'active',nn);
    incumb(nn).initialize(conf.incumbent.radius(nn),database_var,eventQueue(nn));
    
end

eventQueue(conf.deployment.numOfIncumbents+1) =  eventSLS_Cls('operator',netOperatorEventGen.newDemandCnt,'new_demand',1);



%% initialize network generator events
netOperatorEventGen.initialize(eventQueue(conf.deployment.numOfIncumbents+1),[cls.networkOP_ID],conf.database.numOfLSAChannels);

pointNext = 1;

%% initialize queue to store results

[nextTime, minIdx] = getNextEventIdx(eventQueue);
% nextEvent  = eventSLS_Cls(ev_ty, ev_ti, ev_ac, ev_ID);
nextEvent_activeIncumbents = false(1,length(eventQueue));
nextEvent_inactiveIncumbents = false(1,length(eventQueue));
nextEvent_rescheduleOperator = false(1,length(eventQueue));


% resultsLogQueue     = [];
res.out = 0;
res.SINR_vec = 0;
res.SNR_vec = 0;
res.SNR_IM_vec = 0;
res.assignedSpectrum = zeros(1,conf.deployment.numOfNetworkOperators);
res.isMisbehaving = false(1,conf.deployment.numOfNetworkOperators);
res.incumbentActivity = false(1,conf.deployment.numOfIncumbents);
res.overallIncumbentActivity = false(1,conf.deployment.numOfIncumbents);
res.detection.activeBSs_idx = 0;
res.detection.detectedBSs_idx = 0;
nextLogElement      = 0;



%% running event handler

while nextTime < conf.general.simulatedTime 
    
    
    for evNum = 1:length(minIdx)
       
        idxEv = minIdx(evNum);
        
        switch [eventQueue(idxEv).type '_' eventQueue(idxEv).activity]
            
            case 'incumbent_active'
                
                nextEvent_activeIncumbents(idxEv) = true;
                
            case 'incumbent_inactive'
                
                nextEvent_inactiveIncumbents(idxEv) = true;
                
            case 'operator_new_demand'
                
                nextEvent_rescheduleOperator(idxEv) = true;
                 
        end
        
        
    end
    
    for id_Inc = find(nextEvent_activeIncumbents | nextEvent_inactiveIncumbents);
        
        incumb(id_Inc).update();
        
    end
    
%     if any(nextEvent_activeIncumbents)
%        disp('stop' );
%     end
    
    database_var.updateDatabase();
    [cellsPerCarrier, opPerLSAchannel] = LSA_controller.runL1_algorithm([cls.complexPosition],[cls.orientation],conf.cellType.height(1),...
    conf.antenna,conf.carrier.cellType.antennaGain(:,1).',[carrierArray.electricalDowntilt],...
    conf.carrier.cellType.maxTxPw(:,1).' ,cellsPerCarrier,conf.database.numOfLSAChannels,nextTime);
    cellsPerCarrier_no_misbehaviour = cellsPerCarrier;
    if nextLogElement > 0
            resultsLogQueue(nextLogElement).setFinishingTime(nextTime-1);
    end           
            
    %%%% MISBEHAVIOUR AND DETECTION 

    
    if any(nextEvent_activeIncumbents)
        
        [cellsPerCarrier, isMisbehaving] = performMisbehaviour(conf,netOperatorEventGen,cellsPerCarrier,opPerLSAchannel);  
%                  out = runSnapshotSimulation(conf,carrierArray,cellTypeArray,cls,cellsPerCarrier,carrierTypeArray,numOfCarriers);
        out = runSnapshotSimulationQuick(conf,cellsPerCarrier,netOperatorEventGen);
        idx_TX_BS = cellsPerCarrier( conf.deployment.numOfNetworkOperators,netOperatorEventGen.cell_idx_per_NO(opPerLSAchannel).list) ;
        idx_enabled_BSs = cellsPerCarrier_no_misbehaviour( conf.deployment.numOfNetworkOperators,netOperatorEventGen.cell_idx_per_NO(opPerLSAchannel).list );
        N_enabled_BSs = sum( idx_enabled_BSs );
        posInc = compute_incunbent_positions( 5,areaLimVector);     
%         [Est_Loc_active_BS, Pd_table, Pfa_table, Est_energy_active_BS, MT] = Core_function_correct( N_enabled_BSs ,...
%             conf.detection.num_of_sensors, [cls(1:57).complexPosition], posInc, find(idx_TX_BS),nextLogElement);
%         [Est_Loc_active_BS, Pd_table, Pfa_table, Est_energy_active_BS, MT] = Core_function_simulate( N_enabled_BSs ,...
%             conf.detection.num_of_sensors, [cls(1:57).complexPosition], posInc, find(idx_TX_BS),nextLogElement);
        [Est_Loc_active_BS, Pd_table, Pfa_table] = Core_function_dumb(isMisbehaving);        
        
%         if isMisbehaving
%             disp('is misbehaving' );
%         end
        
    %             
        % log the detected misbehaviour
        LSA_controller.logMisbehaviour(Est_Loc_active_BS,idx_enabled_BSs,Pd_table,Pfa_table);
        n_cells_per_operator = cellTypeArray(1).numOfCellsPerOperator;
        incumbentActivity =  false(1,conf.deployment.numOfIncumbents);
        for  id_Inc = find(nextEvent_activeIncumbents );
            [SINR_incumbent, SNR_incumbent, SNR_incumbent_with_interference_margin] = incumb(id_Inc).computeSINR([cls(1:n_cells_per_operator).complexPosition],[cls(1:n_cells_per_operator).orientation],conf.cellType.height(1),...
            conf.antenna,conf.carrier.cellType.antennaGain(:,1).',[carrierArray.electricalDowntilt],...
            conf.carrier.cellType.maxTxPw(:,1).', idx_enabled_BSs);
            if strcmp( incumb(id_Inc).activityState,'active' )
                incumbentActivity(id_Inc) =  true;
            end
        end
        overallIncumbentActivity =  false(1,conf.deployment.numOfIncumbents);
        for id_Inc = 1:conf.deployment.numOfIncumbents
           if strcmp( incumb(id_Inc).activityState,'active' )
                overallIncumbentActivity(id_Inc) =  true;
            end 
        end
        activeBSs_idx = find(idx_TX_BS);
        detectedBSs_idx = Est_Loc_active_BS;
        assignedSpectrum = zeros(1,conf.deployment.numOfNetworkOperators);
        assignedSpectrum(opPerLSAchannel) = N_enabled_BSs/57;
        isMisbeh = false(1,conf.deployment.numOfNetworkOperators);
        isMisbeh(opPerLSAchannel) = isMisbehaving;

    else
        out = runSnapshotSimulationQuick(conf,cellsPerCarrier,netOperatorEventGen);
        idx_enabled_BSs = cellsPerCarrier_no_misbehaviour( conf.deployment.numOfNetworkOperators,netOperatorEventGen.cell_idx_per_NO(opPerLSAchannel).list );
        N_enabled_BSs = sum( idx_enabled_BSs );
        SINR_incumbent = 0;
        SNR_incumbent  = 0;
        SNR_incumbent_with_interference_margin = 0;
        activeBSs_idx = 0;
        assignedSpectrum = zeros(1,conf.deployment.numOfNetworkOperators);
        assignedSpectrum(opPerLSAchannel) = 1;
        detectedBSs_idx = 0;
        isMisbeh = false(1,conf.deployment.numOfNetworkOperators);
        incumbentActivity =  false(1,conf.deployment.numOfIncumbents);
        overallIncumbentActivity =  false(1,conf.deployment.numOfIncumbents);
%             res.detection.detectedBSs_idx = 0;
    end
    
%     cellsPerCarrier = performMisbehaviour(conf,netOperatorEventGen,cellsPerCarrier,opPerLSAchannel);   
% %             out = runSnapshotSimulation(conf,carrierArray,cellTypeArray,cls,cellsPerCarrier,carrierTypeArray,numOfCarriers);
%     out = runSnapshotSimulationQuick(conf,cellsPerCarrier,netOperatorEventGen);
%     % detectMisbehaviour here
%     idx_TX_BS = cellsPerCarrier( 3,netOperatorEventGen.cell_idx_per_NO(opPerLSAchannel).list) ;
%     idx_enabled_BSs = cellsPerCarrier_no_misbehaviour( 3,netOperatorEventGen.cell_idx_per_NO(opPerLSAchannel).list );
%     N_enabled_BSs = sum( idx_enabled_BSs );
%     posInc = compute_incunbent_positions( 5,areaLimVector);     
%     [Est_Loc_active_BS, Pd_table, Pfa_table, Est_energy_active_BS, MT] = Core_function_correct( N_enabled_BSs ,...
%         conf.detection.num_of_sensors, [cls(1:57).complexPosition], posInc, find(idx_TX_BS),nextLogElement);
    %%%% END MISBEHAVIOUR AND DETECTION 
            
    %%%% NO MISBEHAVIOUR AND DETECTION 
%             out = runSnapshotSimulation(conf,carrierArray,cellTypeArray,cls,cellsPerCarrier,carrierTypeArray,numOfCarriers);
%     out = runSnapshotSimulationQuick(conf,cellsPerCarrier,netOperatorEventGen);
%     
%     if any(nextEvent_activeIncumbents)
%         
%         
%         idx_TX_BS = cellsPerCarrier( 3,netOperatorEventGen.cell_idx_per_NO(opPerLSAchannel).list) ;
%         idx_enabled_BSs = cellsPerCarrier_no_misbehaviour( 3,netOperatorEventGen.cell_idx_per_NO(opPerLSAchannel).list );
%         Est_Loc_active_BS = [];
%         Pd_table = 1;
%         Pfa_table = 0;
%         %%%% END  NO MISBEHAVIOUR AND DETECTION 
%     %             
%         % log the detected misbehaviour
%         LSA_controller.logMisbehaviour(Est_Loc_active_BS,idx_enabled_BSs,Pd_table,Pfa_table);
%         n_cells_per_operator = cellTypeArray(1).numOfCellsPerOperator;
%         for  id_Inc = find(nextEvent_activeIncumbents );
%             [SINR_incumbent, SNR_incumbent, SNR_incumbent_with_interference_margin] = incumb(id_Inc).computeSINR([cls(1:n_cells_per_operator).complexPosition],[cls(1:n_cells_per_operator).orientation],conf.cellType.height(1),...
%             conf.antenna,conf.carrier.cellType.antennaGain(:,1).',[carrierArray.electricalDowntilt],...
%             conf.carrier.cellType.maxTxPw(:,1).', idx_enabled_BSs);
%         end
%         activeBSs_idx = find(idx_TX_BS);
%         detectedBSs_idx = Est_Loc_active_BS;
%     else
%         SINR_incumbent = 0;
%         SNR_incumbent  = 0;
%         SNR_incumbent_with_interference_margin = 0;
%         activeBSs_idx = 0;
%         detectedBSs_idx = 0;
% %             res.detection.detectedBSs_idx = 0;
%     end
    
    %%%% END  NO MISBEHAVIOUR AND DETECTION 

%             [ estimated_TX_BS_idx, P_d, P_fa, estimated_TX_pow ] = detectMisbehaviour( N_enabled_BSs ,...
%                 conf.detection.num_of_sensors, [cls(1:n_cells_per_operator).cls.complexPosition], incumb(1).position, idx_TX_BS);
%             resultsLogQueue(nextLogElement) = snapshotResultsLog_Cls(nextEvent.time);         
%             resultsLogQueue(nextLogElement).logData(out);
    nextLogElement = nextLogElement + 1;
    res.out = out;
    res.SINR_vec = SINR_incumbent;
    res.SNR_vec = SNR_incumbent;
    res.SNR_IM_vec = SNR_incumbent_with_interference_margin;
    res.detection.activeBSs_idx = activeBSs_idx;
    res.assignedSpectrum = assignedSpectrum;
    res.detection.detectedBSs_idx = detectedBSs_idx;
    res.isMisbehaving = isMisbeh;
    res.incumbentActivity = incumbentActivity;
    res.overallIncumbentActivity = overallIncumbentActivity;
%     res.assignedSpectrum = 

    resultsLogQueue(nextLogElement) = snapshotResultsLog_Cls(nextTime,res); 

    pointNext = pointNext +1;
    
    

    
%    


%     disp('------');
%     disp('Another event has been processed');
%     disp(['Time ' num2str(nextTime)]);
    [nextTime, minIdx] = getNextEventIdx(eventQueue);
    nextEvent_activeIncumbents = false(1,length(eventQueue));
    nextEvent_inactiveIncumbents = false(1,length(eventQueue));
    nextEvent_rescheduleOperator = false(1,length(eventQueue));
    
    
end



disp('End and saving data');
simFilename = [ sim_label '_n_' num2str(sim_idx)];
save(simFilename);

result = LSA_controller.schedData.penalty;

% % % carrierTypeArray = carrierCellTypeArray;
% % % type2TypeArray = cellType2cellTypeArray;
% % % 
% % % % function mainSim_v2(conf,carrierArray,cellTypeArray,carrierTypeArray,...
% % % %         type2TypeArray)
% % % % 
% % % %     global macroInd;
% % % %     global picoInd;
% % %     
% % % %     varSaveVars = 1;
% % % %     
% % % %     if varSaveVars  == 1
% % % %         varfnamevars = 'lsaSimulatorVars';
% % % %     end
% % % %     
% % % numOfCarriers = conf.general.numOfCarriers;
% % % isEiRP = 0;
% % % %%% Simulator profiling
% % % if conf.results.profiling
% % %     profile on;
% % % end
% % % 
% % % %random generators
% % % seedFileName = strcat(conf.storage.inputFolderPath,conf.storage.seedFileName);
% % % seedNumber = getSeedNumber(conf.randomness.randGeneratorState,isEiRP,...
% % %     seedFileName);
% % % RandStream.setGlobalStream(RandStream('mt19937ar','seed',seedNumber));

%     if varSaveVars  == 1
%         save(varfnamevars,'seedNumber');
%     end
%     

