function out = runSnapshotSimulation(conf,carrierArray,cellTypeArray,cls,...
    cellsPerCarrier,carrierTypeArray,numOfCarriers)

global macroInd;
global picoInd;
global smallInd;

% number of cell types: macro, pico, small cells, ...
% numOfCellTypes = conf.implementation.numOfCellTypes; 

% % number of carriers that will be considered in the current simulation.
% numOfCarriers = conf.general.numOfCarriers;


%% initialization of users
tic;
ues =initializeUsersMultipleFreq_v3(cls,cellTypeArray,...
    conf.general,conf.deployment,conf.user.height);

% if ~isempty(varPoints)
%     plot(real(varPoints),imag(varPoints),'r','lineWidth',3);
% end

% fill(real(varPoints),imag(varPoints),'y');

% hold off

%     if varSaveVars  == 1
%     
%         varuepos = [ues.pos];
%         varuepos_x = [varuepos.x];
%         varuepos_y = [varuepos.y];
% 
%         save(varfnamevars,'varuepos_x','-append');
%         save(varfnamevars,'varuepos_y','-append');
%     end
%     
varAssigUEcarrier = 0;    
if varAssigUEcarrier == 1
    mtxUserVsFreq = assignFreqToUes(length(ues),carrierArray,...
        conf.randomness.usersRandomization,seedNumber);
end

% if isEiRP==1
%     for p=1:numOfCarriers
%         nu=0;
%         for ii=xmin:dx:xmax
%             for j=ymin:dy:ymax
%                 nu=nu+1;
%                 userIndex = carrierArray(p).usersInCarrier(nu);
%                 ues(userIndex).setPosition(ii,j,ues(userIndex).pos.z);
%             end
%         end
%         %N=configV(p).general.numOfUsers;
%         N = length(carrierArray(p).usersInCarrier);
%         for ii=nu+1:N
%             %uesV{p}(ii).setPosition(xmax+3*dx,ymax+3*dy,uesV{p}(i).pos.z);
%             userIndex = carrierArray(p).usersInCarrier(ii);
%             ues(userIndex).setPosition(xmax+3*dx,ymax+3*dy,ues(userIndex).pos.z);
%         end
%     end
% end

tfin = toc;
disp(strcat({'Initialize users: '},num2str(tfin),{' s'}));

% if isEiRP==1
%     seedNumber = sum(100*clock);
%     RandStream.setGlobalStream(RandStream('mt19937ar','seed',seedNumber));
% end;

    tic;      

    % compute channel coefficients
    varDistCellUeNoWA = generateChannelFnc_v2(cls,ues,conf.antenna,conf.channel,...
        conf.user.height,carrierArray,cellTypeArray,...
        cellsPerCarrier,carrierTypeArray,conf.general.useWA,...
        conf.deployment.numOfSites,conf.deployment.numOfCellsPerSite,conf.deployment.numOfNetworkOperators);

    %initializeFastFadingParameters();
    tfin = toc;
    disp(strcat({'Channel matrix computation: '},num2str(tfin),{' s'}));

    totNumOfUes = length(ues);
    totNumOfCells = length(cls);


%%  Attaching users to the cells
% ---------------------------------------------------------------------
tic;
%%% random number generator
if conf.randomness.randGeneratorStateDecoupling
    seedNumber = sum(100*clock);
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',seedNumber));
end

data = data_Cls();

[data.userVsCell_allocationMtx, data.interfererIndexMatrix] = ...
    attachUsersFnc_v2(cls, ues, carrierArray, conf.general,...
    conf.channel.fastfadingFlag,varDistCellUeNoWA);

tfin = toc;
disp(strcat({'Attaching users to cells: '},num2str(tfin),{' s'}));

%     if varSaveVars  == 1    
%         varcellue = [ues.servingCellIndex];
% 
%         save(varfnamevars,'varcellue','-append');
%     end
% ---------------------------------------------------------------------

%% initalization of traffic buffer
% ---------------------------------------------------------------------
tic;
% for each UE initializes the following UE properties related to
% traffic: 1) lambda corresponding to the serving cell (this value is used
% if the traffic model is the burst traffic, 2) ueArrivalTimer, 
% 3) initVarTraffic, 4) ueTXtimer, 5) pckCnt 6) transmittedDataCnt
for ii=1:totNumOfUes
    ues(ii).initializeTraffic(conf.traffic,cls);
end

tfin = toc;
disp(strcat({'Initalize traffic buffer: '},num2str(tfin),{' s'}));
% ---------------------------------------------------------------------

%% Allocating variables for saving long run simulation data

TTIblockDim = conf.results.TTIblockDataLength;
varSimDataArray = longRunSimData_Cls(conf.results,0,0,0,0,0);

for ii=1:numOfCarriers

    % values computed for each UE

    % in each iteration is computed the average capacity (shannon
    % formula) of each UE, as the quotient between: i)sum of capacity 
    % of UE over all assigned PRBs and ii)the number of assigned PRBs 
    % to that UE. This average is added to the same metric computed 
    % in previous iterations. The sum is the value of <ueCumSpectralEfficiency> 
    % in each carrier.
    carrierArray(ii).ueCumSpectralEfficiency = zeros(1,totNumOfUes);

    % For each UE, it is the number of iterations where at least one 
    % PRB has been assigned to that UE.
    carrierArray(ii).ueSpectralEfficiency_cnt = zeros(1,totNumOfUes);       

    varSimDataArray(ii) = longRunSimData_Cls(conf.results,...
        conf.general.numOfTTIperSimulation,TTIblockDim,totNumOfUes,...
        totNumOfCells,carrierArray(ii).numOfAvailablePRB);    
end


% ITERATIONS
totTTI = conf.general.numOfTTIperSimulation;

x = cellTypeArray(macroInd);
varNumOfMacro = x.numOfCells;
fstPicoInd = cellTypeArray(picoInd).firstCellId;
lstPicoInd = fstPicoInd + cellTypeArray(picoInd).numOfCells - 1;
picoAllocType = cellTypeArray(picoInd).allocationType;

simProgress(1:numOfCarriers) = 0;

% load('carloSimulatorVars','varRandomPrb');

for varTTI=1:totTTI

    for ii=1:totNumOfUes            
        ues(ii).dataGenerator(conf.traffic);
    end

    fixedSpecAllocInitVar = zeros(1,numOfCarriers);

    for ii=1:numOfCarriers

        objCarrier = carrierArray(ii);

        % values computed for each cell
        % these values are computed in each iteration. 

        % for each cell, the throughput of their UEs are sorted from 
        % lowest to highest. The throughput that is in the position 
        % round(0.05*numOfUsersOfCell) is the value assigned as <outage> of
        % the cell.
        objCarrier.cellOutage = zeros(1,totNumOfCells);

        % the following property stores, for each cell, the sum of 
        % UE throughputs that are attached to the cell. 
        objCarrier.cellThroughput = zeros(1,totNumOfCells);

        totNumPrbs = objCarrier.numOfAvailablePRB;
        % for each <cell, prb> it stored the SINR corresponding to the pair
        % <ue,prb>, if ue is attached to the cell. Otherwise it stores 0.
        carrierArray(ii).cellVsPRB_SINR_mtx = zeros(totNumOfCells,...
           totNumPrbs);

        objCarrier.allocationMatrix = zeros(totNumOfUes,totNumPrbs);

        numOfCellsInCarrier = length(objCarrier.cellsInCarrier);
        for j = 1:numOfCellsInCarrier

            n_cell = objCarrier.cellsInCarrier(j);
            typeIndex = cls(n_cell).typeIndex;
            objType = cellTypeArray(typeIndex);
            % generate cell traffic

            cls(n_cell).updateBS(conf.traffic.trafficGeneratorType, ...
                [ues(:).UEbuffer]);

            % run spectrum allocation method

            cls(n_cell).runSpectrumAllocation(objCarrier,cellTypeArray,...
                fixedSpecAllocInitVar(ii),conf.specAlloc.softReusePWoffset,...
                conf.specAlloc.RBs_percentage,varTTI,conf.fsu,...
                conf.scheduling.manPRBvector);            

            % run user scheduling method
            % cls(n_cell).random_prbs = varRandomPrb{n_cell,varTTI};


            cls(n_cell).runUsersScheduler([carrierTypeArray(ii,:).schedulingType],...
                objCarrier,[ues(:).UEbuffer],conf.scheduling.UEperTTI,conf.PHY,...
                cellTypeArray,data.interfererIndexMatrix);      

%                 if varSaveVars  == 1    
%                     varSpectVect{n_cell} = [cls(n_cell).ueSchedVars.userSpectVect];
%                     varSpectAlloc{n_cell} = [cls(n_cell).raVars.spectrumAvailable];
%                     save(varfnamevars,'varSpectVect','-append');
%                     save(varfnamevars,'varSpectAlloc','-append');
%                 end

        end

        fixedSpecAllocInitVar(ii) = 1;


        varSchedArray = reshape([cls.ueSchedVars],numOfCarriers, totNumOfCells);
        varSchedArray2 = varSchedArray(ii,:);

        varRaArray = reshape([cls.raVars],numOfCarriers, totNumOfCells);
        varRaArray2 = varRaArray(ii,:);

        compute_SINR_matrix_v2(objCarrier, [ues.servingCellIndex],varSchedArray2,...
            varRaArray2, varTTI,conf.PHYmodelling,conf.l2s.capacityComputation,...
            conf.results.SINR_statistic, totNumOfCells, totNumOfUes);              

%             if varSaveVars==1
%                 varSinr = objCarrier.SINR_mtx;
%                 save(varfnamevars,'varSinr','-append');
%             end
%             
        if conf.results.additionalData

            computeResults(objCarrier,data.userVsCell_allocationMtx,...
                varSchedArray2,varRaArray2,conf.results,totNumOfCells);
        end
    end
        %% Traffic delivery evaluation  

    % the total capacity per UE considering all carriers.
    varSum = zeros(totNumOfUes,1);

    for ii=1:numOfCarriers
        varSum = varSum + sum(carrierArray(ii).C_mtx,2);
    end
    varSum = varSum.';

    varUEbuffer = [ues.UEbuffer];
    transmittedDataVector = conf.PHY.TTIduration .* varSum;    
    effectiveTransmittedData = min(varUEbuffer,transmittedDataVector);


    for ii=1:totNumOfCells
        cls(ii).evaluateTransmittedData(varUEbuffer,...
            conf.traffic.trafficGeneratorType,transmittedDataVector,...
            effectiveTransmittedData);
    end

    if conf.traffic.trafficGeneratorType == 1
        txData = transmittedDataVector;
    else
        txData = effectiveTransmittedData;
    end

    for ii=1:totNumOfUes
        ues(ii).evaluateTransmittedData(conf.traffic.trafficGeneratorType,...
            txData(ii));
    end

    %%% saving data
    for ii=1:numOfCarriers
        varSimDataArray(ii).saveDataOneIteration(conf.results,carrierArray(ii),...
            varTTI,TTIblockDim,cls,ues);
    end

    if rem(varTTI, TTIblockDim) == 0
        save(strcat(conf.storage.outputFolderPath,'longRunDataBlock',num2str( floor( varTTI/ TTIblockDim) )), 'varSimDataArray','conf');
    end

    simProg(ii) = floor(varTTI/conf.general.numOfTTIperSimulation*(100/10))*10;


    % diplaying simulation progress
    if simProg(ii) > simProgress(ii)
        disp(strcat({'Snap-shot simulation progress: '},num2str(simProg(ii)),{' %'}));
        simProgress(ii) = simProg(ii);
    end        
end

vartmp = zeros(1,totNumOfUes);
for ii=1:numOfCarriers
    %% Final results analysis
    dataAnalysisFnc_v2(carrierArray(ii), varSimDataArray(ii),conf.traffic.trafficGeneratorType);       

     vartmp = vartmp + varSimDataArray(ii).numberOfPRBperUE;
end

objFinalResults = finalResults_Cls();

objFinalResults.computeFinalResults(ues, cls, conf.PHY.TTIduration,...
    conf.traffic, data.userVsCell_allocationMtx, vartmp, totTTI,57,...
    cellTypeArray,conf.deployment.numOfUEperCellType,conf.deployment.numOfNetworkOperators);

%      save('lsaSimulatorVars','varSimDataArray','-append');
%     save('lsaSimulatorVars','objFinalResults','-append');
save('lsaSimulatorVars','objFinalResults');
save('lastSim');

disp('Snap-shot simulation loop finished');
disp('');
tfin = toc;
disp(strcat({'Snap-shot simulation: '},num2str(tfin),{' s'}));  

%%% Simulator profiling
if conf.results.profiling
    profile off;
end

% end

%%% save results here
out = objFinalResults; 



end