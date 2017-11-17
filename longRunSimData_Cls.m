classdef longRunSimData_Cls < handle
    
    properties
        % array (num of iterations, num of UEs)
        % longRunSimData_Cls::bandPerUE(ii,j) stores the number of assigned 
        % PRBs to UE j in the iteration ii
        bandPerUE;
        
        % array (num of iterations, num of cells)
        % each element stores the number of UEs that are attached to the
        % cell and that have a buffer size greater than 0.
        cellLoadOverTime;
        
        % array (num of iterations, num of cells). 
        % for each iteration ii and cell j, the throughput of UEs attached
        % to cell j are  sorted from lowest to highest. The throughput in 
        % position round(0.05*numOfUsersOfCell) is the value assigned as 
        % <outage> of cell j in the iteration ii.
        cellOutage;
        
        % array (num of iterations, num of cells, num of PRBs)
        % for each iteration and for each <cell ii, prb j> this property 
        % stores the SINR corresponding to the pair <ue k,prb j>, if k is 
        % attached to cell ii. Otherwise cellSpectEff_mtx(it, ii, j) = 0.
        cellSpectEff_mtx;
        
        % for each iteration, the following property stores, for each 
        % cell, the sum of throughputs corresponding to the UEs that are 
        % attached to the cell. 
        cellThroughput;
        
        % for each iteration, the following property stores the
        % interference signal corresponding to each pair <ue, prb> over the
        % current carrier
        int_mtx;
        
        % array (num of iterations, num of cells, num of PRBs)
        % this property stores the same values as <cellSpectEff_mtx>        
        SINR_mtx;        
        
        % array(num of UEs,num of PRBs,num of Iterations)
        % it stores the SINR per UE and PRB for each iteration
        SINR_UEperPRBarray;
        
        % array (num of iterations, num of UEs)
        % it stores the average of SINR of each UE over the assigned PRBs in
        % the current iteration 
        SINR_UEvector;   
        
        % array (min(TTIblockDim, numTTIperSimu),totNumUes)
        % for each iteration, it is stored the sum of capacity of each UE 
        % considering the capacity obtained for each assigned PRB
        userThOverTTIs_mtx;
        
        % array (1, number of UEs)
        % for each UE is computed once the simulation iterations are done, 
        % the quotient between 1) sum of average capacity of each UE over
        % iterations, and 2) total number of iterations where the UE had at
        % least one assigned PRB.
        UEspectralEfficiecy;
        
        % array (1, number of UEs)
        % for each UE, once the simulation iteratios are done, 
        % is computed the following expression:
        % 10*log10(2 .^ UEspectralEfficiency-1), where
        % <UEspectralEfficiency> is one of the properties of the current
        % object and its description is given in this file.
        
        UE_SINR_dB;
        
        % array (1, number of UEs)
        % it stores the quotient between sum of capacity of each UE, over
        % all carriers and over the assigned PRBs and the time
        % corresponding to the number of simulation iterations
        uth;
        
        % array (1, number of cells)
        % it stores the quotient between sum of capacity of each UE 
        % (attached to the cell, over all carriers and over the assigned 
        % PRBs and the time corresponding to the number of simulation 
        % iterations, where at least one of its attached UEs has a buffer
        % size greater than 0.
        cellTH;
        
        % array (number of iterations, number of UEs)
        % number of assigned PRBs to each UE, per iteration
        numberOfPRBperUE;
        
        % the following properties are assigned when the traffic model is
        % different from full buffer model.
        meanUEtrafficLoad;
        cellLoad;
        meanCellLoad;
        resourceUsagePercentage;
        % --------------------------------------------
    end
    
    methods
        function obj=longRunSimData_Cls(resultCfg,numTTIperSimu,...
                TTIblockDim,totNumUes,totNumCells,totNumPRBs)
            
            obj.UEspectralEfficiecy = zeros(1, totNumUes);
            obj.UE_SINR_dB = zeros(1, totNumUes);
            obj.meanUEtrafficLoad = zeros(1, totNumUes);
            obj.uth = zeros(1, totNumUes);
            obj.cellTH = zeros(1, totNumCells);
            obj.cellLoad = zeros(1, totNumCells);
            obj.meanCellLoad = 0;
            obj.resourceUsagePercentage = 0;
            obj.numberOfPRBperUE = zeros(1,totNumUes);
            
            if resultCfg.thPerUser > 0
                obj.userThOverTTIs_mtx = zeros(min(TTIblockDim, numTTIperSimu),totNumUes);
            end
            
            if resultCfg.thPerCell > 0
                obj.cellThroughput = zeros(numTTIperSimu, totNumCells);          
                obj.cellOutage = zeros(numTTIperSimu, totNumCells);
                obj.cellSpectEff_mtx = zeros(numTTIperSimu, ...
                    totNumCells, totNumPRBs);
            end
            
            if resultCfg.additionalData > 0
                obj.cellLoadOverTime = zeros(numTTIperSimu, totNumCells);
                obj.int_mtx = zeros(totNumUes,numTTIperSimu);
            end
            
            if resultCfg.SINR_UEperPRBmatrix
                obj.SINR_UEperPRBarray = zeros(totNumUes,totNumPRBs,numTTIperSimu);
            end
            
            if resultCfg.SINR_UEvector
                obj.SINR_UEvector(currTTI,:) = zeros(numTTIperSimu,totNumUes);
            end
            
            if resultCfg.bandPerUE
                obj.bandPerUE = zeros(numTTIperSimu,totNumUes);
            end
            
            if resultCfg.SINR_statistic
                obj.SINR_mtx = zeros(numTTIperSimu, ...
                    totNumCells, totNumPRBs);
            end
            
            
        end
    
        
        function saveDataOneIteration(obj,resultCfg,objCarrier,...
            currTTI,TTIblockDim, cls, ues)
            
            vartmp = (objCarrier.allocationMatrix).';
            varsum = sum(vartmp,1);
            obj.numberOfPRBperUE = obj.numberOfPRBperUE + varsum;
        
            if resultCfg.thPerUser
                obj.userThOverTTIs_mtx(rem(currTTI - 1, TTIblockDim) +1,:) = objCarrier.thPerUE;
            end

            if resultCfg.SINR_UEperPRBmatrix
                obj.SINR_UEperPRBarray(:,:,currTTI) = objCarrier.SINR_mtx;
            end
            
            if resultCfg.SINR_UEvector
                obj.SINR_UEvector(currTTI,:) = objCarrier.avSINRperUE;
            end
            
            if resultCfg.bandPerUE
                obj.bandPerUE(currTTI,:) = objCarrier.bandPerUE;
            end

            if resultCfg.thPerCell
                obj.cellOutage(currTTI,:) = objCarrier.cellOutage;
                obj.cellThroughput(currTTI,:) = objCarrier.cellThroughput;
                obj.cellSpectEff_mtx(currTTI,:,:) = objCarrier.cellVsPRB_SINR_mtx;
            end

            if resultCfg.SINR_statistic
                obj.SINR_mtx(currTTI,:,:) = objCarrier.cellVsPRB_SINR_mtx;
            end
            
            if resultCfg.additionalData >0
                varUEbuffer = [ues.UEbuffer];
                totNumCells = length(cls);
                
                for ii=1:totNumCells
                    varUEs = cls(ii).usersList;
                    varUEbufferTmp = varUEbuffer(varUEs);
                    obj.cellLoadOverTime(currTTI,ii) = sum(varUEbufferTmp>0);
                end
                obj.int_mtx(:,currTTI)= objCarrier.int_mtx;
            end            
        end
    end
end

