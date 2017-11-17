function [cls,conf] = initializeCells_noSLS(conf,cellTypeArray,type2TypeArray,...
    carrierTypeArray,arrayAvPRBs,N_sites_extedndedScenario)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                         SIMULATION MODE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % constructor:
        
    global macroInd;
    global picoInd;
    global smallInd;
    
    
    cls = cell_Cls(cellTypeArray(macroInd),0,0,0,[],0,0);       
    

    if strcmp(conf.general.mode,'montecarlo') || strcmp(conf.general.mode,'static')     % MODE: SIMULATION
        
        varCellScenario = sort(conf.deployment.cellScenario,1);
        
        varIntersiteDist = conf.deployment.intersiteDistance; 
        varLen = length(varCellScenario);
        
        varDeploy = conf.deployment;
        varMacroInTheNet = 0;
        
        picoAllocType = cellTypeArray(picoInd).allocationType;
        
        for ii = 1:varLen
            switch varCellScenario(ii)
                case 1
                    
                    maxTxPwArray = [carrierTypeArray(:,macroInd).maxTxPw];
                    
                    cls = initializeMacroHomogScenario_noSLS(cls,cellTypeArray(macroInd),...
                        varIntersiteDist,conf.general.numOfCarriers,maxTxPwArray,varDeploy,N_sites_extedndedScenario);
                    
                otherwise
                    
                    error('Scenario not implemented');
                                                        
            end           
        end   
        
        numOfMacroCells = cellTypeArray(macroInd).numOfCells;
        numOfPicoCells = cellTypeArray(picoInd).numOfCells;
        numOfSmallCells = cellTypeArray(smallInd).numOfCells;
        
        totNumOfCells = numOfMacroCells + numOfPicoCells + numOfSmallCells;
        
        conf.deployment.numOfCells = totNumOfCells;
    end    
end

 