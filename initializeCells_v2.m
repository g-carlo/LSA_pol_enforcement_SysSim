function [cls,conf] = initializeCells_v2(conf,cellTypeArray,type2TypeArray,...
    carrierTypeArray,arrayAvPRBs)
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
                    
                    cls = initializeMacroHomogScenario_v2(cls,cellTypeArray(macroInd),...
                        varIntersiteDist,conf.general.numOfCarriers,picoAllocType,...
                        conf.scheduling.timerQoSEV,conf.general.numOfTTIperSimulation,...
                        arrayAvPRBs,maxTxPwArray,varDeploy);
                    
                case 2
                    maxTxPwArray = [carrierTypeArray(:,macroInd).maxTxPw];
                    
                    cls = initializeMacroHomogScenario_v2(cls,cellTypeArray(macroInd),...
                        varIntersiteDist,conf.general.numOfCarriers,picoAllocType,...
                        conf.scheduling.timerQoSEV,conf.general.numOfTTIperSimulation,...
                        arrayAvPRBs,maxTxPwArray,varDeploy);
                    
                    maxTxPwArray = [carrierTypeArray(:,picoInd).maxTxPw];
                    
                    cls = initializeHetNetScenario_v2(cls,cellTypeArray,...
                        type2TypeArray,varIntersiteDist,conf.general.numOfCarriers,...
                        conf.scheduling.timerQoSEV,conf.general.numOfTTIperSimulation,...
                        arrayAvPRBs,maxTxPwArray);   
                    
                    varMacroInTheNet = 1;
                
                case 3          % SMALL CELL WITH SPPP DISTRIBUTION WITHIN A SQUARE OF MAX. LENGTH, uniform distribution of pico cells per macro
                    if varMacroInTheNet == 0
                        
                        maxTxPwArray = [carrierTypeArray(:,macroInd).maxTxPw];
                        
                        cls = initializeMacroHomogScenario_v2(cls,cellTypeArray(macroInd),...
                            varIntersiteDist,conf.general.numOfCarriers,picoAllocType,...
                            conf.scheduling.timerQoSEV,conf.general.numOfTTIperSimulation,...
                            arrayAvPRBs,maxTxPwArray,varDeploy);
                        
                        varMacroInTheNet = 1;
                    end
                    
                    maxTxPwArray = [carrierTypeArray(:,smallInd).maxTxPw];
                    
                    cls = initializeSmallCellsScenario_SPPP_square_numPerMacro_v2(cls,...
                        cellTypeArray,type2TypeArray,varIntersiteDist,conf.general.numOfCarriers,...
                        conf.scheduling.timerQoSEV,conf.general.numOfTTIperSimulation,arrayAvPRBs,...
                        maxTxPwArray);
                
                % case 4: a random number of small cells are randomly distributed within a square 
                % Side length of the square is conf.deployment.squareLength
                % the random number of small cells is given by: max(poissrnd(cellDensity*(squareLength^2),1))
                % where cellDensity is conf.deployment.cellDensity
                
                case 4       
                    objType = cellTypeArray(smallInd);
                    
                    maxTxPwArray = [carrierTypeArray(:,smallInd).maxTxPw];
                    
                    
                    [numOfCells,cls] = initializeSmallCellsScenario_SPPP_square_v2(cls,...
                        varDeploy,objType,type2TypeArray(smallInd,smallInd).minDistance,...
                        conf.general.numOfCarriers,conf.scheduling.timerQoSEV,...
                        conf.general.numOfTTIperSimulation,arrayAvPRBs,maxTxPwArray);
                    
                    
                    cellTypeArray(smallInd).numOfCells = numOfCells;
                    cellTypeArray(smallInd).numOfUEs = varDeploy.numOfUEperCellType(smallInd)*numOfCells;
                    
                    
                case 5          % SMALL CELL WITH SQUARED GRID DISTRIBUTION WITHIN A SQUARE 
                    objType = cellTypeArray(smallInd);
                    
                    maxTxPwArray = [carrierTypeArray(:,smallInd).maxTxPw];
                    
                    [varDeploy,numOfCells,cls] = initializeSmallCellsScenario_square_v2(cls,...
                        varDeploy,objType,conf.general.numOfCarriers,conf.scheduling.timerQoSEV,...
                        conf.general.numOfTTIperSimulation,arrayAvPRBs,maxTxPwArray);  
                    
                    
                    cellTypeArray(smallInd).numOfCells = numOfCells;
                    cellTypeArray(smallInd).numOfUEs = varDeploy.numOfUEperCellType(smallInd)*numOfCells;
                    
                
                case 6 % SMALL CELL WITH SPPP DISTRIBUTION WITHIN A SQUARE OF MAX. LENGTH
                    if varMacroInTheNet == 0
                        
                        maxTxPwArray = [carrierTypeArray(:,macroInd).maxTxPw];
                        
                        cls = initializeMacroHomogScenario_v2(cls,cellTypeArray(macroInd),...
                            varIntersiteDist,conf.general.numOfCarriers,picoAllocType,...
                            conf.scheduling.timerQoSEV,conf.general.numOfTTIperSimulation,...
                            arrayAvPRBs,maxTxPwArray,varDeploy);
                        
                        varMacroInTheNet = 1;
                    end
                    
                    
                     maxTxPwArray = [carrierTypeArray(:,smallInd).maxTxPw];
                    
                    [numOfCells,cls] = initSmallCellsSPPPsquareMaxLen_v2(cls,...
                        conf.deployment,cellTypeArray(smallInd),...
                        type2TypeArray(smallInd,smallInd).minDistance,...
                        cellTypeArray(macroInd),conf.general.numOfCarriers,...
                        conf.scheduling.timerQoSEV,conf.general.numOfTTIperSimulation,...
                        arrayAvPRBs,maxTxPwArray);    
                    
                                        
                    cellTypeArray(smallInd).numOfCells = numOfCells;
                    cellTypeArray(smallInd).numOfUEs = varDeploy.numOfUEperCellType(smallInd)*numOfCells;
                    
                                    
                                                        
            end           
        end   
        
        numOfMacroCells = cellTypeArray(macroInd).numOfCells;
        numOfPicoCells = cellTypeArray(picoInd).numOfCells;
        numOfSmallCells = cellTypeArray(smallInd).numOfCells;
        
        totNumOfCells = numOfMacroCells + numOfPicoCells + numOfSmallCells;
        
        conf.deployment.numOfCells = totNumOfCells;
    end    
end

 