function [conf,carrierArray,cellTypeArray,carrierCellTypeArray] = ...
    calculateDependentConfParameters_v2(conf,macroInd,picoInd,smallInd,...
        cellTypeArray,carrierArray,cellType2cellTypeArray,carrierCellTypeArray)
           
    
    
    numOfCarriers = conf.general.numOfCarriers;
    numOfCellTypes = conf.implementation.numOfCellTypes;
       
    varCellScenario = conf.deployment.cellScenario;
    
    % scenarios 2, 3, and 6 include macro homogeneous scenario. For that
    % reason, if some of these scenarios are specified in addition to the 
    % scenario 1, which is the macro homogeneous scenario, then scenario 1
    % is removed from the list because it is already included with the 
    % other ones.
    
    if any(ismember([2,3,6], varCellScenario)) && any(varCellScenario==1)
        varIndex = varCellScenario==1;
        varCellScenario(varIndex)=[];
    end
    
    % just one scenario that include small cells is allowed. Scenarios with
    % small cells are 3,4,5 and 6, with different cell distributions. Just 
    % one of these may be specified per simulation.
    
    varIndex = find(varCellScenario>=3);
    varLen = length(varIndex);
    if varLen>1
        varCellScenario(varIndex(1:varLen-1))=[];
    end
    
    conf.deployment.cellScenario = varCellScenario;
    
    % Number of cells per type: macro, pico, small
    
    % if one of the following scenarios will be used during the simulation: 
    % 1, 2, 3, or 6, then a given number of macro cells will be deployed
    % in the current simulation.
    
    if any(ismember([1,2,3,6], varCellScenario))  
        if conf.general.SLS_activated
            cellTypeArray(macroInd).numOfCellsPerOperator = ... %% number of cells
                conf.deployment.numOfCellsPerSite * conf.deployment.numOfSites;
            cellTypeArray(macroInd).numOfCells = cellTypeArray(macroInd).numOfCellsPerOperator ...
                * conf.deployment.numOfNetworkOperators;
        else 
            num_of_sites = conf.deployment.N_sites_extedndedScenario;
            conf.deployment.numOfSites = num_of_sites;
            cellTypeArray(macroInd).numOfCellsPerOperator = ...
                conf.deployment.numOfCellsPerSite * conf.deployment.numOfSites;
            cellTypeArray(macroInd).numOfCells = cellTypeArray(macroInd).numOfCellsPerOperator ...
                * conf.deployment.numOfNetworkOperators;
            if num_of_sites == 11^2;
                conf.deployment.incumbentSouthBound             = -2500;%-4750;
                conf.deployment.incumbentNorthBound             = 2500;%+4750;
                conf.deployment.incumbentEastBound              = 2165;%3897;
                conf.deployment.incumbentWestBound              = -2165;%-3897;
            elseif num_of_sites == 21^2;
                conf.deployment.incumbentSouthBound             = -4750;
                conf.deployment.incumbentNorthBound             = +4750;
                conf.deployment.incumbentEastBound              = 3897;
                conf.deployment.incumbentWestBound              = -3897;
            end    
        end
    else
        cellTypeArray(macroInd).numOfCells = 0;
    end
    
    % pico cells will be deployed if one of the scenarios is 2:
    % heterogeneous network: macro + pico
    if any(varCellScenario==2)
        cellTypeArray(picoInd).numOfCells = cellTypeArray(macroInd).numOfCells * ...
            cellType2cellTypeArray(macroInd,picoInd).numOfOtherCells;          
        if conf.deployment.UEdistributionConfiguration==1
            conf.deployment.numOfUEperCellType(picoInd) = 0; 
        end       
              
    else
        cellTypeArray(picoInd).numOfCells = 0;
        if any(conf.cellType.allocationType == 4) % FSU version II
            error(['Allocation type 4 is not allowed when no pico cell is '...
            'part of the network.']);
        end
    end
    
    % scenario 3: per each macro cell a given number of small cells will be
    % deployed.
    if any(varCellScenario==3)
        cellTypeArray(smallInd).numOfCells = cellTypeArray(macroInd).numOfCells * ...
            cellType2cellTypeArray(macroInd,smallInd).numOfOtherCells;          
    else
        if ~(any(ismember([4,5,6],varCellScenario)))
            cellTypeArray(smallInd).numOfCells = 0;
        end        
    end
    
    varLen = length(varCellScenario);
    
    for ii=1:varLen
        switch varCellScenario(ii)
            case 1 % HOMOGENEOUS NETWORK: ONLY MACRO-CELLS
                % n. of macro-cells in the grid
                if conf.deployment.UEdistributionConfiguration == 2            
                    conf.deployment.UEdistributionConfiguration = 1;
                end                
            case {3,4,5,6}  % SMALL CELLS WITHIN A SQUARE
                cellTypeArray(smallInd).minDistToUE = (conf.deployment.cellDensity)^(-1/2)*0.01;
                for j=1:numOfCarriers
                    carrierCellTypeArray(j,smallInd).antennaType = 3;
                end
        end      
    end
    
    varNumPerType = conf.deployment.numOfUEperCellType;
    
    for ii=1:numOfCellTypes
        cellTypeArray(ii).numOfUEs = varNumPerType(ii)*cellTypeArray(ii).numOfCells;
        
        switch ii
            case macroInd
                cellTypeArray(ii).firstCellId = 1;
            case picoInd
                cellTypeArray(ii).firstCellId = 1+ cellTypeArray(macroInd).numOfCells;
            case smallInd
                cellTypeArray(ii).firstCellId = 1 + cellTypeArray(macroInd).numOfCells + ...
                    cellTypeArray(picoInd).numOfCells;
        end      
    end
        
    for j=1:numOfCellTypes     
        switch j
            case macroInd
                cellTypeArray(j).firstUserId = 1;
            case picoInd
                cellTypeArray(j).firstUserId = cellTypeArray(macroInd).numOfUEs+1;
            case smallInd
                cellTypeArray(j).firstUserId = cellTypeArray(macroInd).numOfUEs+...
                    cellTypeArray(picoInd).numOfUEs + 1;
        end
                
        for ii=1:numOfCarriers    
             % electrical antenna downtilt
            switch carrierArray(ii).scenarioCase
                case 1;
                    carrierArray(ii).electricalDowntilt = 15;
                    
                
                case 2;
                    carrierArray(ii).electricalDowntilt = 6;
                
                otherwise;
                    error('Scenario case must be either 1 or 2');
            end
            
            
        end        
    end
    
    for ii=1:numOfCarriers    
        carrierArray(ii).operatorsInCarrier = find(conf.carrier.operator.carrierPerOperator(ii,:));           
    end        
    
end