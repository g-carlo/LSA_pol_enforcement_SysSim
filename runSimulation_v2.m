function runSimulation_v2()
  
    
    clear;
    clc;
  
    global macroInd;
    global picoInd;
    global smallInd;
    
    conf = getInputParameters();
    
    % set working folder in the path
    
    % conf.storage.workingFolder is the path that defines the current 
    % directory. Folders,subfolders and files that are in the current 
    % directory are added to the matlab path.
    
    cd(conf.storage.workingFolder);
    addpath(genpath('.'));
    
    % number of cell types: macro, pico, small cells, ...
    numOfCellTypes = conf.implementation.numOfCellTypes; 
    
    % number of carriers that will be considered in the current simulation.
    numOfCarriers = conf.general.numOfCarriers;
        
    % for each carrier, an object carrier_Cls is instantiated. Some of its 
    % properties (those that have the same name as the property) 
    % are assigned from the input parameters conf.carrier
    
    carrierArray = instantiateObjectFromInputParam(conf.carrier,...
        'carrier_Cls',1,numOfCarriers);
    
    % cell type properties that are assigned from input parameters
    % conf.cellType
    cellTypeArray = instantiateObjectFromInputParam(conf.cellType,...
        'cellType_Cls', 1,numOfCellTypes);
    
    d = conf.deployment.intersiteDistance; % intersite distance
    l = d * cos(pi/6);                     

    for ii=1:numOfCellTypes
        cellTypeArray(ii).index = ii;
        cellTypeArray(ii).numOfCells = 0;
        cellTypeArray(ii).numOfUEs = 0;
        
        if conf.general.useWA == 1  
            cellTypeArray(ii).shiftArray = [0 -3*l+1i*3.5*d 2*l+1i*4*d ...
                 5*l+1i*0.5*d 3*l-1i*3.5*d -2*l-1i*4*d -5*l-1i*0.5*d];
        else
            cellTypeArray(ii).shiftArray = [0 0 0 0 0 0 0];       
        end
    end
    
    macroInd = find(strcmp(conf.cellType.description,'macro')==1);
    picoInd = find(strcmp(conf.cellType.description,'pico')==1);
    smallInd = find(strcmp(conf.cellType.description,'small')==1);
 
     %%% replica positions for wrap around structure.
        % vector with the position of replicas of the WA structure
        
%     cellTypeArray(macroInd).shiftArray = [0 -3*l+1i*3.5*d 2*l+1i*4*d ...
%             5*l+1i*0.5*d 3*l-1i*3.5*d -2*l-1i*4*d -5*l-1i*0.5*d];
%     
%     cellTypeArray(picoInd).shiftArray = [0 -3*l+1i*3.5*d 2*l+1i*4*d ...
%             5*l+1i*0.5*d 3*l-1i*3.5*d -2*l-1i*4*d -5*l-1i*0.5*d];
%         
%     cellTypeArray(smallInd).shiftArray = [0 -3*l+1i*3.5*d 2*l+1i*4*d ...
%             5*l+1i*0.5*d 3*l-1i*3.5*d -2*l-1i*4*d -5*l-1i*0.5*d];
%     
    carrierCellTypeArray = instantiateObjectFromInputParam(conf.carrier.cellType,...
        'carrierCellType_Cls', numOfCarriers, numOfCellTypes);
   
    carrierCellTypeArray = reshape(carrierCellTypeArray,numOfCarriers,numOfCellTypes);
    
    for ii=1:numOfCarriers
        carrierArray(ii).index = ii;        
        for j = 1:numOfCellTypes
            carrierCellTypeArray(ii,j).carrierIndex = ii;
            carrierCellTypeArray(ii,j).typeIndex = j;
            carrierCellTypeArray(ii,j).determineCoeffPathLoss(...
                carrierArray(ii).frequency,macroInd,picoInd,smallInd);
        end
    end
    
    
    cellType2cellTypeArray = instantiateObjectFromInputParam(conf.cellType.otherCellType,...
        'cellType2cellType_Cls', numOfCellTypes, numOfCellTypes);
   
    cellType2cellTypeArray = reshape(cellType2cellTypeArray,numOfCellTypes,...
        numOfCellTypes);
         
    [conf,carrierArray,cellTypeArray,carrierCellTypeArray] = ...
        calculateDependentConfParameters_v2(conf, macroInd,picoInd,smallInd,...
            cellTypeArray,carrierArray,cellType2cellTypeArray,carrierCellTypeArray);

    varAntennaType = unique(reshape(conf.carrier.cellType.antennaType(1:numOfCarriers,:),1,numOfCarriers*numOfCellTypes));
    antennaArray = cell(3,1);
    for ii = 1:length(varAntennaType)
        if varAntennaType(ii) == 1
            varElectricalDowntilt = [carrierArray(:).electricalDowntilt];
            antennaArray{varAntennaType(ii)} = antennaMkI_Cls(conf.antenna,varElectricalDowntilt);
        else
            antennaArray{varAntennaType(ii)} = omnidirectAntenna_Cls();
        end
    end
        
    mainSim_v2(conf,carrierArray,cellTypeArray,carrierCellTypeArray,...
        cellType2cellTypeArray);
end

function objArray = instantiateObjectFromInputParam(structName,...
    varClassName,dim1,dim2)
    
    numOfElem = dim1*dim2;
    objArray = feval(varClassName);    
    names=fieldnames(structName)';
    noOfFields=length(names);
    
    values = cell(noOfFields,1);
    for m = 1:noOfFields
        if ~isstruct(structName.(names{m}))
            values{m} = getfield(structName,names{m});
            if iscell(values{m})
                vartmp = values{m};
            else
                vartmp = values{m}(1:dim1,1:dim2);
            end
            values{m} = reshape(vartmp,1,numOfElem);          
        end
    end
    
    
    for n = 1:numOfElem
        objArray(n) = feval(varClassName);
        for m = 1:noOfFields
            if ~isstruct(structName.(names{m}))
                objArray(n) = setfield(objArray(n),names{m},...
                    values{m}(n));
            end
        end        
    end    
end

