%%% SON SIMULATOR: cell class
%
%   Author: Carlo Galiotto
%
%   Date:   August 2010
%   Last Modification: March, 2015
%
%   Description: this class implements the cell (both pico and macro) class, defines the
%   properties of a eNB and the methods 
%

classdef cell_Cls < handle
    
    properties        
        % object instantiated from class: position_Cls 
        % position of the cell BS
        pos;                                
        
        % integer that identifies the cell type corresponding to the current cell 
        % it coincides with the position where the object from class cellType_Cls
        % is stored in the array of cellType_Cls objects.
        typeIndex;               
        
        % complex number x+iy that specifies the position in x-y plane
        complexPosition;
        
        % real number. It is the antenna boresight in x-y plane. It is
        % given in radians.
        orientation;            
        
        % integer. It is the cell identifier. It coincides with the position 
        % where the current object is stored in the array of objects 
        % instantiated from the class cell_Cls.
        cellID;
        
        % integer. This is the network operator identifier
        networkOP_ID; 
        
        % identifier of the macro cell to which the current cell is linked.
        macroCellBelongingID;        
        
        % integer. Number of users that are attached to the current cell. 
        loadNOU;
        
        % array of integers. Each element of the array corresponds to 
        % one UE identifier.
        
        usersList;                          
             
        % boolean. It is true, if the current cell is a macro cell.
        % Otherwise, it is false.
        isMacro;                            
        
        % boolean. It is true, if the current cell is a pico cell.
        % Otherwise, it is false.
        isPico;     
        
        % boolean. It is true, if no UE is attached to the current cell. 
        % Otherwise, it is false. 
        isEmptyCell;
        
        % Integer. Index of sector within the site          
        macroCellSectorIndex;                               
              
         % array of integers. UE identifiers that are attached to 
         % the current cell. According to the traffic generator type, 
         % these users must have their buffer size greater than 0.
         
        usersToBeServedVect;
        
        % Integer. Index of the next user to be served by the cell,
        % according to the user scheduling algorithm.
        nextUserIndex;
        
        cumDataTXperCell;
        
        activeTTIcntPerCell;
                        
        % array of objects fsuVariables_Cls. 
        % One element per carrier
        fsuVars; 
        
        % array of resource allocation variables represented by objects
        % from class resourceAllocation_Cls
        % this array has oone element per carrier.
        raVars;
        
        % array of user scheduling variables, which are objects from class
        % userScheduling_Cls.
        % this array has one element per carrier.
        ueSchedVars;
        
        % THE FOLLOWING PROPERTIES HAVE BEED ADDED JUST TO COMPARE RESULTS
        % BETWEEN CARLO SIMULATOR AND THIS ONE.
        % random_prbs;
        % randUesRR;
    end     % end of properties
    
    methods
        
        % constructor
        function obj = cell_Cls(objType, varCellId, numOfCarriers,...
                varMaxTxPw,arrayAvPRBs,timerQoSEV, numOfTTIperSimu)         
            
            if nargin > 4   % if number of arguments in functions are >4
            
                obj.typeIndex = objType.index;
                obj.cellID = varCellId;
                obj.macroCellBelongingID = varCellId;

                obj.pos = position_Cls(0,0,0);
                obj.macroCellSectorIndex = 0;
                obj.complexPosition = 0;
                obj.orientation = 0;           
                obj.loadNOU = 0;                       
                obj.isEmptyCell = true;
                obj.usersList = [];
                obj.usersToBeServedVect = [];
                obj.isMacro = false;
                obj.isPico = false;

                obj.nextUserIndex = 1;

                obj.cumDataTXperCell = 0;

                obj.activeTTIcntPerCell = 0;

                if strcmp(objType.description,'macro')==1
                    obj.isMacro = true;                
                elseif strcmp(objType.description,'pico')==1                
                    obj.isPico = true;
                end            


                obj.fsuVars = fsuVariables_Cls(0,0,0,0);
                obj.raVars = resourceAllocation_Cls(0,0);
                obj.ueSchedVars = userScheduling_Cls(0);


                for ii=1:numOfCarriers
                    obj.fsuVars(ii) = fsuVariables_Cls(timerQoSEV,numOfTTIperSimu,...
                        objType.alphaPF,objType.betaPF);
                    obj.raVars(ii) = resourceAllocation_Cls(varMaxTxPw(ii),arrayAvPRBs(ii));
                    obj.ueSchedVars(ii) = userScheduling_Cls(arrayAvPRBs(ii));
                end
            
            elseif  nargin == 4
                
                obj.typeIndex = objType.index;
                obj.cellID = varCellId;
                obj.macroCellBelongingID = varCellId;

                obj.pos = position_Cls(0,0,0);
                obj.macroCellSectorIndex = 0;
                obj.complexPosition = 0;
                obj.orientation = 0;           
                obj.loadNOU = 0;                       
                obj.isEmptyCell = true;
                obj.isMacro = false;
                obj.isPico = false;

                obj.nextUserIndex = 1;

                if strcmp(objType.description,'macro')==1
                    obj.isMacro = true;                
                elseif strcmp(objType.description,'pico')==1                
                    obj.isPico = true;
                end            
                
            end
            
        end     % end of function
        
        function evaluateTransmittedData(obj,varUEbuffer,trafficGenType,...
                transmittedDataVector,effectiveTransmittedData)
            switch trafficGenType
                case 1  % full buffer                    
                    obj.cumDataTXperCell = obj.cumDataTXperCell + ...
                        sum(transmittedDataVector(obj.usersList));
                    
                    obj.activeTTIcntPerCell = obj.activeTTIcntPerCell + 1;
                    
                case {2,3}  % FTP model 1
                    
                    obj.cumDataTXperCell = obj.cumDataTXperCell + ...
                        sum(effectiveTransmittedData(obj.usersList));
                    
                    if sum(varUEbuffer(obj.usersList)>0)>0
                        obj.activeTTIcntPerCell = obj.activeTTIcntPerCell + 1;
                    end
                otherwise
                    error('No other traffic models available');
            
            end
        end     % end of function
        
        function setPosition(obj,x,y,z)         % set position of BS of the cell
            
            obj.pos.setPosition(x,y,z);
            obj.complexPosition = x+1i*y;
            
            
        end     % end of function
        
        function y = getComplexPosition(obj)    % get the x-y plane position of the BS in complex coordinates
            
             y = obj.pos.getComplexCoordinate();
            
        end     % end of function
        
                
        function y = belongsToCell(obj,x,varIntersiteDist)       % returns 1 if the complex position x belongs to the cell, 0 otherwise
            
            l = abs(x-obj.pos.getComplexCoordinate());
            alpha = angle(x-obj.pos.getComplexCoordinate());
            
            d = varIntersiteDist/2;
            
            if obj.isMacro
                    if angle(exp(1i*(alpha-obj.orientation)))>= -pi/6 && angle(exp(1i*(alpha-obj.orientation)))< pi/6
                        
                        tmp = abs(l*cos(alpha-obj.orientation));
                        if tmp >= 0 && tmp < d
                            y = 1;
                        else 
                            y = 0;
                        end
                        
                    elseif angle(exp(1i*(alpha-obj.orientation))) >= -pi/3 && angle(exp(1i*(alpha-obj.orientation))) < -pi/6
                        
                        tmp = abs(l*cos(-pi/3 + obj.orientation - alpha));
                        if tmp >= 0 && tmp < d
                            y = 1;
                        else
                            y = 0;
                        end
                        
                    elseif  angle(exp(1i*(alpha-obj.orientation))) >= pi/6 && angle(exp(1i*(alpha-obj.orientation))) < pi/3
                        
                        tmp = abs(l*cos(pi/3 + obj.orientation - alpha));
                        if  tmp >= 0 && tmp < d
                            y = 1;
                        else
                            y = 0;
                        end
                    else
                        y = 0;
                    end
                    
            elseif obj.isPico      % pico
                    
                    if l >= 0 && l < d
                        y = 1;
                    else 
                        y = 0;
                    end
                    
            else      % small hexagonal cell
                    
                     tmp = l*cos( alpha - ( pi/6 + fix(alpha./(pi/3)) * pi/3 ) );
                     if tmp >= 0 && tmp < d
                         y = 1;
                     else 
                         y = 0;
                     end
                    
            end
            
            y = y>0;
            
        end     % end of method
        
        
% input parameters:
% ---------------------------------------
% - trafficTypeCarrier: traffic generator type (full buffer or burst traffic),
% which has been specified for the current carrier

% - varUEbuffer: buffer of each UE, considering all UEs in the network.

    function updateBS(obj,trafficTypeCarrier,varUEbuffer)
    
        switch trafficTypeCarrier
            case 1                    
                obj.usersToBeServedVect = obj.usersList;                    
            case {2,3}                    
                active_ues = varUEbuffer(obj.usersList)>0;
                obj.usersToBeServedVect = obj.usersList(active_ues);
        end

        obj.loadNOU = length(obj.usersList);

    end

% input parameters:
% ---------------------------------------
% - trafficTypeCarrier: traffic generator type (full buffer or burst traffic),
% which has been specified for the current carrier

% - varUEbuffer: buffer of each UE, considering all UEs in the network.

% - objCarrier: object carrier_Cls, which represent data of a given
% carrier.

% - macroFSUfeedbackFlag: this variable is specified per cell type. It
% is true if the cell type is macro and if the pico allocation type is 4.
% Otherwise, it is false.

% - currTTI: current TTI of the simulation (iteration number)

% - schedCfg: input parameters related to scheduling process
% (conf.scheduling)

% - fsuCfg: input parameters related to fsu process (conf.fsu)

% - reuseScheme: reuse scheme corresponding to the type of the current cell
    
       
        
    function setAntenna(obj,conf,type)              % set antenna type of the BS

        if strcmp(type, 'Antenna MkI')
            obj.antenna = antennaMkI_Cls(conf);
        elseif strcmp(type, 'Antenna Calib')
            obj.antenna = antennaCalib_Cls(conf);
        elseif strcmp(type, 'Omnidirectional')
            obj.antenna = omnidirectAntenna_Cls;
        else
            obj.antenna = omnidirectAntenna_Cls;
        end

    end       % end of function

        function y = getAntennaGain(obj,x_phi,x_theta)          % get antenna gaiin along x_phi direction on the x-y plane and along the x_theta tilt angle
            
            y = obj.antenna.getGain(angle(exp(1i*(x_phi-obj.orientation))),x_theta);
            
        end       % end of function
        
        
        function attachUsers(obj,x)        % attach user with index x to the cell
            
            if isempty(obj.usersList)
                obj.usersList = x;
            else
                obj.usersList = sort([obj.usersList x]);
            end
        end         % end of function
        
        
        function runSpectrumAllocation(obj,objCarrier,cellTypeArray,...
                fixedSpecAllocInitVar,softReusePWoffset,RBs_percentage,...
                currTTI,fsuCfg,manualPRBvector)     % run spectrum allocation algorithm (FSU)
            
           global objCell;
           
           objCell = obj;
            
            varAllocType = cellTypeArray(obj.typeIndex).allocationType;
            
            fsuObj = obj.fsuVars(objCarrier.index);
            
            switch varAllocType
                
                case 1 
                    
                    if ~fixedSpecAllocInitVar  
                        reuseScheme = cellTypeArray(obj.typeIndex).reuseScheme;
                        
                        obj.fixedSpectrumAllocation(objCarrier,...
                            softReusePWoffset,reuseScheme);                       
                    end
                    
                case 2
                    error('Not implemented')
                    
                case 3
                    obj.manualPRBsetting(objCarrier.index,objCarrier.numOfAvailablePRB,...
                        manualPRBvector);
                    
                case 4
                    error('Not implemented')
                                       
                case 5
                   
                    obj.frequencyALOHA(objCarrier.index,objCarrier.numOfAvailablePRB,...
                        RBs_percentage);
                    
                otherwise
                    
                    error('No any other spectrum allocation schemes implemented')
                    
            end
            
        end         % enf of function
            
        
        function runUsersScheduler(obj,varSchedTypes,objCarrier,...
                varUEbuffer,UEperTTI,PHYcfg,cellTypeArray,interfererIndexMatrix)        
            
            global macroInd;
            global picoInd;
            
            carrierNum = objCarrier.index;
            varSchedAlgorithm = varSchedTypes(obj.typeIndex);
            
            if varSchedAlgorithm == 2 || varSchedAlgorithm == 4
                fsuObj = obj.fsuVars(carrierNum);
                ueSchObj = obj.ueSchedVars(carrierNum);
                
                objType = cellTypeArray(obj.typeIndex);

                currReuseScheme = objType.reuseScheme;
                currAllocType = objType.allocationType;
                gammaPF_type = objType.gammaPF;       
                
                if ~ueSchObj.initPFVar
                    objMacro = cellTypeArray(macroInd);
                    objPico = cellTypeArray(picoInd);

                    ueSchObj.initPFVar = 1;
                    obj.RoundRobinScheduler(objCarrier,varUEbuffer,UEperTTI);
                
                    if objPico.allocationType == 4
                        fsuObj.initializeSchedulerAndFSUParamters(objCarrier,objMacro,...
                            objPico,interfererIndexMatrix);
                    end
                    ueSchObj.lastActiveUEindex = varUEbuffer(obj.usersList)>0;
                    ueSchObj.presThVect = zeros(size(obj.usersList));
                    ueSchObj.pastThVect = ueSchObj.presThVect;                
                else  
                    if currAllocType == 1 && currReuseScheme == 4                        
                        obj.ProportionalFairScheduler_SFR(objCarrier,UEperTTI,...
                            PHYcfg.TTIduration,fsuObj.alphaPF,fsuObj.betaPF,gammaPF_type,...
                            varUEbuffer);
                    else    
                        if varSchedAlgorithm == 2
                            obj.ProportionalFairScheduler(varUEbuffer,objCarrier,...
                                UEperTTI,PHYcfg.TTIduration,fsuObj.alphaPF,fsuObj.betaPF,gammaPF_type);                       
                        else
                            obj.PropFairSchedMultiCarriers(varUEbuffer,objCarrier,...
                                UEperTTI,gammaPF_type,PHYcfg.TTIduration);
                        end
                    end
                end
            else
                switch  varSchedAlgorithm
                    case 1
                        obj.RoundRobinScheduler(objCarrier,varUEbuffer,UEperTTI);
                    case 3
                        obj.FrequencySelectiveScheduler(objCarrier,UEperTTI,varUEbuffer);
                end                
            end
            
            if ~obj.isEmptyCell
               raObj = obj.raVars(carrierNum);
               ueSchObj = obj.ueSchedVars(carrierNum);
            
               % number of UEs attached to the current cell that have been
               % assigned with some PRB in the current carrier.
               
               % the array ueSchedObj.userSpectVect stores the UE index,
               % that have been allocated to the PRBs allocated to the
               % current cell.
               raObj.currentResourceUsage = sum(ueSchObj.userSpectVect>0);
               raObj.currentResourceAvailability = length(raObj.spectrumAvailable);
            end
                
        end     % end of function
        
    end     % end of methods
 

    methods (Access=private)        % private methods

        function ProportionalFairScheduler(obj,varUEbuffer,objCarrier,...
                UEperTTI,TTIduration,alphaPF,betaPF,gammaPF)
            
            carrierNum = objCarrier.index;
            ueSchObj = obj.ueSchedVars(carrierNum);
            raObj = obj.raVars(carrierNum);
            % varNumOfAvPRB = objCarrier.numOfAvailablePRB;
            
            
            if ueSchObj.initPFVar >0
            
                av_sp = raObj.spectrumAvailable;
                
                ueSchObj.userSpectVect = zeros(1,length(av_sp));
                lou = obj.usersToBeServedVect;
                
                
                if ~isempty(lou) && ~isempty(av_sp>0)
                    
                    %%% STEP 1: update past throughput
                    
                    nou = length(lou);
                    n_uePerTTI = UEperTTI;
                    cellUEindex = obj.usersList;
                    av_sp = av_sp(av_sp>0);
                    n_prb = length(av_sp);                    
                    min_nOfUEs = min(n_uePerTTI,nou);
                    
%                     
%                     SINR_mtx_cellUEs = obj.dataVar.SINR_mtx_estim_fsu(cellUEindex,av_sp);
%                     experiencedTHperPRB_mtx = shannonMkI(SINR_mtx_cellUEs);
                    experiencedTHperPRB_mtx = objCarrier.C_mtx_estim_fsu_current(cellUEindex,av_sp);   % current SINR (or capacity) is for computation of delivered data (or throughput), NOT for estimation
                    estimatedTHperPRB_mtx = objCarrier.C_mtx_estim_fsu(cellUEindex,av_sp);
                    alloc_mtx = objCarrier.allocationMatrix(cellUEindex,av_sp);
                    d_n_vect = sum((alloc_mtx.*experiencedTHperPRB_mtx).');
                    % obj.dataVar.UEbuffer(obj.usersList)>0
                    ueSchObj.presThVect(ueSchObj.lastActiveUEindex) = gammaPF .* ueSchObj.pastThVect(ueSchObj.lastActiveUEindex) ...
                        + (1 - gammaPF).* d_n_vect(ueSchObj.lastActiveUEindex);
                    
                    activeUEset = varUEbuffer(obj.usersList)>0;
                    
                    r_n_mtx = repmat(ueSchObj.presThVect(activeUEset).',1,n_prb);
                    r_n_mtx(r_n_mtx<=0) = 2^-1023;
                    den = (r_n_mtx).^ betaPF;
                    den(den<=1e-6) = 1e-6;
                    num = estimatedTHperPRB_mtx(activeUEset,:);
                    M = (num) .^ alphaPF ./ den;
                    M(num<=0) = 0;
                                                           
                    %%% STEP 2: find the best UE for each PRB 
                    
                    M_reshaped = M(:);
                    
%                     stopAllocFlag = false;
%                     SINR_mtx = obj.dataVar.SINR_mtx_estim_fsu(lou,av_sp);
                    C_mtx = TTIduration * objCarrier.PRBbw * objCarrier.C_mtx_estim_fsu(lou,av_sp);
                    
                    [~,M_sorted_ind] = sort(M_reshaped,'descend');
                    PRBallocationVect = -ones(1,n_prb);
                    UEallocationVect = false(1,nou);
                    availableUEvect = true(1,nou);
                    UEdataTxVect = zeros(1,nou); 
                    dataToTx_activeUE = varUEbuffer(lou);
                    UEperPRBallocMtx = false(nou,n_prb);
%                     C_mtx(UE_ind,PRB_ind);
                    
                    % comment here when MEX file is used (BEGINNING)

                    nAllocatedUEs = 0;
                    nAllocatedPRBs = 0;
                    nNA_UEs = 0;
                    next_ind = 1;
                    PF_loop_initFlag = true;
                    stopAllocFlag = false;
                    
                    while ~stopAllocFlag

                        m_ind = M_sorted_ind(next_ind);
                        UE_ind = rem(m_ind-1,nou)+1;
                        PRB_ind = ceil(m_ind/nou);

                        if PRBallocationVect(PRB_ind) < 0
                            % allocate on this PRB
                            if availableUEvect(UE_ind) == true

                                UEperPRBallocMtx(UE_ind,PRB_ind) = true;
                                UEdataTxVect(UE_ind) = UEdataTxVect(UE_ind) + C_mtx(UE_ind,PRB_ind);
                                if UEdataTxVect(UE_ind) >= dataToTx_activeUE(UE_ind)
                                    availableUEvect(UE_ind) = false;
                                    nNA_UEs = nNA_UEs + 1;
                                end
                                PRBallocationVect(PRB_ind) = UE_ind;
                                nAllocatedPRBs = nAllocatedPRBs + 1;
                                if UEallocationVect(UE_ind) == false
                                    nAllocatedUEs = nAllocatedUEs + 1;                
                                    UEallocationVect(UE_ind) = true;
                                end

                            end     

                        end

                        next_ind = next_ind + 1;
                            % don't allocate on this PRB
                            % don't allocate on this PRB

                        tmpPRB_flag = false;
%                         tmpUE_flag = false;
                        tmpMaxNumUE_flag = false;

                        if nAllocatedPRBs >= n_prb 
                            % all PRBs have been allocated
                            tmpPRB_flag = true;
                        end
                        if nAllocatedUEs == min_nOfUEs && PF_loop_initFlag
                            PF_loop_initFlag = false;
                            for n = 1:nou
                                if UEallocationVect(n) == false
                                    availableUEvect(n) = false;
                                end
                            end
                        end
                        if nNA_UEs >= min_nOfUEs
                            % all PRBs have been allocated
                            tmpMaxNumUE_flag = true;
                        end
                        stopAllocFlag = tmpPRB_flag || tmpMaxNumUE_flag;

                    end
                    
                    % comment here when MEX file is used (END)
                    
%                     [out1 out2] =  subRoutinePFmex(min_nOfUEs,nou,PRBallocationVect,availableUEvect+0,UEdataTxVect,C_mtx,dataToTx_activeUE,M_sorted_ind);

                    objCarrier.allocationMatrix(obj.usersList,:) = 0;   
                    
%                     UEallocationVect = out2>0;
%                     UEperPRBallocMtx = out1>0;
                    
                    
                    UE_to_beAllocated = 1:length(UEallocationVect);
                    UE_to_beAllocated = UE_to_beAllocated(UEallocationVect>0);
                    %%%%THIS IS THE OUTPUT OF THE FUNCTION - BEGIN
                    for n = 1:length(UE_to_beAllocated) 
                        userToAllocate = lou(UE_to_beAllocated(n)); 
                        prbIndex = UEperPRBallocMtx(UE_to_beAllocated(n),:);
                        PRBtoAllocate = av_sp( prbIndex); 
                        
                        % ueSchObj.userSpectVect(PRBtoAllocate) = userToAllocate; 
                        ueSchObj.userSpectVect(prbIndex) = userToAllocate; 
                        objCarrier.allocationMatrix(userToAllocate,PRBtoAllocate) = 1; 
                    end 
                    
                    %%%%THIS IS THE OUTPUT OF THE FUNCTION - END
                    
                    ueSchObj.pastThVect = ueSchObj.presThVect;                    
                    ueSchObj.lastActiveUEindex = varUEbuffer(obj.usersList)>0;

                else
                     objCarrier.allocationMatrix(obj.usersList,:) = 0;

                end
            end
        end         % end of function
        
        % PF scheduling
        function [ueSpectVect, allocMtx] = PropFairSchedMultiCarriers(obj,...
                varUEbuffer,objCarrier,UEperTTI,gammaPF,TTIduration,...
                objMacro,objPico,interfererIndexMatrix)
            
            numAvailPRB = objCarrier.numOfAvailablePRB;
            carrierNum = objCarrier.index;
            C_mtx_extim_fsu_current = objCarrier.C_mtx_estim_fsu_current;
            C_mtx_estim_fsu = objCarrier.C_mtx_estim_fsu;
            allocMtx = objCarrier.allocationMatrix;
            PRBbw = objCarrier.PRBbw;
            
            usObj = obj.ueSchedVars(carrierNum);
            raObj = obj.resAllocVariables(carrierNum);
            fsuObj = obj.fsuVars(carrierNum);
            
            alphaPF = fsuObj.alphaPF;
            betaPF = fsuObj.betaPF;
            
            if usObj.initPFVar > 0
                
                % the array usSpectVect stores UE index. Element
                % ueSpectVect(k) stores the UE index that has been assigned
                % the PRB k.
                
                %no UE has been scheduled yet
                % numAvailPRB is the total number of PRB corresponding to
                % the carrier <carrierNum>
                ueSpectVect = zeros(1,numAvailPRB);
                
                % array with the UE index that has to be served by the
                % current cell
                lou = obj.usersToBeServedVect;
                
                % PRBs of the carrier <carrierNum> that have been assigned
                % to the current cell
                
                availSpect = raObj.spectrumAvailable;
                
                % if there is some UE attached to the current cell and
                % there are spectrum allocated to the current cell
                if ~isempty(lou) && ~isempty(availSpect>0)
                    
                    %%% STEP 1: update past throughput
                    
                    % number of UEs that must be served 
                    nou = length(lou);
                    cellUEindex = obj.usersList;
                    
                    % available PRBs
                    av_sp = availSpect(availSpect>0);
                    n_prb = length(av_sp);
                    
                    % the number of UEs that must be scheduled is the
                    % minimum between the number of UEs that must be served
                    % by the cell and the input parameter UEperTTI that is
                    % the maximum number of UEs that must be scheduled per
                    % TTI
                
                    min_nOfUEs = min(UEperTTI,nou);
                    
                    % capacity that has been computed in the previous
                    % iteration. It is obtained the values for each pair 
                    % <UE that is attached to the cell; PRB that is 
                    % available to be assigned to the UE>
                    
                    % current SINR (or capacity) is for computation 
                    % of delivered data (or throughput), NOT for estimation
                    
                    experiencedTHperPRB_mtx = C_mtx_extim_fsu_current(cellUEindex,av_sp);   
                    
                    % estimatedTHperPRB_mtx = C_mtx_extim_fsu(cellUEindex,av_sp);
                    
                    % UE scheduling of previous iteration
                    alloc_mtx = allocMtx(cellUEindex,av_sp);
                    
                    % for each UE is obtained the sum of capacity
                    % considering the available PRBs. The capacity UE,PRB
                    % is added to the sum if PRB has been assigned to UE in
                    % the previous iteration
                    d_n_vect = sum((alloc_mtx .* experiencedTHperPRB_mtx).');
                    
                    usObj.presThVect(usObj.lastActiveUEindex) = gammaPF .* ...
                        usObj.pastThVect(usObj.lastActiveUEindex) + (1 - gammaPF) ...
                        .* d_n_vect(usObj.lastActiveUEindex);
                    
                    % UEs that are attached to the current cell and that
                    % have a buffer greater than zero.
                    
                    activeUEset = varUEbuffer(obj.usersList)>0;
                    
                    r_n_mtx = repmat(usObj.presThVect(activeUEset).',1,n_prb);
                    r_n_mtx(r_n_mtx<=0) = 2^-1023;
                    den = (r_n_mtx) .^ betaPF;
                    den(den<=1e-6) = 1e-6;
                    
                    % this is the main difference between single and
                    % multiple user scheduler
                    
                    
                    num = sumOfTHperPRB_v2(cls,carrierNum,numOfCarriers,...
                        mtxUserVsFreq, carrierArray);
                                       
                    %---------
                    
                    M = (num) .^ alphaPF ./ den;
                    M(num<=0) = 0;
                                                           
                    %%% STEP 2: find the best UE for each PRB 
                    
                    M_reshaped = M(:);
                    
                    C_mtx = TTIduration * PRBbw .* C_mtx_estim_fsu(lou,av_sp);
                    
                    [~,M_sorted_ind] = sort(M_reshaped,'descend');
                    PRBallocationVect = -ones(1,n_prb);
                    UEallocationVect = false(1,nou);
                    availableUEvect = true(1,nou);
                    UEdataTxVect = zeros(1,nou); 
                    dataToTx_activeUE = [ues(lou).UEbuffer];
                     
                    UEperPRBallocMtx = false(nou,n_prb);

                    nAllocatedUEs = 0;
                    nAllocatedPRBs = 0;
                    nNA_UEs = 0;
                    next_ind = 1;
                    PF_loop_initFlag = true;
                    stopAllocFlag = false;
                    
                    while ~stopAllocFlag

                        m_ind = M_sorted_ind(next_ind);
                        UE_ind = rem(m_ind-1,nou)+1;
                        PRB_ind = ceil(m_ind/nou);

                        if PRBallocationVect(PRB_ind) < 0
                            % allocate on this PRB
                            if availableUEvect(UE_ind) == true

                                UEperPRBallocMtx(UE_ind,PRB_ind) = true;
                                UEdataTxVect(UE_ind) = UEdataTxVect(UE_ind) + C_mtx(UE_ind,PRB_ind);
                                if UEdataTxVect(UE_ind) >= dataToTx_activeUE(UE_ind)
                                    availableUEvect(UE_ind) = false;
                                    nNA_UEs = nNA_UEs + 1;
                                end
                                PRBallocationVect(PRB_ind) = UE_ind;
                                nAllocatedPRBs = nAllocatedPRBs + 1;
                                if UEallocationVect(UE_ind) == false
                                    nAllocatedUEs = nAllocatedUEs + 1;                
                                    UEallocationVect(UE_ind) = true;
                                end

                            end     

                        end

                        next_ind = next_ind + 1;
                            % don't allocate on this PRB
                            % don't allocate on this PRB

                        tmpPRB_flag = false;
%                         tmpUE_flag = false;
                        tmpMaxNumUE_flag = false;

                        if nAllocatedPRBs >= n_prb 
                            % all PRBs have been allocated
                            tmpPRB_flag = true;
                        end
                        if nAllocatedUEs == min_nOfUEs && PF_loop_initFlag
                            PF_loop_initFlag = false;
                            for n = 1:nou
                                if UEallocationVect(n) == false
                                    availableUEvect(n) = false;
                                end
                            end
                        end
                        if nNA_UEs >= min_nOfUEs
                            % all PRBs have been allocated
                            tmpMaxNumUE_flag = true;
                        end
                        stopAllocFlag = tmpPRB_flag || tmpMaxNumUE_flag;

                    end
                    

                    allocMtx(obj.usersList,:) = 0;   
                    
                    UE_to_beAllocated = 1:length(UEallocationVect);
                    UE_to_beAllocated = UE_to_beAllocated(UEallocationVect>0);
                    %%%%THIS IS THE OUTPUT OF THE FUNCTION - BEGIN
                    for n = 1:length(UE_to_beAllocated) 
                        userToAllocate = lou(UE_to_beAllocated(n)); 
                        PRBtoAllocate = av_sp( UEperPRBallocMtx(UE_to_beAllocated(n),:) ); 
                        ueSpectVect(PRBtoAllocate) = userToAllocate; 
                        allocMtx(userToAllocate,PRBtoAllocate) = 1; 
                    end 
                    
                    %%%%THIS IS THE OUTPUT OF THE FUNCTION - END
                    
                    usObj.pastThVect = usObj.presThVect;
                    usObj.lastActiveUEindex = varUEbuffer(obj.usersList)>0;
                else
                    
                     allocMtx(obj.usersList,:) = 0;
                end
                
            else

                usObj.initPFVar = 1;
                obj.RoundRobinScheduler(objCarrier,varUEbuffer,UEperTTI);
                if obj.conf.pico.allocationType == 4
                    
                    fsuObj = obj.fsuVars(objCarrier.index);
                    fsuObj.initializeSchedulerAndFSUParamters(objCarrier,...
                        objMacro, objPico,interfererIndexMatrix);
        
                end
                usObj.lastActiveUEindex = varUEbuffer(obj.usersList)>0;                
                usObj.presThVect = zeros(size(obj.usersList));
                usObj.pastThVect = usObj.presThVect;
                
            end
        end         % end of function
        
        % PF scheduling
        function ProportionalFairScheduler_SFR(obj, objCarrier,UEperTTI,...
                TTIduration,alphaPF,betaPF,gammaPF,varUEbuffer)         
            
            carrierNum = objCarrier.index;
            ueSchObj = obj.ueSchedVars(carrierNum);
            
            if ueSchObj.initPFVar > 0
                raObj = obj.raVars(carrierNum);
                av_sp = raObj.spectrumAvailable;
               % varNumOfAvPRB = objCarrier.numOfAvailablePRB;
                
                % ueSchObj.userSpectVect = zeros(1,varNumOfAvPRB);
                ueSchObj.userSpectVect = zeros(1,length(av_sp));
                lou = obj.usersToBeServedVect;
                
                
                if ~isempty(lou) && ~isempty(av_sp>0)
                    
                    %%% STEP 1: update past throughput
                    nou = length(lou);                   
                    cellUEindex = obj.usersList;
                    av_sp = av_sp(av_sp>0);
                    PRBbandType_vect = raObj.PRBbandTypeVector(av_sp);
                    n_prb = length(av_sp);
                    min_nOfUEs = min(UEperTTI,nou);                  
                    
                    % current SINR (or capacity) is for computation of delivered data (or throughput), not for estimation                                        
                    experiencedTHperPRB_mtx = objCarrier.C_mtx_estim_fsu_current(cellUEindex,av_sp); 
                    
                    % estimated SINR (or capacity) is for estimating 
                    estimatedTHperPRB_mtx = objCarrier.C_mtx_estim_fsu(cellUEindex,av_sp);                  
                    alloc_mtx = objCarrier.allocationMatrix(cellUEindex,av_sp);
                    d_n_vect = sum((alloc_mtx.*experiencedTHperPRB_mtx).');
                    
                    ueSchObj.presThVect(ueSchObj.lastActiveUEindex) = gammaPF .* ueSchObj.pastThVect(ueSchObj.lastActiveUEindex) + ...
                        (1- gammaPF).* d_n_vect(ueSchObj.lastActiveUEindex);
                    
                    activeUEset = varUEbuffer(obj.usersList)>0;
                    r_n_mtx = repmat(ueSchObj.presThVect(activeUEset).',1,n_prb);
                    r_n_mtx(r_n_mtx<=0) = 2^-1023;
                    den = (r_n_mtx) .^ betaPF;
                    den(den<=1e-6) = 1e-6;
                    num = estimatedTHperPRB_mtx(activeUEset,:);
                    M = (num) .^ alphaPF ./ den;
                    M(num<=0) = 0;         % UEs x PRBs
                                        
                    
                    if sum(activeUEset) > 1
                        
                        %%% STEP 2: Allocation of UEs in the high-power band
                    
                        M_highPWband = M(:,PRBbandType_vect == 2);
                        M_reshaped = M_highPWband(:);                  

                        stopAllocFlag = false;
    %                     SINR_mtx = obj.dataVar.SINR_mtx_estim_fsu(lou,av_sp);
                        C_mtx_hiPW = TTIduration * objCarrier.PRBbw ...
                            .* objCarrier.C_mtx_estim_fsu(lou,av_sp(PRBbandType_vect == 2));
                        n_prb_hiPW = sum(PRBbandType_vect == 2);

                        [~,M_sorted_ind] = sort(M_reshaped,'descend');
                        PRBallocationVect_hiPW = -ones(1,n_prb_hiPW);
                        UEallocationVect = false(1,nou);
                        availableUEvect = true(1,nou);
                        UEdataTxVect = zeros(1,nou); 
                        
                        dataToTx_activeUE = [varUEbuffer(lou)];
                        
                        UEperPRBallocMtx_hiPW = false(nou,n_prb_hiPW);
    %                     C_mtx(UE_ind,PRB_ind);
    
                        nAllocatedUEs = 0;
                        nAllocatedPRBs_hiPW = 0;
                        nNA_UEs = 0;
                        next_ind = 1;
                        PF_loop_initFlag = true;
                        max_nOfUEs_highPWband = ceil(min_nOfUEs * 0.3);

                        while ~stopAllocFlag

                            m_ind = M_sorted_ind(next_ind);
                            UE_ind = rem(m_ind-1,nou)+1;
                            PRB_ind = ceil(m_ind/nou);

                            if PRBallocationVect_hiPW(PRB_ind) < 0
                                % allocate on this PRB
                                if availableUEvect(UE_ind) == true

                                    UEperPRBallocMtx_hiPW(UE_ind,PRB_ind) = true;
                                    UEdataTxVect(UE_ind) = UEdataTxVect(UE_ind) + C_mtx_hiPW(UE_ind,PRB_ind);
                                    if UEdataTxVect(UE_ind) >= dataToTx_activeUE(UE_ind)
                                        availableUEvect(UE_ind) = false;
                                        nNA_UEs = nNA_UEs + 1;
                                    end
                                    PRBallocationVect_hiPW(PRB_ind) = UE_ind;
                                    nAllocatedPRBs_hiPW = nAllocatedPRBs_hiPW + 1;
                                    if UEallocationVect(UE_ind) == false
                                        nAllocatedUEs = nAllocatedUEs + 1;                
                                        UEallocationVect(UE_ind) = true;
                                    end

                                end     

                            end

                            next_ind = next_ind + 1;
                                % don't allocate on this PRB
                                % don't allocate on this PRB

                            tmpPRB_flag = false;
    %                         tmpUE_flag = false;
                            tmpMaxNumUE_flag = false;

                            if nAllocatedPRBs_hiPW >= n_prb_hiPW 
                                % all PRBs have been allocated
                                tmpPRB_flag = true;
                            end
                            if nAllocatedUEs == max_nOfUEs_highPWband && PF_loop_initFlag
                                PF_loop_initFlag = false;
                                for n = 1:nou
                                    if UEallocationVect(n) == false
                                        availableUEvect(n) = false;
                                    end
                                end
                            end
                            if nNA_UEs >= max_nOfUEs_highPWband
                                % all PRBs have been allocated
                                tmpMaxNumUE_flag = true;
                            end
                            stopAllocFlag = tmpPRB_flag || tmpMaxNumUE_flag;

                        end

                        %%% STEP 3: Allocation of UEs in the low-power band

                        M_lowPWband = M(:,PRBbandType_vect == 1);
                        M_reshaped = M_lowPWband(:);                  

                        stopAllocFlag = false;
                        C_mtx_lowPW = TTIduration * objCarrier.PRBbw ...
                            .* objCarrier.C_mtx_estim_fsu(lou,av_sp(PRBbandType_vect == 1));
                        n_prb_lowPW = sum(PRBbandType_vect == 1);

                        [~,M_sorted_ind] = sort(M_reshaped,'descend');
                        PRBallocationVect_lowPW = -ones(1,n_prb_lowPW);
                        availableUEvect = true(1,nou);
                        availableUEvect(UEallocationVect) = false;
                        nNA_UEs = nAllocatedUEs;
                        UEperPRBallocMtx_lowPW = false(nou,n_prb_lowPW);
                        nAllocatedPRBs_lowPW = 0;
                        next_ind = 1;
                        PF_loop_initFlag = true;

                        while ~stopAllocFlag

                            m_ind = M_sorted_ind(next_ind);
                            UE_ind = rem(m_ind-1,nou)+1;
                            PRB_ind = ceil(m_ind/nou);

                            if PRBallocationVect_lowPW(PRB_ind) < 0
                                % allocate on this PRB
                                if availableUEvect(UE_ind) == true

                                    UEperPRBallocMtx_lowPW(UE_ind,PRB_ind) = true;
                                    UEdataTxVect(UE_ind) = UEdataTxVect(UE_ind) + C_mtx_lowPW(UE_ind,PRB_ind);
                                    if UEdataTxVect(UE_ind) >= dataToTx_activeUE(UE_ind)
                                        availableUEvect(UE_ind) = false;
                                        nNA_UEs = nNA_UEs + 1;
                                    end
                                    PRBallocationVect_lowPW(PRB_ind) = UE_ind;
                                    nAllocatedPRBs_lowPW = nAllocatedPRBs_lowPW + 1;
                                    if UEallocationVect(UE_ind) == false
                                        nAllocatedUEs = nAllocatedUEs + 1;                
                                        UEallocationVect(UE_ind) = true;
                                    end

                                end     

                            end

                            next_ind = next_ind + 1;
                                % don't allocate on this PRB
                                % don't allocate on this PRB

                            tmpPRB_flag = false;
    %                         tmpUE_flag = false;
                            tmpMaxNumUE_flag = false;

                            if nAllocatedPRBs_lowPW >= n_prb_lowPW 
                                % all PRBs have been allocated
                                tmpPRB_flag = true;
                            end
                            if nAllocatedUEs == min_nOfUEs && PF_loop_initFlag
                                PF_loop_initFlag = false;
                                for n = 1:nou
                                    if UEallocationVect(n) == false
                                        availableUEvect(n) = false;
                                    end
                                end
                            end
                            if nNA_UEs >= min_nOfUEs
                                % all PRBs have been allocated
                                tmpMaxNumUE_flag = true;
                            end
                            stopAllocFlag = tmpPRB_flag || tmpMaxNumUE_flag;

                        end     % end of while

                        UEperPRBallocMtx = false(nou,n_prb);
                        UEperPRBallocMtx(:,PRBbandType_vect == 2) = UEperPRBallocMtx_hiPW;
                        UEperPRBallocMtx(:,PRBbandType_vect == 1) = UEperPRBallocMtx_lowPW;   
                        
                    else    % only 1 UE to be allocated 
                        
                        %%% STEP 2: Allocation of UEs in the high-power band
                    
                        M_highPWband = M(:,PRBbandType_vect == 2);
                        M_reshaped = M_highPWband(:);                  

                        stopAllocFlag = false;
                        C_mtx_hiPW = TTIduration * objCarrier.PRBbw ...
                            .* C_mtx_estim_fsu(lou,av_sp(PRBbandType_vect == 2));
                        n_prb_hiPW = sum(PRBbandType_vect == 2);

                        [~,M_sorted_ind] = sort(M_reshaped,'descend');
                        PRBallocationVect_hiPW = -ones(1,n_prb_hiPW);
                        UEallocationVect = false(1,nou);
                        availableUEvect = true(1,nou);
                        UEdataTxVect = zeros(1,nou); 
                        dataToTx_activeUE = [varUEbuffer(lou)];
                        UEperPRBallocMtx_hiPW = false(nou,n_prb_hiPW);
    %                     C_mtx(UE_ind,PRB_ind);
                        nAllocatedUEs = 0;
                        nAllocatedPRBs_hiPW = 0;
                        nNA_UEs = 0;
                        next_ind = 1;
                        PF_loop_initFlag = true;

                        while ~stopAllocFlag

                            m_ind = M_sorted_ind(next_ind);
                            UE_ind = rem(m_ind-1,nou)+1;
                            PRB_ind = ceil(m_ind/nou);

                            if PRBallocationVect_hiPW(PRB_ind) < 0
                                % allocate on this PRB
                                if availableUEvect(nou) == true

                                    UEperPRBallocMtx_hiPW(UE_ind,PRB_ind) = true;
                                    UEdataTxVect(UE_ind) = UEdataTxVect(UE_ind) + C_mtx_hiPW(UE_ind,PRB_ind);
                                    if UEdataTxVect(UE_ind) >= dataToTx_activeUE(UE_ind)
                                        availableUEvect(UE_ind) = false;
                                        nNA_UEs = nNA_UEs + 1;
                                    end
                                    PRBallocationVect_hiPW(PRB_ind) = UE_ind;
                                    nAllocatedPRBs_hiPW = nAllocatedPRBs_hiPW + 1;
                                    if UEallocationVect(UE_ind) == false
                                        nAllocatedUEs = nAllocatedUEs + 1;                
                                        UEallocationVect(UE_ind) = true;
                                    end

                                end     

                            end

                            next_ind = next_ind + 1;
                                % don't allocate on this PRB
                                % don't allocate on this PRB

                            tmpPRB_flag = false;
    %                         tmpUE_flag = false;
                            tmpMaxNumUE_flag = false;

                            if nAllocatedPRBs_hiPW >= n_prb_hiPW 
                                % all PRBs have been allocated
                                tmpPRB_flag = true;
                            end
                            if nAllocatedUEs == min_nOfUEs && PF_loop_initFlag
                                PF_loop_initFlag = false;
                                for n = 1:nou
                                    if UEallocationVect(n) == false
                                        availableUEvect(n) = false;
                                    end
                                end
                            end
                            if nNA_UEs >= min_nOfUEs
                                % all PRBs have been allocated
                                tmpMaxNumUE_flag = true;
                            end
                            stopAllocFlag = tmpPRB_flag || tmpMaxNumUE_flag;

                        end

                        UEperPRBallocMtx = false(nou,n_prb);
                        UEperPRBallocMtx(:,PRBbandType_vect == 2) = UEperPRBallocMtx_hiPW;

                        
                    end
                    
                    allocationMtx(obj.usersList,:) = 0;   

                    UE_to_beAllocated = 1:length(UEallocationVect);
                    UE_to_beAllocated = UE_to_beAllocated(UEallocationVect>0);
                    for n = 1:length(UE_to_beAllocated) 
                        userToAllocate = lou(UE_to_beAllocated(n)); 
                        prbInd = UEperPRBallocMtx(UE_to_beAllocated(n),:);
                        PRBtoAllocate = av_sp( prbInd); 
                        % ueSchObj.userSpectVect(PRBtoAllocate) = userToAllocate; 
                        ueSchObj.userSpectVect(prbInd) = userToAllocate; 
                        allocationMtx(userToAllocate,PRBtoAllocate) = 1; 
                    end 

                   
                    ueSchObj.pastThVect = ueSchObj.presThVect;
                    ueSchObj.lastActiveUEindex = varUEbuffer(obj.usersList)>0;

                else
                    objCarrier.allocationMatrix(obj.usersList,:) = 0;
                end
            end
        end         % end of function

        % RRschduling
        function RoundRobinScheduler(obj,objCarrier,varUEbuffer,UEperTTI)

            carrierNum = objCarrier.index;
            
            raObj = obj.raVars(carrierNum);
            ueSchObj = obj.ueSchedVars(carrierNum);
            
            numOfAvPRBs = objCarrier.numOfAvailablePRB;
            
            av_sp = raObj.spectrumAvailable;
            
            % userSpectVect is a property of carrier_Cls
            % size = 1 x Num of PRBs in the carrier
            
            % ueSchObj.userSpectVect = zeros(1,numOfAvPRBs);
            ueSchObj.userSpectVect = zeros(1,length(av_sp));
            
            % allocationMtx is a property of carrier_Cls
            % size = total number of UEs in the network x Num of PRBs in
            % the carrier
            
            objCarrier.allocationMatrix(obj.usersList,1:numOfAvPRBs) = 0;
            
            varNextUserInd = obj.nextUserIndex;
            
            
            
            lou = obj.usersToBeServedVect;
            
            if ~isempty(lou) && ~isempty(av_sp>0)
               
                nou = length(lou);
                
                av_sp = av_sp(av_sp>0);
                n_prb = length(av_sp);
                
                if isempty(objCarrier.SINR_mtx_estim)
                    SINR_mtx = 0.1 * ones(nou,n_prb);
                else
                    SINR_mtx = objCarrier.SINR_mtx_estim(lou,av_sp);
                end
                
                
                n_served_ue = min(nou,UEperTTI);
                
                UE_index_vect = 1:nou;
                UEcnt = 0;
                PRBcnt = 0;
                
                if varUEbuffer(varNextUserInd) > 0 
                    nextUEindexVect = UE_index_vect(lou >= varNextUserInd);
                    nextUEindex = nextUEindexVect(1);
                else
                    nextUEindexVect = UE_index_vect(lou > varNextUserInd);
                    if ~isempty(nextUEindexVect)
                        nextUEindex = nextUEindexVect(1);
                    else
                        nextUEindex = 1;
                    end
                    
                end
                
                max_nPRBs_perUE = ceil(n_prb/n_served_ue);
                
                % COMMENT: carlo SIMULATOR COMPARISON
                PRB_scrambledSeq = randSequenceFromBin(n_prb,n_prb);
                
                % PRB_scrambledSeq = obj.random_prbs;
                 
                    
                while ( UEcnt< UEperTTI && PRBcnt <n_prb)
                    
                    currentUEindex_reindexing = rem(UEcnt-1+nextUEindex,nou)+1;
                    candidatePRBs = PRB_scrambledSeq( min(PRBcnt+1,n_prb) : min(PRBcnt+max_nPRBs_perUE,n_prb));
                    estimatedDeliveredTraffic = shannonMkI(SINR_mtx(currentUEindex_reindexing,candidatePRBs));
                    cumulativeDeliveredTraffic = cumsum(estimatedDeliveredTraffic);
                    tmp_vect = 1:max_nPRBs_perUE;
                    
                    varDataBuffer = varUEbuffer(lou(currentUEindex_reindexing));
                                       
                    tmp_neededPRBs = tmp_vect(cumulativeDeliveredTraffic > varDataBuffer);
                    
                    if isempty(tmp_neededPRBs)
                        prbInd = candidatePRBs;
                        % PRBtoAllocate = av_sp( candidatePRBs );
                    else
                        PRBs_reindex = 1:tmp_neededPRBs(1);
                        prbInd = candidatePRBs(PRBs_reindex);
                        % PRBtoAllocate = av_sp(prbInd);
                    end
                    
                    PRBtoAllocate = av_sp(prbInd);
                    UEcnt = UEcnt + 1;
                    PRBcnt = PRBcnt + length(PRBtoAllocate);
                    % ueSchObj.userSpectVect( PRBtoAllocate ) = lou(currentUEindex_reindexing);
                    ueSchObj.userSpectVect( prbInd ) = lou(currentUEindex_reindexing);
                    objCarrier.allocationMatrix( lou(currentUEindex_reindexing) , PRBtoAllocate ) = 1;
                end
                %%% next UE to be served (at the next scheduling round)   
                
                obj.nextUserIndex = lou(rem(UEcnt-1+nextUEindex,nou)+1);
               
           
            end
        end     % end of function
        
        function FrequencySelectiveScheduler(obj,objCarrier,UEperTTI,varUEbuffer)               % RRschduling
        
            carrierNum = objCarrier.index;
            
            ueSchObj = obj.ueSchedVars(carrierNum);
            raObj = obj.raVars(carrierNum);
            
            if ueSchObj.initPFVar >0

                av_sp = raObj.spectrumAvailable;
                
                n_tot_PRBs = objCarrier.numOfAvailablePRB;
                
                % ueSchObj.userSpectVect = zeros(1,n_tot_PRBs);
                ueSchObj.userSpectVect = zeros(1,length(av_sp));
                
                objCarrier.allocationMatrix(obj.usersList,:) = 0;

                n_init = obj.nextUserIndex;
                
                lou = obj.usersToBeServedVect;
                
                
                if ~isempty(lou) && ~isempty(av_sp > 0)

                    lou = obj.usersToBeServedVect;
                    nou = length(lou);
                    n_uePerTTI = UEperTTI;
                    av_sp = av_sp(av_sp>0);
                    n_prb = length(av_sp);
                    min_n_ues = min(nou, n_uePerTTI);
                    n_served_ue = min([min_n_ues n_prb]);
                    %%% use data.SINR_mtx_estim for estimated SINR, since
                    %%% this has SINR only for PRBs in use by the cell
                    SINR_mtx = objCarrier.SINR_mtx_estim(lou,av_sp);

                    users_to_be_served = rem(n_init+(1:n_served_ue)-1,nou)+1;
                    
                    varRandUes = randSequenceFromBin(n_served_ue,n_served_ue);
                    % varRandUes = obj.randUesRR;
                    scrambledUEs = users_to_be_served(varRandUes);
                    prbPerUE = fix(n_prb/min_n_ues); 
                    addPRB = rem(n_prb,min_n_ues);
                    SINR_mtx_small = SINR_mtx(scrambledUEs,:);

                    %%%%THIS IS THE OUTPUT OF THE FUNCTION - BEGIN
                    for n = 1:length(scrambledUEs)    

                        if n <= addPRB
                            totPRBperUE = prbPerUE + 1;
                        else
                            totPRBperUE = prbPerUE;
                        end    
                        userToAllocate = lou(scrambledUEs(n));
                        [~, PRBind] = sort(SINR_mtx_small(n,:),2,'descend');
                        prbInd = PRBind(1:totPRBperUE);
                        
                        PRBtoAllocate = av_sp(prbInd);
                        
                        % ueSchObj.userSpectVect(PRBtoAllocate) = userToAllocate;
                        ueSchObj.userSpectVect(prbInd) = userToAllocate;
                        objCarrier.allocationMatrix(userToAllocate,PRBtoAllocate) = 1;
                        SINR_mtx_small(:,PRBind(1:totPRBperUE)) = -1;                                            
                    end
                    obj.nextUserIndex = rem( obj.nextUserIndex + n_served_ue - 1,nou)+1;       
                    raObj.pwPerPRB = raObj.maxTxPw / n_tot_PRBs;
                    %%%%THIS IS THE OUTPUT OF THE FUNCTION - END
                else                    
                     raObj.pwPerPRBvector = zeros(1,n_tot_PRBs);
                end
            else
                    ueSchObj.initPFVar = 1;
                    obj.RoundRobinScheduler(objCarrier,varUEbuffer,UEperTTI);
            end
        end     % end of function        
        
        
        function fixedSpectrumAllocation(obj, objCarrier, ...
                softReusePWoffset,reuseScheme)
            
            carrierNum = objCarrier.index;
            
            raObj = obj.raVars(carrierNum);
            
            n_tot_PRBs = objCarrier.numOfAvailablePRB;
            raObj.cellsVsPRB_allocationMatrix = zeros(1,n_tot_PRBs);
            
            varPwPerPRB = raObj.maxTxPw / n_tot_PRBs;
            pwNotUsedPRBs = zeros(1,n_tot_PRBs);
            defaultType = 1;
            defaultTypeNotUsed = zeros(1,n_tot_PRBs);
            availableSpectrum =  1:n_tot_PRBs;
            
            if reuseScheme == 4 
                pwNotUsedPRBs = varPwPerPRB * softReusePWoffset;
                defaultType = 2;
                defaultTypeNotUsed = ones(1,n_tot_PRBs);
            end
                        
            if (reuseScheme == 2 || reuseScheme==4) && obj.isMacro
                varPortion = ceil(n_tot_PRBs/3);
                varSector = obj.macroCellSectorIndex;
                availableSpectrum = (varSector-1)*varPortion+1:...
                    min(n_tot_PRBs,varSector*varPortion);
            elseif reuseScheme==3 
                varPortion = floor(n_tot_PRBs/2);
                if obj.isMacro
                    availableSpectrum = 1:varPortion;
                else
                    availableSpectrum = varPortion+1:n_tot_PRBs;
                end
            end
            
            raObj.spectrumAvailable = availableSpectrum;
            raObj.pwPerPRBvector = pwNotUsedPRBs;
            raObj.pwPerPRBvector(availableSpectrum) = varPwPerPRB;
            raObj.PRBbandTypeVector = defaultTypeNotUsed;
            raObj.PRBbandTypeVector(availableSpectrum) = defaultType;
            raObj.pwPerPRBAllowedVector = raObj.pwPerPRBvector;
            
            if reuseScheme == 4
                raObj.spectrumAvailable = 1:n_tot_PRBs; 
            end
                      
            raObj.cellsVsPRB_allocationMatrix(raObj.spectrumAvailable) = 1;
            
        end       % end of function
        
        
        function manualPRBsetting(obj,carrierNum,numOfAvailablePRB,manualPRBvector)
           
            raObj = obj.raVars(carrierNum);
            
            varPwPerPRB = raObj.maxTxPw / numOfAvailablePRB;
            raObj.spectrumAvailable = manualPRBvector;
            raObj.pwPerPRBvector = zeros(1,numOfAvailablePRB);
            raObj.pwPerPRBvector(manualPRBvector) = varPwPerPRB .* ones(size(manualPRBvector));                                
            raObj.pwPerPRBAllowedVector = varPwPerPRB .* ones(1,numOfAvailablePRB);
            
        end         % end of function
        
        function frequencyALOHA(obj,carrierNum,numOfAvailablePRB,...
                RBs_percentage)
            
            raObj = obj.raVars(carrierNum);
            
            varPwPerPRB = raObj.maxTxPw/numOfAvailablePRB;
            
            raObj.spectrumAvailable = sort(randSequenceFromBin(ceil(RBs_percentage*numOfAvailablePRB/100),numOfAvailablePRB));
            raObj.pwPerPRBvector = zeros(1,numOfAvailablePRB);
            raObj.pwPerPRBvector(raObj.spectrumAvailable) = varPwPerPRB .* ones(size(raObj.spectrumAvailable));                                
            raObj.pwPerPRBAllowedVector = varPwPerPRB .* ones(1,numOfAvailablePRB);
            
        end
        
    end     % end of methods

end         % end of class


