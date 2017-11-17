classdef fsuVariables_Cls < handle
    
    properties
        alphaPF; % scalar
        AV_PRBvector;
        betaPF; % scalar
        BSstate;
        
        C_esteem_cls_mtx;
        
        cumUserTh; % cumulative throughput for each user that is 
        
        % attached to the cell.
             
        HIIboolean; % boolean variable: true/false
        HIIcounterVector;
        HIItimeStamp;
        
        indexParamPF; % integer
        initFSUvar; % scalar
        
        lastSignalledPRBs; % array of PRBs
        
        lastTH_measure; % complex number 
        macroBS_HIItx_indexVector;
        macroHIIsignalStruct; % struct with two fields: list and timeStamp
        macroToPicoInterfererIndexMatrix;
        minSpectrum_pBS;
        oldNofUsedPRBs; % scalar
        
        pico_NA_PRB_timer;
        picoInterfererVector;
        
        picoTimerNew;        
        PRB_timeStamp;
        PRBsHIIStateVector; % array of numbers. The number of elements 
        % corresponds to the number of PRBs that form the available 
        % spectrum for the cell.
        prevQoSevState;
        
        QoSevState; % integer
        resendHIItimer;
        syncFSUtimeSlot; % scalar
        
        timerQoSEV; % integer
        
    end
    
    methods
        function obj = fsuVariables_Cls(varTimerQoSEV,varNumOfTTIperSimu,...
                alphaPF_type,betaPF_type)
            
            obj.HIIboolean = false;
            obj.macroHIIsignalStruct.timestamp = 0;
            obj.macroHIIsignalStruct.list = [];
            obj.initFSUvar = false;
            obj.QoSevState = 0;
            obj.C_esteem_cls_mtx = [];
            obj.prevQoSevState = 0;
            obj.timerQoSEV = varTimerQoSEV;
            obj.HIItimeStamp = 0;
            obj.BSstate = [];
            if varTimerQoSEV > 0
                obj.BSstate = zeros(5,ceil(varNumOfTTIperSimu/varTimerQoSEV));
            end
            obj.minSpectrum_pBS = 15;
            
            obj.alphaPF = alphaPF_type;
            obj.betaPF = betaPF_type;
        end
        
        
        function updateAndSendHII(obj,av_sp,schedCfg,resendHIItimer,currTTI)
            
            global objCell;
            
            obj.HIIboolean = false;
            obj.macroHIIsignalStruct.list = [];
                       
            if obj.timerQoSEV <= 0
                
                estimatedUTHvector = obj.cumUserTh/schedCfg.timerQoSEV;
                AUTH = mean(estimatedUTHvector);
                QoS = CDFfncMkII(estimatedUTHvector,5);
                obj.timerQoSEV = schedCfg.timerQoSEV;
                obj.cumUserTh = zeros(1,length(objCell.usersList));
                 
                %%% state evaluation
                
                if AUTH > schedCfg.alphaMinAUTH 
                    
                    if QoS > schedCfg.alphaMinQoS
                        obj.QoSevState = 1;
                    else
                        obj.QoSevState = 4;
                    end
                else
                    
                    if QoS > schedCfg.alphaMinQoS
                        obj.QoSevState = 2;
                    else
                        deltaTH_perc = abs(obj.lastTH_measure - (AUTH+1i*QoS)) / abs(obj.lastTH_measure);
                        if obj.QoSevState == 3 && deltaTH_perc < 0.05
                            obj.QoSevState = 5;
                        else
                            obj.QoSevState = 3;
                        end
                    end
                end
                
                %%% change PF paramter   
                
                switch obj.QoSevState
                    
                    case 1
                       obj.PRBsHIIStateVector(av_sp) = 0;
                    case 2
                        %%% decrease PF constant
                        obj.indexParamPF = min( obj.indexParamPF+1 , length(schedCfg.PFalphaSet));
                        obj.alphaPF = schedCfg.PFalphaSet( obj.indexParamPF );
                        obj.betaPF = schedCfg.PFbetaSet( obj.indexParamPF );
                    case 3
                        %%% send HII 
                        PRBvector_1stHII = obj.PRBsHIIStateVector == 0;
                        PRBvector_noMoreHII = obj.PRBsHIIStateVector(av_sp) == 2;
                        if any(PRBvector_1stHII)
                            %%% find UEs
                            [~, UEind] = sort(estimatedUTHvector());
                            worstUEsInd = UEind(1:ceil(objCell.loadNOU*0.05));
                            %%% find PRBs

                            UEperPRBmat = obj.C_esteem_cls_mtx(worstUEsInd,av_sp);
                            UTHmetric = mean(UEperPRBmat,1);
                            PRBvector_atLeastOneHIIindex = obj.PRBsHIIStateVector(av_sp) > 0;
                            UTHmetric(PRBvector_atLeastOneHIIindex) = -1;
                            %%% the BEST PRBs of the WORST USERs are chosen
                            %%% for sending the HII
                            [PRBval, PRBind] = sort(UTHmetric,'descend');
                            PRBforHII = PRBind( 1:4 );
                            PRBforHIIval = PRBval(1:4);
                            PRBforHII = PRBforHII(PRBforHIIval>-1);
                            PRBforHII = av_sp( PRBforHII );
                            % send HII for PRB with index PRBforHII
                            obj.PRBsHIIStateVector(PRBforHII) = 1;
                            obj.lastSignalledPRBs = PRBforHII;
                            obj.HIIboolean  = true;
                            obj.macroHIIsignalStruct.list = PRBforHII;
                            obj.HIItimeStamp = rem(obj.HIItimeStamp +1,1000);
                            obj.macroHIIsignalStruct.timestamp = objCell.cellID * 10000 + obj.HIItimeStamp;
                        elseif all(PRBvector_noMoreHII)
                            if obj.resendHIItimer < 0 
                                obj.resendHIItimer = resendHIItimer;
                            end
                        end
                        
                    case 4
                        %%% increase PF constant
                        obj.indexParamPF = max( obj.indexParamPF-1 , 1);
                        obj.alphaPF = schedCfg.PFalphaSet( obj.indexParamPF );
                        obj.betaPF = schedCfg.PFbetaSet( obj.indexParamPF );
                        
                    case 5
                         
                        PRBvector_2ndHII = obj.PRBsHIIStateVector == 1;
                        if any(PRBvector_2ndHII)
                            % sendHII for PRBs
                            obj.HIIboolean = true;
                            obj.macroHIIsignalStruct.list = obj.lastSignalledPRBs;
                            obj.macroHIIsignalStruct.timestamp = objCell.cellID*10000 + obj.HIItimeStamp;
                            obj.PRBsHIIStateVector(obj.lastSignalledPRBs) = 2;
                        end
                end
                
                obj.BSstate(:,sum(obj.BSstate(1,:)>0)+1) = [obj.QoSevState currTTI AUTH QoS obj.indexParamPF].';    
                obj.prevQoSevState = obj.QoSevState;
                obj.lastTH_measure = AUTH+1i*QoS;
            
            end
            
            if obj.resendHIItimer == 0
                
                obj.PRBsHIIStateVector(av_sp) = 0;
                
            end
            
        end         % end of function

        
        
        % varThPerUE: throughput per UE in the current carrier
        % usersList: list of UEs that are attached to the current cell
        function measureUserTH(obj,varThPerUE)
            global objCell;
            
            if isempty(objCell.usersList)
                obj.cumUserTh = [];
            else
                if isempty(varThPerUE)
                    obj.cumUserTh = zeros(1,length(objCell.usersList));
                else
                    obj.cumUserTh = obj.cumUserTh + varThPerUE(objCell.usersList);
                end             
            end         
        end 
        
        function FSUalgorithm_MkI(obj,objCarrier,currTTI,fsuCfg,FSU_SINT_th)
            
            global objCell;
            
            carrierNum = objCarrier.index;
            
            numOfAvPRB = objCarrier.numOfAvailablePRB;
            
            raObj = objCell.raVars(carrierNum);
            
            t = rem(currTTI - 1, fsuCfg.numSyncSlots)+1;

            if ~isempty(objCell.usersToBeServedVect)

                if obj.initFSUvar <= 0

                    raObj.cellsVsPRB_allocationMatrix = zeros(1, numOfAvPRB);
                    
                    varMin = min(fsuCfg.initPRBoccupancy, numOfAvPRB);
                    
                    raObj.spectrumAvailable = randSequenceFromBin(varMin, numOfAvPRB);
                    obj.oldNofUsedPRBs = length(raObj.spectrumAvailable);
                    obj.initFSUvar = 1;
                    
                    raObj.pwPerPRBvector = zeros(1, numOfAvPRB);
                    raObj.pwPerPRBvector(raObj.spectrumAvailable) = (raObj.maxTxPw / numOfAvPRB);
                    raObj.pwPerPRBAllowedVector = (raObj.maxTxPw / numOfAvPRB) .* ones(1, numOfAvPRB);


                else 
                    
                    if  t == obj.syncFSUtimeSlot

                        raObj.cellsVsPRB_allocationMatrix = zeros(1, numOfAvPRB);
                        lou = objCell.usersToBeServedVect;
                        SINR_mtx = (length(raObj.spectrumAvailable)/obj.oldNofUsedPRBs) .* ...
                            objCarrier.SINR_mtx_estim_fsu(lou,:);
                        
                        SINR_mtx_dB = 10*log10(SINR_mtx);

                        if size(SINR_mtx_dB,1) == 1
                            meanSINR = SINR_mtx_dB;
                        else
                            meanSINR = mean(SINR_mtx_dB);
                        end

                        vectPRB = 1:numOfAvPRB ;
                        tmpSpectrum = vectPRB(meanSINR > FSU_SINT_th);

                        newLen = length(tmpSpectrum);
                        deltaLength = newLen - obj.oldNofUsedPRBs;

                        if deltaLength > fsuCfg.deltaN_PRB 

                            newLength = obj.oldNofUsedPRBs + fsuCfg.deltaN_PRB  ;
                            if newLength < fsuCfg.minPRBperCell
                                newLength = fsuCfg.minPRBperCell;
                            end
                            unsortedSINR = meanSINR;
                            [~, sortedPRBindex] = sort(unsortedSINR,2,'descend');
                            newSpectrum = sort(sortedPRBindex(1:newLength));


                        elseif deltaLength < -fsuCfg.deltaN_PRB  

                            newLength = obj.oldNofUsedPRBs - fsuCfg.deltaN_PRB;
                            if newLength < fsuCfg.minPRBperCell
                                newLength = fsuCfg.minPRBperCell;
                            end
                            unsortedSINR = meanSINR;
                            [~, sortedPRBindex] = sort(unsortedSINR,2,'descend');
                            if newLength > length(sortedPRBindex)
                               disp('ciao o bel'); 
                            end
                            newSpectrum = sort(sortedPRBindex(1:newLength));

                        else

                            newLength = length(tmpSpectrum);
                            if newLength < fsuCfg.minPRBperCell
                                newLength = fsuCfg.minPRBperCell;
                                unsortedSINR = meanSINR;
                                [~, sortedPRBindex] = sort(unsortedSINR,2,'descend');
                                newSpectrum = sort(sortedPRBindex(1:newLength));
                            else
                                newSpectrum = tmpSpectrum;
                            end


                        end

                        obj.oldNofUsedPRBs = length(raObj.spectrumAvailable);
                        raObj.spectrumAvailable = newSpectrum;
                        raObj.cellsVsPRB_allocationMatrix(raObj.spectrumAvailable) = 1;

                    end     %   end of if t == obj.syncFSUtimeSlot

                end         % end of if obj.initFSUvar <= 0

            end     %  end of ~isempty(obj.usersToBeServedVect)        
            
        end         % end of function

        
        function FSUalgorithm_MkII(obj,objCarrier,fsuCfg)
            
            global objCell;
            global cls;

            carrierNum = objCarrier.index;
            
            raObj = objCell.raVars(carrierNum);
            
            numOfAvPRBs = objCarrier.numOfAvailablePRB;
            
            if ~isempty(objCell.usersToBeServedVect)

                if ~obj.initFSUvar      %%% intialization of FSU algorithm

                    raObj.cellsVsPRB_allocationMatrix  = zeros(1,numOfAvPRBs);
                    varMin = min(fsuCfg.initPRBoccupancy, numOfAvPRBs);
                    
                    raObj.spectrumAvailable = sort(randSequenceFromBin(varMin, numOfAvPRBs));
                    obj.initFSUvar = true;
                    obj.picoTimerNew = fsuCfg.picoTimerNew;
                    obj.AV_PRBvector = ones(1, numOfAvPRBs) > 0;
                    raObj.pwPerPRBvector = zeros(1, numOfAvPRBs);
                    raObj.pwPerPRBvector(raObj.spectrumAvailable) = (raObj.maxTxPw / numOfAvPRBs);
                    raObj.pwPerPRBAllowedVector = (raObj.maxTxPw / numOfAvPRBs) .* ones(1, numOfAvPRBs);

                else                %%% main body of FSU algorithm
                    
                    %%% update input signals
                    
                    HIIflag = [cls(obj.macroBS_HIItx_indexVector).HIIboolean];
                    HIIsenderIndex = obj.macroBS_HIItx_indexVector(HIIflag==1);
                    
                    
                    %%% increase number of used PRBs
                    
                    if obj.picoTimerNew == 0
                        
                        raObj.cellsVsPRB_allocationMatrix = zeros(1, numOfAvPRBs);
                        listOfUsedPRBs = raObj.spectrumAvailable;
                        SINR_mtx_dB = obj.C_esteem_cls_mtx;

                        meanSINR = mean(SINR_mtx_dB,1);

                        vectPRB = 1:numOfAvPRBs;
                        
                        NA_PRBlist = vectPRB(~obj.AV_PRBvector);
                        meanSINR(NA_PRBlist) = -1000;
                        meanSINR(listOfUsedPRBs) = -1000;
                        [sortedPRB_SINR, sortedPRBindex] = sort(meanSINR,2,'descend');
                        AV_PRBvect = sortedPRBindex(sortedPRB_SINR > -1000);
                        if length(AV_PRBvect) >= fsuCfg.deltaN_PRB
                            newPRBsIndex = AV_PRBvect(1:fsuCfg.deltaN_PRB);
                        else
                            newPRBsIndex = AV_PRBvect;
                        end
                        raObj.spectrumAvailable = sort([raObj.spectrumAvailable newPRBsIndex]);
                        raObj.cellsVsPRB_allocationMatrix(raObj.spectrumAvailable) = true;
                        obj.picoTimerNew = fsuCfg.picoTimerNew;
                        
                        obj.BSstate(:,sum(obj.BSstate(1,:)>0)+1) = [1 currTTI length(raObj.spectrumAvailable) sum(obj.AV_PRBvector>0) 0].';  
                        
                    end
                    
                    if any(obj.pico_NA_PRB_timer == 0)
                        
                        % set PRBs with > 0  HII to AV
                        obj.AV_PRBvector(obj.pico_NA_PRB_timer == 0) = true;
                        PRBvec = 1:numOfAvPRBs;
                        % reset timestamp list
                        for n = PRBvec(obj.pico_NA_PRB_timer == 0)
                            obj.PRB_timeStamp(n).list = [];
                            obj.HIIcounterVector(n) = 0;
                        end
                        
                        obj.BSstate(:,sum(obj.BSstate(1,:)>0)+1) = [2 currTTI length(raObj.spectrumAvailable) sum(obj.AV_PRBvector>0) 0].'; 
                        
                    end
                    
                    if any(HIIflag)
                        
                        %%% check if 
                        for senderCnt = HIIsenderIndex
                            
                            HII_prblist = cls(senderCnt).fsuVars(carrierNum).macroHIIsignalStruct.list;
                            HII_timeStamp = cls(senderCnt).fsuVars(carrierNum).macroHIIsignalStruct.timestamp;                            
                            for n = HII_prblist
                                if obj.HIIcounterVector(n) == 0
                                    if   length(raObj.spectrumAvailable) > obj.minSpectrum_pBS
                                        %%% increment HII counter
                                        obj.HIIcounterVector(n) = obj.HIIcounterVector(n) + 1;
                                        %%% set LONG NA TIMER
                                        obj.pico_NA_PRB_timer(n) = fsuCfg.picoTimerTL;
                                        %%% add timestamp to the list 
                                        obj.PRB_timeStamp(n).list = HII_timeStamp;
                                        %%% remove PRB required in HII
                                        raObj.spectrumAvailable(raObj.spectrumAvailable == n) = [];
                                        raObj.cellsVsPRB_allocationMatrix(n) = false;
                                        %%% set PRBs as NON AVAILABLE
                                        obj.AV_PRBvector(n) = false;
                                    end
                                elseif obj.HIIcounterVector(n) == 1
                                    [row, col] = find(obj.PRB_timeStamp(n).list == HII_timeStamp);
                                    if ~isempty(row) 
                                        %%% this is a dulicate HII =>
                                        %%% remove timestamp from the list
                                        obj.PRB_timeStamp(n).list(col) = [];
                                    else
                                        %%% this is a new HII
                                        %%% a) find if MBS sent a previsous
                                        %%% HII
                                        [BSrow, BScol] = find( fix( obj.PRB_timeStamp(n).list /10000 )  == senderCnt);
                                        if ~isempty(BSrow) 
                                            obj.PRB_timeStamp(n).list(BScol) = HII_timeStamp;   % replace old time stame with new one
                                        else
                                            obj.pico_NA_PRB_timer(n) = fsuCfg.picoTimerTL;    
                                            obj.PRB_timeStamp(n).list = [obj.PRB_timeStamp(n).list HII_timeStamp];  % insert new timestamp
                                        end
                                    end
                                    if isempty(obj.PRB_timeStamp(n).list) 
                                        obj.pico_NA_PRB_timer(n) = fsuCfg.picoTimerTs;
                                    end
                                end
                            end
                        end
                        
                        obj.BSstate(:,sum(obj.BSstate(1,:)>0)+1) = [3 currTTI length(raObj.spectrumAvailable) sum(obj.AV_PRBvector>0) 0].'; 
   
                    end
                    
                    obj.pwPerPRBvector = zeros(1, numOfAvPRBs);
                    obj.pwPerPRBvector(raObj.spectrumAvailable) = (raObj.maxTxPw / numOfAvPRBs);
                end         % end of if obj.initFSUvar <= 0
            else
                
                raObj.cellsVsPRB_allocationMatrix = zeros(1, numOfAvPRBs);            
                raObj.pwPerPRBvector = zeros(1, numOfAvPRBs);
                raObj.pwPerPRBvector(raObj.spectrumAvailable) = (raObj.maxTxPw / numOfAvPRBs);
                raObj.pwPerPRBAllowedVector = (raObj.maxTxPw /  numOfAvPRBs) ...
                    .* ones(1, numOfAvPRBs);

            end     %  end of ~isempty(obj.usersToBeServedVect)
            
        end         % end of function
        
        
        function initializeSchedulerAndFSUParamters(obj,objCarrier,objMacro,...
                objPico,interfererIndexMatrix)
            
            global objCell;
            global cls;
            
            numOfAvPRBs = objCarrier.numOfAvailablePRB;
            av_sp = objCell.raVars(objCarrier.index).spectrumAvailable;
            
            %%% initialize MeNB algorithm
            if objCell.isMacro
                
                if isempty(obj.cumUserTh)
                    obj.cumUserTh = zeros(1,length(objCell.usersList));
                end
                
                obj.indexParamPF = 3;
                obj.PRBsHIIStateVector = -ones(1, numOfAvPRBs);
                obj.PRBsHIIStateVector(av_sp) = 0;
                obj.resendHIItimer = -1;
                
                usersList = objCell.usersList;
                
                lastPicoId = objPico.firstCellId + objPico.numOfCells - 1;
                
                interfererList = [];
                for n = usersList 
                    picoIndVector = interfererIndexMatrix(n,:);
                    picoInterferers = picoIndVector(picoIndVector >= objPico.firstCellId);
                    picoInterferers = picoInterferers(picoInterferers <= lastPicoId);
                    if ~isempty(picoInterferers)
                        interfererList = [interfererList picoInterferers(1)];
                    end
                end
                interfererList = unique(interfererList);
                obj.picoInterfererVector = interfererList;
                for n = interfererList
                    obj.macroToPicoInterfererIndexMatrix(n - objPico.firstCellId + 1) = true;
                end
            end
            
            if objCell.isPico 
                fstMacroInd = objMacro.firstCellId;
                lstMacroInd = fstMacroInd + objMacro.numOfCells - 1;
                
                macroIndVector = fstMacroInd:lstMacroInd;
                
                arrayFlagTmp = cls(macroIndVector).macroToPicoInterfererIndexMatrix(objCell.cellId - objPico.firstCellId + 1);
                
                obj.macroBS_HIItx_indexVector = macroIndVector(arrayFlagTmp==1);
                obj.picoTimerNew = fsuCfg.picoTimerNew;
                obj.pico_NA_PRB_timer = zeros(1, numOfAvPRBs);
                obj.HIIcounterVector = zeros(1, numOfAvPRBs);
                for n = 1:numOfAvPRBs
                    obj.PRB_timeStamp(n).list = [];
                end                
            end 
            
            obj.C_esteem_cls_mtx = zeros(length(objCell.usersList),...
                numOfAvPRBs);
        end         % end of function

        
        
        
        

    end
    
end

