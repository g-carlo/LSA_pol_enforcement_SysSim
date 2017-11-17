%%% SON SIMULATOR: update channel function
%
%   Author: Carlo Galiotto
%
%   Date:   July 2010
%   Last Modification: April, 2015
%
%   Description: this function compute and update the matrix of the channel 
%   coefficients

function varDistNoWA = generateChannelFnc_v2(cls, ues, strAnt, strChannel,...
    varHeightRx, carrierArray, cellTypeArray, carrierPerCell, carrierCellType,...
    useWA,numOfSitesMacroCells,numOfMacroCellsPerSite,numOfOperators)

    global macroInd;
        
    numOfCarriers = length(carrierArray);
    numOfCellTypes = length(cellTypeArray);
    
    % total number of cells in the network
    totCells = length(cls);
    
    % total number of UEs in the network
    nUes = length(ues);
       
    LOSflag = strChannel.LOS;
    varPenetrationLoss = strChannel.penetrationLoss;
    shadowingFlag = strChannel.shadowingFlag;
    
    cellPos = [cls.complexPosition];
    boresightMtx = repmat([cls.orientation].',1,nUes);
    
    varRa = [cls.raVars];
    varRa = reshape(varRa,numOfCarriers,totCells);
        
    ue_posVect = [ues.complexPosition];
          
    % DISTANCE WRAP AROUND MODEL, ANGLES BETWEEN CELLS AND UEs
    varDistWA = zeros(totCells,nUes);
    varPosCellDistWA = zeros(totCells,nUes);
    varPhiAngle = zeros(totCells,nUes);
    varThetaAngle = zeros(totCells,nUes);
    varDistNoWA = zeros(totCells,nUes);
    
    
        
    for ii=1:numOfCellTypes
        objType = cellTypeArray(ii);  
        nCells = objType.numOfCells;
        if nCells > 0
            fstInd = objType.firstCellId;
            lstInd = fstInd + nCells - 1;

            [wa_cell_UE_mtx,wa_pos_mtx,phi_angle,theta_angle,distanceNoWA] = ...
                computeDistanceToUes(objType, cellPos(fstInd:lstInd),...
                ue_posVect, varHeightRx,useWA);

            varDistWA(fstInd:lstInd,:) = wa_cell_UE_mtx;
            varPosCellDistWA(fstInd:lstInd,:) = wa_pos_mtx;
            varPhiAngle(fstInd:lstInd,:) = phi_angle;
            varThetaAngle(fstInd:lstInd,:) = theta_angle;
            varDistNoWA(fstInd:lstInd,:) = distanceNoWA;
        end
    end
        
    %  ---------------------------------------------------------
        
    % CHANNEL COEFFICIENTS PER CARRIER
    for carrierNum = 1:numOfCarriers
        objCarrier = carrierArray(carrierNum);
        
        varAreaType = objCarrier.scenarioCase;
        varFreq = objCarrier.frequency;
        elecDowntiltCarrier = objCarrier.electricalDowntilt;
        
        % LOS PROPERTY (CELL, UE)
        varIsLOS = ones(totCells,nUes);
        fsPloss = zeros(totCells,nUes);
        PLmtx = zeros(totCells,nUes);
        sfCoeffMtx = zeros(totCells,nUes);
        antGainMtx = ones(totCells,nUes);
        
        for ii = 1:numOfCellTypes
            objType = cellTypeArray(ii);
            nCells = objType.numOfCells;
            if nCells>0
                
                objCarrierType = carrierCellType(carrierNum,ii);
                
                fstInd = objType.firstCellId;
                lstInd = fstInd + nCells - 1;
                varHeightTx = objType.height;
                
                                
                phi_angle = varPhiAngle(fstInd:lstInd,:);                
                varOrientation = boresightMtx(fstInd:lstInd,:);

                antGainMtx(fstInd:lstInd,:) = objCarrierType.antennaGain;
                
                if ii == macroInd
                    theta_angle = varThetaAngle(fstInd:lstInd,:);
                    antGainMtx(fstInd:lstInd,:) = antGainMtx(fstInd:lstInd,:) ...
                        .* directiveAnt_v2(strAnt, varOrientation,...
                        phi_angle, theta_angle, elecDowntiltCarrier);
                    
                end
                
                varDist = varDistWA(fstInd:lstInd,:);
                
                if LOSflag
                    varIsLOS(fstInd:lstInd,:) = setLOSproperty(ii,varAreaType,...
                        varDist, phi_angle,varOrientation,strAnt);
                end
                
                varPlModel = objCarrierType.pathLossModel;
                C = objCarrierType.coeff;

                [PLmtx(fstInd:lstInd,:), fsPloss(fstInd:lstInd,:)] = ...
                    computePathloss(varPlModel, C, varFreq, varIsLOS(fstInd:lstInd,:),...
                    varHeightRx, varHeightTx, varAreaType, varDist, nCells, nUes);
                                
                if shadowingFlag
                    
                    constShadowDev = objType.shadowingStdDev;
                    sfCoeffMtx(fstInd:lstInd,:) = computeShadowing(...
                        constShadowDev, ii, nCells, nUes, numOfSitesMacroCells,...
                        numOfMacroCellsPerSite, numOfOperators, varAreaType, varFreq,...
                        varDist, PLmtx(fstInd:lstInd,:), fsPloss(fstInd:lstInd,:),...
                        varPlModel);                  
                end              
            end
        end
        
        chCoeff = -PLmtx - varPenetrationLoss + sfCoeffMtx;
        objCarrier.mtx = 10 .^ (chCoeff ./ 10);
           
        % computation of antenna gain

        vartmp = (carrierPerCell(carrierNum,:)).';
        vartmp = repmat(vartmp,1,nUes);
        
        objCarrier.mtx(~vartmp) = 0;
        objCarrier.antennaGainMatrix = antGainMtx;
        
        vartmp = [varRa(carrierNum,:).maxTxPw];
        varMaxTxPw = repmat(vartmp.',1,nUes);
        
        objCarrier.chPlusAntenna = objCarrier.mtx .* antGainMtx;   
        
        objCarrier.rxPw = objCarrier.chPlusAntenna .* varMaxTxPw;
        
%         varSaveVars = 1;
%         if varSaveVars == 1
%             
%             varFnameVars = 'lsaSimulatorVars';
%             
%             vardistance = varDistWA;
%             varphiangle = varPhiAngle;
%             varthetaangle = varThetaAngle;
%             varpathloss = PLmtx;
%             varchcoeff = objCarrier.mtx;
%             varantgain = antGainMtx;
%             varchplusant = objCarrier.chPlusAntenna;
%             varisloss = varIsLOS;
%             varrxpw = objCarrier.rxPw;
%             save(varFnameVars,'vardistance','varphiangle','varthetaangle',...
%                 'varpathloss','varchcoeff','varantgain','varchplusant',...
%                 'varisloss', 'varrxpw','-append');
%         end
    
    end     
end



        