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

classdef LSA_controller_Cls < handle
    
    properties        
        
        %num of Channels
        numOfLSAChannels;     
        
        channelFrequencies;
        
        %num of incumbents
        numOfIncumbents;

        %num of incumbents
        numOfOperators;
        
        % list of incumbents
        listOfIncumbents;
        
        % operator obj
        netOpEventGen;
        
        % memory of 
        schedData;
                
        % parameter variable
        paramVar;

    end     % end of properties
    
    methods
        
        % constructor
        function obj = LSA_controller_Cls(numOfChannels,numOfIncumbents,numOfOperators,param,freqVector)         
            
            
            obj.numOfLSAChannels    = numOfChannels;
            obj.channelFrequencies  = freqVector;
            obj.numOfIncumbents     = numOfIncumbents;
            obj.numOfOperators      = numOfOperators;
            obj.listOfIncumbents    = cell(1,numOfIncumbents);
            obj.paramVar            = param;
%             obj.scheduleIncumbents();  

            
        end     % end of function
        
        
        function initializeLSAcontroller(obj,incList,netOpGenEvent,penaltyVector)         % set position of BS of the cell
            
            % initialize incumbent list
            for n = 1:obj.numOfIncumbents
               obj.listOfIncumbents{n} = incList(n);
            end
            
            obj.schedData.currentOperator       = [];
            obj.schedData.nextEstSchedulingTime = [];
            obj.schedData.leftScheduledTime     = [];
            obj.schedData.lastSchedulingTime    = [];
            obj.schedData.lastSumOfEnabledCells = [];
            obj.schedData.penalty               = penaltyVector;
%             obj.schedData.penalty               = zeros(1,obj.numOfOperators);
            obj.schedData.misBehFlag            = false(1,obj.numOfOperators);
            obj.schedData.P_detection           = 1;
            obj.schedData.P_falseAlarm          = 0;
            obj.schedData.lastMisbehaveTime     = zeros(1,obj.numOfOperators); %Change: Added last misbehave time parameter
            
            obj.netOpEventGen   = netOpGenEvent;
            
            % schedule incumbents
            obj.scheduleIncumbents(); 
            
            % sc
            obj.schedData.lastScheduledOp = [];
            
        end     % end of function

        
% input parameters:
% ---------------------------------------
% - trafficTypeCarrier: traffic generator type (full buffer or burst traffic),
% which has been specified for the current carrier

% - varUEbuffer: buffer of each UE, considering all UEs in the network.

        function [cellsPerCarrier, opPerLSAchannel] = runL1_algorithm(obj,clsPositionVector,clsBoresightVector,...
                clsAntennaHeight,antennaStruct,antennaGain,elecDowntiltCarrier,maxTxPower,cellsPerCarrier,numOfLSAchannels,currentTime)

            
            %%% R1 algorithm for 1 channel/ 2 operator / 1 channel
            
            if numOfLSAchannels == 1
                
                %%%%%%% ONLY 1 Incumbent 
                
                % choose network operator;
                
                if isempty(obj.schedData.currentOperator)
                    
                    obj.scheduleOperator(currentTime);
                    obj.schedData.leftScheduledTime = obj.paramVar.L1_scheduling.schedulingSlotLength;
                    
                elseif currentTime >= obj.schedData.nextEstSchedulingTime
                    
                    obj.scheduleOperator(currentTime);
                    obj.schedData.leftScheduledTime = obj.paramVar.L1_scheduling.schedulingSlotLength;
                    
                elseif currentTime < obj.schedData.nextEstSchedulingTime
                    
                    n_cell_per_operator = obj.paramVar.deployment.numOfSites*obj.paramVar.deployment.numOfCellsPerSite;
                    if obj.schedData.lastSumOfEnabledCells > 0
                        obj.schedData.leftScheduledTime = obj.schedData.leftScheduledTime - ...
                            (currentTime - obj.schedData.lastSchedulingTime)*(obj.schedData.lastSumOfEnabledCells/ n_cell_per_operator);
                    end
                    
                    
                end   
                
                clearToTxMat = true(obj.numOfIncumbents,length(clsPositionVector));
                
                channelID = [];
                
                for inc_idx = 1:obj.numOfIncumbents

                    if isempty(channelID)
                        channelID = obj.listOfIncumbents{inc_idx}.channelID;
                    else 
                        if channelID ~= obj.listOfIncumbents{inc_idx}.channelID;
                           error('error'); 
                        end
                            
                    end 
                    
                    if strcmp(obj.listOfIncumbents{inc_idx}.activityState,'active')

                            % i) check which cell of operator inc_idx on channel
                        

    %                     NO_idx  = obj.netOpEventGen.operatorPerCarrier(channelID); 
                        if obj.paramVar.L1_scheduling.excZoneMethod == 1
                            clearToTxMat( inc_idx ,:) = obj.computeChannelMatrix(inc_idx,clsPositionVector,clsBoresightVector,...
                               clsAntennaHeight,antennaStruct,antennaGain(channelID),elecDowntiltCarrier(channelID),...
                               maxTxPower(channelID),obj.channelFrequencies(channelID));
                        elseif obj.paramVar.L1_scheduling.excZoneMethod == 2
                            clearToTxMat( inc_idx ,:) = obj.computeDistanceMatrix(inc_idx,clsPositionVector);
                        end
                       
                        

                        % reset 
% % %                         operatorToBeSched = obj.schedData.currentOperator;
% % %                         cell_idx_bool = false(size(clearToTx ) );
% % %                         cell_idx_bool( obj.netOpEventGen.cell_idx_per_NO(operatorToBeSched).list )  = true;
% % %                         cellsPerCarrier(channelID, : ) = 0;
% % %                         cellsPerCarrier(channelID, cell_idx_bool & clearToTx ) = 1;
% % %                         sumOfEnabledCells = sum(cell_idx_bool & clearToTx );
% % %                         opPerLSAchannel = operatorToBeSched; 
                        % inc_idx are allowed to TX
                        % ii) allocate such cell on channel inc_idx 
% 
%                     else
                        

% % %                         channelID = obj.listOfIncumbents{1}.channelID;
% % % 
% % %                         operatorToBeSched = obj.schedData.currentOperator;
% % %     %                     NO_idx  = obj.netOpEventGen.operatorPerCarrier(channelID); 
% % % 
% % %                         % reset 
% % %                         cellsPerCarrier(channelID, : ) = 0;
% % %                         cellsPerCarrier(channelID, obj.netOpEventGen.cell_idx_per_NO(operatorToBeSched).list) = 1;
% % %                         sumOfEnabledCells = obj.paramVar.deployment.numOfSites*obj.paramVar.deployment.numOfCellsPerSite;
% % %                         opPerLSAchannel = operatorToBeSched; 
                        % i) allocate operator inc_idx on channel inc_idx
                        % ii) allocate all cells of operator inc_idx on channel inc_idx

                    end


                end
                
                clearToTxTot = all( clearToTxMat , 1); 
                operatorToBeSched = obj.schedData.currentOperator;
                cell_idx_bool = false(size(clearToTxTot ) );
                cell_idx_bool( obj.netOpEventGen.cell_idx_per_NO(operatorToBeSched).list )  = true;
                cellsPerCarrier(channelID, : ) = 0;
                cellsPerCarrier(channelID, cell_idx_bool & clearToTxTot ) = 1;
                sumOfEnabledCells = sum(cell_idx_bool & clearToTxTot );
                opPerLSAchannel = operatorToBeSched;
                
                % update variables before next event
                
                n_cell_per_operator = obj.paramVar.deployment.numOfSites*obj.paramVar.deployment.numOfCellsPerSite;
                obj.schedData.nextEstSchedulingTime = ceil( ( obj.schedData.leftScheduledTime)/ ...
                    (sumOfEnabledCells/n_cell_per_operator )) + currentTime;
                obj.schedData.lastSchedulingTime    = currentTime;
                obj.schedData.lastSumOfEnabledCells = sumOfEnabledCells;
                               
                % update event
                obj.netOpEventGen.updateEvent(obj.schedData.nextEstSchedulingTime);           
                
                
                %%%% MORE THAN 1 INCUMBENT
                
                
                
            elseif numOfLSAchannels == 2
            
            %%% R1 algorithm for 2 channel / 2 operators / 2 incumbents
                        
                for inc_idx = 1:obj.numOfIncumbents

                    if strcmp(obj.listOfIncumbents{inc_idx}.activityState,'active')

                        % i) check which cell of operator inc_idx on channel
                        channelID = obj.listOfIncumbents{inc_idx}.channelID;

                        NO_idx  = obj.netOpEventGen.operatorPerCarrier(channelID); 
%                         if obj.paramVar.L1_scheduling.excZoneMethod == 1
%                             clearToTxMat( inc_idx ,:) = obj.computeChannelMatrix(obj,inc_idx,clsPositionVector,clsBoresightVector,...
%                                clsAntennaHeight,antennaStruct,antennaGain(channelID),elecDowntiltCarrier(channelID),...
%                                maxTxPower(channelID),obj.channelFrequencies(channelID));
%                         elseif obj.paramVar.L1_scheduling.excZoneMethod == 2
%                             clearToTxMat( inc_idx ,:) = obj.computedistanceMatrix(obj,inc_idx,clsPositionVector);
%                         end
                        clearToTx = obj.computeChannelMatrix(obj,inc_idx,clsPositionVector,clsBoresightVector,...
                           clsAntennaHeight,antennaStruct,antennaGain(channelID),elecDowntiltCarrier(channelID),...
                           maxTxPower(channelID),obj.channelFrequencies(channelID));

                        % reset 
                        cell_idx_bool = false(size(clearToTx ) );
                        cell_idx_bool( obj.netOpEventGen.cell_idx_per_NO(NO_idx).list )  = true;
                        cellsPerCarrier(channelID, : ) = 0;
                        cellsPerCarrier(channelID, cell_idx_bool & clearToTx ) = 1;
                        % inc_idx are allowed to TX
                        % ii) allocate such cell on channel inc_idx 

                    else
                        channelID = obj.listOfIncumbents{inc_idx}.channelID;

                        NO_idx  = obj.netOpEventGen.operatorPerCarrier(channelID); 

                        % reset 
                        cellsPerCarrier(channelID, : ) = 0;
                        cellsPerCarrier(channelID, obj.netOpEventGen.cell_idx_per_NO(NO_idx).list) = 1;
                        % i) allocate operator inc_idx on channel inc_idx
                        % ii) allocate all cells of operator inc_idx on channel inc_idx



                    end
                    % check if incumbent(inc_idx) is active

                    % if not 
                        % i) allocate operator inc_idx on channel inc_idx
                        % ii) allocate all cells of operator inc_idx on channel inc_idx

                    % if so
                        % i) check which cell of operator inc_idx on channel
                        % inc_idx are allowed to TX
                        % ii) allocate such cell on channel inc_idx 

                end
            
            else 
                error('WRong number of LSA channels');
                
            end
            
        end
        
        function logMisbehaviour(obj,estim_active_BSs,idx_enabled_BSs,P_d,P_fa)
            
            estim_active_BSs_bool = false(1,length(idx_enabled_BSs));                 
            estim_active_BSs_bool( estim_active_BSs ) = true;
            
            misBehFlag = any(estim_active_BSs_bool(~idx_enabled_BSs ) );
%             obj.schedData.misBehFlag = false(1,obj.numOfOperators);
            obj.schedData.misBehFlag( obj.schedData.currentOperator) = obj.schedData.misBehFlag( obj.schedData.currentOperator) | misBehFlag;
            
            obj.schedData.P_detection           = P_d;
            obj.schedData.P_falseAlarm          = P_fa;
            
        end

          
    end

 

    methods (Access=private)        % private methods

        function emptyFunction(obj)
            
        end    

        function scheduleIncumbents(obj)
             
            for n = 1:obj.numOfIncumbents
                obj.listOfIncumbents{n}.channelID   = obj.paramVar.general.numOfCarriers;
                obj.listOfIncumbents{n}.carrierFreq = obj.channelFrequencies( obj.listOfIncumbents{n}.channelID ) ;
            end
             
        end
        
        function clearToTx = computeChannelMatrix(obj,inc_idx,clsPositionVector,clsBoresightVector,...
                clsAntennaHeight,antennaStruct,antennaGain,elecDowntiltCarrier,maxTxPower,freq_Hz)
               
            BS_to_inc_dist = abs(clsPositionVector - obj.listOfIncumbents{inc_idx}.position); 
            radius_to_dist = abs(obj.listOfIncumbents{inc_idx}.radius ./ BS_to_inc_dist);

            posIncumbentReceiver = radius_to_dist.*clsPositionVector ...
                + (1-radius_to_dist) .* obj.listOfIncumbents{inc_idx}.position; 

            BS_to_incReceiver_dist = abs( clsPositionVector - posIncumbentReceiver)/1000;
            freq_MHz = freq_Hz/1e6;
            [PLmatrix, ~] = ExtendedHata(diag(BS_to_incReceiver_dist), freq_MHz, obj.listOfIncumbents{inc_idx}.receiverHeight, ...
                clsAntennaHeight, obj.paramVar.incumbent.propagationClutter(inc_idx) );

            pathlossCoeff = diag(PLmatrix) - obj.paramVar.database.ShadFadingSafetyMargin;

            phi_angleVector = angle( posIncumbentReceiver - clsPositionVector);
            theta_angleVector =  atan( (clsAntennaHeight - obj.listOfIncumbents{inc_idx}.receiverHeight)...
                ./ BS_to_incReceiver_dist);   


            antennaGainCoeff = antennaGain .* directiveAnt_v2(antennaStruct, clsBoresightVector, phi_angleVector,...
            theta_angleVector, elecDowntiltCarrier);

%             chCoeff_dB = -pathlossCoeff - obj.paramVar.channel.penetrationLoss;
            chCoeff_dB = -pathlossCoeff;
            chCoeff_lin = 10 .^ (chCoeff_dB ./ 10);

            chPlusAntenna = chCoeff_lin.' .* antennaGainCoeff;

            rxPowerVector = maxTxPower .* chPlusAntenna;

            noisePower = obj.paramVar.PHY.noiseFloor * obj.paramVar.incumbent.noiseFigure(inc_idx)  ...
                * obj.paramVar.incumbent.channelBandwidth(inc_idx);

            clearToTx = rxPowerVector/noisePower < obj.paramVar.incumbent.IR_threshold(inc_idx) & ...
                abs(clsPositionVector - obj.listOfIncumbents{inc_idx}.position) > obj.listOfIncumbents{inc_idx}.radius;


            % check if incumbent 


        end
        
        function clearToTx = computeDistanceMatrix(obj,inc_idx,clsPositionVector)
               
            BS_to_inc_dist = abs(clsPositionVector - obj.listOfIncumbents{inc_idx}.position); 
            radius_to_dist = abs(obj.listOfIncumbents{inc_idx}.radius ./ BS_to_inc_dist);

            posIncumbentReceiver = radius_to_dist.*clsPositionVector ...
                + (1-radius_to_dist) .* obj.listOfIncumbents{inc_idx}.position; 

            BS_to_incReceiver_dist = abs( clsPositionVector - posIncumbentReceiver);

            clearToTx = BS_to_incReceiver_dist > obj.paramVar.incumbent.exclusionZoneDistance(inc_idx);


            % check if incumbent 


        end
        
        function scheduleOperator(obj,currentTime)
            
            KV_penalty_threshold = obj.paramVar.L1_scheduling.penaltyThreshold;
            Dummy_op = obj.numOfOperators;

%             % round robin
                    if isempty(obj.schedData.lastSchedulingTime)
                        operatorToBeSched = randsample(1:obj.numOfOperators-1,1);
                    else

                        obj.LSA_penalty(currentTime);   %Calculate penalty
                        %reset misbehavior flag flag
                        obj.schedData.misBehFlag  = false(1,obj.numOfOperators);
                        if all( obj.schedData.penalty == 0)
                            operatorToBeSched = randi([1 obj.numOfOperators-1],1,1);
        %                     operatorToBeSched = rem(obj.schedData.currentOperator,obj.numOfOperators-1)+1;
                        else 
                            [KV_penalty, opIdx]  = min( obj.schedData.penalty(1:obj.numOfOperators-1) );
                            if KV_penalty>KV_penalty_threshold
                                opIdx = Dummy_op;
                            end
                            operatorToBeSched = opIdx;
                        end
                    end
            obj.schedData.currentOperator = operatorToBeSched;

            % check if incumbent 
%             disp(operatorToBeSched);


        end
        
        function Penalty = LSA_penalty(obj, currentTime)

            a = 100; %magnitude of penalty without considering Pd and Pfa
            
%             opIdx = obj.schedData.currentOperator;

            weighted_penalty = a * max(0,( obj.schedData.P_detection - obj.schedData.P_falseAlarm)) .* obj.schedData.misBehFlag;
%             pen_val = obj.schedData.penalty( opIdx);
%             Delta_time = time - obj.schedData.lastSchedulingTime;
            Delta_time = currentTime - obj.schedData.lastMisbehaveTime;    %Change: Calculating elapsed time based on the last misbehavour time 
            if (strcmp(obj.paramVar.L1_scheduling.penaltyCalculation,'exponential'))
                Elapsed_value = floor(obj.schedData.penalty.* exp(-(Delta_time)./obj.paramVar.L1_scheduling.penaltyLambda));
            elseif (strcmp(obj.paramVar.L1_scheduling.penaltyCalculation,'linear'))
                Elapsed_value = max(0,obj.schedData.penalty+obj.paramVar.L1_scheduling.penaltyBeta*Delta_time);
            else
                disp('check Penalty Function');
                pause
            end
            
            obj.schedData.lastMisbehaveTime = obj.schedData.misBehFlag .* currentTime + ...
                obj.schedData.lastMisbehaveTime .* (~obj.schedData.misBehFlag);  %Change: added track of the last misbehaviour time
            obj.schedData.penalty = Elapsed_value + weighted_penalty;
%             obj.schedData.penalty
%             disp(obj.schedData.penalty);
        end
        
        
    end     % end of methods

end         % end of class


