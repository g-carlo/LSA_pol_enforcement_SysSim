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

classdef incumbent_Cls < handle
    
    properties 
        % ID
        ID;
        
        % double
        radius;
        
        % double
        position;
        
        % area of interest
        areaLimits;
        
        % activity
        activityState;
        
        % double
        activityCnt;
        
        % incumbent receiver height
        receiverHeight;
        
        % incumbent receiver height
        transmitterHeight;
        
        % double 
        lambdaIdle;
        
        % double 
        lambdaActive;       
        
        % channel ID
        channelID;
                        
        %carrier frequency
        carrierFreq;
        
        % incumbent RX vector
        incumbent_RX_vector;
        
        % param variable
        paramVar;
        
    end
        
        properties (Access=private)
        % event  pointer
        eventPointer;
        
        % database pointer
        databasePointer;
        

    end     % end of properties
    
    methods
        
        % constructor
        function obj = incumbent_Cls(ID,paramVarInput,areaLimVector)         
            
            
            obj.ID               = ID;
            obj.carrierFreq      = [];
            obj.paramVar         = paramVarInput;
            obj.lambdaActive     = paramVarInput.incumbent.lambdaActive(ID);
            obj.lambdaIdle       = paramVarInput.incumbent.lambdaIdle(ID);
            obj.activityCnt      = obj.generateIdleTime(); 
            obj.areaLimits       = areaLimVector;      % [east west south north]
            obj.position         = 0;
            

             
        end     % end of function                      
        
        % constructor
        function initialize(obj,radius,database,event)         
            
            obj.databasePointer = database;
            obj.radius          = radius;
            obj.activityState    = [];     % incumbent is inizialized as inactive
%             obj.setNewRandomPosition();
            obj.eventPointer    = event;
            obj.receiverHeight  = obj.paramVar.incumbent.receiverHeight(obj.ID);
            obj.transmitterHeight  = obj.paramVar.incumbent.transmistterHeight(obj.ID);
            obj.incumbent_RX_vector.numOfReceivers = obj.paramVar.incumbent.numOfReceivers;              
            obj.update();
            
        end     % end of function  

        function update(obj)
            
            if isempty(obj.activityState)
                obj.activityState = 'active';
                obj.activityCnt = 1; 
%                 obj.updateDatabaseMatrix(false,obj.ID);
                obj.eventPointer.setProperties('incumbent',obj.activityCnt,'inactive',obj.ID);             
            
            else 
                
                if strcmp(obj.activityState ,'inactive')  % switch from inactive to active

                    obj.activityState = 'active';
                    obj.activityCnt = obj.activityCnt + obj.generateActiveTime( );
                    obj.setNewRandomPosition();
                    obj.incumbent_RX_vector.position = obj.generatePositionRX();
    %                 obj.updateDatabaseMatrix(true,obj.ID);
                    obj.eventPointer.setProperties('incumbent',obj.activityCnt,'inactive',obj.ID); 

                elseif strcmp(obj.activityState ,'active') % switch from active to inactive

                    obj.activityState = 'inactive';
                    obj.activityCnt = obj.activityCnt + obj.generateIdleTime( ); 
    %                 obj.updateDatabaseMatrix(false,obj.ID);
                    obj.eventPointer.setProperties('incumbent',obj.activityCnt,'active',obj.ID); 

                else

                    error('Wrong activity state');

                end
            end
            
%             switch obj.stateVar
%                 case 0  % INITIALIZATION state
%                     if obj.activityCnt <= 0
%                         obj.stateVar = 1;
%                     end
%                
%                 case 1  % SWITCH TO ACTIVE  state
%                     obj.activityState = true;
%                     obj.activityCnt = obj.generateActiveTime( ); 
%                     obj.updateDatabaseMatrix( true , obj.ID); 
%                     obj.stateVar = 2;
% 
%                 case 2  % ACTIVE-WAIT state
%                     if obj.activityCnt <= 0
%                         obj.stateVar = 3;
%                     end
%                     
%                 case 3  % SWITCH TO IDLE state
%                     obj.activityState = false;
%                     obj.activityCnt = obj.generateIdleTime( ); 
%                     obj.updateDatabaseMatrix( false ,obj.ID );
%                     obj.stateVar = 4;
%                     
%                 case 4  % IDLE-WAIT state
%                     if obj.activityCnt <= 0
%                         obj.stateVar = 1;
%                     end
% 
%             end
%             obj.activityCnt = obj.activityCnt -1;

        end
        
        function [SINR_vector, SNR_vector, SNR_vector_with_interference_margin] = computeSINR(obj,clsPositionVector,clsBoresightVector,...
                clsAntennaHeight,antennaStruct,antennaGain,elecDowntiltCarrier,maxTxPower,idx_enabled_BSs)

               
            %%%% BSs to incumbent receivers
            n_of_BSs = length(clsPositionVector); 
            n_of_RXs = obj.incumbent_RX_vector.numOfReceivers(obj.ID);

            posIncRX_matrix = repmat( obj.incumbent_RX_vector.position,[n_of_BSs,1]); 
            posBSs_matrix   = repmat( clsPositionVector.',[1,n_of_RXs]); 

            BS_to_incReceiver_dist = abs( posIncRX_matrix - posBSs_matrix)/1000;
            freq_MHz = obj.carrierFreq/1e6;
            [PLmatrix, ~] = ExtendedHata(BS_to_incReceiver_dist, freq_MHz, obj.receiverHeight, ...
                clsAntennaHeight, obj.paramVar.incumbent.propagationClutter(obj.ID) );
            
            shadowFading        = 10.^((8*randn(n_of_BSs,n_of_RXs)/10));
            rayleighFading      = exprnd(1,[n_of_BSs n_of_RXs]);

            phi_angleVector = angle( posIncRX_matrix - posBSs_matrix);
            theta_angleVector =  atan( (clsAntennaHeight - obj.receiverHeight)...
                ./ BS_to_incReceiver_dist);   

            antennaGainCoeff = antennaGain(obj.channelID) .* directiveAnt_v2(antennaStruct, repmat( clsBoresightVector.', [1 n_of_RXs]), phi_angleVector,...
            theta_angleVector, elecDowntiltCarrier(obj.channelID));

%             chCoeff_dB = -pathlossCoeff - obj.paramVar.channel.penetrationLoss;
            chCoeff_dB = -PLmatrix;
            chCoeff_lin = 10 .^ (chCoeff_dB ./ 10);
            
            chPlusAntenna = chCoeff_lin .* antennaGainCoeff .* shadowFading .* rayleighFading;
            
            interferenceVector = maxTxPower(obj.channelID) * repmat(idx_enabled_BSs.',1,n_of_RXs) .* chPlusAntenna;
            
            %%%% incumbent transmitters to incumbent receivers
            
            incTX_to_incRX_dist = abs( obj.position - obj.incumbent_RX_vector.position ) / 1000;
            [PLmatrix_incumb, ~] = ExtendedHata(incTX_to_incRX_dist, freq_MHz, obj.receiverHeight, ...
                obj.transmitterHeight, obj.paramVar.incumbent.propagationClutter(obj.ID) );
            
            chCoeff_inc_dB = -PLmatrix_incumb;
            chCoeff_inc_lin = 10 .^ (chCoeff_inc_dB ./ 10);
            
            shadowFading_inc        = 10.^((8*randn(1,n_of_RXs)/10));
            rayleighFading_inc      = exprnd(1,[1 n_of_RXs]);
            
            chTotal = chCoeff_inc_lin .* shadowFading_inc .* rayleighFading_inc;
            rxPowerVector = obj.paramVar.incumbent.txPower(obj.ID) .* chTotal;   
            
            noisePower = obj.paramVar.PHY.noiseFloor * obj.paramVar.incumbent.noiseFigure(obj.ID)  ...
                * obj.paramVar.incumbent.channelBandwidth(obj.ID);
            noisePower_with_interference_margin = noisePower*( 1 + obj.paramVar.incumbent.interferenceMargin(obj.ID) ) ;
            SINR_vector   = rxPowerVector ./ ( sum( interferenceVector )  + noisePower);
            SNR_vector   = rxPowerVector ./ noisePower;
            SNR_vector_with_interference_margin   = rxPowerVector ./ noisePower_with_interference_margin;
            

            % check if incumbent 


        end

        
    end     % end of methods
 

    methods (Access=private)        % private methods

         function out = generateIdleTime(obj)
            
             out = ceil(-log(1-rand(1,1))./obj.lambdaIdle);

         end
         
         function out = generateActiveTime(obj)
            
             out = ceil(-log(1-rand(1,1))./obj.lambdaActive);

         end
         
         function updateDatabaseMatrix(obj,activityState_to_set,incumbentID)
             
             if activityState_to_set
                 obj.databasePointer.updateMatrix_active(incumbentID);
             else
                 obj.databasePointer.updateMatrix_idle(incumbentID);
             end
             
         end
         
         function setNewRandomPosition(obj)
            
             
             if obj.ID == 1
                 
                 obj.position = 800;
                 
             else
                 
                 %%% position in a square
                 X_len = obj.areaLimits(2) - obj.areaLimits(1); 
                 Y_len = obj.areaLimits(4) - obj.areaLimits(3); 

                 posIncumbents = zeros(1,obj.databasePointer.numOfIncumbents);
                 for n = 1:obj.databasePointer.numOfIncumbents
                     posIncumbents(n) = obj.databasePointer.listOfIncumbents{n}.position;
                 end

                 new_pos = (rand(1)*X_len  - X_len/2)  + 1i*( rand(1)*Y_len - Y_len/2) ;
                 dist_to_other_incumbents = abs(new_pos-posIncumbents);
                 isTooClose = dist_to_other_incumbents <= 4*obj.radius;
                 isTooClose(obj.ID) = false;
                 isin = inpolygon( real( new_pos ) , imag( new_pos ), ...
                     [ obj.areaLimits(1) obj.areaLimits(2) obj.areaLimits(2) obj.areaLimits(1)  ] , [ obj.areaLimits(3) obj.areaLimits(3) obj.areaLimits(4) obj.areaLimits(4)  ] );
                 while ~isin || any(isTooClose)
                     new_pos = (rand(1)*X_len  - X_len/2)  + 1i*( rand(1)*Y_len - Y_len/2) ;
                     dist_to_other_incumbents = abs(new_pos-posIncumbents);
                     isTooClose = dist_to_other_incumbents <= 4*obj.radius;
                     isTooClose(obj.ID) = false;
                     isin = inpolygon( real( new_pos ) , imag( new_pos ), ...
                     [ obj.areaLimits(1) obj.areaLimits(2) obj.areaLimits(2) obj.areaLimits(1)  ] , [ obj.areaLimits(3) obj.areaLimits(3) obj.areaLimits(4) obj.areaLimits(4)  ] );
                 end
                 obj.position = new_pos;
                 
             end
             


% % % %              %%% position in a portion of annulus
% %              R1 = 5000; %3270;
% %              R2 = 0;     %630;
% %              rho = sqrt(rand(1)*(R1^2-R2^2)+R2^2);
% %              theta = pi/6*(rand(1)-0.5);
% % % %              
% % %              R1 = 10;
% % %              rho = R1;
% % %              theta = 0;
% % 
% %              X_len = rho .* cos(theta);
% %              Y_len = rho .* sin(theta);
% % 
% %              obj.position = X_len + 1i*Y_len;

         end
         
         function posVector = generatePositionRX(obj)

             
             R1 = obj.radius; %3270;
             R2 = obj.radius*0.05;     %630;
             rho = sqrt(rand(1,obj.incumbent_RX_vector.numOfReceivers(obj.ID))*(R1^2-R2^2)+R2^2);
             theta = 2*pi*(rand(1,obj.incumbent_RX_vector.numOfReceivers(obj.ID)));

             X_len = rho .* cos(theta);
             Y_len = rho .* sin(theta);

             posVector = obj.position + ( X_len + 1i*Y_len );
             
         end
         

        
    end

end         % end of class


