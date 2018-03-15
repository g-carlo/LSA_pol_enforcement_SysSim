%%% SON SIMULATOR: UserEquipment class
%
%   Author: Carlo Galiotto
%
%   Date:   July 2010
%
%   Description: this class implements the UE class, defines the
%   properties of a UE and the methods 
%

classdef ue_Cls < handle

    properties
       
        pos;                % position [object] of the UE 
        complexPosition;
        antenna;            % antenna [antenna] of the UE
        networkID;          % network ID (not used yet)
        % integer. This is the network operator identifier
        networkOP_ID; 
        servingCellIndex;   % ID of the serving cell of this UE
        interfererCells;
        UEbuffer;
        ueArrivalTimer;
        pckCnt;  
        lambda;
        initVarTraffic;
        ueTXtimer;
        transmittedDataCnt;
        % activeUEset;
    end     % end of properties
    
    methods
        
        function obj = ue_Cls(x,y,z)     % class constructor
            
            obj.pos = position_Cls(x,y,z);
            obj.servingCellIndex = 0;
        end     % end of function
        
        function setPosition(obj,x,y,z)         % set position of BS of the UE
            
            obj.pos.setPosition(x,y,z);
            obj.complexPosition = x + 1i*y;
            
        end     % end of function
        
        function y = getComplexPosition(obj)    % get the x-y plane position of the UE in complex coordinates
            
             y = obj.pos.getComplexCoordinate();
            
        end     % end of function        
        
        function attachUser(obj,cell_ind)       % attach the UE to the cell with index cell_ind
            
            obj.servingCellIndex = cell_ind;
            
        end         % end of method
        
        function dataGenerator(obj, trafficCfg)                 % generates data
            
            obj.ueArrivalTimer = obj.ueArrivalTimer-1;
            
            switch trafficCfg.trafficGeneratorType
                
                case 1                  % full buffer                    
                    if obj.initVarTraffic==false                        
                        obj.initVarTraffic = true;
                        obj.UEbuffer = trafficCfg.FTP_model1_pckSize;
                        % obj.activeUEset = obj.UEbuffer >0;                       
                    end 
                    
                    
                case {2,3}          % FTP model 1
                    
                    %%% decrease counter for new arrivals                    
                    % generate new traffic                     
                    if obj.initVarTraffic                        
                        if obj.ueArrivalTimer<=0
                            obj.pckCnt = obj.pckCnt + 1;
                            obj.UEbuffer = obj.UEbuffer + trafficCfg.FTP_model1_pckSize;
                            newArrivalTimes = -log(1-rand(1,1))/obj.lambda * 1000; 
                            obj.ueArrivalTimer = newArrivalTimes;                            
                        end
                        % obj.activeUEset = obj.UEbuffer >0;                       
                    else                       
                        obj.ueArrivalTimer = -log(1-rand(1,1))./obj.lambda * 1000;
                        obj.initVarTraffic = true;   
                        obj.UEbuffer = trafficCfg.FTP_model1_pckSize;
                    end                     
                otherwise
                    
                    disp('No other traffic models implemented');                
            end
            
        end     % end of function
        
        function initializeTraffic(obj,trafficCfg,cls)         % constructor
            
            obj.ueArrivalTimer = -1;
            obj.initVarTraffic = false;
            
            obj.ueTXtimer = 0;
            
            obj.pckCnt = 0;
            obj.transmittedDataCnt = 0;           
            
            switch trafficCfg.trafficGeneratorType                
                case 2
                    totNumOfCells = length(cls);
                    nUe_perCell = zeros(1,totNumOfCells);
                    for n = 1:totNumOfCells
                        nUe_perCell(n) = length(cls(n).usersList);
                    end
                    lambdaCell = trafficCfg.FTP_model1_TrafficOfferedPerCell ./ ...
                        (nUe_perCell .* trafficCfg.FTP_model1_pckSize);                    
    
                    obj.lambda = lambdaCell(obj.servingCellIndex);
                case 3                    
                    obj.lambda = trafficCfg.FTP_model1_lambda_perUE;                
            end
        end     % end of function
        
        
        function evaluateTransmittedData(obj,trafficGenType,transmittedData)
            switch trafficGenType
                case 1              % full buffer
                    obj.ueTXtimer = obj.ueTXtimer + 1;
                    
                case {2,3}          % FTP model 1
                    
                    if obj.UEbuffer>0
                        obj.ueTXtimer = obj.ueTXtimer + 1;
                    end 
                    
                    obj.UEbuffer = obj.UEbuffer - transmittedData;                                       
                otherwise
                    error('No other traffic models available');            
            end
            obj.transmittedDataCnt = obj.transmittedDataCnt + transmittedData;            
        end
    
    end     % end of methods
    
    

end         % end of class