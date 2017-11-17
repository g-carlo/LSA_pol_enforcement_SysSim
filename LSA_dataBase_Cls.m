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

classdef LSA_dataBase_Cls < handle
    
    properties        
        
        %num of Channels
        numOfLSAChannels;                                
        
        %num of incumbents
        numOfIncumbents;
        
        % pixel position
        pixelPosition;
        
        % cell of N_R x N_C matrices
        matrices;
        
        % list of incumbents
        listOfIncumbents;
        
        % operator obj
        netOpEventGen;
        
        % polygon of Interest
        polygonVertices;
        
        % dimension of the area of interest
        X_length; 
        Y_length;
        
        % parameter variable
        paramVar;

    end     % end of properties
    
    methods
        
        % constructor
        function obj = LSA_dataBase_Cls(numOfChannels,numOfIncumbents,param)         
            
            
            obj.numOfLSAChannels    = numOfChannels;
            obj.numOfIncumbents     = numOfIncumbents;
            obj.listOfIncumbents    = cell(1,numOfIncumbents);
            obj.paramVar            = param;
%             obj.scheduleIncumbents();  

            
        end     % end of function
        
        
        function initializeDataBase(obj,incList,netOpGenEvent)         % set position of BS of the cell

            obj.matrices = cell(1,obj.numOfLSAChannels);
            
%             X_half_length = 2.5 * obj.paramVar.deployment.intersiteDistance * cos(pi/6);
%             Y_half_length = 2.5 * obj.paramVar.deployment.intersiteDistance;
%             obj.X_length = 2*X_half_length;
%             obj.Y_length = 2*Y_half_length;
%             N_R = 2 *ceil( X_half_length / obj.paramVar.database.spaceResolution );
%             N_C = 2 *ceil( Y_half_length / obj.paramVar.database.spaceResolution );
            obj.X_length = ( obj.paramVar.deployment.incumbentWestBound -obj.paramVar.deployment.incumbentEastBound ) ;
            obj.Y_length = ( obj.paramVar.deployment.incumbentNorthBound - obj.paramVar.deployment.incumbentSouthBound);
            N_R = ceil( obj.X_length / obj.paramVar.database.spaceResolution );
            N_C = ceil( obj.Y_length / obj.paramVar.database.spaceResolution );
            
            % initialize matrices
            switch obj.paramVar.deployment.cellScenario
                
                case 1
                    
                    for n = 1:obj.numOfLSAChannels; 
                        obj.matrices{n} = false( N_R , N_C);
                    end
                                        
                otherwise
                    error('Not implemented');
                    
            end
            
            % initialize incumbent list
            for n = 1:obj.numOfIncumbents
               obj.listOfIncumbents{n} = incList(n);
            end
            
            % initialize pixel
            x_pos = (( -N_C/2:(N_C/2-1 ) ) + 0.5 )*obj.paramVar.database.spaceResolution;
            y_pos = flipud((( -N_R/2:(N_R/2-1 ) ) + 0.5 ).')*obj.paramVar.database.spaceResolution;
            
            x_pos_rep = x_pos(ones(1,N_R),:);
            y_pos_rep = y_pos(:,ones(1,N_C));
            
            obj.pixelPosition = x_pos_rep + 1i * y_pos_rep;
            
            % initialize polygon vertices
            obj.polygonVertices = 2.5 * obj.paramVar.deployment.intersiteDistance .* exp( 1i*( pi/6:pi/3:2*pi )  );    
            
            obj.netOpEventGen   = netOpGenEvent;
            
            % schedule incumbents
%             obj.scheduleIncumbents(); 
            
            
        end     % end of function

        
% input parameters:
% ---------------------------------------
% - trafficTypeCarrier: traffic generator type (full buffer or burst traffic),
% which has been specified for the current carrier

% - varUEbuffer: buffer of each UE, considering all UEs in the network.

        function updateDatabase(obj)

            for n =  1:obj.numOfIncumbents
                obj.matrices{obj.listOfIncumbents{n}.channelID}(:) = false;
                
                if strcmp(obj.listOfIncumbents{n}.activityState ,'active')
                    obj.updateMatrix_active(n);
                end
               
                
            end


        end

        function updateMatrix_active(obj,incumbent_ID)

            obj.matrices{obj.listOfIncumbents{incumbent_ID}.channelID}(:) = false;
            pos = obj.listOfIncumbents{incumbent_ID}.position;
            dist = abs( pos - obj.pixelPosition );
            pixel_bool = dist <= obj.listOfIncumbents{incumbent_ID}.radius; 
            obj.matrices{obj.listOfIncumbents{incumbent_ID}.channelID}(pixel_bool) = true;

        end
        
        function updateMatrix_idle(obj,incumbent_ID)

            obj.matrices{obj.listOfIncumbents{incumbent_ID}.channelID}(:) = false;

        end
  
    end

 

    methods (Access=private)        % private methods

        function emptyFunction(obj)
            
        end    

        function scheduleIncumbents(obj)
             
            % FIXED ALLOCATION
            
            for n = 1:obj.numOfIncumbents
                obj.listOfIncumbents{n}.channelID = obj.paramVar.incumbent.frequencyCarrierIndex(n);
            end
             
        end
        

        
    end     % end of methods

end         % end of class


