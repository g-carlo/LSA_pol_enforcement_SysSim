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

classdef networkOperator_eventGenerator_Cls < handle
    
    properties 
        
        % double
        numOfNOs;
        
        % double
        demandVector;       
       
        % activity
        newDemandCnt;
                
        % double 
        lambdaNewDemand;
                
        % database pointer
        xxPointer;        
        
        % param variable
        var;
        
        % event pointer
        eventPointer;       
        
        % operator per carrier
        operatorPerCarrier;
        
        % cell idx per operator
        cell_idx_per_NO;

    end     % end of properties
    
    methods
        
        % constructor
        function obj = networkOperator_eventGenerator_Cls(numOfNOs,lam_newDem)         
            
            obj.numOfNOs            = numOfNOs;     % incumbent is inizialized as inactive
            obj.lambdaNewDemand     = lam_newDem;
            obj.newDemandCnt        = obj.generateNewDemandTime();
%             obj.demandVector        = zeros(1,numOfNOs);
            obj.demandVector        = obj.generateNewDemandVector();
            
        end     % end of function                      
        
        % constructor
        function initialize(obj,event,operatorPerCellVector,numOfLSAChannels)         
            
            if numOfLSAChannels == 1
                
                
                obj.eventPointer        = event;
                if obj.numOfNOs == 2
                    obj.operatorPerCarrier  = [1 2 2];
                end

                for n = 1:obj.numOfNOs
                    obj.cell_idx_per_NO(n).list = find(operatorPerCellVector == n) ;
                end
                
            elseif numOfLSAChannels == 2
                obj.eventPointer        = event;
                obj.operatorPerCarrier  = [1 2 1 2];

                for n = 1:obj.numOfNOs
                    obj.cell_idx_per_NO(n).list = find(operatorPerCellVector == n) ;
                end
            else
               error('Wrong number of LSA channels'); 
            end
            
        end     % end of function  

        function update(obj)           
            
            obj.newDemandCnt        = obj.newDemandCnt + obj.generateNewDemandTime(); 
            obj.demandVector        = obj.generateNewDemandVector();
            obj.eventPointer.setProperties('operator',obj.newDemandCnt,'new_demand',1); 
            
        end
        
        function updateEvent(obj,timeEvent)           
            
            obj.newDemandCnt        = obj.newDemandCnt + obj.generateNewDemandTime(); 
            obj.demandVector        = obj.generateNewDemandVector();
            obj.eventPointer.setProperties('operator',timeEvent,'new_demand',1); 
            
        end
        
    end     % end of methods
 

    methods (Access=private)        % private methods

         function out = generateNewDemandTime(obj)
            
             % static 
             out = 60;
             
             % exponential 
%              out = ceil(-log(1-rand(1,1))./obj.lambdaNewDemand);

         end                     
         
          function out = generateNewDemandVector(obj)
            
             % generates a vector of demands, multiple of 0.25
             out = 0.25*(randi(5,1,obj.numOfNOs)-1);
             while all( out == 0)
                 out = 0.25*(randi(5,1,obj.numOfNOs)-1);
             end

         end 
        
    end

end         % end of class


