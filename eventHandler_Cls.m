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

classdef eventHandler_Cls < handle
    
    properties 
      
        % double
        eventQueue;
        numOfIncumbents;

    end     % end of properties
    
    methods
        
        % constructor
        function obj = eventHandler_Cls(numOfIncumbents)         
            
            obj.numOfIncumbents     = numOfIncumbents;
            obj.initializeQueue(numOfIncumbents,numOfNOs);

        end     % end of function                      
        
        % constructor

        function update(obj)

        end
        
    end     % end of methods
 

    methods (Access=private)        % private methods

         function emptyFunction(obj)          

         end         
         
         function initializeQueue(obj,numOfIncumbents,numOfNOs)
             
             field1 = 'type';  
             value1 = cell(1,numOfIncumbents+numOfNOs);
             field2 = 'time';  
             value2 = zeros(1,numOfIncumbents+numOfNOs);
             obj.eventQueue = struct(field1,value1,field2,value2); 
             
         end
        
    end

end         % end of class


