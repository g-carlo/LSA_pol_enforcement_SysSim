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

classdef eventSLS_Cls < handle
    
    properties 
      
        % double
        type;
        time;
        ID;
        activity;

    end     % end of properties
    
    methods
        
        % constructor
        function obj = eventSLS_Cls(type,time,activity,ID)         
            
            if ~ischar(type) || numel(time) ~= 1  ||~ischar(activity)
                error('Wrong inputs of class eventsSLS_cls');
            end
            obj.type        = type;
            obj.time        = time;
            obj.activity    = activity;
            obj.ID          = ID;

        end     % end of function                      
        
        % constructor
        function [type, time, activity, ID] = getProperties(obj)         
            
            type        = obj.type;
            time        = obj.time;
            activity    = obj.activity;
            ID          = obj.ID;
            
        end     % end of function 
        
        function setProperties(obj,type,time,activity,ID)         
            
            if ~ischar(type) || numel(time) ~= 1  ||~ischar(activity)
                error('Wrong inputs of class eventsSLS_cls');
            end
            obj.type        = type;
            obj.time        = time;
            obj.activity    = activity;
            obj.ID          = ID;

        end     % end of function 
        
        function out = isEmpty(obj)         
            
            if obj.time < 0 
                out = true;
            else 
                out = false;
            end

        end     % end of function 
        
% 
%         function update(obj)
% 
%         end
        
    end     % end of methods
 

%     methods (Access=private)        % private methods
% 
%          function emptyFunction(obj)          
% 
%          end         
%          
%          function initializeQueue(obj,numOfIncumbents,numOfNOs)
%              
%              field1 = 'type';  
%              value1 = cell(1,numOfIncumbents+numOfNOs);
%              field2 = 'time';  
%              value2 = zeros(1,numOfIncumbents+numOfNOs);
%              obj.eventQueue = struct(field1,value1,field2,value2); 
%              
%          end
%         
%     end

end         % end of class


