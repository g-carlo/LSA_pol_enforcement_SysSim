%%% SON SIMULATOR: position class
%
%   Author: Carlo Galiotto
%
%   Date:   July 2010
%
%   Description: this class implements the position class 

classdef position_Cls < handle

    properties
       
        x;         % x in the x-y-z coordinate system
        y;         % y in the x-y-z coordinate system
        z;         % z in the x-y-z coordinate system
        
    end     % end of properties
    
    methods
        
        % constructor
        
        function obj = position_Cls(x,y,z)
            
            obj.x = x;
            obj.y = y;
            obj.z = z;
            
        end
        
        % methods
        
        function setPosition(obj,x,y,z)     % sets the position
            
            obj.x = x;
            obj.y = y;
            obj.z = z;
            
        end
        
        function [x,y,z] = getPosition(obj)     % returns the vector of coordinates
            
            x = obj.x;
            y = obj.y;
            z = obj.z;
            
        end
        
        function a = getAbstractValue(obj)       % returns the abs coordinate vector
            
            a = sqrt( obj.x.^2 + obj.y.^2 + obj.z.^2);
            
        end
        
        function a = getDistance(obj,x)         % returns the abs of projection of the postion vector on the x-y plane
            
            a = sqrt( (obj.x - x.x).^2 + (obj.y-x.y).^2 + (obj.z-x.z).^2);
            
        end
        
        function a = getComplexCoordinate(obj)  % returns the x-y coordinates as a complex number
            
            a = obj.x + j*obj.y;
            
        end
    
    end     % end of methods

end         % end of class