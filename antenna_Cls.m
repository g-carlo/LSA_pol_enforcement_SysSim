%%% SON SIMULATOR: antenna class
%
%   Author: Carlo Galiotto
%
%   Date:   July 2010
%
%   Description: this class implements the abstract antenna class


classdef antenna_Cls < handle

    properties
        
        antenna_type;
       
    end     % end of properties
    
    methods (Abstract)
        % Returns antenna gain as a function of theta
        antenna_gain = getGain(obj,theta,phi)

    end

end         % end of class