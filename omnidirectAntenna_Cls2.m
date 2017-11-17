%%% SON SIMULATOR: omnidirectional antenna class
%
%   Author: Carlo Galiotto
%
%   Date:   July 2010
%
%   Description: this class implements the omnidirectional antenna class


classdef omnidirectAntenna_Cls2 < antenna_Cls    % extends antenna_Cls class

    
    methods
        
        function obj = omnidirectAntenna_Cls2()      % class constructor
            obj.antenna_type = 'Omnidirectional';
        end

        function antennaGain = getGain(obj,theta,phi)   % it returns the antenna gain
            
            antennaGain = 1;
            
        end       
    end   

end         % end of class