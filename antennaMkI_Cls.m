%%% SON SIMULATOR: antenna MkI class
%
%   Author: Carlo Galiotto
%
%   Date:   July 2010
%
%   Description: this class implements the antenna class with the directive antenna
%   pattern given by 3GPP TR 36.814


classdef antennaMkI_Cls < antenna_Cls

    properties

        elet_downtilt;          % electrical downtilt
        phi3db;                 % 3dB beam-width of the horizontal pattern of the directional antenna
        theta3db;               % 3dB beam-width of the vertical pattern of the directional antenna
        gainMinTheta;           % minimum gain of the vertical pattern of directional antenna 
        gainMinPhi;             % minimum gain of the horizontal pattern of directional antenna 
       
    end     % end of properties
    
    methods
        function obj = antennaMkI_Cls(strAntenna,electricalDowntilt)     % antenna constructor
            
            obj.antenna_type = 'Antenna MkI';
            obj.elet_downtilt = electricalDowntilt;
            obj.phi3db = strAntenna.horizontalPhi3dB;
            obj.theta3db = strAntenna.verticalTheta3dB;
            obj.gainMinTheta = strAntenna.verticalMinGain;
            obj.gainMinPhi = strAntenna.horizontalMinGain;

        end

        function antennaGain = getGain(obj,phi,theta)       % return the antenna gain over the direction phi (x-y plane) and theta (tilt angle)
            
            phi = phi/pi*180;
            theta = theta/pi*180;          
            gainV = -min(12.*(( theta - obj.elet_downtilt )./obj.theta3db ).^2,obj.gainMinTheta);
%             gainV = 0;
            gainH = -min(12.*(( phi )./obj.phi3db ).^2,obj.gainMinPhi);
            antennaGain = 10.^( (-min(-(gainV+gainH),obj.gainMinPhi) /10 ));
            
        end
        
    end
    
    

end         % end of class



