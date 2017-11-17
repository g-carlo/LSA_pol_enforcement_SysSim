%%% SON SIMULATOR: cell class
%
%   Author: Carlo Galiotto
%
%   Date:   August 2010
%
%   Description: this class implements the cell (both pico and macro) class, defines the
%   properties of a eNB and the methods 
%

classdef data_Cls < handle
    
    properties
        userVsCell_allocationMtx;       % matrix of allocation; rows = users, columns = cells
       
        %(UEs x number of interferers per UE)
        interfererIndexMatrix; 
    end     % end of properties
    
    methods
        
        function obj = data_Cls()         % constructor
           
            obj.userVsCell_allocationMtx = [];    % zeros(obj.conf.general.numOfCells,obj.conf.general.numOfUsers);
                        
        end     % end of function
    
    end     % end of methods

end         % end of class


