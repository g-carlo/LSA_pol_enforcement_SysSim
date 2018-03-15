classdef userScheduling_Cls < handle
    
    properties
        % scheduling variables
        lastActiveUEindex;
        pastThVect;         % past throughput vector for PF scheduling
        presThVect;         % current throughput vector for PF scheduling
        initPFVar;          % initialization flag variable for PF scheduling 
        
        % array, each column corresponds to one PRB that has been 
        % assigned to the current cell, on the current carrier,
        % Each element of this array stores a UE index,
        % the identifier of the UE that is assigned the PRB
        userSpectVect;         
        
    end
    
    methods
        
        function obj = userScheduling_Cls(numOfAvPRBs)     
            obj.lastActiveUEindex = 0;
            obj.pastThVect = [];
            obj.presThVect = [];
            obj.initPFVar = 0;  
            obj.userSpectVect = zeros(1,numOfAvPRBs);            
        end
    end
    
end

