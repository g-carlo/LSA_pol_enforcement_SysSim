classdef resourceAllocation_Cls < handle
    properties
        maxTxPw; % array where the number of elements is the number of 
                 % carriers
        spectrumAvailable; 
        
        
        pwPerPRBvector;
        pwPerPRBAllowedVector;        
        PRBbandTypeVector;                  % (for SFR) vector reporting PRB with high/low power
              
        currentResourceUsage;
        currentResourceAvailability;  
        cellsVsPRB_allocationMatrix;
    end
    
    methods
        
        function obj = resourceAllocation_Cls(varMaxTxPw,numOfAvailablePRB)
            obj.maxTxPw = varMaxTxPw;
            obj.spectrumAvailable = [];
            obj.pwPerPRBvector = zeros(1,numOfAvailablePRB);
            obj.pwPerPRBAllowedVector = zeros(1,numOfAvailablePRB);
            obj.PRBbandTypeVector = zeros(1,numOfAvailablePRB);            
            obj.cellsVsPRB_allocationMatrix = zeros(1,numOfAvailablePRB);
            
            obj.currentResourceUsage = 0;
            obj.currentResourceAvailability = 0;
        end
    end
    
end

