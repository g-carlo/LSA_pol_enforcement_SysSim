classdef carrier_Cls < handle
    
    properties
        
        index;
        
        % the following properties are assigned from the inputParameters.m
        availableBW;
        frequency;
        loadPtxData;
        numOfAvailablePRB;
        PRBbw;
        scenarioCase;
        uePercentage;  
        user_Pn_PRB;
        
        % the following properties are assigned according to other
        % parameter values.
        electricalDowntilt;
                 
        
        % index of cells that can use the current carrier
        cellsInCarrier;
        
        % index of cells that can use the current carrier
        operatorsInCarrier;
        
        % index of UEs that can use the current carrier
        usersInCarrier;
        
        % channel coefficients (shadowing,penetrationloss,pathloss)
        mtx; 
        antennaGainMatrix;
        % mtx times antenna gain
        chPlusAntenna; 
        RSRPmatrix;
        rxPw;
        
        % SINR and capacity
        
        % total number of UEs in the network X total number of PRBs in the
        % current carrier
        C_mtx_estim_fsu_current;
        C_mtx_estim_fsu;
        lastEstimated_C_mtx;
        
        % SINR_mtx_estim_fsu_current;
        SINR_mtx_estim_fsu;
        SINR_mtx_estim;
        lastEstimated_SINR_mtx;
        SINR_mtx;
        C_mtx;
        int_mtx;
        
        %other measures
        
        thPerUE;
        ueCumSpectralEfficiency;
        ueSpectralEfficiency_cnt;
        avSINRperUE;
        cumResourceUsage;
        cumResourceAvailability;
        bandPerUE;
        cellOutage;
        cellThroughput;
                
        % matrix with SINR of each cell per PRB; 
        % the SINR is the one for the user allocated by the cell over that PRB
        cellVsPRB_SINR_mtx;             
        % resource allocation
                
        % total number of UEs x number of PRBs in the current carrier
        % it stores 1 if PRB j is assigned to UE ii. Otherwise it stores 0.
        allocationMatrix;
        
    end     % end of properties
    
    methods
        
        function obj = carrier_Cls()         % constructor
            obj.thPerUE = [];
            obj.cumResourceUsage = 0;
            obj.cumResourceAvailability = 0;
        end     % end of function
    
    end     % end of methods

end         % end of class


