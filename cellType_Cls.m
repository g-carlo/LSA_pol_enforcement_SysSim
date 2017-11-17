classdef cellType_Cls<handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties  
        index;
        description;
        height;
        alphaPF;
        betaPF;
        gammaPF;
        allocationType;
        reuseScheme;
        radius;
        shadowingStdDev;
        firstCellId;
        firstUserId;
        
        numOfCells;
        minDistToUE;
        numOfUEs;
        numOfUEsPerOperator;        
        numOfCellsPerOperator;
        
        shiftArray;
                
        % FSU
        FSU_SINT_th;
        macroFSUfeedbackFlag;
    end
    
    methods
       function obj = cellType_Cls()         % constructor 
           obj.firstCellId = 1;
           obj.description = '';
           obj.radius = 0;
           obj.numOfCells = 1;
           obj.numOfCellsPerOperator = 0;           
           obj.minDistToUE = 0;
           obj.alphaPF = 0;
           obj.betaPF = 0;
           obj.gammaPF = 0;
           obj.macroFSUfeedbackFlag = false;
       end
       
       function [numOfCells, numOfUsers] = getNumOfCellsAndUsers(obj)
           numOfCells = obj.numOfCells;
           numOfUsers = obj.numOfUEs;           
       end
       
       function [firstCellId,lastCellId] = getFirstLastCellId(obj)
           firstCellId = obj.firstCellId;
           lastCellId = firstCellId + obj.numOfCells - 1;
       end
       
       function [firstUserId,lastUserId] = getFirstLastUserId(obj)
           firstUserId = obj.firstUserId;
           lastUserId = firstUserId + obj.numOfUEs - 1;
       end      
    end
end


