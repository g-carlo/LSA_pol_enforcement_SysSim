%%% SON SIMULATOR: initializeUsersMultipleFreq
%
%   Author: Almir Maric
%
%   Date:   July 2013
%
%   Description: function initialize users for each frequency. Users are assigned to frequency randomly or quasi stationary.

function ues = initializeUsersMultipleFreq_v2(cls,cellTypeArray,genConfig,...
    deployConfig,UEheight)
    
    global macroInd;
    global picoInd;
    global smallInd;
       
    varMode = genConfig.mode;
    
    if strcmp(varMode,'montecarlo') || strcmp(varMode,'static')     % MODE: SIMULATION
        ues = ue_Cls(0,0,0);

        totalNumOfUsers = sum([cellTypeArray(:).numOfUEs]);
        
        for ii=1:totalNumOfUsers
            ues(ii) = ue_Cls(0,0,0);
        end
        
        if cellTypeArray(macroInd).numOfCells > 0
            
             placeUserFunction_3GPPscenario_v2(ues,deployConfig,UEheight,...
                 cls,cellTypeArray);
             
             firstIndex = cellTypeArray(macroInd).firstUserId;
             varNumOfUes = cellTypeArray(macroInd).numOfUEs;
             lastIndex = firstIndex + varNumOfUes - 1;
             
             drawPositions(ues,varNumOfUes,macroInd,firstIndex,lastIndex,...
                 0,0,'userPos.txt',0);
             
             if cellTypeArray(picoInd).numOfUEs > 0
                 firstIndex = cellTypeArray(picoInd).firstUserId;
                 varNumOfUes = cellTypeArray(picoInd).numOfUEs;
                 lastIndex = firstIndex + varNumOfUes - 1;

                 drawPositions(ues,varNumOfUes,picoInd,firstIndex,lastIndex,...
                     0,0,'userPos.txt',0);
             end
             
                              
             for ii = 1: ( cellTypeArray(macroInd).numOfCellsPerOperator * deployConfig.numOfUEperCellType( macroInd ))
                 for jj = 1:deployConfig.numOfNetworkOperators   
                     ue_ind = cellTypeArray(macroInd).numOfCellsPerOperator * (jj  - 1) + ii;
                     ues(ue_ind).networkOP_ID = jj;
                     
                 end
             end
             
             
             
        end        
        
        if cellTypeArray(smallInd).numOfCells > 0
            sqLen = deployConfig.squareLength;
                    
            firstIndex = cellTypeArray(smallInd).firstUserId;
            varNumOfUes = cellTypeArray(smallInd).numOfUEs;
            lastIndex = firstIndex + varNumOfUes - 1;

            firstCellId = cellTypeArray(smallInd).firstCellId;
            lastCellId = firstCellId + ...
                cellTypeArray(smallInd).numOfCells - 1;

            posOfCells = [cls(firstCellId:lastCellId).complexPosition];
            minDistance = cellTypeArray(smallInd).minDistToUE;

            placeUE_SmallCell_v2(posOfCells,ues,sqLen,minDistance,...
                firstIndex,lastIndex,UEheight);

            drawPositions(ues,varNumOfUes,smallInd,firstIndex,...
                lastIndex,0,0,'userPos.txt',0);
        end
    end                                   
 end