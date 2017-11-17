%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%          GENERATE SMALL CELL SCENARIO SPPP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [numOfCells, cls] = initSmallCellsSPPPsquareMaxLen_v2(cls,...
    varDeploy,objType,minDistSmallCells,objMacro,numOfCarriers,...
    timerQoSEV,numOfTTIperSimu,arrayAvPRBs,maxTxPwArray)
    
    global smallInd;
    
    firstIndex = objType.firstCellId;
    txHeight = objType.height;
    
    firstIndexMacro = objMacro.firstCellId;
    lastIndexMacro = firstIndexMacro + objMacro.numOfCells - 1;
    
    squareLength = varDeploy.squareLength;
    halfSqLen = squareLength/2;
    
    cell_density = varDeploy.cellDensity;
    varIntersiteDist = varDeploy.intersiteDistance;
        
    % numOfCells = max(poissrnd(cell_density*(squareLength^2)),1);    
    numOfCells = 20;
    lastIndex = firstIndex + numOfCells - 1;
    posOfNewCells = zeros(1,numOfCells);
    ii = 1;
    for cell_id = firstIndex:lastIndex
        
        cls(cell_id) = cell_Cls(objType,cell_id,numOfCarriers,...
            maxTxPwArray,arrayAvPRBs,timerQoSEV,numOfTTIperSimu);
        
        tmp = 1;

        if ii==1
            while tmp
                xy = (squareLength*rand(1)-halfSqLen) + 1i*(squareLength*rand(1)-halfSqLen);
                j = firstIndexMacro;
                while j<=lastIndexMacro && cls(j).belongsToCell(xy,varIntersiteDist)==0
                    j = j + 1;
                end
                if j<=lastIndexMacro
                    tmp = 0;
                    cls(cell_id).macroCellBelongingID = j;
                end
            end
        else
            while tmp
                xy = (squareLength*rand(1)-halfSqLen) + 1i*(squareLength*rand(1)-halfSqLen);
                                
                if min(abs( xy - posOfNewCells)) >= minDistSmallCells
                    j=firstIndexMacro;
                    while j<=lastIndexMacro && cls(j).belongsToCell(xy,varIntersiteDist)==0
                        j = j + 1;
                    end
                    if j<=lastIndexMacro
                        tmp = 0;
                        cls(cell_id).macroCellBelongingID = j;
                    end
                end
            end
        end 

        posOfNewCells(ii) = xy;
        ii = ii + 1;
        cls(cell_id).setPosition( real(xy) , imag(xy) , txHeight);        
    end
    objType.numOfCells = numOfCells;
    
   drawPositions(cls,numOfCells,smallInd,firstIndex,lastIndex,...
       0,0,'smallCellPositions.txt',1);   
end
