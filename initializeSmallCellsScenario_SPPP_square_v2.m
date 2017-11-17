%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%          GENERATE SMALL CELL SCENARIO SPPP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [numOfCells,cls] = initializeSmallCellsScenario_SPPP_square_v2(cls,...
    varDeploy,objType,minDistSmallCells,numOfCarriers,timerQoSEV,numOfTTIperSimu,...
    arrayAvPRBs,maxTxPwArray)
    
    global smallInd;
    
    squareLength = varDeploy.squareLength;
    halfSqLen = squareLength/2;
    
    cellDensity = varDeploy.cellDensity;
    txHeight = objType.height;
    firstIndex = objType.firstCellId;
    
    
   % numOfCells = max(poissrnd(cellDensity*(squareLength^2),1));
   
   numOfCells = 20;
   
    lastIndex = firstIndex + numOfCells - 1;
    posOfNewCells = zeros(1,numOfCells);
    
    ii = 1;
    for cell_ind = firstIndex:lastIndex

        cls(cell_ind) = cell_Cls(objType,cell_ind,numOfCarriers,...
            maxTxPwArray,arrayAvPRBs,timerQoSEV,numOfTTIperSimu);
        
        xy = (squareLength*rand(1)-halfSqLen) + 1i*(squareLength*rand(1)-halfSqLen);
        if ii>1
            while min(abs( xy - posOfNewCells)) < minDistSmallCells 
                xy = (squareLength*rand(1)-halfSqLen) + 1i*(squareLength*rand(1)-halfSqLen);
            end
        end

        posOfNewCells(ii) = xy;
        ii = ii + 1;             
        cls(cell_ind).setPosition( real(xy) , imag(xy) , txHeight);             
    end 
    
    objType.numOfCells = numOfCells;
    
    drawPositions(cls,numOfCells,smallInd,firstIndex,lastIndex,...
        0,0,'smallCellPositions.txt',1);   
end
