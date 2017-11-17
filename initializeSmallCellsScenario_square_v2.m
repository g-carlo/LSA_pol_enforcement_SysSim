%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%          GENERATE SMALL CELL SCENARIO SPPP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varDeploy,numOfCells, cls] = initializeSmallCellsScenario_square_v2(cls,...
    varDeploy,objType,numOfCarriers,timerQoSEV,numOfTTIperSimu,arrayAvPRBs,...
    maxTxPwArray)
    
    global smallInd;
    
    squareLength = varDeploy.squareLength;
    cell_density = varDeploy.cellDensity;
    txHeight = objType.height;
    firstIndex = objType.firstCellId;
    
    numOfCell_1 = cell_density*(squareLength^2);
    numOfCell_perEdge = ceil( ( sqrt(numOfCell_1) -1) / 2 ) * 2+1;
    n_perEdge = (numOfCell_perEdge - 1)/2;
    squareArea = numOfCell_perEdge^2 / numOfCell_1 * squareLength^2;
    intersiteDistance = sqrt( squareArea ) / numOfCell_perEdge ;
    numOfCells = numOfCell_perEdge^2;
    
    
    outSquareLength = sqrt( squareArea );
   
    pos_grid = intersiteDistance * repmat( (-n_perEdge : n_perEdge), numOfCell_perEdge,1 ) + ...
        1i * intersiteDistance * repmat(flipud((-n_perEdge : n_perEdge).'), 1 , numOfCell_perEdge );

    lastIndex = firstIndex + numOfCells - 1;
    ii = 1;
    for cell_id = firstIndex:lastIndex
        
        cls(cell_id) = cell_Cls(objType,cell_id,numOfCarriers,...
            maxTxPwArray,arrayAvPRBs,timerQoSEV,numOfTTIperSimu);

        xy = pos_grid(ii);     

        %%% THIS IS THE OUTPUT OF THE FUNCTION - BEGIN
        cls(cell_id).setPosition( real(xy) , imag(xy) , txHeight);
        %%% THIS IS THE OUTPUT OF THE FUNCTION - END
        
        ii = ii + 1;
    end
    varDeploy.squareLength = outSquareLength;
    
    objType.numOfCells = numOfCells;
    
    drawPositions(cls,numOfCells,smallInd,firstIndex,lastIndex,...
       0,0,'smallCellPositions.txt',1);   
 end

