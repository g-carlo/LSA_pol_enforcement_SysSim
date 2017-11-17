%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%          GENERATE SMALL CELL SCENARIO SPPP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cls = initializeSmallCellsScenario_SPPP_square_numPerMacro_v2(cls,...
    objType,type2type,varIntersiteDist,numOfCarriers,timerQoSEV,numOfTTIperSimu,...
    arrayAvPRBs,maxTxPwArray)

    global macroInd;
    global smallInd;
    
    firstIndex = objType(smallInd).firstCellId;
    smallCellsHeight = objType(smallInd).height;
    description = objType(smallInd).description;
    numOfSmallCells = objType(smallInd).numOfCells;
    
    numOfMacroCells = objType(macroInd).numOfCells;
    
    minDistSmallCells = type2type(smallInd,smallInd).minDistance;
    smallCellsPerMacro = type2type(macroInd,smallInd).numOfOtherCells;
    
    posOfNewCells = zeros(1,numOfSmallCells);

    l = varIntersiteDist/2;     % radius of inscribed circle
    b = l/(cos(pi/6));          % radius of circumscribed circle

    ii = 1;
    varTmpDist = 0;
    cell_id = firstIndex;
    for p=1:smallCellsPerMacro
        for j=1:numOfMacroCells
            
            macroPosition = cls(j).pos.getComplexCoordinate();
            sectorIndex = cls(j).macroCellSectorIndex;
            
            xy = get_smallCellPosInSectorOfMacro(macroPosition,...
                sectorIndex,b,l);            
            
           while min(abs( xy - posOfNewCells)) < varTmpDist || ...
            cls(j).belongsToCell(xy,varIntersiteDist)==0
                xy = get_smallCellPosInSectorOfMacro(macroPosition,...
                    sectorIndex,b,l);
            end
            

            posOfNewCells(ii) = xy;
            ii = ii + 1;
            
            cls(cell_id) = cell_Cls(objType(smallInd),cell_id,numOfCarriers,...
                maxTxPwArray,arrayAvPRBs,timerQoSEV,numOfTTIperSimu);
            cls(cell_id).setPosition( real(xy) , imag(xy) , smallCellsHeight);
            cls(cell_id).macroCellBelongingID = j;
            varTmpDist = minDistSmallCells;
            cell_id = cell_id+1;
        end
    end
    
    % drawPositions(cls,numOfSmallCells,smallInd,firstIndex,lastIndex,...
      %  0,0,'smallCellPositions.txt',1);   
end

function xy = get_smallCellPosInSectorOfMacro(macroPosition,...
                sectorIndex,b,l)
    switch sectorIndex
        case 1
            xy = macroPosition + (b*rand(1)) + 1i*(2*l*rand(1)-l);
        case 2
            xy = macroPosition + (b*rand(1)-b) + 1i*(2*l*rand(1)-l);
        case 3
            xy = macroPosition + (2*b*rand(1)-b) + 1i*(l*rand(1)-l);
    end;
end
