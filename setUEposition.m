function setUEposition(ues,UEheight,cellTypes,cellTypeArray,initC,cls,...
    firstUEindex,currCellType,deploy)

    numOfUsersPerCell = deploy.numOfUEperCellType(currCellType);
    cellFirstIndex = cellTypeArray(currCellType).firstCellId;
    
    numOfCellTypes = length(cellTypes);
    cellPos = cell(numOfCellTypes,1);
    minDist = zeros(1,numOfCellTypes);
    varIntersiteDist = deploy.intersiteDistance;
    
    l = varIntersiteDist/2;
    b = l/(cos(pi/6));       
    
    ind = 1;
       
    for ii=1:numOfCellTypes        
        typeIndex = cellTypes(ii);
        
        % number of a given cell type
        numOfCells = cellTypeArray(typeIndex).numOfCells;
        
        if numOfCells>0
        
        % minimum distance between a given cell type and the UE
            minDist(ii) = cellTypeArray(typeIndex).minDistToUE;

            % positions of a given cell type
            posArray = [cls(ind:(ind+numOfCells-1)).pos];
            cellPos{ii} = [posArray(:).x] + 1i*[posArray(:).y];        

            ind = ind + numOfCells;  
        end
    end
    
    
    cellLastIndex = cellFirstIndex + cellTypeArray(currCellType).numOfCells - 1;
    for n = cellFirstIndex:cellLastIndex
        varCellPos = cls(n).pos.getComplexCoordinate();
        sectorIndex = cls(n).macroCellSectorIndex;
        C(1:4) = initC(1:4);
        switch sectorIndex
            case 1
                C(2) = 0;                        
            case 3
                C(1) = 2*b;                        
                C(3) = l;
        end
        for k = 1:numOfUsersPerCell      
            tmp = 1;
            while tmp    
                
                varRand1 = rand(1);
                varRand2 = rand(1);
                xy = varCellPos + (C(1)*varRand1-C(2)) + 1i*(C(3)*varRand2-C(4));
                
                if cls(n).belongsToCell(xy,varIntersiteDist)                    
                    tmp = 0;
                    ii = 1;
                    varbreak = 0;
                    while ii<=numOfCellTypes && varbreak==0
                        if ~isempty(cellPos{ii}) 
                            minVarDist = min(abs(xy - cellPos{ii}));
                            if  minVarDist >= minDist(ii) 
                                ii = ii + 1;
                            else
                                varbreak = 1;
                            end
                        else
                            ii = ii + 1;
                        end
                    end                    
                    if ii<=numOfCellTypes
                        tmp = 1;
                    end
                end
            end
                     
            ind = (firstUEindex-1) + (n - cellFirstIndex) * numOfUsersPerCell + k;

            %%% THIS IS THE OUTPUT OF THE FUNCTION - BEGIN
            ues(ind).setPosition( real(xy) ,imag(xy) , UEheight);
            %%% THIS IS THE OUTPUT OF THE FUNCTION - END
        end    % end of for
    end     %%%%%%%%%%%%%%%%%%%%%%%
end