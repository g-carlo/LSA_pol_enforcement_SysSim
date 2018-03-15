function [carrierArray,cellsPerCarrier,varPoints] = setCellsPerCarrier(...
    carrierArray, varCellPosArray, varRndAssigCellsToCarriers,varCellRatio)

    numOfCarriers = length(carrierArray);
    totNumCells = length(varCellPosArray);
    
    cellsPerCarrier = zeros(numOfCarriers,totNumCells);
    cellsPerCarrier(1,:) = 1;
    varPoints = [];
    
    if varRndAssigCellsToCarriers == 0
        varPoints = [-600-1i*450 600-1i*450 600+1i*450 ...
          -600+1i*450 -600-1i*450];            
        for ii = 2:numOfCarriers
            cellsPerCarrier(ii,:) = determineCellsPerCarrier(...
                varCellPosArray,varPoints,1);          
        end
    else
        varArrayRatio = varCellRatio;
        for ii = 2:numOfCarriers
            varRatio = varArrayRatio(ii);
            varNumCells = floor(varRatio * totNumCells);
            
            varSelectedCells = randSequenceFromBin(varNumCells,totNumCells);
            cellsPerCarrier(ii,varSelectedCells) = 1;      
        end
    end
    
    carrierArray = assignCellsPerCarrier(carrierArray,cellsPerCarrier);

end

