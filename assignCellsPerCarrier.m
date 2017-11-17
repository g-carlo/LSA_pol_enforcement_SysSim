function carrierArray = assignCellsPerCarrier( carrierArray,cellsPerCarrier)
    
    numOfCarriers = length(carrierArray);
    
    for ii =  1:numOfCarriers
        carrierArray(ii).cellsInCarrier = find(cellsPerCarrier(ii,:)==1);            
    end
end

