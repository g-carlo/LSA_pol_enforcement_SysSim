function [carrierArray,cellsPerCarrier,varPoints] = setCellsPerCarrier_v1(...
    carrierArray, varCellPosArray, varRndAssigCellsToCarriers,varCellRatio,cls_networkOP_ID)

    numOfCarriers = length(carrierArray);
    totNumCells = length(varCellPosArray);
    
    cellsPerCarrier = zeros(numOfCarriers,totNumCells);
    
       
    for kk = 1:numOfCarriers
        for op_jj = carrierArray(kk).operatorsInCarrier
            cellsPerCarrier(kk,cls_networkOP_ID == op_jj) = 1;
        end
    end
    
        
%     cellsPerCarrier(1,:) = 1;
%     varPoints = [];
%     
%     if varRndAssigCellsToCarriers == 0
%         varPoints = [-600-1i*450 600-1i*450 600+1i*450 ...
%           -600+1i*450 -600-1i*450];            
%         for ii = 2:numOfCarriers
%             cellsPerCarrier(ii,:) = determineCellsPerCarrier(...
%                 varCellPosArray,varPoints,1);          
%         end
%     else
%         varArrayRatio = varCellRatio;
%         for ii = 2:numOfCarriers
%             varRatio = varArrayRatio(ii);
%             varNumCells = floor(varRatio * totNumCells);
%             
%             varSelectedCells = randSequenceFromBin(varNumCells,totNumCells);
%             cellsPerCarrier(ii,varSelectedCells) = 1;      
%         end
%     end
%     
    carrierArray = assignCellsPerCarrier(carrierArray,cellsPerCarrier);

end

