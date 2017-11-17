% this function receives:
% <varCellPosArray>: it is an array of cell positions. 
% Each element of this array is a complex number, where the real and 
% imaginary part indicate the coordinates x and y, respectively. These 
% coordinates correspond to the position of each cell.

% <varPoints>: It is an array that contains points that define the polygon.
% Each point should be given as a complex number. The number of points must
% be greater than 2.

% <incCellsInPoly>: this variable must be 1, to indicate that cells inside
% the polygon can use the carrier. In order to indicate the opposite: cells
% that are inside the polygon cannot use the carrier, this input parameter
% must be 0.
function cellsPerCarrier = determineCellsPerCarrier(varCellPosArray, ...
    varPoints, incCellsInPoly)
    
    if length(varPoints)>2
        
        varX = real(varPoints);
        varY = imag(varPoints);

        nCells = length(varCellPosArray);

        cellsPerCarrier = ones(1,nCells) .* (1-incCellsInPoly);

        for j=1:nCells
            varCellPos = varCellPosArray(j);
            varAns = inpolygon(real(varCellPos),imag(varCellPos),varX,varY);
            if varAns == incCellsInPoly
                cellsPerCarrier(j) = incCellsInPoly;
            end
        end
    else
        error('At least 3 points must be given to define the work area.');
    end
end

