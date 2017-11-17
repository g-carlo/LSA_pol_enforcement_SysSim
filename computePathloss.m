function [PLmatrix, fsPLoss] = computePathloss(varPlModel,C,freq,varIsLOS,height_rx,...
    height_tx, varAreaType, varDist, nCells,nUes)
    
    PLmatrix = zeros(nCells,nUes);
    fsPLoss = zeros(nCells,nUes);
    
    switch varPlModel        
        case {1,3,4}
            PLmatrix = C(1,1) + C(1,2) .* log10(varDist);                        
        case 2
            PLLOS = C(1,1) + C(1,2) .* log10(varDist);
            PLmatrix(varIsLOS==1) = PLLOS(varIsLOS==1);
            PLNLOS = C(2,1) + C(2,2) .* log10(varDist);
            PLmatrix(varIsLOS==0) = PLNLOS(varIsLOS==0);
        case {5,6}
            
            if varAreaType == 2
                varAreaType = 4;
            else
                varAreaType = 6;
            end
            
            freq = freq/1e6;
            [PLmatrix, fsPLoss] = ExtendedHata(varDist, nCells, nUes, freq, height_rx, ...
                height_tx, varAreaType);
    end          
end