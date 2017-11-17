function sfCoeffMtx = computeShadowing(constShadowDev, typeInd,...
    nCells, nUes, numOfSites, numOfMacroPerSite, numOfOperators,varAreaType, freq,...
    varDist, pathloss, fsPathloss, varPlModel)
    
    global macroInd
    global smallInd
    
    stdDev = constShadowDev;
    
    if varPlModel == 5
        stdDev = 3.5  .* ones(nCells,nUes);
        
        vartmp = varDist>0.04;
        stdDev(vartmp) = 3.5 + ((12-3.5)/(0.1-0.04) .* (varDist(vartmp)-0.04));
        
        vartmp = varDist>0.1;
        stdDev(vartmp) = 12;
        
        vartmp = varDist>0.2;
        stdDev(vartmp) = 12+(9-12)/(0.6-0.2)*(varDist(vartmp)-0.2);
        
        vartmp = varDist>0.6;
        stdDev(vartmp) = 9;        
            
    elseif varPlModel == 6
        
        if varAreaType == 2
            varAreaType = 4;
        else
            varAreaType = 6;
        end
        
        switch varAreaType
            case {7, 6}
                A = 5.2;
            otherwise
                A = 6.6;
        end

        stdDev = 0.65 * (log10(freq)^2) - 1.3*log10(freq) + A; 
        stdDev = stdDev .* ones(nCells,nUes);
        vartmp = (pathloss - fsPathloss)./2.326348;
        stdDev = min(stdDev, vartmp);    
    end
    
    
    if typeInd == smallInd
        sfCoeff = randn(nCells, nUes);
        sfCoeffMtx = stdDev .* sfCoeff;        
    else
        if typeInd == macroInd
            randVals = randn(1,numOfSites);
            sfCoeffCells = repmat( reshape( repmat( randVals, ...
                numOfMacroPerSite * numOfOperators, 1), numOfSites * numOfMacroPerSite * numOfOperators,1) ,1, nUes);
        else
            randVals = randn(nCells,1);
            sfCoeffCells = repmat( randVals, 1, nUes);
        end            
        sfCoeffUes = repmat(randn( 1 , nUes),nCells,1); 
        
        % COMMENT: JUST TO COMPARE TO CARLO SIMULATOR
%         if typeInd == macroInd
%             load('carloSimulatorVars','varRandShadowMacro','varRandUeShadow');
%             sfCoeffCells = varRandShadowMacro;
%         end
%         
%         if typeInd == picoInd
%             load('carloSimulatorVars','varRandShadowPico','varRandUeShadow');
%             sfCoeffCells = varRandShadowPico;
%         end
%         
%         sfCoeffUes = repmat(varRandUeShadow,nCells,1);
        sfCoeffMtx = sqrt(0.5) * stdDev .* (sfCoeffCells + sfCoeffUes);
    end
    
end 
