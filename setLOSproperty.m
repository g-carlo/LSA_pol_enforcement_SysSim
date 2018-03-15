function varIsLOS = setLOSproperty(typeInd,varAreaType,varDist,phi_angle,...
    boresightMtx,strAnt)
    
    global macroInd
    global picoInd
    
    nCells = size(varDist,1);
    nUes = size(varDist,2);
    
    if typeInd == macroInd || typeInd == picoInd
        
        
        arrayOne = ones(nCells,nUes);
        
        if varAreaType == 1           
            if typeInd == macroInd
                varExp = exp(-varDist ./ 0.063);
                probLOS = min(0.018./varDist, arrayOne) .* (1 - varExp) + ...
                    varExp;
            else
                probLOS = 0.5 - min(0.5 .* arrayOne, 5.*exp(-0.156./varDist))...
                    + min(0.5 .* arrayOne, 5 .* exp(-varDist./0.03));       
            end        
        elseif varAreaType == 2
            if typeInd == macroInd           
                probLOS = exp(-(varDist-0.01));
            else
                probLOS = 0.5 - min(0.5 .* arrayOne, 3* exp(-0.3./varDist)) + ...
                    min(0.5 .* arrayOne, 3 .* exp(-varDist./0.095));
            end            
        end
        
%        load('carloSimulatorVars','varRandomLos');
        % COMMENT: JUST TO COMPARE TO CARLO SIMULATOR
         rndVarMatrix = rand(nCells,nUes);
        
%         if typeInd==macroInd
%             rndVarMatrix = varRandomLos(1:57,:);
%         else
%             rndVarMatrix = varRandomLos(58:285,:);
%         end
%         % 1 is LOS, 0 is NLOS

        varIsLOS = rndVarMatrix <= probLOS;    
        if typeInd == macroInd
            LOS_angle_matrix = abs(angle(exp(1i*(phi_angle-boresightMtx)))...
                ./pi .*180);
            varIsLOS(LOS_angle_matrix >= strAnt.horizontalPhi3dB * ...
                sqrt(strAnt.horizontalMinGain/12) ) = false;        
        end
    else
        varIsLOS = zeros(nCells,nUes);
    end
end