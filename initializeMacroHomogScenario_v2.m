%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%          GENERATE MACRO HOMOGENEOUS SCENARIO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cls = initializeMacroHomogScenario_v2(cls,objCellType,...
    varIntersiteDist,numOfCarriers,picoAllocType,timerQoSEV,numOfTTIperSimu,...
    arrayAvPRBs,maxTxPwArray,varDeploy)
    
    
    varNumOfCellsPerSite = varDeploy.numOfCellsPerSite;
    numOfCells = objCellType.numOfCells;
    firstAngle = varDeploy.firstSectorAngle;
    angleWidth = 2*pi/varDeploy.numOfSectors;
    firstIndex = objCellType.firstCellId;
    lastIndex = firstIndex + numOfCells - 1;
    
    % generation of macro cells;
    for n = firstIndex:lastIndex   
        
        cls(n) = cell_Cls(objCellType,n,numOfCarriers,maxTxPwArray,...
            arrayAvPRBs,timerQoSEV,numOfTTIperSimu);          
        
        sectorId = rem(n-1,varDeploy.numOfSectors);
        cls(n).macroCellSectorIndex = sectorId+1;
        cls(n).orientation = rem(firstAngle + sectorId*angleWidth,2*pi);    
        for ii = 1:numOfCarriers
            cls(n).fsuVars(ii).syncFSUtimeSlot = 1;
        end
        
    end
   
    % placement of macro cells;
%     placeMacroCellsFnc_v2(cls,objCellType,varIntersiteDist,varNumOfCellsPerSite);    
    placeMacroCells_andAssigntToNO_Fnc_v1(cls,objCellType,varIntersiteDist,varNumOfCellsPerSite);
    
%     drawPositions(cls,numOfCells,objCellType.index,firstIndex,lastIndex,...
%         varIntersiteDist,varNumOfCellsPerSite,'cellPositions.txt',...
%         1);
    
    % assign cell to network operator
        
%     switch varDeploy.numOfNetworkOperators
% 
% 
%         case 1
%             for n = firstIndex:lastIndex   
%                     cls(n).networkOP_ID = 1;          
%             end
%             
%         case 2
%             op_1_idx_vector = [1:12 22:36];
% %                 op_2_idx_vector = [13:21 37:57];
%             for n = firstIndex:lastIndex   
%                 if any( op_1_idx_vector  == n)  
%                     cls(n).networkOP_ID = 1;
%                 else 
%                     cls(n).networkOP_ID = 2;
%                 end
%             end
% 
%         case 3
%             op_1_idx_vector = [1 4:9 22:30 55:57];
%             op_2_idx_vector = [2 10:15 31:42 ];
%             for n = firstIndex:lastIndex   
%                 if any( op_1_idx_vector  == n)  
%                     cls(n).networkOP_ID = 1;
%                 elseif any( op_2_idx_vector  == n)  
%                     cls(n).networkOP_ID = 2;
%                 else
%                     cls(n).networkOP_ID = 3;
%                 end
%             end
% 
%         otherwise
%             error('Scenario with more than 3 NOs not implemented');
% 
%     end   
              
end

function placeMacroCells_andAssigntToNO_Fnc_v1(cls,objCellType,varIntersiteDist,numCellsPerSite)

% global cls;
% global config;
txHeight = objCellType.height;
numOfCells = objCellType.numOfCells;
numOfCellsPerOperator = objCellType.numOfCellsPerOperator;
numOfOperators = numOfCells/numOfCellsPerOperator;
if rem( numOfCells,numOfCellsPerOperator) ~= 0
    error(' numOfCells must be a multiple of numOfCellsPerOperator');
end
l = varIntersiteDist;

% placement of macro cells in 19 sites, hexagonal shape


   
cell_angle = [zeros(1,3) reshape(repmat(30:60:330,numCellsPerSite,1),1 ,length(30:60:330)*numCellsPerSite )...
    reshape(repmat(30:30:360,numCellsPerSite,1),1 ,length(30:30:360)*numCellsPerSite )]*pi/180;

cell_radius = [zeros(1,3) l*ones(1,6*numCellsPerSite) ...
    reshape(repmat(reshape( repmat([2*l (l/cos(pi/6)+l*tan(pi/6))].',1,6),1,12),3,1),1,12*3)];

pos = cell_radius.*exp(1i*cell_angle);

for k = 1: numOfOperators
    
    for n = 1:numOfCellsPerOperator
            
        cell_idx = n + (k - 1)*numOfCellsPerOperator;

        cls(cell_idx).setPosition( real(pos(n)) ,imag(pos(n)) , txHeight);
        cls(cell_idx).networkOP_ID = k;
    
    end

end


        
%         cell_angle = [zeros(1,3) reshape(repmat(30:60:330,config.general.macroPerSite,1),1 ,length(30:60:330)*config.general.macroPerSite )]*pi/180;
% 
%         cell_radius = [zeros(1,3) l*ones(1,6*config.general.macroPerSite) ];
% 
%         pos = cell_radius.*exp(1i*cell_angle);
% 
%         for n = 1:config.general.numOfMacroCells
% 
%             cls(n).setPosition( real(pos(n)) ,imag(pos(n)) , config.macro.height);
% 
%         end

end

  
function placeMacroCellsFnc_v2(cls,objCellType,varIntersiteDist,numCellsPerSite)

    numOfCells = objCellType.numOfCells;
    txHeight = objCellType.height;
   
    firstIndex = objCellType.firstCellId;
    
    % placement of macro cells in N sites, hexagonal shape

    %varIntersiteDist = config.general.intersiteDistance;
    %numCellsPerSite = config.general.macroPerSite;
    %txHeight = config.macro.height;

    %numOfCells = 111;
    %varIntersiteDist = 500;
    %numCellsPerSite = 3;
    %cls = zeros(1,numOfCells);

    firstAngle = pi/6;
    tierNumber = 1;
    incX(1:6) = [-1 -1 0 1 1 0];
    varY(1:6) = [0 1 0 0 -1 0];
    varConst = cos(pi/6);
    varRatio = varIntersiteDist*varConst.*incX;

    varInd = firstIndex + numCellsPerSite - 1;
    lastIndex = firstIndex + numOfCells - 1;
    
    
    
    while varInd<lastIndex
        numOfSites = 6*tierNumber;
        varAngleResolution = 2*pi/numOfSites; %radians
        lastAngle = firstAngle + varAngleResolution*(numOfSites-1);
        varAngleTan = firstAngle:varAngleResolution:lastAngle;
        varTmp = mod(varAngleTan,pi/2)==0;
        varAngleTan(~varTmp) = tan(varAngleTan(~varTmp));
        varAngleTan(varTmp) = 0;

        varTierSide = tierNumber*varIntersiteDist;
        varXtierHex(1:6) = varTierSide*varConst.*[1 0 -1 -1 0 1];
        varTmp = varY.*varTierSide;

        varSite = 1;
        for k=1:6        
            for s2 = 0:(tierNumber-1)
                varX = varXtierHex(k) + varRatio(k)*s2;            
                %cls(varInd+1:varInd+numCellsPerSite)= complex(varX,varAngleTan(varSite)*varX + varTmp(k));
                for n=1:numCellsPerSite
                    cls(varInd+n).setPosition(varX ,varAngleTan(varSite)*varX + varTmp(k), txHeight)
                end
                varInd = varInd + numCellsPerSite;
                varSite = varSite + 1;
                varTmp(k) = 0;
            end
        end   
        tierNumber = tierNumber + 1;
    end
end