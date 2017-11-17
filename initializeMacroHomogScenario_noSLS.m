%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%          GENERATE MACRO HOMOGENEOUS SCENARIO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cls = initializeMacroHomogScenario_noSLS(cls,objCellType,...
    varIntersiteDist,numOfCarriers,maxTxPwArray,varDeploy,N_sites_extedndedScenario)

    varNumOfCellsPerSite = varDeploy.numOfCellsPerSite;
    numOfCells = objCellType.numOfCells;
    firstAngle = varDeploy.firstSectorAngle;
    angleWidth = 2*pi/varDeploy.numOfSectors;
    firstIndex = objCellType.firstCellId;
    lastIndex = firstIndex + numOfCells - 1;
    
    % generation of macro cells;
    for n = firstIndex:lastIndex   
        
        cls(n) = cell_Cls(objCellType,n,numOfCarriers,maxTxPwArray);          
        
        sectorId = rem(n-1,varDeploy.numOfSectors);
        cls(n).macroCellSectorIndex = sectorId+1;
        cls(n).orientation = rem(firstAngle + sectorId*angleWidth,2*pi);    
        
    end
   
    % placement of macro cells;
    placeMacroCells_andAssigntToNO_Fnc_noSLS(cls,objCellType,varIntersiteDist,varNumOfCellsPerSite,N_sites_extedndedScenario);
    
 
end

function placeMacroCells_andAssigntToNO_Fnc_noSLS(cls,objCellType,varIntersiteDist,numCellsPerSite,N_sites_extedndedScenario)

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

if N_sites_extedndedScenario == 11^2

    pos_vc_y = l*1i*(5:-1:-5).';
    pos_vc_y_coordinate = [ repmat( [pos_vc_y l/2*1i+pos_vc_y] ,1 ,5 ) pos_vc_y];

    pos_vc_x = l*sqrt(3)/2*(-5:5);
    pos_vc_x_coordinate = repmat( pos_vc_x ,11 ,1 );       

    pos_sites = pos_vc_x_coordinate + pos_vc_y_coordinate;

    pos_cells = repmat( pos_sites(:).', [numCellsPerSite 1 ] );
    pos  = pos_cells(:);

elseif N_sites_extedndedScenario == 21^2
    
    pos_vc_y = l*1i*(10:-1:-10).';
    pos_vc_y_coordinate = [ repmat( [pos_vc_y l/2*1i+pos_vc_y] ,1 ,10 ) pos_vc_y];
    
    pos_vc_x = l*sqrt(3)/2*(-10:10);
    pos_vc_x_coordinate = repmat( pos_vc_x ,21 ,1 );       
    
    % this is the matrix with the positions of the sites
    pos_sites = pos_vc_x_coordinate + pos_vc_y_coordinate;

    % we need to generate the position of the cells
    pos_cells = repmat( pos_sites(:).', [numCellsPerSite 1 ] );
    pos  = pos_cells(:);
    
end

   
% cell_angle = [zeros(1,3) reshape(repmat(30:60:330,numCellsPerSite,1),1 ,length(30:60:330)*numCellsPerSite )...
%     reshape(repmat(30:30:360,numCellsPerSite,1),1 ,length(30:30:360)*numCellsPerSite )]*pi/180;
% 
% cell_radius = [zeros(1,3) l*ones(1,6*numCellsPerSite) ...
%     reshape(repmat(reshape( repmat([2*l (l/cos(pi/6)+l*tan(pi/6))].',1,6),1,12),3,1),1,12*3)];



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