%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%          GENERATE MACRO HETEROGENEOUS SCENARIO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cls = initializeHetNetScenario_v2(cls,objType,objType2Type,...
    varIntersiteDist,numOfCarriers,timerQoSEV,numOfTTIperSimu,arrayAvPRBs,...
    maxTxPwArray)

    global picoInd;
    global macroInd;
    
    firstIndex = objType(picoInd).firstCellId;
    numOfCells = objType(picoInd).numOfCells;
    lastIndex = firstIndex + numOfCells - 1;
    
       
    % generation of pico cells;
    for n = firstIndex:lastIndex    
        
        
        cls(n) = cell_Cls(objType(picoInd),n,numOfCarriers,...
            maxTxPwArray,arrayAvPRBs,timerQoSEV,numOfTTIperSimu);        
    end
   
    % placement of pico cells; 

    placePicoCellsFnc_v2(cls,objType,objType2Type,varIntersiteDist);    
    
    
    numOfMacro = objType(macroInd).numOfCells;
    
    syncInitVect = 2 .* ones(numOfCarriers,numOfMacro);

    for n = firstIndex:lastIndex

        macInd = cls(n).macroCellBelongingID;
        
        for ii = 1:numOfCarriers
            cls(n).fsuVars(ii).syncFSUtimeSlot = syncInitVect(ii,macInd);
            syncInitVect(ii,macInd) = syncInitVect(ii,macInd) +1;
        end        
    end

    %drawPositions(cls,numOfCells,picoInd,...
     %   firstIndex,lastIndex,varIntersiteDist,0,'picoCellPositions.txt',1);    
end


function placePicoCellsFnc_v2(cls,objType,objType2Type,varIntersiteDist)

global macroInd;
global picoInd;

numOfMacroCells = objType(macroInd).numOfCells;
firstMacroIndex = objType(macroInd).firstCellId;
lastMacroIndex = firstMacroIndex + numOfMacroCells - 1;

numOfPicoCells = objType(picoInd).numOfCells;
picoHeight = objType(picoInd).height;
firstPicoIndex = objType(picoInd).firstCellId;

numOfPicoPerMacro = objType2Type(macroInd,picoInd).numOfOtherCells;
minDistMacroPico = objType2Type(macroInd,picoInd).minDistance;
minDistPicoPico = objType2Type(picoInd,picoInd).minDistance;

b = varIntersiteDist/2/(cos(pi/6));

posOfPico = zeros(1,numOfPicoCells);
ii = 1;
for n = firstMacroIndex:lastMacroIndex
    macroPos = cls(n).pos.getComplexCoordinate();
    macroOrient = cls(n).orientation;
        
    for k = 1:numOfPicoPerMacro        
        xy = macroPos + ((b-minDistMacroPico)*rand(1)+minDistMacroPico)...
            *exp( 1i*( 2/3*pi*rand(1)-pi/3  + macroOrient) );

        if ii==1         
            while ( ~cls(n).belongsToCell(xy,varIntersiteDist) || abs(xy - macroPos) < minDistMacroPico)
                xy = macroPos + ((b-minDistMacroPico)*rand(1)+minDistMacroPico)...
            *exp( 1i*( 2/3*pi*rand(1)-pi/3  + macroOrient) );
            end
        else
            while ( ~cls(n).belongsToCell(xy, varIntersiteDist) || ...
                    abs(xy - macroPos) < minDistMacroPico || ...
                    min(abs( xy - posOfPico(1:(ii-1)))) < minDistPicoPico)

                xy = macroPos + ((b-minDistMacroPico)*rand(1)+minDistMacroPico)...
            *exp( 1i*( 2/3*pi*rand(1)-pi/3  + macroOrient) );
            end
        end     % end of if

        posOfPico(ii) = xy;
        ii = ii + 1;
        ind = (firstPicoIndex - 1) + (n-1) * numOfPicoPerMacro + k;

        %%% THIS IS THE OUTPUT OF THE FUNCTION - BEGIN
        cls(ind).setPosition( real(xy) ,imag(xy) , picoHeight);
        cls(ind).macroCellBelongingID = n;
        %%% THIS IS THE OUTPUT OF THE FUNCTION - END
    end      
end
end