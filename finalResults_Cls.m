classdef finalResults_Cls < handle
   
    properties
        meanUEtrafficLoad;
        cellLoad;
        meanCellLoad;
        uth;
        cellTH;
    end
    
    methods
        function obj = finalResults_Cls()
        end
        
        function computeFinalResults(obj,ues,cls, TTIduration,trafficCfg,...
                userVsCell_allocationMtx, numOfPrbPerUe,totNumIt,nMacroTh,...
                cellTypeArray, varNumPerType,numOfNOs)

            global macroInd;
            global picoInd;
            global smallInd;            
            
            totNumCells = length(cls);
                       
            totNumMacro = cellTypeArray(macroInd).numOfCells;
            totNumPico = cellTypeArray(picoInd).numOfCells;
            
            varCellsId = [];
            varUesId = [];
            
            if totNumMacro>0 
                varCellsId = 1:nMacroTh;
                varUesId = 1:(nMacroTh*varNumPerType(macroInd)*numOfNOs);
            end
            
            if totNumPico>0
                nPicoPerMacro = totNumPico/totNumMacro; 
                nPicoTh = nMacroTh * nPicoPerMacro;
                lstPicoTh = totNumMacro + nPicoTh;
                varCellsId = [varCellsId totNumMacro+1:lstPicoTh];
                
                fstPicoUe = cellTypeArray(picoInd).firstUserId;
                lstPicoUe = fstPicoUe + nPicoTh*varNumPerType(picoInd)-1;
                varUesId = [varUesId fstPicoUe:lstPicoUe];
            end
            
            totNumSmall = cellTypeArray(smallInd).numOfCells;
            
            if totNumSmall >0
                
                if totNumMacro>0
                    nSmallPerMacro = totNumSmall/totNumMacro; 
                    nSmallTh = nMacroTh * nSmallPerMacro;
                else
                    nSmallTh = totNumSmall;
                end
                lstSmallTh = totNumMacro + totNumPico + nSmallTh;
                varCellsId = [varCellsId totNumMacro+totNumPico+1:lstSmallTh];
                
                fstSmallUe = cellTypeArray(smallInd).firstUserId;
                lstSmallUe = fstSmallUe + nSmallTh*varNumPerType(smallInd)-1;
                varUesId = [varUesId fstSmallUe:lstSmallUe];
            end
            
            numUesTh = length(varUesId);
            numCellsTh = length(varCellsId);
            
            %%% UE throughput
            dataTX = [ues(varUesId).transmittedDataCnt]; 
            timeElapsed = [ues(varUesId).ueTXtimer] .* TTIduration;
            obj.uth = dataTX ./ timeElapsed;
            obj.uth(dataTX <= 0) = 0;

            %%% Cell throughput 
            cellDataTX = [cls(varCellsId).cumDataTXperCell];
            activeTTIperCell = [cls(varCellsId).activeTTIcntPerCell] .* TTIduration;
            obj.cellTH = cellDataTX ./ activeTTIperCell;
            obj.cellTH(cellDataTX<=0) = 0;

            avgNumPrbPerUe = numOfPrbPerUe(varUesId)./totNumIt;
            varAvgPrbAllUes = sum(avgNumPrbPerUe)/numUesTh;
            varAvgThAllUes = sum(obj.uth)/numUesTh;
            varAvgThAllCells = sum(obj.cellTH)/numCellsTh;
                      
            %%% network load
            if trafficCfg.trafficGeneratorType > 1

                trafficGenrationRate = [ues.pckCnt] .* trafficCfg.FTP_model1_pckSize/(numOfTTIperSimu .* TTIduration);
                meanTrafficGenerationRate = mean(trafficGenrationRate);
                trafficDeliveryRate = mean(obj.uth);
                obj.meanUEtrafficLoad = meanTrafficGenerationRate./trafficDeliveryRate;
                
                cellTrafficGenerationRate = sum(repmat(trafficGenrationRate,totNumCells,1) .* userVsCell_allocationMtx,2).';
                obj.cellLoad = cellTrafficGenerationRate ./ obj.cellTH;
                obj.cellLoad(cellTrafficGenerationRate <= 0 ) = 0;
                
                obj.meanCellLoad = mean(obj.cellLoad);
                
            end 
            
% % %             val1 = min(obj.uth);
% % %             val2 = max(obj.uth);
% % %             val3 = min(obj.cellTH);
% % %             val4 = max(obj.cellTH);
% % %             val5 = min(avgNumPrbPerUe);
% % %             val6 = max(avgNumPrbPerUe);
% % %             
% % %             figure 
% % % 
% % %             % varX = size(obj.uth,2);
% % %             hist(obj.uth);
% % %             % plot(1:varX, obj.uth);
% % %             title({'User throughput (Multiple carriers - Simulator)';...
% % %                 ['Average Th over all UEs: ', num2str(varAvgThAllUes)]});
% % %             % xlabel('UE Identifier');
% % %             % ylabel('Throughput');
% % %             
% % %             xlabel('Average UE Throughput');
% % %             ylabel('Number of UEs');
% % %                         
% % %             figure 
% % % 
% % %             % varX = size(obj.cellTH,2);
% % %             hist(obj.cellTH);
% % %             % plot(1:varX,obj.cellTH);
% % %             title({'Cell Throughput (Multiple Carriers - Simulator)';...
% % %                 ['Average Th over all Cells: ', num2str(varAvgThAllCells)]});
% % %             
% % %             xlabel('Average Cell Throughput');
% % %             ylabel('Number of Cells');
% % %             
% % %             % xlabel('Cell Identifier');
% % %             % ylabel('Throughput');
% % %             
% % %             figure 
% % %             % varX = size(avgNumPrbPerUe,2);
% % %             hist(avgNumPrbPerUe);
% % %             % plot(1:varX, avgNumPrbPerUe);
% % %             title({'Average number of PRBs per UE (Multiple Carriers - Simulator)';...
% % %                 ['Avg number of assigned PRBs (over all UEs): ',num2str(varAvgPrbAllUes)]});
% % %             
% % %             % xlabel('UE Identifier');
% % %             % ylabel('Average number of assigned PRBs');
% % %             
% % %             xlabel('Average number of assigned PRBs to each UE');
% % %             ylabel('Number of UEs');
% % %             
        end
    end
    
end

