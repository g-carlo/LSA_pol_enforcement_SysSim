function computeResults(objCarrier,userVsCell_allocationMtx,...
    varSchedArray, varRaArray, resultCfg,nOfCel)

    if resultCfg.additionalData
        nOfPRB = objCarrier.numOfAvailablePRB;
        
        mt = zeros(nOfCel,nOfPRB);
        
        for n = 1:nOfCel
            
            % UEs that are attached to the cell n and that uses some of the
            % PRBs in the current carrier
            tmp = varSchedArray(n).userSpectVect;
            
            % no UE that is attached to the cell n has been scheduled
            if sum(tmp>0)<=0
                mt(n,:) = zeros(1,nOfPRB);
            else
                % some UE that is attached to the cell n has been scheduled
                % PRBs allocated to cell n
                av_sp = varRaArray(n).spectrumAvailable;
                varLen = length(av_sp);
                for prbIndex = 1:varLen
                    prbID = av_sp(prbIndex);
                    if tmp(prbIndex) <= 0
                        mt(n,prbID) = 0;
                    else
                        mt(n,prbID) = objCarrier.SINR_mtx(tmp(prbIndex),prbID);
                    end
                end
            end 

        end

        objCarrier.cellVsPRB_SINR_mtx = mt;

       
        if resultCfg.thPerCell > 0 || resultCfg.SINR_statistic > 0
            
%             numAvPrb = objCarrier.numOfAvailablePRB;
            
            for n = 1:nOfCel

                if resultCfg.thPerCell > 0
                    tmp = objCarrier.thPerUE(userVsCell_allocationMtx(n,:)>0);
                    if isempty(tmp)
                        objCarrier.cellOutage(n) = 0;
                        objCarrier.cellThroughput(n) = 0;        
                    else
                        [xTarget, ~] = CDFfncMkII(tmp,5);
                        objCarrier.cellOutage(n) = xTarget;
                        objCarrier.cellThroughput(n) = sum(tmp);
                    end
                end

%                 if resultCfg.SINR_statistic > 0
% 
%                     tmp1 = 1:numAvPrb;
%                     prb_list = tmp1(varSchedArray(n).userSpectVect>0);
% 
%                     if isempty(prb_list)
%                         objCarrier.cellVsPRB_SINR_mtx(n,:) = zeros(1,numAvPrb);
%                     else
%                         for prb_index = prb_list
%                             objCarrier.cellVsPRB_SINR_mtx(n,prb_index) = objCarrier.SINR_mtx(varSchedArray(n).userSpectVect(prb_index),prb_index);
%                         end
%                     end   
%                 end

            end     % end of for

        end      % end of IF (config.results.thPerCell > 0 ...
    end
 end

