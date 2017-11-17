%%% SON SIMULATOR: attach user to the cell function
%
%   Author: Carlo Galiotto
%
%   Date:   August 2010
%   Last modification: April, 2015
%
%   Description: this function attach the users to the cell

%conf.general.Prx_threshold
%conf.general.attachUserMethod

function [outUeCellAlloc, outInterfIndMtx] = attachUsersFnc_v2(cls, ues, ...
    carrierArray,strGen, inFFflag, varDistCellUeNoWA)

    global macroInd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         SIMULATION MODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nCells = length(cls);
    nUes = length(ues);
    
    inPrxThreshold = strGen.Prx_threshold;
    inAttachUeMethod = strGen.attachUsersMethod;
    inInterfPerUe = strGen.interferersPerUE;
    
    outUeCellAlloc = false(nCells,nUes);
    
    nCarriers = length(carrierArray);
    
    Prx_mtx = zeros(nCells,nUes,nCarriers);
    
    for carrierNum=1:nCarriers
        Prx_mtx(:,:,carrierNum) = [carrierArray(carrierNum).rxPw];        
        carrierArray(carrierNum).RSRPmatrix = Prx_mtx(:,:,carrierNum);
    end
    
    switch inAttachUeMethod
        case 1      % attach user to the cell based on the maximum power recieved from the BSs
            
            Prx_mtx = sum(Prx_mtx,3);
            
            % save('lsaSimulatorVars','Prx_mtx','-append');
            
            [cell_val, cellsIndexNumber] = sort(Prx_mtx,1,'descend');

             for n = 1:nUes
                if cell_val(1,n) >= inPrxThreshold
                   cellId = cellsIndexNumber(1,n); 
                   ues(n).attachUser(cellId);
                   cls(cellId).attachUsers(n);
                   outUeCellAlloc(cellId,n) = 1;                   
                else 
                    for carrierNum = 1:nCarriers
                        carrierArray(carrierNum).RSRPmatrix(:,n) = 0;
                        carrierArray(carrierNum).mtx(:,n) = 0;
                        carrierArray(carrierNum).chPlusAntenna(:,n) = 0;                        
                    end                
                end
            end

        case 2  % attaching based on distance
            distanceMtx = varDistCellUeNoWA;

            [cell_val, cellsIndexNumber] = sort(distanceMtx,1);
            %----------------


            %%% THIS IS THE OUTPUT OF THE FUNCTION - BEGIN

            for n = 1:nUes
                %Check is center of BS is center of site, than 3 macro cells
                %have the same center. We need to decide which one of the cells is closer to user.
                cellId = cellsIndexNumber(1,n);

                if cls(cellId).typeIndex == macroInd

                    [macroCellsIndexNumber, temp_index] = sort(cellsIndexNumber(1:3,n));
                    temp_cell_val = cell_val(1:3,n);
                    macroCell_val = temp_cell_val(temp_index);

                    centerOfMacro=cls(cellId).complexPosition;
                    positionOfUser = ues(n).complexPosition;
                    angleOfDistanceVector = angle(positionOfUser-centerOfMacro)*180/pi;
                    if (angleOfDistanceVector<=-30) && (angleOfDistanceVector>-150)
                        % it is third one
                        varInd = 3;                        
                    else
                        if (angleOfDistanceVector<=90) && (angleOfDistanceVector>-30)
                            %it is first one
                            varInd = 1;
                        else
                            %it is second one
                            varInd = 2;
                        end
                    end

                    maxVal = macroCell_val(varInd);
                    cellId = macroCellsIndexNumber(varInd);                   

                else
                    %cell is not macro
                    maxVal = cell_val(1,n);
                    cellId = cellsIndexNumber(1,n);
                end

                if maxVal >= inPrxThreshold                    
                    ues(n).attachUser(cellId);
                    cls(cellId).attachUsers(n);
                    outUeCellAlloc(cellId,n) = 1;
                else
                    for carrierNum = 1:nCarriers
                        carrierArray(carrierNum).RSRPmatrix(:,n) = 0;
                        carrierArray(carrierNum).mtx(:,n) = 0;
                        carrierArray(carrierNum).chPlusAntenna(:,n) = 0;                         
                    end                
                end                
            end
            
        case 3  % based on network operator
            
            Prx_mtx = sum(Prx_mtx,3);
            
            % save('lsaSimulatorVars','Prx_mtx','-append');
            
            [cell_val, cellsIndexNumber] = sort(Prx_mtx,1,'descend');

            for n = 1:nUes
                nn_tmp = 1;
                candidate_cell_idx = cellsIndexNumber(nn_tmp,n);  
                while ues(n).networkOP_ID ~= cls(candidate_cell_idx).networkOP_ID
                    nn_tmp = nn_tmp+1;
                    if nn_tmp > length(cellsIndexNumber);
                       error('Index out of bound: fix it!!'); 
                    end
                    candidate_cell_idx = cellsIndexNumber(nn_tmp,n);  
                end
                if cell_val(candidate_cell_idx,n) >= inPrxThreshold 
                    cellId = candidate_cell_idx; 
                    ues(n).attachUser(cellId);
                    cls(cellId).attachUsers(n);
                    outUeCellAlloc(cellId,n) = 1;                   
                else 
                    for carrierNum = 1:nCarriers
                        carrierArray(carrierNum).RSRPmatrix(:,n) = 0;
                        carrierArray(carrierNum).mtx(:,n) = 0;
                        carrierArray(carrierNum).chPlusAntenna(:,n) = 0;                        
                    end                
                end
            end
            
        otherwise        
    end

    if inFFflag<=0
        for n = 1:nCells
           cls(n).isEmptyCell = isempty(cls(n).usersList);
        end
    end

    outInterfIndMtx = cellsIndexNumber(1:1+inInterfPerUe,:).';
end