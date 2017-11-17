%%% SON SIMULATOR: UserEquipment class
%
%   Author: Carlo Galiotto
%
%   Date:   July 2010
%
%   Description: this class implements the eNB class, defines the
%   properties of a eNB and the methods 


% <cell_ind_vector> is an array that stores for each UE, the index 
% corresponding to the serving cell of that UE.

function compute_SINR_matrix_v2(objCarrier,cell_ind_vector,varSchedArray,...
    varRaArray, currTTI, PHYmodelling, capComputation, SINR_statistic,nOfCel,nOfUsers)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          AWGN CHANNEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
% without fast-fading

    % noise power per PRB of the current carrier. It is computed according
    % to the following expression:
    
    % noisePSD * noiseFigure * PRBbandwidth
    
    % where:
    % noisePSD is the power spectral density (dB/Hz). It is an input
    % parameter (conf.PHY.noiseFloor).
    
    % noiseFigure is the user noise figure. It is an input parameter
    % (conf.user.noiseFigure)
    
    % PRBbw is the bandwidth per PRB, which depends on each carrier. It is
    % computed according to the following expression: (available bandwidth 
    % of current carrier/ total number of PRBs in the current carrier)
    
    varUE_Pn_PRB = objCarrier.user_Pn_PRB;
    
    % number of Physical Resource Blocks in the current carrier
    nOfPRB = objCarrier.numOfAvailablePRB;
   
    % channel coefficients between each pair of cell and UE in the current
    % carrier. If the cell does not use the current carrier, this
    % coefficient is 0.
    
    % <CELLS X UEs>
    ch_mtx = objCarrier.chPlusAntenna;
    
    
    ch_int_mtx = ch_mtx;
    ch_Prx_mtx = zeros(nOfCel,nOfUsers);

    % carrier_Cls::allocationMatrix(ii,j) stores 1, if the PRB j has been
    % assigned the UE ii. Otherwise it stores 0. The UE scheduling has 
    % been performed over the current carrier.
    
    % size: <number of UEs x number of PRBs>
    allMtx = objCarrier.allocationMatrix;
    
    
    % P_tx_mtx(ii,j) stores the power corresponding to PRB j in the cell ii, 
    % if the PRB j has been assigned to some UE attached to cell ii. Otherwise, 
    % it stores 0.
    
      
    % P_tx_mtx(ii,j) stores 1 if PRB j has been assigned to some UE that 
    % is attached to cell ii. Otherwise, it stores 0.
    
    % P_tx_mtx(ii,j) stores the power corresponding to PRB j in the cell ii, 
    % if the PRB j has been assigned to some UE attached to cell ii. Otherwise, 
    % it stores 0.
        
    P_tx_mtx =  zeros(nOfCel, nOfPRB);
    
    % P_tx_mtx_fsu(ii,j) stores the power corresponding to PRB j in cell ii.
    P_tx_mtx_fsu =  zeros(nOfCel, nOfPRB);
    
    for ii=1:nOfCel
        raCell = varRaArray(ii);
        av_sp = raCell.spectrumAvailable;
        varLen = length(av_sp);
        varUes = varSchedArray(ii).userSpectVect;
        for j=1:varLen
            prbID = av_sp(j);
            if varUes(j)>0 && prbID>0
                P_tx_mtx(ii,prbID) = raCell.pwPerPRBvector(prbID);
            end
        end
        P_tx_mtx_fsu(ii,:) = [raCell.pwPerPRBAllowedVector];
    end
    
    % array of UE index, considering all UEs in the network (1, 2, ..
    % nOfUsers)
    temp_vector=1:nOfUsers;
    
    % UE index that are not attached to any cell are removed from 
    % <temp_vector>
    temp_vector(cell_ind_vector==0)=[];
    
    % elements 0 are removed from <cell_ind_vector>
    cell_ind_vector(cell_ind_vector==0)=[];
    
    % Each element ch_int_mtx(ii, j) where ii is the serving cell of UE j,
    % is assigned 0. Thus, these pairs do not contribute to the interference
    % signal 
    
    lin_indices = sub2ind(size(ch_int_mtx),cell_ind_vector , temp_vector);
    ch_int_mtx(lin_indices) = 0;
    
    % channel coefficient between each UE and its serving cell are set to
    % the channel coefficient corresponding to the pair. The other pairs
    % are assigned 0.
    ch_Prx_mtx(lin_indices) = ch_mtx(lin_indices);
    
    P_u_mask = allMtx >0;
    
    % n_UEs x nPRBs matrix
    
    % interference signal per UE and PRB.
    int_mtx =  (ch_int_mtx.') * P_tx_mtx ;    
    
    % ch_Prx_mtx.' is a 2D array (M, N), where M is the number of UEs, and
    % N is the number of cells in the network. Each element ch_Prx_mtx(ii,j)
    % stores the channel coefficient between the UE ii and cell j, if the
    % cell j is the serving cell of UE ii, else it stores 0.

    
    % P_tx_mtx(ii, j) stores the power per PRB j in the cell ii, if the PRB j
    % has been assigned to some UE that is attached to cell ii, else it
    % stores 0.
    
    % P_u_mtx is a 2D array (M, N), where M is the total number of UEs in 
    % the network, and N is the total number of PRBs in the current carrier.
    
    % Each element of P_u_mtx(ii,j) stores 0, if the PRB j has not been
    % assigned to UE ii. Otherwise, it stores the product <channel
    % coefficient between UE ii and its serving cell> * <power
    % corresponding to PRB j in the serving cell of UE ii>
    
    P_u_mtx = P_u_mask .*( (ch_Prx_mtx .')* P_tx_mtx );     % n_UEs x nPRBs matrix
    
    P_u_mtx_estim = ( (ch_Prx_mtx .') * P_tx_mtx );    
    
    P_u_mtx_fsu_estim = ((ch_Prx_mtx .') * P_tx_mtx_fsu);     % n_UEs x nPRBs matrix

    SINR_mtx = P_u_mtx./(int_mtx + varUE_Pn_PRB);   % n_UEs x nPRBs matrix
    
    
    %%% SINR estimate for FSU and scheduling
    
    switch PHYmodelling.CQIreportType
        
        case 1
            
            %%% SINR_mtx_estim indicates estimated SINR only for PRBs in use by the cell
            
            SINR_mtx_estim = P_u_mtx_estim ./ (int_mtx + varUE_Pn_PRB);   % n_UEs x nPRBs matrix
            
            %%% SINR_mtx_estim_fsu indicates estimated SINR for all PRBs, also for the ones not in use by the cell
            SINR_mtx_estim_fsu = P_u_mtx_fsu_estim./(int_mtx + varUE_Pn_PRB);       % n_UEs x nPRBs matrix
            if capComputation == 0
                C_mtx_estim_fsu = shannonMk0(SINR_mtx_estim_fsu);               % n_UEs x nPRBs matrix
            elseif capComputation == 1
                C_mtx_estim_fsu = shannonMkI(SINR_mtx_estim_fsu);               % n_UEs x nPRBs matrix
            end
            
            % the following MTXs have dim: n_UEs x nPRBs
            objCarrier.C_mtx_estim_fsu_current = C_mtx_estim_fsu;                % estimation Capacity matrix (all over the spectrum) - Current value (for delivered data evaluation)      
            % objCarrier.SINR_mtx_estim_fsu_current = SINR_mtx_estim_fsu;          % estimation Capacity matrix  
            objCarrier.C_mtx_estim_fsu = C_mtx_estim_fsu;                         
            objCarrier.SINR_mtx_estim_fsu = SINR_mtx_estim_fsu;   
            objCarrier.SINR_mtx_estim = SINR_mtx_estim;
            
        case 2
            
            error('No implemented yet');
            
        case 3
           
            SINR_mtx_estim_fsu = P_u_mtx_fsu_estim./(int_mtx + varUE_Pn_PRB);       % n_UEs x nPRBs matrix
            if capComputation == 0
                C_mtx_estim_fsu = shannonMk0(SINR_mtx_estim_fsu);               % n_UEs x nPRBs matrix
            elseif capComputation == 1
                C_mtx_estim_fsu = shannonMkI(SINR_mtx_estim_fsu);               % n_UEs x nPRBs matrix
            end 
            objCarrier.C_mtx_estim_fsu_current = C_mtx_estim_fsu;
            % objCarrier.SINR_mtx_estim_fsu_current = SINR_mtx_estim_fsu;
           
            
            if currTTI < PHYmodelling.CQIreportDelay
                objCarrier.C_mtx_estim_fsu = 0.1*ones(nOfUsers,nOfPRB);
                objCarrier.SINR_mtx_estim = zeros(nOfUsers,nOfPRB);                   % n_UEs x nPRBs matrix
                objCarrier.SINR_mtx_estim_fsu = zeros(nOfUsers,nOfPRB);           % n_UEs x nPRBs matrix 
            elseif rem(currTTI - 1, PHYmodelling.CQIreportDelay) == 0
                objCarrier.C_mtx_estim_fsu = objCarrier.lastEstimated_C_mtx;
                objCarrier.SINR_mtx_estim_fsu = objCarrier.lastEstimated_SINR_mtx;
            end
            
            if rem(currTTI - 1, PHYmodelling.CQIreportSamplingTime) == 0
                %%% SINR_mtx_estim indicates estimated SINR only for PRBs in use by the cell
                %%% SINR_mtx_estim_fsu indicates estimated SINR for all
                %%% PRBs, also for the ones not in use by the cell
                objCarrier.lastEstimated_C_mtx = C_mtx_estim_fsu;
                objCarrier.lastEstimated_SINR_mtx = SINR_mtx_estim_fsu;                
            end
            
        otherwise
            error('No other CQI report types available');
            
    end

    SINR_mtx(P_u_mtx<=0) = 0;                               % n_UEs x nPRBs matrix
    
    objCarrier.SINR_mtx = SINR_mtx;                               % n_UEs x nPRBs matrix
    objCarrier.int_mtx = int_mtx;
    
    if SINR_statistic
        % number of PRBs assigned to each UE
        avSINRperUE_den = sum(objCarrier.allocationMatrix.');
        
        % each element of this array (ii, j) stores the sum of
        % (ch_coeff(ii, cell k)*pwPRB(cell k, PRB j), from k = 1 to the
        % total number of cells in the network. Consider that ch_coeff(ii,
        % cell k) = 0, if cell k is the serving cell of UE ii.
        
        % mean of int_mtx per UE
        objCarrier.int_mtx = mean(int_mtx,2);
        
        % for each UE is obtained the sum of powers corresponding to those
        % PRBs that have been assigned to the UE. 
        objCarrier.avSINRperUE = sum((P_u_mtx.*objCarrier.allocationMatrix).')./(sum( ((int_mtx + varUE_Pn_PRB) .* objCarrier.allocationMatrix).' ));
        objCarrier.avSINRperUE(avSINRperUE_den<=0) = 0; 
        ueSpectralEfficiency = (sum(shannonMk0(SINR_mtx),2).')./avSINRperUE_den ;
        objCarrier.ueCumSpectralEfficiency(avSINRperUE_den>0) = objCarrier.ueCumSpectralEfficiency(avSINRperUE_den>0) + ueSpectralEfficiency(avSINRperUE_den>0);
        objCarrier.ueSpectralEfficiency_cnt(avSINRperUE_den>0) = objCarrier.ueSpectralEfficiency_cnt(avSINRperUE_den>0) + 1;
    end
    
    %%% CAPACITY
    
    if capComputation == 0
        C = objCarrier.PRBbw .* shannonMk0(SINR_mtx);               % n_UEs x nPRBs matrix
    elseif capComputation == 1
        C = objCarrier.PRBbw .* shannonMkI(SINR_mtx);               % n_UEs x nPRBs matrix
    end
    
    C(SINR_mtx<=0) = 0;
    objCarrier.thPerUE = sum(C.');                                % 1 x UEs vector
    objCarrier.C_mtx = C;
    
    %%% RESOURCES USAGE
    
    % for each cell  and each carrier, the sum of UEs that have been
    % assigned with some PRB is stored in the property
    % cell_Cls::resourceAllocation_Cls(carrierNum)::currentResourceUsage
    
    % carrier_Cls::cumResourceUsage stores the sum of UEs that have been 
    % assigned per cell in the current carrier, added to the sum computed in the previous iterations.
    objCarrier.cumResourceUsage = objCarrier.cumResourceUsage + ...
        sum([varRaArray.currentResourceUsage]);
    
    % currentResourceAvailability stores the number of PRBs that have been
    % assigned to each cell in a given carrier in the current iteration
    % these values of each cell are summed and are added to the cumulative
    % sum of previous iterations. 
    objCarrier.cumResourceAvailability = objCarrier.cumResourceAvailability + ...
        sum([varRaArray.currentResourceAvailability]);   
    % carrier_Cls::allocationMatrix is 2D array with size 
    % <Number of UEs in all the network x Number of PRBs in the current 
    % carrier> It stores 1 if the UE ii has been assigned PRB j. Otherwise,
    % it stores 0.
    
    % carrier_Cls::bandPerUE stores the number of PRBs that have been
    % assigned to each UE in the current iteration
    objCarrier.bandPerUE = sum(objCarrier.allocationMatrix,2).';
    %%%  
end