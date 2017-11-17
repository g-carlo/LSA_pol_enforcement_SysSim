function [wa_cell_UE_mtx,wa_pos_mtx,phi_angle,theta_angle,distanceNoWA] = ...
    computeDistanceToUes(objType,cel_posVect_noWA,ue_posVect,varHeightRx,...
    useWA)
   
    nUe = length(ue_posVect);
    nCel = objType.numOfCells;
    
    varHeightTx = objType.height;
    
    if useWA==0
    
        ue_posVect = repmat(ue_posVect,nCel,1);
        cel_posVect_noWA = repmat(cel_posVect_noWA.',1,nUe);
    
        wa_cell_UE_mtx = abs( cel_posVect_noWA - ue_posVect)./1000;
        wa_pos_mtx = cel_posVect_noWA;
       
        phi_angle = angle( ue_posVect - cel_posVect_noWA);
        distanceNoWA = wa_cell_UE_mtx .* 1000;
        theta_angle = atan( (varHeightTx - varHeightRx) ./ distanceNoWA);   
    else
    
        shiftVector = objType.shiftArray;


        nShift = length(shiftVector);

        % UE pos array repeated for each cell and each shiftArray element
        pos_ue_array = repmat(ue_posVect,[nCel 1 nShift]);

        pos_cel_array_noWA = repmat(cel_posVect_noWA.',[1 nUe nShift]);
        shift_array = repmat(shiftVector.',[1 nCel nUe]);
        shift_array = permute(shift_array,[2 3 1]);
        pos_cel_array_WA = pos_cel_array_noWA+shift_array;
        dist_array_WA = abs(pos_ue_array-pos_cel_array_WA);
        [minVal_mtx, indexMinVal_mtx] = min(dist_array_WA,[],3);
        wa_cell_UE_mtx = minVal_mtx./1000;
        index_dim1 = reshape(repmat((1:nCel).',1,nUe),1,nCel*nUe);
        index_dim2 = reshape(repmat((1:nUe),nCel,1),1,nCel*nUe);
        index_dim3 = reshape(indexMinVal_mtx,1,nCel*nUe);
        wa_pos_mtx_notShaped = pos_cel_array_WA(sub2ind(size(pos_cel_array_WA),...
            index_dim1,index_dim2,index_dim3));
        wa_pos_mtx = reshape(wa_pos_mtx_notShaped,nCel,nUe);

        phi_angle = angle( repmat(ue_posVect,nCel,1) - wa_pos_mtx);
        theta_angle = atan( (varHeightTx - varHeightRx) ./ minVal_mtx);   

        ue_posVect = repmat(ue_posVect,nCel,1);
        cel_posVect_noWA = repmat(cel_posVect_noWA.',1,nUe);
        distanceNoWA = abs( cel_posVect_noWA - ue_posVect);
    end
end