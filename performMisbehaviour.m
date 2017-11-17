function [cellsPerCarrier, isMisbehaving] = performMisbehaviour(conf,netOperatorEventGen,cellsPerCarrier,op_idx)  
    
    n_misbehaving_cells = 3;
    isMisbehaving = rand(1) <=conf.operator.misBehProbability(op_idx);
    if isMisbehaving
        numCarriers    = length(conf.carrier.frequency);
        numLSAcarriers = conf.database.numOfLSAChannels;
        index_LSA_channels = (numCarriers - numLSAcarriers + 1):numCarriers;
        list_of_cells = netOperatorEventGen.cell_idx_per_NO(op_idx).list;
    
        bool_TX_cells = cellsPerCarrier( index_LSA_channels ,list_of_cells);
        idx_non_TX_cells = find(~bool_TX_cells);
        if ~isempty(idx_non_TX_cells)
            misbehavingCells_idx_toRemap = randSequenceFromBin(n_misbehaving_cells ,length(idx_non_TX_cells));    
            misbehavingCells_idx = idx_non_TX_cells( misbehavingCells_idx_toRemap);
            bool_cells = false(size( list_of_cells) );
            bool_cells(misbehavingCells_idx) = true;
            cellsPerCarrier( index_LSA_channels ,list_of_cells )  = bool_cells | cellsPerCarrier( index_LSA_channels ,list_of_cells );    
        end
                
    end
    

end