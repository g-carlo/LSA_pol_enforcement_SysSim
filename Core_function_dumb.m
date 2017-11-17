function [Est_Loc_active_BS, Pd_table, Pfa_table]  = Core_function_dumb(isMisbehaving)


if isMisbehaving 
    
    Est_Loc_active_BS = [1 2]; 
    Pd_table = 1;
    Pfa_table = 0; 
        
else
    
    Est_Loc_active_BS = []; 
    Pd_table = 1;
    Pfa_table = 0;
    
end