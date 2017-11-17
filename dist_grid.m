function  dist_mat = dist_grid(Sens_Loc, BS_Loc)
% Computes the distances between the sensors and the BS
N = length(Sens_Loc);
M = length(BS_Loc); 

dist_mat = zeros(N,M);

for i=1:N
    
    dist_mat(i,:) = abs(Sens_Loc(i) - BS_Loc);
    
end

end

