function  H_pathloss = pathloss_grid_1( Network_len,alpha, unit)
%pathloss between nodes on a 2D geographical area
%   Detailed explanation goes here
Network_grid = Network_len^2;

% Compute the distances between cells fo the geographical area
D1 = combvec(1:Network_len, 1:Network_len);
Dif = D1(2,:)- D1(1,:);
D_mid1 = combvec(Dif.^2,Dif.^2);
D_mid = flipud(D_mid1);
Dist = sqrt(sum(D_mid,1));
D = reshape(Dist,Network_len,Network_grid,[]);
D2 = permute(D,[1 3 2]);
D3 = reshape(D2,Network_grid,Network_grid);
D4 = unit*D3;
D4(D4==0)=0.8;
%alpha = 1.8 + 0.4*rand(size(D4));
H_pathloss = D4.^(-alpha);

end

