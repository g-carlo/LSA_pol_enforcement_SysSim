N = 9;

M = 10;

d = 250;

x = -250:250:250;
y = -250:250:250;

[x_pos,y_pos] = meshgrid( x,fliplr(y));

% H_pathloss = pathloss_grid_1( 3,4, 250);
% 
% H_pathloss = 

z_pos = x_pos+1i*y_pos;
z_pos = z_pos(:);
z_pos(10) = z_pos(9);


sens_pos_1 = 1i*500;
sens_pos_2 = -500;
sens_pos_3 = -1i*500;
sens_pos_4 = 500;

dist_sens_1_to_BSs = abs(sens_pos_1-z_pos);	
dist_sens_2_to_BSs = abs(sens_pos_2-z_pos);
dist_sens_3_to_BSs = abs(sens_pos_3-z_pos);
dist_sens_4_to_BSs = abs(sens_pos_4-z_pos);

for n = 1:M
   PL_1(n) = -ExtendedHata(dist_sens_1_to_BSs(n)/1000, 1,1, 2300, 32, 1.5, 6); 
   PL_2(n) = -ExtendedHata(dist_sens_2_to_BSs(n)/1000, 1,1, 2300, 32, 1.5, 6); 
   PL_3(n) = -ExtendedHata(dist_sens_3_to_BSs(n)/1000, 1,1, 2300, 32, 1.5, 6); 
   PL_4(n) = -ExtendedHata(dist_sens_4_to_BSs(n)/1000, 1,1, 2300, 32, 1.5, 6); 
end

PL_1 = PL_1 + 8*randn(1,M);
PL_2 = PL_2 + 8*randn(1,M);
PL_3 = PL_3 + 8*randn(1,M);
PL_4 = PL_4 + 8*randn(1,M);

PL_1 = 10.^( PL_1./10);
PL_2 = 10.^( PL_2./10);
PL_3 = 10.^( PL_3./10);
PL_4 = 10.^( PL_4./10);



PL_1 = PL_1 .* ( raylrnd(1,1,M) + 1i*raylrnd(1,1,M));
PL_2 = PL_2 .* ( raylrnd(1,1,M) + 1i*raylrnd(1,1,M));
PL_3 = PL_3 .* ( raylrnd(1,1,M) + 1i*raylrnd(1,1,M));
PL_4 = PL_4 .* ( raylrnd(1,1,M) + 1i*raylrnd(1,1,M));


G = [PL_1; PL_2; PL_3; PL_4];
BS_to_BS_channel = zeros(1,N);

