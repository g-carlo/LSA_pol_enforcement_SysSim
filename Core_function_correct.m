% MISBEHAVIOR DETECTION FOR TCD LSA SIMULATOR Cofe Function 
% 5/10/2016 

function [Est_Loc_active_BS, Pd_table, Pfa_table, Est_energy_active_BS, Comp_time] = Core_function_correct( N_white_BS, NSens, BS_positions, IU_Loc, Loc_active_BS,file_label)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_S = 10;

% Define a topology for the LSA network
Cells = 57; % number of cells (3*19 sites = 57 cells)

BS_P = BS_positions;
BS_Loc = transpose(BS_P(:,1:Cells));
N_IUs = length(IU_Loc);

% parameters to set
height_BS = 32;
height_Sens = 1;
antenna_gain = 10.^(14/10);
frequency = 2300;
noise_thresh = 3e-3; % 3e-3 5dBm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Power level 
LU_power = 0.72;  % 43 dBm over 5MHz so (180/5000)*20 W = 0.72 W per channel
IU_power = 9e-4;  % 20 dBm over 20MHz so (180/20000)*20 W = 9e-4 W per channel

% Define the location of the IU
% Unit = 500;
% Sen_range = 3;
%IU_Loc = Sen_range*Unit + Unit*rand*sign(randn) + 1i*(Sen_range*Unit + Unit*rand*sign(randn));

% Define the Sensing Network and its topology

% Dis_dif_max = zeros(N_IUs,1);
% for n=1:N_IUs
% Dis_dif = abs(IU_Loc(n) - BS_Loc);
% Dis_dif_max(n) = max(Dis_dif);
% end

Dis_dif = abs(IU_Loc(1) - BS_Loc);
[~,m_ind] = max(Dis_dif);
Center_of_the_SN = 0+1i*0; % (IU_Loc + BS_Loc(m_ind))/2;
R_of_the_SN = abs(Center_of_the_SN-BS_Loc(m_ind));

% [~,m_ind] = max(Dis_dif_max);
% Center_of_the_SN = 0+1i*0; %(IU_Loc(m_ind) + BS_Loc(m_ind))/2;
% R_of_the_SN = abs(Center_of_the_SN-IU_Loc(m_ind));
% Generate random events for the computation of the probabilities
Sens_Loc = Center_of_the_SN + rand(NSens,1)*R_of_the_SN.*sign(randn(NSens,1)) + 1i* rand(NSens,1)*R_of_the_SN.*sign(randn(NSens,1));

% Define the Sensing Network and its topology
dist_mat_BS = dist_grid(Sens_Loc, BS_Loc);
dist_mat_IU = dist_grid(Sens_Loc, IU_Loc');

% Uncomment to display the topology of the LSA Network
% figure(1);
% plot(BS_Loc,'.','MarkerSize',8)
% hold on
% plot(Sens_Loc,'o','MarkerSize',8)
% plot(IU_Loc,'*')
% %plot(Center_of_the_SN,'v','MarkerSize',8)
% legend('LU Net.','Sensor Net.', 'IU');
% title('LSA Network topology');
% hold off
% % 
% %Compute the distances between sensors and BSs
% dist_mat = dist_grid(Sens_Loc, BS_Loc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path loss according to the Hata model
Path_loss_mat = [-ExtendedHata(dist_mat_BS/1000, frequency, height_Sens, height_BS, 6) -ExtendedHata(dist_mat_IU/1000, frequency, height_Sens, height_BS, 6)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add Shadow fading 
PL_SH = Path_loss_mat + 8*randn(NSens, Cells + N_IUs);
PL_SH_lin = 10.^( PL_SH./10);

% Add the Antenna gains
load('BS_antenna_orientation_file.mat');

% matrix of antenna orientations in a matrix format, N_Sens x N_BS matrix
boresightMtx = repmat(BS_antenna_orientation,NSens,1);

% position of base stations in a matrix format, N_Sens x N_BS matrix
BS_Loc_matrix = repmat( BS_Loc.',NSens,1 ) ;

% position of sensors in a matrix format, N_Sens x N_BS matrix
Sens_Loc_matrix = repmat( Sens_Loc,1,Cells ) ;

% horizontal angle between sensors and BSs, N_Sens x N_BS matrix
phi_angle = angle( Sens_Loc_matrix - BS_Loc_matrix );

% vertical angle between sensors and BSs, N_Sens x N_BS matrix
theta_angle = atan( (height_BS-height_Sens) ./ dist_mat_BS);    

% antenna gains from BSs to Sensors, N_Sens x N_BS matrix 
antennaGainMatrix1 = antenna_gain.*directiveAntExternalFnc(phi_angle,theta_angle,boresightMtx);
antennaGainMatrix = [antennaGainMatrix1 ones(NSens,N_IUs)];

% Compute the channels 
G1 = sqrt(antennaGainMatrix.*PL_SH_lin);

% Add Rayleigh fading and compute the channels
% Variance of the normal Gaussian distribution for the simulation of the Rayleigh channel
G = G1.*( randn(NSens,Cells + N_IUs) + 1i*randn(NSens,Cells + N_IUs) )/sqrt(2);

% Dimension of the noise vector
m = 1;

% Number of actiive BSs
LU_Transmitters = length(Loc_active_BS);
IU_Transmitters = length(IU_Loc);

% Noise st.deviation
sig_n = sqrt(5.68*1e-15); % (180KHz / 5000 KHZ)*20 W, SNR = 32;

% LU's that each transmits with equal power over the detection period
% L = sort(randsample(setdiff(1:cm,S),Transmitters));
% L = sort(randsample(1:Cells,LU_Transmitters));
L = [Loc_active_BS Cells+1:1:Cells+N_IUs];
E_star = zeros(Cells + N_IUs,1);
  % for (180/5000)*20 W
E_star(L(1:LU_Transmitters)) = LU_power;
E_star(L(LU_Transmitters+1:end)) = IU_power;
En_irwl = zeros(Cells+N_IUs,1);

% Initialization 
En_sum = 0;
T = 0;

y_clean_tot = zeros(NSens,T_S);
y_tot       = zeros(NSens,T_S);

for t=1:T_S

     % AWGN for the linear equation
     Gamma_n = sig_n^2;
     XY_n = mvnrnd(zeros(2*m,1),[Gamma_n zeros(m,m); zeros(m,m) Gamma_n]/2, NSens);
     X_n = XY_n(:,1);
     Y_n = XY_n(:,2);
     v = complex(X_n,Y_n);  
    
%     % Transmitted signal
%     x_star = zeros(Cells,1);
%     n = length(L); %number of symbols
%     x = 2*(randi([0 1],1,n)-0.5) + 1i*2*(randi([0 1],1,n)-0.5); %Random QPSK symbols
%     x_star(L) = sqrt(LU_power*1/2).*x;
    
     % Transmitted signal
    x_star = zeros(Cells+N_IUs,1);
    n = length(L); %number of symbols
    x1 = 2*(randi([0 1],1,n-N_IUs)-0.5) + 1i*2*(randi([0 1],1,n-N_IUs)-0.5); %Random QPSK symbols
    x2 = 2*(randi([0 1],1,N_IUs)-0.5) + 1i*2*(randi([0 1],1,N_IUs)-0.5); %Random QPSK symbols
    x_star(L(1:LU_Transmitters)) = sqrt(LU_power*1/2).*x1;
    x_star(L(LU_Transmitters+1:end)) = sqrt(IU_power*1/2).*x2;

    % Compute the clean signal
    y_clean = G*x_star;
    
    % Measure the SNR
    SNR=10*log10(sum(abs(y_clean).^2)/sum(abs(v).^2));
    y = y_clean + v;
    LL = x_star~=0; 
    
    y_clean_tot(:,t) = y_clean;
    y_tot(:,t) = y;
      
%    if N_white_BS == 0
%  
%        Pd_table = 0.7;
%        Pfa_table = 0.5;
%     
%     % 2 
%     % Test STELA for LASSO
%     mu = 3e-13;
%     tic;
%     [x_stela, ~, ~] = stela(G, y, mu);
%     time = toc;
%     %X2 = [x_star x_stela];
%     
%     en_est = abs(x_stela).^2;
%     En_sum = En_sum + en_est;
%     
%     T = T + time;
%         
%    elseif N_white_BS == 15
       
    Pd_table = 0.95;
    Pfa_table = 0.05;
    
    % 5 
    % IRW-LASSO 
    mu_stela = 7.51e-14;
    mu5 = 4.51e-14;
    tic;
    [x_st, ~, ~] = stela(G, y, mu_stela);
    [x_irwl, ~, ~, ~, ~] = irwl_w_stela(  mu5, G, y, x_st,1e-3,1 );
    time=toc;
    
    en_est = abs(x_irwl).^2;
    En_sum = En_sum + en_est;    
    T = T + time;
    
%    else
%        
%        Pd_table = 0.8;
%        Pfa_table = 0.3;
%     
%     % 5 
%     % IRW-LASSO 
%     mu5 = 3e-13;
%     tic;
%     [x_st, ~, ~] = stela(G, y, mu);
%     [x_irwl, ~, ~, ~, ~] = irwl_w_stela(  mu5, G, y, x_st,1e-3,1 );
%     time=toc;
%     
%     en_est = abs(x_irwl).^2;
%     En_sum = En_sum + en_est;    
%     T = T + time;
%     
%    end
   
end
 
    Est_energy_active_BS = En_sum/T_S; 
    Est_energy_active_BS(Cells+1:end) = 0;
    Est_Loc_active_BS = find(Est_energy_active_BS>noise_thresh);
    Comp_time = T;
    
    if length(Est_Loc_active_BS) > 10
        
        % matrix G the output y (from the linear equation) and also the locations of IUs and the Sensors. And the locations of the misbehaving BSs of the LU.
        fileName = ['item_n_'  num2str(file_label) '_n_active_BS_' num2str(length(Est_Loc_active_BS))];
        save( fileName, 'NSens' , 'BS_positions', 'IU_Loc', 'Loc_active_BS' , 'G', ...
            'y_clean_tot', 'y_tot','Sens_Loc' );
                
    end
    
end


% % function [ x, t, S] = stela(A, b, mu, Max_N)
% % 
% % if nargin<4,
% %     Max_N=20000;
% % end;
% % 
% % AA = A'*A;
% % M = size(A,2);
% % x = zeros(M,1);
% % c=0;
% % precision = 1e-6;
% % 
% % % Set the fast soft thresholding function
% % fast_sthresh = @(x,th) sign(x).*max(abs(x) - th,0);
% % 
% % % gamma opereator
% % gamma_op = @(x) min(max(x,0),1);
% % 
% % 
% % 
% % for t=1:Max_N
% %     
% %     x_old=x;
% %  
% %     r = diag(AA).*x - AA*x + A'*b;
% %     Bxt = fast_sthresh(r+c*x,mu)./(diag(AA)+c);
% %  
% %     g_p = (A*(Bxt-x))'*(A*(Bxt-x));
% %     gamma = gamma_op(- ( (A*x-b)'*A*(Bxt-x) + mu*(norm(Bxt,1)-norm(x,1)) )/g_p);
% % 
% %     x = x + gamma*(Bxt-x);
% % 
% %     S = sort(find(x));
% %     
% %     err = norm(x_old-x);
% %     if err< precision
% %         break;
% %     end
% %     
% % end
% % 
% % 
% % end
% % 
% % 
% % function [ x_hat, i, S2, g_vec,error ] = irwl_w_stela(  mu, A, y, x_stela,  eps, gamma)
% % %UNTITLED2 Summary of this function goes here
% % %   Detailed explanation goes here
% % 
% % PRECISION =1e-12;
% % Max_N2 = 20;
% % 
% % for i=1:Max_N2
% %    
% %     x_old = x_stela;
% %  
% %     [x_hat, ~, S2,g_vec] = stela3(mu, A, y, x_stela, eps, gamma);
% %     
% %     error = norm(x_old - x_hat);
% %     x_stela = x_hat;
% %     if error<PRECISION
% %         return;
% %     end
% %     
% % end
% % 
% % 
% % end
