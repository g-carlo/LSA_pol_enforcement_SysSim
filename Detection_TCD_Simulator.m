% MISBEHAVIOR DETECTION FOR TCD LSA SIMULATOR
% 5/10/2016 

clear all

rseed=1;
rand('seed', rseed);
randn('seed', rseed);

% Parallel processing option for matlab 
%  matlabpool('open',4);

% Number of independent realizations 
It =100;
% Sampling period
T = 10; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define a topology for the LSA network
Cells = 57; % number of cells (3*19 sites = 57 cells)

load('BS_position_file.mat');
BS_P = BS_positions;
BS_Loc = transpose(BS_P(:,1:Cells));

% parameters to set
height_BS = 32;
height_Sens = 1;
antenna_gain = 10.^(14/10);
frequency = 2300;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the location of the IU
Unit = 500;
Sen_range = 3;
IU_Loc = Sen_range*Unit + Unit*rand*sign(randn) + 1i*(Sen_range*Unit + Unit*rand*sign(randn));

% Define the Sensing Network and its topology
NSens = 40;
Dis_dif = abs(IU_Loc - BS_Loc);
[mv,m_ind] = max(Dis_dif);
Center_of_the_SN = 0+1i*0; % (IU_Loc + BS_Loc(m_ind))/2;
R_of_the_SN = abs(Center_of_the_SN-BS_Loc(m_ind));

load('BS_antenna_orientation_file.mat');

% matrix of antenna orientations in a matrix format, N_Sens x N_BS matrix
boresightMtx = repmat(BS_antenna_orientation,NSens,1);

% Dimension of the noise vector
m = 1;

% Licensee transmitters
% Percentage of sparsity
per_x = 0.1;
% Number of actiive BSs
LU_Transimtters = round(per_x*Cells);

% Noise st.deviation
sig_n = sqrt(5.68*1e-15); % (180KHz / 5000 KHZ)*20 W, SNR = 32;
noise_thresh = 1e-5; % -41dBm

val = 0.72;  % for (180/5000)*20 W


% Initialization 
MT1 = zeros(1,It);
MT2 = zeros(1,It);
MT5 = zeros(1,It);

R_ols = zeros(1,It);
R_stela = zeros(1,It);
R_irwl = zeros(1,It);

Av_r_err_ols = zeros(1,It);
Av_r_err_stela = zeros(1,It);
Av_r_err_irwl = zeros(1,It);


Total_sup_cor_ols = 0;
Total_sup_err_ols = 0;
Pdc_ols = 0;
Pfa_ols = 0;

Total_sup_cor_stela = 0;
Total_sup_err_stela = 0;
Pdc_stela = 0;
Pfa_stela = 0;

Total_sup_cor_irwl = 0;
Total_sup_err_irwl = 0;
Pdc_irwl = 0;
Pfa_irwl = 0;

% hb = waitbar(0,'Please wait...0 %');

for i=1:It

% Initialization 

T1 =0; 
T2 =0; 
T5 =0; 

Rat_ols = 0;
Rat_stela = 0;
Rat_irwl = 0;

En_stela = zeros(Cells+1,1);
En_ols = zeros(Cells+1,1);
En_irwl = zeros(Cells+1,1);

% Generate random events for the computation of the probabilities
Sens_Loc = Center_of_the_SN + rand(NSens,1)*R_of_the_SN.*sign(randn(NSens,1)) + 1i* rand(NSens,1)*R_of_the_SN.*sign(randn(NSens,1));

% % Uncomment to see the topology of the LSA Network
% figure(1);
% plot(BS_Loc,'.','MarkerSize',8)
% hold on
% plot(Sens_Loc,'o','MarkerSize',8)
% plot(IU_Loc,'*')
% %plot(Center_of_the_SN,'v','MarkerSize',8)
% legend('LU Net.','Sensor Net.', 'IU');
% title('LSA Network topology');
% hold off

%Compute the distances between sensors and BSs
dist_mat_BS = dist_grid(Sens_Loc, BS_Loc);
dist_mat_IU = dist_grid(Sens_Loc, IU_Loc);

% Path loss according to the Hata model
Path_loss_mat = [-ExtendedHata(dist_mat_BS/1000, frequency, height_Sens, height_BS, 6) -ExtendedHata(dist_mat_IU/1000, frequency, height_Sens, height_BS, 6)];

% Add Shadow fading 
PL_SH = Path_loss_mat + 8*randn(NSens,Cells +1);
PL_SH_lin = 10.^( PL_SH./10);

% Add the Antenna gains
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
antennaGainMatrix = [antennaGainMatrix1 ones(NSens,1)];

% Compute the channels 
G1 = sqrt(antennaGainMatrix.*PL_SH_lin);

% Add Rayleigh fading and compute the channels
% Variance of the normal Gaussian distribution for the simulation of the Rayleigh channel
G = G1.*( randn(NSens,Cells + 1) + 1i*randn(NSens,Cells + 1) )/sqrt(2);

% LU's that each transmits with equal power over the detection period
% L = sort(randsample(setdiff(1:cm,S),Transimtters));
L = [sort(randsample(1:Cells,LU_Transimtters)) Cells+1];
E_star = zeros(Cells+1,1);
E_star(L) = val;

for t=1:T

     % AWGN for the linear equation
     Gamma_n = sig_n^2;
     XY_n = mvnrnd(zeros(2*m,1),[Gamma_n zeros(m,m); zeros(m,m) Gamma_n]/2, NSens);
     X_n = XY_n(:,1);
     Y_n = XY_n(:,2);
     v = complex(X_n,Y_n);  
    
    % Transmitted signal
    x_star = zeros(Cells+1,1);
    n = length(L); %number of symbols
    x = 2*(randi([0 1],1,n)-0.5) + 1i*2*(randi([0 1],1,n)-0.5); %Random QPSK symbols
    x_star(L) = sqrt(val*1/2).*x;
    
    % Compute the clean signal
    y_clean = G*x_star;
    
    % Measure the SNR
    SNR=10*log10(sum(abs(y_clean).^2)/sum(abs(v).^2));
    y = y_clean + v;
    LL = x_star~=0; 
      
    % 1
    % Least Squares Solution
    tic;
    x_ols = pinv(G)*y;
    t1=toc;
    
    % Perform Energy detection
    en_ls = abs(x_ols).^2;
    En_ols = En_ols + en_ls;
    
    ratio_ols = norm(x_star)^2/norm(x_star - x_ols)^2;
    Rat_ols = Rat_ols + ratio_ols;
    
    % 2 
    % Test STELA for LASSO
    mu = 5e-13;
    tic;
    [x_stela, n_it, S] = stela(G, y, mu);
    t2 = toc;
    %X2 = [x_star x_stela];
    
    en_st = abs(x_stela).^2;
    En_stela = En_stela + en_st;
    
    ratio_stela = norm(x_star)^2/norm(x_star - x_stela)^2;
    Rat_stela = Rat_stela + ratio_stela;
  
    
%     % 3
%     % Adaptive LASSO solution - one-step-approach with LS weight
%     tic;
%     x_ls = pinv(G)*y;
%     [x_wls, nw_it1, S1] = stela3(mu3, G, y, x_ls, 0 , 3);
%     t4=toc;
%     
%      en_wls = abs(x_wls).^2;
%      En_wls = En_wls + en_wls;
%      sup_wls = find(En_wls>sig_n);
%     
%     ratio_wls = norm(x_star)^2/norm(x_star - x_wls)^2;
%     Rat_wls = Rat_wls + ratio_wls; 
%     
%     % 4
%     % Adaptive LASSO solution - one-step-approach with stela weight
%     tic;
%     [x_st, ~, ~] = stela(G, y, 1);
%     [x_wstela, nw_it2, S2] = stela3( mu3, G, y, x_st, 1e-3 , 1);
%     t5=toc;
%     
%      en_wstela = abs(x_wstela).^2;
%      En_wstela = En_wstela + en_wstela;
%      sup_wstela = find(En_wstela>sig_n);
%     
%      ratio_wstela = norm(x_star)^2/norm(x_star - x_wstela)^2;
%      Rat_wstela = Rat_wstela + ratio_wstela;
    
    % 5 
    % IRW-LASSO 
    mu5 = 5e-13;
    tic;
    [x_st, ~, ~] = stela(G, y, mu);
    [x_irwl, nw_it4, S4, gvec2, err_pr ] = irwl_w_stela(  mu5, G, y, x_st,1e-3,1 );
    t5=toc;
    
    en_irwl = abs(x_irwl).^2;
    En_irwl = En_irwl + en_irwl;
    
    ratio_irwl = norm(x_star)^2/norm(x_star - x_irwl)^2;
    Rat_irwl = Rat_irwl + ratio_irwl;
    
    %X = [x_star x x_stela x_ws x_wls x_wstela x_irwl]; 
   
    % Compute the time required
    T1 = T1 + t1;
    T2 = T2 + t2;
    T5 = T5 + t5;
    
end

    MT1(i) = T1/T;
    MT2(i) = T2/T;
    MT5(i) = T5/T;
    
    % Average the energy and check the support
    Av_En_ols = En_ols/T;
    sup_ols = find(Av_En_ols>sig_n^2);
    
    Av_En_stela = En_stela/T;
    sup_stela = find(Av_En_stela>sig_n^2);
    sup_stela2 = find(Av_En_stela>noise_thresh);
      
    Av_En_irwl = En_irwl/T;
    sup_irwl = find(Av_En_irwl>sig_n^2);    
    
    % Average the SER  
    R_ols(i) = Rat_ols/T;
    R_stela(i) = Rat_stela/T;
    R_irwl(i) = Rat_irwl/T;
    
    % Error    
    Av_r_err_ols(i)= norm(E_star - Av_En_ols,2)*100/norm(E_star,2);
    Av_r_err_stela(i)= norm(E_star - Av_En_stela,2)*100/norm(E_star,2);
    Av_r_err_irwl(i)= norm(E_star - Av_En_irwl,2)*100/norm(E_star,2);
 
    
   %%%%%%%%%%%%%%%%%%%%%%%%%% 1. OLS prob. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Soft decision probability OLS
    E_es_ls = Av_En_ols > noise_thresh;
    
    % Identified transmitions
    sup_cor_ls = sum(LL&E_es_ls)*100/sum(LL);
    Total_sup_cor_ols = Total_sup_cor_ols + sup_cor_ls;
    
    % False alarm transmitions
    if sum(E_es_ls)~=0;
    sup_err_ls = (sum(E_es_ls)-sum(LL&E_es_ls))*100/sum(E_es_ls);
    else
        sup_err_ls=0;
    end
    Total_sup_err_ols = Total_sup_err_ols + sup_err_ls;

    % Hard decision probability LS
    Sup_ls = find(Av_En_ols > noise_thresh);
    Sdc_ls = setdiff(L,Sup_ls);
    Sfa_ls = setdiff(Sup_ls,L);
    
    if isempty(Sdc_ls)==1
       Pdc_ols = Pdc_ols + 1;
    end
    
    if isempty(Sfa_ls)==0
       Pfa_ols = Pfa_ols + 1;
    end
       
    %%%%%%%%%%%%%%%%%%%%%% 2. STELA prob.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Soft decision probability STELA
    E_es_stela = Av_En_stela > noise_thresh;
    
    % Identified transmitions
    sup_cor_stela = sum(LL&E_es_stela)*100/sum(LL);
    Total_sup_cor_stela = Total_sup_cor_stela + sup_cor_stela;
    
    % False alarm transmitions
    if sum(E_es_stela)~=0;
    sup_err_stela = (sum(E_es_stela)-sum(LL&E_es_stela))*100/sum(E_es_stela);
    else
        sup_err_stela=0;
    end
    Total_sup_err_stela = Total_sup_err_stela + sup_err_stela;
    
    % Hard decision probability STELA
    Sup_STELA = find(Av_En_stela > noise_thresh);
    S2dc = setdiff(L,Sup_STELA);
    S2fa = setdiff(Sup_STELA,L);
    
    if isempty(S2dc)==1
       Pdc_stela = Pdc_stela + 1;
    end
    
    if isempty(S2fa)==0
       Pfa_stela = Pfa_stela + 1;
    end
        
   %%%%%%%%%%%%%%%%%%%%%%  5. IRWL-LASSO prob.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % Soft decision probability irwLASSO
    E_es_irwl = Av_En_irwl > noise_thresh;
    
    % Identified transmitions
    sup_cor_irwl = sum(LL&E_es_irwl)*100/sum(LL);
    Total_sup_cor_irwl = Total_sup_cor_irwl + sup_cor_irwl;
    
    % False alarm transmitions
    if sum(E_es_irwl)~=0;
    sup_err_irwl = (sum(E_es_irwl)-sum(LL&E_es_irwl))*100/sum(E_es_irwl);
    else
        sup_err_irwl=0;
    end
    Total_sup_err_irwl = Total_sup_err_irwl + sup_err_irwl;

    % Hard decision probability wLS
    Sup_irwl = find(Av_En_irwl > noise_thresh);
    Sdc_irwl = setdiff(L,Sup_irwl);
    Sfa_irwl = setdiff(Sup_irwl,L);
    
    if isempty(Sdc_irwl)==1
       Pdc_irwl = Pdc_irwl + 1;
    end
    
    if isempty(Sfa_irwl)==0
       Pfa_irwl = Pfa_irwl + 1;
    end
  
%  fi=round(i*1000/It)/10;
%  formatSpec = ' %1$3.1f %2$c';
%  waitbar(fi/100,hb,sprintf(formatSpec,fi,'%'));

end

%  close(hb);


 %%%%%%%%%%%%%%%%%%%%%%%%%% 1. OLS computations %%%%%%%%%%%%%%%%%%%%%%%%%%%
  Time_ols = mean(MT1);
  Relative_error_ols = mean(Av_r_err_ols);
  MSER_ols = 10*log10(mean(R_ols));
  
  % Soft Prob
  PDS_LS = Total_sup_cor_ols/It; 
  PFS_LS = Total_sup_err_ols/It;
  
  % Hard Prob
  PD_LS = 100*Pdc_ols/It;
  PF_LS = 100*Pfa_ols/It;
  
 
 %%%%%%%%%%%%%%%%%%%%%%%%%% 2. STELA computations %%%%%%%%%%%%%%%%%%%%%%%%%
  Time_stela = mean(MT2);
  Relative_error_stela = mean(Av_r_err_stela);
  MSER_stela = 10*log10(mean(R_stela));
  
  % Soft Prob
  PDS_STELA = Total_sup_cor_stela/It; 
  PFS_STELA = Total_sup_err_stela/It;
  
  % Hard Prob
  PD_STELA = 100*Pdc_stela/It;
  PF_STELA = 100*Pfa_stela/It;


%%%%%%%%%%%%%%%%%%%%%%%%%% 5. STELA computations %%%%%%%%%%%%%%%%%%%%%%%%%
  Time_irwl = mean(MT5);
  Relative_error_irwl = mean(Av_r_err_irwl);
  MSER_irwl = 10*log10(mean(R_irwl));
  
  % Soft Prob
  PDS_IRWL = Total_sup_cor_irwl/It; 
  PFS_IRWL = Total_sup_err_irwl/It;
  
  % Hard Prob
  PD_IRWL = 100*Pdc_irwl/It;
  PF_IRWL = 100*Pfa_irwl/It;

  matlabpool('close');
  
  mu = 5e-13;
  mu5 = mu;
  Result_table = [0 PDS_LS PFS_LS PD_LS PF_LS Relative_error_ols MSER_ols Time_ols; mu PDS_STELA PFS_STELA PD_STELA PF_STELA Relative_error_stela MSER_stela Time_stela; mu5 PDS_IRWL PFS_IRWL PD_IRWL PF_IRWL Relative_error_irwl MSER_irwl Time_irwl];
  printmat(Result_table, 'Results', 'OLS STELA IRWL-STELA', 'Regul_par Soft_dec_d Soft_dec_fa Hard_dec_d Hard_dec_fa Rel_error MSER Comp_time');

% % Probability of detection and false alarm in signal level
% figure(1);
% subplot(2,1,1);
% %plot(Sensors_vec, Sup_cor0 ,':','Markersize',3,'LineWidth',1.6); % LS
% hold on
% %plot(Sensors_vec, Sup_cor1 ,'o--','Markersize',3,'LineWidth',1.6); % ADMM
% plot(Sensors_vec, Sup_cor2 ,'x-','Markersize',9,'LineWidth',1.6); % STELA
% %plot(Sensors_vec, Sup_cor3 ,'--or','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, Sup_cor_ls , '-.+black' ,'Markersize',9,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, Sup_cor_wls ,'^:cyan','Markersize',9,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, Sup_cor_wstela, '*--yellow','Markersize',9,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, Sup_cor_irwl ,'s-.green','Markersize',9,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% xlabel('number of sensors');
% ylabel('percentage %');
% title('Probability of Detection');
% legend( 'STELA','OLS','AD-LASSO_{OLS}','AD-LASSO_{LASSO}','IRWL-STELA');
% hold off
% grid on
% % Probability of false alarm
% subplot(2,1,2);
% hold on
% %plot(Sensors_vec,Sup_err0,':','Markersize',3,'LineWidth',1.6); % LS
% %plot(Sensors_vec,Sup_err1,'o--','Markersize',3,'LineWidth',1.6); % ADMM
% plot(Sensors_vec,Sup_err2,'x-','Markersize',9,'LineWidth',1.6); % STELA
% %plot(Sensors_vec,Sup_err3,'--or','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, Sup_err_ls ,'-.+black' ,'Markersize',9,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, Sup_err_wls ,'^:cyan','Markersize',9,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, Sup_err_wstela ,'*--yellow','Markersize',9,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, Sup_err_irwl ,'s-.green','Markersize',9,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% xlabel('number of sensors');
% ylabel('percentage %');
% title('Probability of False Alarm');
% %legend( 'STELA','Proposed','OLS','AD-LASSO_{OLS}','AD-LASSO_{LASSO}','IRWL-STELA');
% hold off
% grid on
% 
% % Relative error
% figure(2);
% %plot(Sensors_vec, Relative_error0,'.:','Markersize',3,'LineWidth',1.6); % LS
% hold on
% %plot(Sensors_vec, Relative_error1,'o--','Markersize',3,'LineWidth',1.6); % ADMM
% plot(Sensors_vec, Relative_error2,'x-','Markersize',9,'LineWidth',1.6); % STELA
% %plot(Sensors_vec, Relative_error3,'--or','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, Relative_error_ols ,'-.+black' ,'Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, Relative_error_wls ,'^:cyan','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, Relative_error_wstela ,'*--yellow','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, Relative_error_irwl ,'s-.green', 'Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% xlabel('number of sensors');
% ylabel('percentage %');
% title('Mean Relative Error');
% legend( 'STELA','OLS','AD-LASSO_{OLS}','AD-LASSO_{LASSO}','IRWL-STELA');
% hold off
% grid on
% 
% % Probability of detection and probability of false alarm hard decision
% figure(3);
% subplot(2,1,1);
% %plot(Sup_err0/100, Sup_cor0/100 ,':','Markersize',3,'LineWidth',1.6); % LS
% hold on
% %plot(Sup_err1/100, Sup_cor1/100 ,'o--','Markersize',3,'LineWidth',1.6); % ADMM
% plot(Sensors_vec, PD2*100 ,'x-','Markersize',9,'LineWidth',1.6); % STELA
% %plot(Sensors_vec, PD3*100   ,'--or','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, PD_LS*100 ,'-.+black' ,'Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, PD_WLS*100  ,'^:cyan','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, PD_WSTELA*100   ,'*--yellow','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, PD_IRWL*100   ,'s-.green','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% xlabel('number of sensors');
% ylabel('percentage %');
% title('Probability of Detection - Hard dec.');
% legend( 'STELA','OLS','AD-LASSO_{OLS}','AD-LASSO_{LASSO}','IRWL-STELA');
% hold off
% grid on
% subplot(2,1,2); 
% %plot(Sup_err0/100, Sup_cor0/100 ,':','Markersize',3,'LineWidth',1.6); % LS
% hold on
% %plot(Sup_err1/100, Sup_cor1/100 ,'o--','Markersize',3,'LineWidth',1.6); % ADMM
% plot(Sensors_vec, PF2*100 ,'x-','Markersize',9,'LineWidth',1.6); % STELA
% %plot(Sensors_vec, PF3*100  ,'--or','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, PF_LS*100  ,'-.+black' ,'Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, PF_WLS*100  ,'^:cyan','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, PF_WSTELA*100  ,'*--yellow','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, PF_IRWL*100  ,'s-.green','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% xlabel('number of sensors');
% ylabel('percentage %');
% title('Probability of False Alarm - Hard dec.');
% %legend( 'STELA','Proposed','OLS','AD-LASSO_{OLS}','AD-LASSO_{LASSO}','IRWL-STELA');
% hold off
% grid on
% 
% % Probability of detection and probability of false alarm soft decision
% figure(4);
% subplot(2,1,1);
% %plot(Sup_err0/100, Sup_cor0/100 ,':','Markersize',3,'LineWidth',1.6); % LS
% hold on
% %plot(Sup_err1/100, Sup_cor1/100 ,'o--','Markersize',3,'LineWidth',1.6); % ADMM
% plot(Sensors_vec, PD2S,'x-','Markersize',9,'LineWidth',1.6); % STELA
% %plot(Sensors_vec, PD3S  ,'--or','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, PDS_LS ,'-.+black' ,'Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, PDS_WLS  ,'^:cyan','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, PDS_WSTELA  ,'*--yellow','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, PDS_IRWL ,'s-.green','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% xlabel('number of sensors');
% ylabel('percentage %');
% title('Probability of Detection - Soft dec.');
% legend( 'STELA','OLS','AD-LASSO_{OLS}','AD-LASSO_{LASSO}','IRWL-STELA');
% hold off
% grid on
% subplot(2,1,2);
% %plot(Sup_err0/100, Sup_cor0/100 ,':','Markersize',3,'LineWidth',1.6); % LS
% hold on
% %plot(Sup_err1/100, Sup_cor1/100 ,'o--','Markersize',3,'LineWidth',1.6); % ADMM
% plot(Sensors_vec, PF2S ,'x-','Markersize',9,'LineWidth',1.6); % STELA
% %plot(Sensors_vec, PF3S  ,'--or','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, PFS_LS ,'-.+black' ,'Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, PFS_WLS ,'^:cyan','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, PFS_WSTELA  ,'*--yellow','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, PFS_IRWL  ,'s-.green','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% xlabel('number of sensors');
% ylabel('percentage %');
% title('Probability of False Alarm - Soft dec.');
% %legend( 'STELA','Proposed','OLS','AD-LASSO_{OLS}','AD-LASSO_{LASSO}','IRWL-STELA');
% hold off
% grid on
% 
% % % Stem plots per nember of sensors
% % figure(5);
% % %subplot(2,1,1);
% % stem(MAv_dif2(:,4),'x-','Markersize',7,'LineWidth',1.6); % STELA
% % hold on
% % stem(MAv_dif3(:,4),'--or','Markersize',6,'LineWidth',1.6); % RW-STELA
% % ylabel(' Errors in Energy');
% % title('STELA');
% % %subplot(2,1,2);
% % % ylabel(' Errors in Energy');
% % % title('RW-STELA');
% % legend( 'L-STELA','IRWL-STELA');
% % hold off
% % grid off
% 
% % Running Time evaluation
% figure(5);
% %subplot(2,1,1);
% hold on
% %plot(Sensors_vec, Time0,':','Markersize',3,'LineWidth',1.6); % LS
% %plot(Sensors_vec, Time1,'o--','Markersize',3,'LineWidth',1.6); % ADMM
% plot(Sensors_vec, Time2,'x-','Markersize',9,'LineWidth',1.6); % STELA
% plot(Sensors_vec, Time_ols ,'-.+black' ,'Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, Time_wls  ,'^:cyan','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, Time_wstela  ,'*--yellow','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec, Time_irwl  ,'s-.green','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% xlabel('number of sensors');
% ylabel('time (sec)');
% title('Running time');
% legend( 'STELA','OLS','AD-LASSO_{OLS}','AD-LASSO_{LASSO}','IRWL-STELA');
% hold off 
% grid on
% % subplot(2,1,2);
% % hold on
% % %plot(Sensors_vec, Time0,':','Markersize',3,'LineWidth',1.6); % LS
% % %plot(Sensors_vec, Time1,'o--','Markersize',3,'LineWidth',1.6); % ADMM
% % plot(Sensors_vec, Mnit2,'x-','Markersize',9,'LineWidth',1.6); % STELA
% % plot(Sensors_vec, Mnit3,'--or','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% % xlabel('number of sensors');
% % ylabel('number of iterations');
% % title('Iterations');
% % legend('STELA','RW-STELA')
% % hold off   
% % grid on
% 
% 
% % Mean Signal to error ratio - MSER
% figure(6);
% %plot(Sensors_vec, Relative_error0,'.:','Markersize',3,'LineWidth',1.6); % LS
% hold on
% %plot(Sensors_vec, Relative_error1,'o--','Markersize',3,'LineWidth',1.6); % ADMM
% plot(Sensors_vec, MSER1,'x-','Markersize',9,'LineWidth',1.6); % STELA
% %plot(Sensors_vec, MSER2,'--or','Markersize',5,'LineWidth',1.6, 'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec,  MSER_ols ,'-.+black' ,'Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec,  MSER_wls  ,'^:cyan','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec,  MSER_wstela  ,'*--yellow','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% plot(Sensors_vec,  MSER_irwl  ,'s-.green','Markersize',5,'LineWidth',1.6,'MarkerFaceColor','red'); % Reweighted Stela
% xlabel('number of sensors');
% ylabel('dB');
% title('Mean Signal to Error Ratio - MSER');
% legend( 'STELA','OLS','AD-LASSO_{OLS}','AD-LASSO_{LASSO}','IRWL-STELA');
% 
% hold off
% grid on



