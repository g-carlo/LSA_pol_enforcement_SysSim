% draw Extended Hata model

clear all;

d = 1:1:2000;
Hm = 1.5;
Hb =  10;

PL_hata = zeros(1,length(d));

PL_hata = ExtendedHata(d/1000, 2300, Hm, Hb, 6);

% for n = 1:length(d)
% 	PL_hata(n) = ExtendedHata(d(n)/1000, 1,1, 2300, Hm, Hb, 6);
%     PL_h
% end

figure(1)
plot(d,PL_hata,'b');hold on;
xlabel('Distance [km]')
ylabel('Attenuation [dB]')
grid on;

% % Noise power
% BW = 5e6;
% NoiseFigure = 9;
% N_pow_dB = -174 + 10*log10(BW) + NoiseFigure;
% N_pow    = 10.^(N_pow_dB/10);
% N_to_I_th_dB =  3;
% N_to_I_th =  10.^(N_to_I_th_dB/10);
% P_TX = 43;              % in dBm for 5MHz
% SW_margin = 0;
% 
% Gant_max = 14;
% Gant_min = 14-25;
% 
% P_rx_min = P_TX + Gant_min - PL_hata-SW_margin;
% P_rx_max = P_TX + Gant_max - PL_hata+SW_margin;
% 
% % 
% 
% figure(2)
% plot(d/1000,P_rx_max,'b'); hold on;
% plot(d/1000,P_rx_min,'r');
% plot(d/1000,(N_pow_dB - N_to_I_th_dB)*ones(1,length(d)),'g');