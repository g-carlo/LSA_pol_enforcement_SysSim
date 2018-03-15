% function LTE shannon formula

function y = shannonMkI(x)

thH = 10^(23.5/10);
BW_eff = 0.75;
% BW_eff = 0.75*8;
SINR_eff = 10^(1.25/10);
% SINR_eff = 10^(10/10);
logConst = log(2);

% sinr_table = thH*ones(size(x));

ind1 = x>thH;
ind2 = x<=thH;

% th_1 = BW_eff*log(1 + sinr_table(ind1)/SINR_eff)/logConst;
th_1 = BW_eff*log(1 + thH/SINR_eff)/logConst;
th_2 = BW_eff*log(1 + x(ind2)/SINR_eff)/logConst;

y = zeros(size(x));
y(ind1) = th_1;
y(ind2) = th_2;

