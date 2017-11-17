% MISBEHAVIOR DETECTION FOR TCD LSA SIMULATOR Core Function 
% 5/12/2016 
% Kostas Voulgaris
% This function simulates the detection algorithm with the use of
% pre-compiled tables. 

function [Est_Loc_active_BS, Pd_table, Pfa_table, Est_energy_active_BS, Comp_time] = Core_function_simulate(N_white_BS, NSens, BS_positions, IU_Loc, Loc_active_BS,file_label)
tic
N = numel(Loc_active_BS); %Number of transmitting LUs --> determines the probabilities Pd and Pfa

Probabilities = [
0.950000000	0.008750000	0.000223375;
0.955000000	0.049993175	0.006136325;
0.947083325	0.023922700	0.001967475;
0.947500000	0.031464000	0.003820850;
0.943250000	0.037368450	0.005047950;
0.943333500	0.041212250	0.006421425;
0.946071425	0.045883900	0.011075000;
0.941562500	0.065631350	0.018826400;
0.940555575	0.060901675	0.017005300;
0.942625000	0.062381350	0.019893550;
0.935340900	0.064701050	0.020434800;
0.935833425	0.066571375	0.023833300;
0.931057600	0.080163050	0.033863650;
0.933481950	0.081239350	0.035784975;
0.933416675	0.078813050	0.037380950;];


Pd = Probabilities(N,1); %Pd from probabilities table 
Pfa = Probabilities(N,3); %Pfa from probabilities table

Loc_inactive_BS = setdiff(1:57,Loc_active_BS);

Est_active_BS = zeros(57,1);
Est_active_BS(Loc_active_BS) = rand(N,1)<=Pd;
Est_active_BS(Loc_inactive_BS) = rand(57-N,1)<=Pfa;

Est_Loc_active_BS = find(Est_active_BS==1);

LU_power = 0.72;

Pd_table = Pd;
Pfa_table = Pfa;
Comp_time = toc;
Est_energy_active_BS = LU_power * Est_active_BS;

end

