function conf = getInputParameters()

    conf.implementation.numOfCellTypes = 3;

    % physical layer parameters
    conf.PHY.TTIduration  = 1e-3;                      % TTI duration [s]
    conf.PHY.noiseFloor   = 10^(-204/10);              % noise PSD [dB/Hz]

    % UE parameters
    conf.user.noiseFigure       = 10^(9/10);    % UE receiver figure noise
    conf.user.height            =   1.5;

    % used for fast fading implementation
    conf.carrier.availableBW                  = [10e6 10e6 10e6 5e6];                          % available BW 
    conf.carrier.frequency                    = [2e9 2e9  2e9 2.31e9];                         % carrier frequency [GHz]
    conf.carrier.scenarioCase                 = [1 1 1 1];  % area type: 1: urban; 2: suburban
    conf.carrier.loadPtxData                  = [0 0 0 0];                                % load PtxData from file: 0 - no     1 - yes
    conf.carrier.uePercentage                 = [1 1 1 1];
    conf.carrier.numOfAvailablePRB            = [50 50 50 25];                   % number of available PRBs
%     conf.carrier.availableBW                  = [10e6 10e6 10e6 10e6];                          % available BW 
%     conf.carrier.frequency                    = [2e9 2e9 2.31e9 2.31e9];                         % carrier frequency [GHz]
%     conf.carrier.scenarioCase                 = [1 1 1 1];  % area type: 1: urban; 2: suburban
%     conf.carrier.loadPtxData                  = [0 0 0 0];                                % load PtxData from file: 0 - no     1 - yes
%     conf.carrier.uePercentage                 = [1 1 1 1];
%     conf.carrier.numOfAvailablePRB            = [50 50 50 50];                   % number of available PRBs
    conf.carrier.PRBbw                        = conf.carrier.availableBW ...
                                                    ./conf.carrier.numOfAvailablePRB;          % bandwidth of a PRB in each carrier
    
    conf.carrier.user_Pn_PRB = conf.PHY.noiseFloor .* conf.carrier.PRBbw .* conf.user.noiseFigure;      % noise power per PRB

    
    conf.general.attachUsersMethod            = 3;       % method of attaching users to the cells.  1 = max received power  2 = min distance
    conf.general.numOfCarriers                = 4;       % number of carriers
    conf.general.numOfSnapshots               = 1;       % number of snapshots per simulation
    conf.general.numOfTTIperSimulation        = 20;     % simulation length in TTIs
    conf.general.simulatedTime                = 600*24*50;     % simulation length in TTIs   
    conf.general.mode                         = 'montecarlo';  % 'montecarlo', 'static', 'calibration'
    conf.general.SLS_activated                = true;
    conf.general.interferersPerUE             = 15;            % number of interferers for each UE
    conf.general.Prx_threshold                = -Inf;          % received power threshold
    conf.general.randomAssigCellsPerCarrier   = 0;       % 0: determine cells per carrier according to the area specified as a polygon with a sequence of points. 
                                                         % 1: determine cells per carrier randomly.
    conf.general.useWA                        = 1;       % 0: don't compute replica positions to compute distances between UEs and cell. Use the position of each cell.
    conf.general.cellRatio                    = [1 1 1]; % ratio between the cells that can use each carrier and the total number of cells.
    
    
    %%% PARAMAMETERS FOR CALIBRATION MODE ONLY
    conf.general.calibrationCase              = 1;                                % [1 2] 1 = 3GPP calibration 2D      2 = calibration 3D
    
    % deployment parameters
    conf.deployment.cellScenario                = 1;    % 1: macro cells, 2: macro + pico cells, 
                                                        % 3: macro + small cells, 4: a random number of small cells are randomly distributed within a square, 
                                                        % 5: small cells WITH SQUARED GRID DISTRIBUTION WITHIN A SQUARE, 
                                                        % 6: small cells WITH SPPP DISTRIBUTION WITHIN A SQUARE OF MAX. LENGTH
    conf.deployment.UEdistributionConfiguration     = 1;                                % [1 2] kind of user distribution per macro-pico. 1 = uniform within macro / 2 = uniform within pico and macro (see Sim. Req. Doc.)
    conf.deployment.cellDensity                     = 10 / 500^2;                        % (ONLY FOR conf.general.cellScenario == [4,5]) number of cells per square km
    conf.deployment.squareLength                    = 1500;                                % length [in m] of square edge in small scell scenario
    conf.deployment.numOfUEperCellType              = [20 10 3];
    conf.deployment.intersiteDistance               = 500;
    % number of sites for homogeneous scenario with macro cells.
    conf.deployment.numOfSites                      = 19;
    % number of macro cells per site
    conf.deployment.numOfCellsPerSite               = 3;
    conf.deployment.firstSectorAngle                = pi/6;  % first angle that is used to compute positions of macro cells in a hexagonal grid.
    conf.deployment.numOfSectors                    = 3; % number of sectors in each site (macro cells only)
    conf.deployment.numOfNetworkOperators           = 3;  
    conf.deployment.numOfIncumbents                 = 1;
    conf.deployment.incumbentSouthBound             = -2500;%-4750;
    conf.deployment.incumbentNorthBound             = 2500;%+4750;
    conf.deployment.incumbentEastBound              = 2165;%3897;
    conf.deployment.incumbentWestBound              = -2165;%-3897;
    conf.deployment.N_sites_extedndedScenario       = 11^2;
    
    conf.incumbent.lambdaIdle                   = 1/45 * ones( 1 , conf.deployment.numOfIncumbents );
    conf.incumbent.lambdaActive                 = 1/5 * ones( 1 , conf.deployment.numOfIncumbents );
    conf.incumbent.radius                       = 50 * ones( 1 , conf.deployment.numOfIncumbents );
    conf.incumbent.receiverHeight               = 1.5 * ones( 1 , conf.deployment.numOfIncumbents );        % height in meters
    conf.incumbent.frequencyCarrierIndex        = 3 * ones( 1 , conf.deployment.numOfIncumbents );        % height in meters
    conf.incumbent.transmistterHeight           = 10 * ones( 1 , conf.deployment.numOfIncumbents );        % height in meters
    conf.incumbent.propagationClutter           = 6 * ones( 1 , conf.deployment.numOfIncumbents );            % 7 = Dense Urban; 6 = Urban;  5 = Open in urban;  4 = Suburban; 3 = Forest; 2 = Rural;  1 = Water
    conf.incumbent.IR_threshold                 = 10.^(-2./10)* ones( 1 , conf.deployment.numOfIncumbents );
    conf.incumbent.noiseFigure                  = 10.^(9/10) * ones( 1 , conf.deployment.numOfIncumbents );
    conf.incumbent.channelBandwidth             = 5e6* ones( 1 , conf.deployment.numOfIncumbents );
    conf.incumbent.numOfReceivers               = 10 * ones( 1 , conf.deployment.numOfIncumbents );
    conf.incumbent.txPower                      = 10.^((20 -30)./10) / 4 * ones( 1 , conf.deployment.numOfIncumbents );   % TX power is 20dBm over 20MHz, we are interested in the power of 5MHz
    conf.incumbent.interferenceMargin           = 10.^((8)./10)* ones( 1 , conf.deployment.numOfIncumbents );
    conf.incumbent.exclusionZoneDistance        = 900 * ones( 1 , conf.deployment.numOfIncumbents );

    
%     conf.incumbent.lambdaIdle                   = [ 1/10 1/36 ];
%     conf.incumbent.lambdaActive                 = [ 1/10 1/10 ];
%     conf.incumbent.radius                       = [ 50 500 ];
%     conf.incumbent.receiverHeight               = [ 1.5 1.5];        % height in meters
%     conf.incumbent.frequencyCarrierIndex        = [ 3 3];        % height in meters
%     conf.incumbent.transmistterHeight           = [ 10 10];        % height in meters
%     conf.incumbent.propagationClutter           = [ 6 6];            % 7 = Dense Urban; 6 = Urban;  5 = Open in urban;  4 = Suburban; 3 = Forest; 2 = Rural;  1 = Water
%     conf.incumbent.IR_threshold                 = 10.^([-2 -2]./10);
%     conf.incumbent.noiseFigure                  = [10.^(9/10) 10.^(9/10)];
%     conf.incumbent.channelBandwidth             = 5e6;
%     conf.incumbent.numOfReceivers               = [10 10];
%     conf.incumbent.txPower                      = 10.^(([20 20]-30)./10) / 4;   % TX power is 20dBm over 20MHz, we are interested in the power of 5MHz
%     conf.incumbent.interferenceMargin           = 10.^(([8 8])./10);

    
    conf.database.spaceResolution               = 50;
    conf.database.numOfLSAChannels              = 1;
    conf.database.timeRes_T_frame               = 3600;
    conf.database.subframePerFrame              = 6;
    conf.database.ShadFadingSafetyMargin        = 16;

    conf.operator.misBehProbability             = [0.8 0.01 0.0];
    conf.operator.initialPenalty                = [0 0 0];
    conf.operator.cellSpectralEfficiency        = 2.23;  
    
    conf.L1_scheduling.schedulingSlotLength     = 24;
    conf.L1_scheduling.penaltyThreshold         = 25;
    conf.L1_scheduling.penaltyCalculation       = 'exponential'; %Can be either exponential or linear
    conf.L1_scheduling.penaltyLambda            = 100*conf.L1_scheduling.schedulingSlotLength/log(2); %? value for the exponential decay          % with X/log(2), penalty halves at X     
    conf.L1_scheduling.penaltyBeta              = -0.18;   %Decay coefficient for the linear decay
    conf.L1_scheduling.excZoneMethod            = 1;                % 1: based on power;  2: based on distance
    
    conf.detection.num_of_sensors               = 50;
    
    % file parameters
%     conf.storage.workingFolder             = 'C:\Repository\ADEL_sim_v2_rep';
    conf.storage.workingFolder             = './';
%     conf.storage.workingFolder             = 'C:\Users\Carlo\Workspace\Matlab\ADEL_simulator_bitBucket';
    conf.storage.inputFolderPath           = 'inputData/';
    conf.storage.seedFileName              = 'randState';
    conf.storage.mainFolderPath                   = 'C:/Workspace/Matlab/SONsim/';    % SONsim folder path
    conf.storage.outputFolderPath                 = './results/';                       % results folder path
    conf.storage.outputFile                       = 'AS_snap_';                        % output file name

    % results
    conf.results.profiling                 = false;                            % [true/false] simulator profiling: true/false

    % randomness
    conf.randomness.randGeneratorState        = 1;                          % possible values = {1,2,other}. 
                                                                            % 1: computes a new seed value and saves it if EiRP is not performed. 
                                                                            % 2: uses the seed number that is stored in a given file.
                                                                            % other: computes a  new seed number, but it is not saved.
    conf.randomness.randGeneratorStateDecoupling     = true;                            % [true/false]   "true": random generator seed is reset after network deployement
    conf.randomness.usersRandomization               = 2;                                % [1 2] :   1 = change of number of users (for fixed num. of cells) just adds new users and keeps positions of the old ones    2 = change of number of users (for fixed num. of cells) creates completely new users in simulation


    % parameters of antenna

    conf.antenna.horizontalPhi0       = 0;                                   % direction of maximum gain in the horizontal pattern of directional antenna 
    conf.antenna.horizontalPhi3dB     = 70;                                  % 3dB beam-width of the directional antenna horizontal pattern [�]
    conf.antenna.verticalTheta3dB     = 10;                                  % 3dB beam-width of the directional antenna vertical pattern [�]
    conf.antenna.horizontalMinGain    = 25;                                  % minimum gain of the directional antenna horizontal pattern [dB]
    conf.antenna.verticalMinGain      = 20;                                  % minimum gain of the directional antenna vertical pattern [dB]


    % parameters of the channel;
    conf.channel.LOS                          = true;                        % [true/false] LOS/NLOS propagation: true/false
    conf.channel.shadowingFlag                = true;                        % [true/false]shadow fading: true/false  
    conf.channel.fastfadingFlag               = false;                       % [true/false] fast fading
    conf.channel.penetrationLoss              = 20;                          % penetration loss  [dB]

   
    % results parameters

    conf.results.thPerUser                = true;                           % [true/false] throughput per user per TTI:  true/false
    conf.results.additionalData           = true;                           % [true/false] additional data recording:    true/false
    conf.results.TTIblockDataLength       = 10000;                           % TH data block length (if thPerUser is active)

    conf.results.SINR_UEperPRBmatrix      = false;                           % [true/false] SINR statistics  -   1 = true;   0 = false
    conf.results.SINR_UEvector            = false;                           % [true/false]
    conf.results.SINR_statistic           = true;                           % [true/false]
    conf.results.bandPerUE                = false;                           % [true/false]
    conf.results.thPerCell                = false;                           % [true/false] throughput per cell statistics  -   1 = true;   0 = false

    % spectrum allocation param parameters
    conf.specAlloc.softReusePWoffset   = 10^(-9/10);                          % power offset for secondary bands in ICIC with soft frequency reuse 
    conf.specAlloc.RBs_percentage      = 50;                                  

    % FSU parameters

    conf.fsu.resendHIItimer           = 4000;                                % only for MBS - reset MBS algorithm timer
    conf.fsu.picoTimerNew             = 1200;                                % only for PBS -  increase n. of PRBs timer 
    conf.fsu.picoTimerTL              = 2000;                                % only for PBS -  long non-avaiable PRB timer 
    conf.fsu.picoTimerTs              = 500;                                 % only for PBS -  short non-available PRB timer
    conf.fsu.interfListUpdatingTime   = 100;                                 % update interfering PBS list timer       
    conf.fsu.macro_FSU_SINT_th        = 11;                                  % Absolute SINR threshold for FSU algorithm
    conf.fsu.pico_FSU_SINT_th         = 11;
    conf.fsu.deltaN_PRB               = 5;          
    conf.fsu.minPRBperCell            = 5;
    conf.fsu.initPRBoccupancy         = 50; 
    conf.fsu.numSyncSlots             = 5;
    conf.fsu.gammaQoS                 = 0.99;


    % scheduling parameters

    conf.scheduling.UEperTTI          = 10;                                       % maximum number of UEs which can be allocated per TTI
    conf.scheduling.PFalphaSet        = [ 1   1  1 1.2  1.4  1.6  1.8  2  4 ];    % set of alpha paremeters for PF scheduler
    conf.scheduling.PFbetaSet         = [ 2  1.2  1  1    1    1    1   1  1 ];   % set of beta paremeters for PF scheduler
    conf.scheduling.alphaMinQoS       = 1.2e5;                                    % alpha Min QoS parameter for FSU
    conf.scheduling.alphaMinAUTH      = 4.5e5;                                    % alpha Min AUTH parameter for FSU
    conf.scheduling.timerQoSEV        = 500;                                      % timer QoS evaluation timer
    conf.scheduling.manPRBvector      =  5:50;

    % Channel implementation parameters

    conf.ff.TTIblock                     = 30;                          % TTI block length for fast-fading computation
    conf.ff.numOfTaps                    = 6;                           % num of taps for fast-fading implementation
    conf.ff.numOfSinusoids               = 8;                           % num of sinusoids for fast-fading implementation    
    conf.ff.freqResolution                 = 2e5;                         % frequency resolution for FF implementation
    conf.ff.PDPvect                      = [0 -0.9 -4.9 -8.0 -7.8 -23.9];             % taps power delay profile power vector [dB]
    conf.ff.delayVect                    = [0 200 800 1200 2300 3700].*1e-9;          % taps power delay profile delay vector [ns]    
    conf.ff.MTspeedUE                         = 3/3.6;                                     % UE speed for FF implementation [m/s]    
    conf.ff.uniqueFFseqLength            = 3420*15*100;                 % length of pre-loaded FF vector 
    conf.ff.nOfFFtables                  = 20;                          % number of pre-loaded FF tables available
    conf.ff.fastfadingImplementationType = 1;                           % [1 2 3 4] fast-fading implementation scheme: 1-3 = real time computation; 4 = pre-loaded data


    conf.l2s.capacityComputation              = 1;                           % [0 1] link-to-system interface capacity computation approx. : 0 = Shannon's formula;  1 = upper-bounded  Shannon formula
    conf.PHYmodelling.CQIreportType           = 1;                           % [1 3] CQI report implementation scheme:   1 = no delay/no time sampling;   3 = delay and time sampling
    conf.PHYmodelling.CQIreportDelay          = 9;                           % CQI report delay [TTI]
    conf.PHYmodelling.CQIreportSamplingTime   = 10;                          % CQI sampling period [TTI]

    % parameters of traffic generator;

    conf.traffic.trafficGeneratorType             = 1;                       %  [1 2] traffic model type:  1 = full buffer  -  2 = burst traffic
    conf.traffic.FTP_model1_TrafficOfferedPerCell = 10e6;                    % traffic offered per cell          
    conf.traffic.FTP_model1_lambda_perUE          = 0.0417*10;               % average arrival rate per UE [ n� packets per second ]
    conf.traffic.FTP_model1_pckSize               = 0.5*8*1e6;               % packet size [bit]
    
    % to redefine
    conf.traffic.lam_newDem                       = 1/400;


    
    % parameters per cell type

    % each row represents one carrier, each column one cell type
    % antenna gain of:  macro, pico, and small cells

    conf.carrier.cellType.antennaGain       =  [10^(14/10) 10^(5/10) 10^(5/10);     
                                                10^(14/10) 10^(5/10) 10^(5/10);
                                                10^(14/10) 10^(5/10) 10^(5/10);
                                                10^(14/10) 10^(5/10) 10^(5/10)];
%     conf.carrier.cellType.antennaGain       =  [10^(14/10) 10^(5/10) 10^(5/10);     
%                                                 10^(14/10) 10^(5/10) 10^(5/10);
%                                                 10^(14/10) 10^(5/10) 10^(5/10);];                                            


    % 1 = 'Antenna MkI'  - 2 = 'Antenna Calib' 3 = 'Omnidirectional'
    conf.carrier.cellType.antennaType            = [1 1 1;
                                                    1 1 1;
                                                    1 1 1;
                                                    1 1 1];
%     conf.carrier.cellType.antennaType            = [1 1 1;
%                                                     1 1 1;
%                                                     1 1 1];

    % maximum transmission power per carrier and cell type
    conf.carrier.cellType.maxTxPw       = [10^((46-30)/10) 10^((30-30)/10) 10^((27-30)/10);
                                           10^((46-30)/10) 10^((30-30)/10) 10^((5.5-30)/10);
                                           10^((43-30)/10) 10^((30-30)/10) 10^((36.99-30)/10);
                                           10^((43-30)/10) 10^((30-30)/10) 10^((36.99-30)/10)]; 
%     conf.carrier.cellType.maxTxPw       = [10^((46-30)/10) 10^((30-30)/10) 10^((27-30)/10);
%                                            10^((46-30)/10) 10^((30-30)/10) 10^((5.5-30)/10);
%                                            10^((43-30)/10) 10^((30-30)/10) 10^((36.99-30)/10)]; 



    % [1 2] scheduler type   1 = Round-Robin  /  2 = Proportional Fair 
    % scheduling type per carrier and cell type
    conf.carrier.cellType.schedulingType         = [1 1 1;
                                                    1 1 1;
                                                    1 1 1;
                                                    1 1 1];    
%     conf.carrier.cellType.schedulingType         = [1 1 1;
%                                                     1 1 1;
%                                                     1 1 1];    


   
    % [1 2 3 4 5 6] path loss model: 1-2 = (ref. to 3GPP TR 36.814)  
    % 3, 4 small cells       
    % 5 = Modified Hata    6 = Extended Hata (OFCOM)                                                                                            
    conf.carrier.cellType.pathLossModel  = [2 2 6; 
                                            2 2 6; 
                                            2 2 6;
                                            2 2 6;];   
%     conf.carrier.cellType.pathLossModel  = [2 2 6; 
%                                             2 2 6; 
%                                             2 2 6;];      
    
    % carrier per operator: carriers has row index, operators has column
    % index
%     conf.carrier.operator.carrierPerOperator  = [1 0; 0 1; 1 0; 0 1;];   
    conf.carrier.operator.carrierPerOperator  = [1 0 0; 0 1 0; 0 0 1; 1 1 1; ];   
                                                
                                                
    % antenna height per cell type[m]
    conf.cellType.height                 = [32 32 5];      

    % shadowing standard deviation
    conf.cellType.shadowingStdDev        = [8 10 8];

    % scheduling parameters per cell type
    % alpha parameter for PF scheduling
    conf.cellType.alphaPF                = [1 1 1];                                   
    % forgetting factor for PF scheduling
    conf.cellType.gammaPF                = [0.95 0.95 0.95];                                
    % beta parameter for PF scheduling
    conf.cellType.betaPF                 = [1 1 1];                                   

    % spectrum allocation algorithm per cell type
    % [1,2,3,4] spectrum allocation type: 1 = fixed reuse;     2 = FSU;
    % 3 = manual, 4: FSU version II, just when pico cells are part of the network;
    conf.cellType.allocationType         = [1 1 1];                                  

    % [1 2 3] frequency reuse scheme :   1 = FULL BAND reuse (REUSE 1) 
    %/ 2 = REUSE 3 / 3 = FULL REUSE OVER HALF BAND                                           
    conf.cellType.reuseScheme            = [1 1 1];                                   

    % minDistance between two cell types
    % element(j,k) indicates the minimun distance that there must be
    % between cell types j and k
    conf.cellType.otherCellType.minDistance =           [0 75 0;
                                                        75 40 0;
                                                         0  0 0];

    % cell radius per cell type
    % serving area radius [m]
    conf.cellType.radius                  = [0 40 0];                                  

    % minimum distance to UE according to the cell type
    conf.cellType.minDistToUE             = [35 10 10];

    
    % element(j,k) indicates the number of cells type k that are 
    % in each cell type j
    conf.cellType.otherCellType.numOfOtherCells        = [0 4 4;
                                                          0 0 0;
                                                          0 0 0];
                                                      
    
    

    conf.cellType.description               = {'macro';'pico';'small'};

    % FSU parameters per cell type
    conf.cellType.FSU_SINT_th        = [11 11 0];        % Absolute SINR threshold for FSU algorithm
end
