% process results

% close all
% clear all
% 
% load('data_adel_1.mat');
isGain = false; 
isMisb = ~isGain;
timeScale = 6;

if isGain 

    num_of_operators = conf.deployment.numOfNetworkOperators;
    num_of_incumbents = conf.deployment.numOfIncumbents;

    lastElemIdx = length(resultsLogQueue);
    if isempty(resultsLogQueue(lastElemIdx).finishingTime)
        lastElemIdx = lastElemIdx - 1;
    end

    timeLineVector = zeros(1,2*lastElemIdx);
    throughputOverTime = zeros(num_of_operators,2*lastElemIdx);
    assignedSpectrumOverTime = zeros(num_of_operators,2*lastElemIdx);
    isMisbehavingOverTime = zeros(2,2*lastElemIdx);
    incActivity = zeros(num_of_incumbents,2*lastElemIdx);
    overallIncActivity = zeros(num_of_incumbents,2*lastElemIdx);
    inc_SINR_vector = [];
    inc_SNR_vector = [];
    inc_SNR_IM_vector = [];
    prob_detection = [];
    num_active_BS = [];
    prob_false_alarm = [];
    all_BS = false(1,57);

    probVector = struct( 'P_d' , cell(1,57),'P_fa' , cell(1,57),'P_d_mean' , cell(1,57),'P_fa_mean',cell(1,57));

    for nn = 1:lastElemIdx

        timeLineVector(nn*2-1) = resultsLogQueue(nn).startingTime;
        timeLineVector(nn*2) = resultsLogQueue(nn).finishingTime;
        resObject = resultsLogQueue(nn).resObj;
        tmp = sum(resObject.out);
        throughputOverTime(:,nn*2-1) = tmp.';
        throughputOverTime(:,nn*2) = tmp.';
        incActivity(:,nn*2-1) = resObject.incumbentActivity.';
        incActivity(:,nn*2) = resObject.incumbentActivity.';
        overallIncActivity(:,nn*2) = resObject.overallIncumbentActivity.';
        overallIncActivity(:,nn*2-1) = resObject.overallIncumbentActivity.';
        if resObject.SNR_vec ~= 0;
            inc_SINR_vector     = [ inc_SINR_vector 10*log10(resObject.SINR_vec) ];
            inc_SNR_vector      = [ inc_SNR_vector 10*log10(resObject.SNR_vec) ];
            inc_SNR_IM_vector      = [ inc_SNR_IM_vector 10*log10(resObject.SNR_IM_vec) ];
        end
    %     if any(resObject.detection.activeBSs_idx) ~= 0 
    %         detected_BS = all_BS; 
    %         detected_BS(resObject.detection.detectedBSs_idx(resObject.detection.detectedBSs_idx<=57) ) = true;
    %         active_BS   = all_BS; 
    %         active_BS(resObject.detection.activeBSs_idx) = true;
    %         
    %         D_set = detected_BS & active_BS;
    %         P_d_snap = sum(D_set) / sum(active_BS);
    %         
    %         F_set = detected_BS & ~(D_set);
    %         P_fa_snap = sum(F_set) / sum(detected_BS); 
    %         prob_detection = [ prob_detection P_d_snap];
    %         num_active_BS = [ num_active_BS sum(active_BS) ];
    %         prob_false_alarm = [ prob_false_alarm P_fa_snap];
    %         
    %         if length(setdiff(resObject.detection.detectedBSs_idx,resObject.detection.activeBSs_idx)) > 10
    %            disp('248stop' )
    %         end
    %     end


    end

    % compute average throughput
    cum_th_over_time = zeros(num_of_operators,1);
    cum_spectrum_over_time = zeros(num_of_operators,1);
    cum_th_over_time_inc = 0;
    cum_th_over_time_inc_1 = 0;

    cellCapacity = conf.operator.cellSpectralEfficiency;
    num_cells = conf.deployment.numOfSites*conf.deployment.numOfCellsPerSite;
    th_owned_band = cellCapacity * num_cells * 10e6;

    for nn = 1:(lastElemIdx-1)

        cum_th_over_time = cum_th_over_time + (throughputOverTime(:,nn*2-1) - th_owned_band)*...
            ( timeLineVector(nn*2) - timeLineVector(nn*2-1)+1);
        cum_spectrum_over_time = cum_spectrum_over_time + assignedSpectrumOverTime(:,nn*2-1)*...
            ( timeLineVector(nn*2) - timeLineVector(nn*2-1)+1);
        cum_th_over_time_inc = cum_th_over_time_inc + (sum(overallIncActivity(:,nn*2-1)) )*...
            ( timeLineVector(nn*2) - timeLineVector(nn*2-1)+1) * cellCapacity * 5e6;
        cum_th_over_time_inc_1 = cum_th_over_time_inc_1 + ( sum(overallIncActivity(:,nn*2-1))>0 )* conf.deployment.numOfIncumbents *...
            ( timeLineVector(nn*2) - timeLineVector(nn*2-1)+1) * cellCapacity * 5e6;

    end

    av_th = cum_th_over_time/timeLineVector(end);

    av_spectrum = cum_spectrum_over_time/timeLineVector(end);

    misBehProb = sum(isMisbehavingOverTime(:,1:2:end),2)/lastElemIdx

    th_incumb = cum_th_over_time_inc/timeLineVector(end);
    
    th_incumb_1 = cum_th_over_time_inc_1/timeLineVector(end);
    th_incumb_2 = conf.deployment.numOfIncumbents *  cellCapacity * 5e6;

    gainASE = (sum(av_th)+th_incumb)/th_incumb
    
    gainASE_1 = (sum(av_th)+th_incumb_1)/th_incumb_1


    figure(1)

    plot(timeLineVector/timeScale,throughputOverTime(1,:)/1e9);hold on;
    plot(timeLineVector/timeScale,throughputOverTime(2,:)/1e9,'r');
    axis([0 max(timeLineVector/timeScale) 0 1.1*max(max(throughputOverTime(1,:)/1e9),max(throughputOverTime(2,:)/1e9))] );
    xlabel('Time [hours]');
    ylabel('Throughput per Network Operator [Gbps]');


    figure(2)
    drawCDF( inc_SNR_vector , 2 , [1 0 0]  ); 
    drawCDF( inc_SNR_IM_vector , 2 , [0 1 0]  ); 
    drawCDF( inc_SINR_vector , 2 , [0 0 1]  ); %hold on;
%     drawCDF( inc_SINR_vector , 2 , [0 1 1]  ); %hold on;
    xlabel('SINR [dB]');
    legend('SNR (no interference)','SNR (with 8dB int. margin)','SINR (std. LSA)','SINR (ADEL LSA)');
    axis([-10 60 0 100]);



% figure(2)
% subplot(1,2,1)
% plot(timeLineVector,throughputOverTime(1,:)/1e9);
% axis([0 max(timeLineVector) 0 2.5]);
% subplot(1,2,2)
% plot(timeLineVector,throughputOverTime(2,:)/1e9);
% axis([0 max(timeLineVector) 0 2.5]);

else 
    
    num_of_operators = conf.deployment.numOfNetworkOperators;

    lastElemIdx = length(resultsLogQueue);
    if isempty(resultsLogQueue(lastElemIdx).finishingTime)
        lastElemIdx = lastElemIdx - 1;
    end

    timeLineVector = zeros(1,2*lastElemIdx);
    throughputOverTime = zeros(num_of_operators,2*lastElemIdx);
    assignedSpectrumOverTime = zeros(num_of_operators,2*lastElemIdx);
    isMisbehavingOverTime = zeros(num_of_operators,2*lastElemIdx);
    incActivity = zeros(1,2*lastElemIdx);
    inc_SINR_vector = [];
    inc_SNR_vector = [];
    inc_SNR_IM_vector = [];
    prob_detection = [];
    num_active_BS = [];
    prob_false_alarm = [];
    all_BS = false(1,57);

    probVector = struct( 'P_d' , cell(1,57),'P_fa' , cell(1,57),'P_d_mean' , cell(1,57),'P_fa_mean',cell(1,57));

    for nn = 1:lastElemIdx

        timeLineVector(nn*2-1) = resultsLogQueue(nn).startingTime;
        timeLineVector(nn*2) = resultsLogQueue(nn).finishingTime;
        resObject = resultsLogQueue(nn).resObj;
        tmp = sum(resObject.out);
        throughputOverTime(:,nn*2-1) = tmp.';
        throughputOverTime(:,nn*2) = tmp.';
        isMisbehavingOverTime(:,nn*2-1) = resObject.isMisbehaving.';
        isMisbehavingOverTime(:,nn*2) = resObject.isMisbehaving.';
        assignedSpectrumOverTime(:,nn*2-1) = resObject.assignedSpectrum.';
        assignedSpectrumOverTime(:,nn*2) = resObject.assignedSpectrum.';
        incActivity(:,nn*2-1) = resObject.incumbentActivity;
        incActivity(:,nn*2) = resObject.incumbentActivity;
        if resObject.SNR_vec ~= 0;
            inc_SINR_vector     = [ inc_SINR_vector 10*log10(resObject.SINR_vec) ];
            inc_SNR_vector      = [ inc_SNR_vector 10*log10(resObject.SNR_vec) ];
            inc_SNR_IM_vector      = [ inc_SNR_IM_vector 10*log10(resObject.SNR_IM_vec) ];
        end
    %     if any(resObject.detection.activeBSs_idx) ~= 0 
    %         detected_BS = all_BS; 
    %         detected_BS(resObject.detection.detectedBSs_idx(resObject.detection.detectedBSs_idx<=57) ) = true;
    %         active_BS   = all_BS; 
    %         active_BS(resObject.detection.activeBSs_idx) = true;
    %         
    %         D_set = detected_BS & active_BS;
    %         P_d_snap = sum(D_set) / sum(active_BS);
    %         
    %         F_set = detected_BS & ~(D_set);
    %         P_fa_snap = sum(F_set) / sum(detected_BS); 
    %         prob_detection = [ prob_detection P_d_snap];
    %         num_active_BS = [ num_active_BS sum(active_BS) ];
    %         prob_false_alarm = [ prob_false_alarm P_fa_snap];
    %         
    %         if length(setdiff(resObject.detection.detectedBSs_idx,resObject.detection.activeBSs_idx)) > 10
    %            disp('248stop' )
    %         end
    %     end


    end

    % compute average throughput
    cum_th_over_time = zeros(num_of_operators,1);
    cum_spectrum_over_time = zeros(num_of_operators,1);

    cellCapacity = conf.operator.cellSpectralEfficiency;
    num_cells = conf.deployment.numOfSites*conf.deployment.numOfCellsPerSite;
    th_owned_band = cellCapacity * num_cells * 10e6;

    for nn = 1:(lastElemIdx-1)

        cum_th_over_time = cum_th_over_time + (throughputOverTime(:,nn*2-1) - th_owned_band)*...
            ( timeLineVector(nn*2) - timeLineVector(nn*2-1)+1);
        cum_spectrum_over_time = cum_spectrum_over_time + assignedSpectrumOverTime(:,nn*2-1)*...
            ( timeLineVector(nn*2) - timeLineVector(nn*2-1)+1);

    end

    av_th = cum_th_over_time/timeLineVector(end);

    av_spectrum = cum_spectrum_over_time/timeLineVector(end);

    misBehProb = sum(isMisbehavingOverTime(:,1:2:end),2)/lastElemIdx

    th_incumb = conf.deployment.numOfIncumbents *  cellCapacity * 5e6;

    gainASE = (sum(av_th)+th_incumb)/th_incumb


    % for bb = 1:length(num_active_BS)
    %     
    %     probVector(num_active_BS(bb)).P_d = [ probVector(num_active_BS(bb)).P_d  prob_detection(bb ) ];
    %     probVector(num_active_BS(bb)).P_fa = [ probVector(num_active_BS(bb)).P_fa  prob_false_alarm(bb ) ];  
    %     
    % end
    % 
    % for bb = 1:57
    %     
    %     if ~isempty(probVector(num_active_BS(bb)).P_d)
    %         probVector(bb).P_d_mean = mean(probVector(bb).P_d);  
    %         probVector(bb).P_fa_mean = mean(probVector(bb).P_fa);  
    %     end
    %     
    % end

    figure(1)
    % <<<<<<< working copy
%     subplot(211)
    plot(timeLineVector/timeScale,throughputOverTime(1,:)/1e9);hold on;
    plot(timeLineVector/timeScale,throughputOverTime(2,:)/1e9,'r');
    plot(timeLineVector/timeScale,throughputOverTime(3,:)/1e9,'r');
%     axis([0 max(timeLineVector/timeScale) 0 1.1*max(max(throughputOverTime(1,:)/1e9),max(throughputOverTime(2,:)/1e9))] );
    xlabel('Time [hours]');
    ylabel('Throughput per Network Operator [Gbps]');
    axis([0 150 0 1.1*max(max(throughputOverTime(1,:)/1e9),max(throughputOverTime(2,:)/1e9))] );
    legend('MNO 1','MNO 2');
    

%     figure(1);
%     subplot(212)
%     plot(timeLineVector/timeScale,incActivity);hold on;
%     axis([0 max(timeLineVector/timeScale) 0 2] );
%     xlabel('Time [hours]');
%     ylabel('Incumbent activity');

    figure();
    plot(timeLineVector/timeScale,isMisbehavingOverTime(1,:));hold on;
    plot(timeLineVector/timeScale,isMisbehavingOverTime(2,:),'r');hold on;
    axis([0 max(timeLineVector/timeScale) 0 2] );
    xlabel('Time [hours]');
    ylabel('Incumbent activity');

    % figure()
%     drawCDF( inc_SINR_vector , 2 , [0 0 1]  ); %hold on;
%     drawCDF( inc_SNR_vector , 2 , [1 0 0]  ); 
%     drawCDF( inc_SNR_IM_vector , 2 , [0 1 0]  ); 
%     xlabel('SINR [dB]');
%     legend('SINR','SNR (no interference)','SNR (with 8dB int. margin');



    % figure(2)
    % subplot(1,2,1)
    % plot(timeLineVector,throughputOverTime(1,:)/1e9);
    % axis([0 max(timeLineVector) 0 2.5]);
    % subplot(1,2,2)
    % plot(timeLineVector,throughputOverTime(2,:)/1e9);
    % axis([0 max(timeLineVector) 0 2.5]);

    
end

