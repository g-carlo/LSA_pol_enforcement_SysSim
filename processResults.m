% process results

close all
% clear all
% 
% load('data_adel_1.mat');
isGain = false; 
isMisb = ~isGain;
timeScale = 1;
n_samples    = 10;
filelabelRes = 'MISBEH_09_Dec_p_0-3_0-3_0-0_T10';
% filelabelRes = 'GAIN_10_incumbents_dist_';
num_of_operators = 3;

if isGain 

    % variables
    av_th_vector       = zeros(num_of_operators,n_samples);
    av_spectrum_vector = zeros(num_of_operators,n_samples);
    misBehProb_vector  = zeros(num_of_operators,n_samples);
    th_incumb_vector   = zeros(num_of_operators,n_samples);
    gainASE_vector     = zeros(num_of_operators,n_samples);



    for n_s = 1:n_samples

        fileNameProcessRes = [ filelabelRes '_n_' num2str(n_s)];
        load(fileNameProcessRes);

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
        incActivity = zeros(1,2*lastElemIdx);
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
            isMisbehavingOverTime(:,nn*2-1) = resObject.isMisbehaving.';
            isMisbehavingOverTime(:,nn*2) = resObject.isMisbehaving.';
            assignedSpectrumOverTime(:,nn*2-1) = resObject.assignedSpectrum.';
            assignedSpectrumOverTime(:,nn*2) = resObject.assignedSpectrum.';
            overallIncActivity(:,nn*2) = resObject.overallIncumbentActivity.';
            overallIncActivity(:,nn*2-1) = resObject.overallIncumbentActivity.';
            if resObject.SNR_vec ~= 0;
                inc_SINR_vector     = [ inc_SINR_vector 10*log10(resObject.SINR_vec) ];
                inc_SNR_vector      = [ inc_SNR_vector 10*log10(resObject.SNR_vec) ];
                inc_SNR_IM_vector      = [ inc_SNR_IM_vector 10*log10(resObject.SNR_IM_vec) ];
            end

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
        misBehProb = sum(isMisbehavingOverTime(:,1:2:end),2)/lastElemIdx;
        th_incumb = cum_th_over_time_inc/timeLineVector(end);
        th_incumb_1 = cum_th_over_time_inc_1/timeLineVector(end);
        th_incumb_2 = conf.deployment.numOfIncumbents *  cellCapacity * 5e6;

        gainASE = (sum(av_th)+th_incumb)/th_incumb;

        av_th_vector(:,n_s)       = av_th;
        av_spectrum_vector(:,n_s) = av_spectrum;
        misBehProb_vector(:,n_s)  = misBehProb;
        th_incumb_vector(n_s)     = th_incumb;
        gainASE_vector(n_s)       = gainASE;

        dispSimTime( n_s,n_samples ) 

    end
   
    
    
else

    % variables
    av_th_vector       = zeros(num_of_operators,n_samples);
    av_spectrum_vector = zeros(num_of_operators,n_samples);
    misBehProb_vector  = zeros(num_of_operators,n_samples);
    th_incumb_vector   = zeros(num_of_operators,n_samples);
    gainASE_vector     = zeros(num_of_operators,n_samples);




    for n_s = 1:n_samples

        fileNameProcessRes = [ filelabelRes '_n_' num2str(n_s)];
        load(fileNameProcessRes);

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

        end

        % compute average throughput
        cum_th_over_time = zeros(num_of_operators,1);
        cum_spectrum_over_time = zeros(num_of_operators,1);

        cellCapacity = conf.operator.cellSpectralEfficiency;
        num_cells = conf.deployment.numOfSites*conf.deployment.numOfCellsPerSite;
        th_owned_band = cellCapacity * num_cells * 10e6;

        for nn = 1:(lastElemIdx-1)

            cum_th_over_time = cum_th_over_time + (throughputOverTime(:,nn*2-1) - th_owned_band)*...
                ( timeLineVector(nn*2) - timeLineVector(nn*2-1));
            cum_spectrum_over_time = cum_spectrum_over_time + assignedSpectrumOverTime(:,nn*2-1)*...
                ( timeLineVector(nn*2) - timeLineVector(nn*2-1));

        end

        av_th = cum_th_over_time/timeLineVector(end);
        av_spectrum = cum_spectrum_over_time/timeLineVector(end);
        misBehProb = sum(isMisbehavingOverTime(:,1:2:end),2)/lastElemIdx;
        th_incumb = conf.deployment.numOfIncumbents *  cellCapacity * 5e6;
        gainASE = (sum(av_th)+th_incumb)/th_incumb;

        av_th_vector(:,n_s)       = av_th;
        av_spectrum_vector(:,n_s) = av_spectrum;
        misBehProb_vector(:,n_s)  = misBehProb;
        th_incumb_vector(n_s)     = th_incumb;
        gainASE_vector(n_s)       = gainASE;

        dispSimTime( n_s,n_samples ) 

    end


    figure(1);
    plot(cumsum(av_spectrum_vector(1,:)) ./ (1:n_samples) );
    hold on;
    plot(cumsum(av_spectrum_vector(2,:)) ./ (1:n_samples) );
    hold on;
    plot(cumsum(av_spectrum_vector(3,:)) ./ (1:n_samples) );
    axis([1 n_samples 0 1])

end





