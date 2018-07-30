function fig1=processResults_fncPenThresh(n_samples,filelabelRes,Thresholds)
% process results

if nargin==2
    Thresholds = [];
end

close all
% clear all
% 
% load('data_adel_1.mat');
isGain = false; 
isMisb = ~isGain;
timeScale = 1;
% n_samples    = 10;
% filelabelRes = 'MISBEH_09_Dec_p_0-3_0-3_0-0_T20';
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

    av_spectrum_vector = av_spectrum_vector./repmat(sum(av_spectrum_vector,1),3,1);
    
    fig1 = figure(1);
    if(~isempty(Thresholds))
        size_Thresholds = size(Thresholds);
        size_samples = size(av_spectrum_vector);
        n = size_samples(1,2)/size_Thresholds(1,2);
        c = zeros(num_of_operators,size_Thresholds(1,2));
        for i = 1:size_Thresholds(1,2)
            for j = 1:num_of_operators
            c(j,i) = mean(av_spectrum_vector(j,(i-1)*n+1:(i)*n-1));
%             c(2,i) = mean(av_spectrum_vector(2,(i-1)*n+1:(i)*n-1));
%             c(3,i) = mean(av_spectrum_vector(3,(i-1)*n+1:(i)*n-1));
            end
        end
        plot(Thresholds,c(1,:),'-vb');
        hold on;
        plot(Thresholds,c(2,:),'--^r');
        plot(Thresholds,c(3,:),'-.og');
        hold off;
        grid on;
        ylim([0 0.8])
        legend('Location','northwest')
        legend('MNO1','MNO2','Unallocated resources')   
        legend('Location','northwest')
        title('LSA SLS Resources vs Penalty')
        xlabel('Penalty Threshold')
        ylabel('Share of LSA resources')
        
    else
        plot(Thresholds,av_spectrum_vector(1,:));
        hold on;
        plot(Thresholds,av_spectrum_vector(2,:));
        plot(Thresholds,av_spectrum_vector(3,:));
    end

end


    filename = ['Res_' filelabelRes];
    save(filename,'n_samples','av_spectrum_vector');
    savefig(fig1,[filename '_1']);
end






