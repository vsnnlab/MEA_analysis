clc;clear;close all;

load('MEA_sample_data_for_LFPevent.mat'); %% 5min sample data, for 4 sample electrode

SF = 1/(tstamp(2) - tstamp(1))*1000;
intv = tstamp(2) - tstamp(1); %% ms
tstart = tstamp(1);
time_of_interest = 0.20*SF; %% half size_of_LFP_signal?

wo = 60/(SF/2);
bw = wo/35;
[b,a] = iirnotch(wo,bw, -25); %% Notch filter

%%% pre-defining the LFP shape parameters
LFP_duration = nan(1, size(vm,2)); %% duration
LFP_tstart = nan(1, size(vm,2)); %% LFP start
LFP_tend = nan(1, size(vm,2)); %% LFP end
LFP_amplitude = nan(1, size(vm,2)); %% LFP amplitude
LFP_minamp = nan(1, size(vm,2)); %% LFP maxamp
LFP_maxamp = nan(1, size(vm,2)); %% LFP minamp

LFP_tstamp = cell(1, size(vm,2));

L = time_of_interest*2+1;
ds_L = (L-1) / 20 + 1; %% size of LFP when downsampled into 1ms scale
LFP_average_signal = nan(ds_L, size(vm,2));

fig_bool = true;
for elec_iter = 1:size(vm,2)
    if mod(elec_iter, size(vm,2)/4) == 0
        elec_iter
    end
    sample_vm = vm(:,elec_iter);
    vm_stop = filter(b,a,sample_vm); %% Notch filter applied
    %% Low-Pass Filter signal
    lpFilt = designfilt('lowpassiir','FilterOrder',4,'PassbandFrequency',150,'DesignMethod','cheby1','SampleRate',SF);
    vm_stop_low = filter(lpFilt, vm_stop);
    mean_stop_low = mean(vm_stop_low);
    std_stop_low = std(vm_stop_low);
    %% LFP from threshold value
    negval_threshold = mean_stop_low-4*std_stop_low;
    
    temp = (imregionalmin(vm_stop_low).*(vm_stop_low<negval_threshold));
    idx_of_LFP_th = find(temp);
    LFP_tstamp(1, elec_iter) = {idx_of_LFP_th};
    
    disp([num2str(numel(idx_of_LFP_th)), ' LFP detected']);
    temp_LFP_w_threshold = [];
    for th_iter = 1:size(idx_of_LFP_th,1)
        tidx = (idx_of_LFP_th(th_iter) - time_of_interest):(idx_of_LFP_th(th_iter)+time_of_interest);
        if tidx(1) <= 0
            size_signal = tidx(end);
            signal = [zeros(length(tidx) - size_signal,1); vm_stop_low(1:tidx(end),1)];
        elseif tidx(end) > size(vm_stop_low,1)
            size_signal = length(tidx(1):size(vm_stop_low,1));
            signal = [vm_stop_low(tidx(1):end,1); zeros(length(tidx) - size_signal,1)];
        else
            signal = vm_stop_low(tidx,1);
        end
        temp_LFP_w_threshold = cat(2, temp_LFP_w_threshold, signal);
    end
    if fig_bool & ~ isempty(temp_LFP_w_threshold)
        figure; hold on;
        plot([-time_of_interest:time_of_interest]/SF * 1000, temp_LFP_w_threshold*1e5, '-');
        plot([-time_of_interest:time_of_interest]/SF * 1000, mean(temp_LFP_w_threshold*1e5, 2), 'k-', 'LineWidth', 2);
        title(['Sample raw trace of LFP events, electrode #', num2str(elec_iter)]);
        xlim([-75 75]);
        xlabel('Time (ms)'); ylabel('uV');
    end
    %% LFP shape detection
    ds_ratio = SF * 0.001; %% downsample to 1ms timescale
    expectedL = (size(temp_LFP_w_threshold,1)-1) / ds_ratio+1;
    downsampled_LFP = nan(expectedL, size(temp_LFP_w_threshold, 2));
    for temp_LFP_iter = 1:size(temp_LFP_w_threshold, 2)
        tempLFP = temp_LFP_w_threshold(1:end-1, temp_LFP_iter);
        reshapedLFP = reshape(tempLFP, [], expectedL-1);
        downsampled_LFP(:, temp_LFP_iter) = [nanmean(reshapedLFP)'; tempLFP(end)];
    end
    
    %%% To determine the LFP duration; use bootstrap method
    downsampled_signal = nanmean(reshape(vm_stop_low, [], size(vm_stop_low, 1) / ds_ratio)', 2);
    bootstrapN = 3000;
    surrogateVal = nan(1, bootstrapN);
    disp('Making surrogate data');
    for ii=1:bootstrapN
        random_index = randsample(1:size(downsampled_signal,1), size(downsampled_LFP,2));
        surrogateVal(1,ii) = mean(downsampled_signal(random_index,1));
    end
    pval = 0.01; %% Originally it was built as 0.001 when 20min recording was used.
    surrogateVal = sort(surrogateVal);
    uppth = surrogateVal(end-pval*bootstrapN);
    lowth = surrogateVal(pval*bootstrapN);
    mean_signal = mean(vm_stop_low);
    avg_LFP = mean(downsampled_LFP, 2);
    LFP_average_signal(:, elec_iter) = avg_LFP;
    ds_tstamp = [(-time_of_interest/SF):(ds_ratio/SF):(time_of_interest/SF)]*1000; %% ms    
    
    %%% If signal does not exceeds threshold more than 5ms, it is considered as a stable signal
    stable_T = (5/1000) * (SF / ds_ratio);
    middle_index = find(ds_tstamp==0);
    stable_index = find(avg_LFP > lowth & avg_LFP < uppth);
    exceed_index = find(avg_LFP < lowth | avg_LFP > uppth);
    
    negative_time = stable_index (stable_index < middle_index);
    temp1 = find(diff(negative_time) ~= 1);
    temp2 = find(diff(temp1) > stable_T);
    if ~isempty(temp2)
        temp2 = temp2(end);
        idx1 = temp1(temp2+1);
    elseif ~isempty(temp1)
        idx1 = temp1(1);
    else
        idx1 = numel(negative_time);
    end
    LFP_tstart(1, elec_iter) = ds_tstamp(negative_time(idx1));
    
    positive_time = stable_index ( stable_index > middle_index );
    temp1 = find(diff(positive_time) ~= 1);
    temp2 = find(diff(temp1) > stable_T);
    if ~isempty(temp2)
        temp2 = temp2(1);
        idx1 = temp1(temp2);
    elseif ~isempty(temp1)
        idx1 = temp1(end);
    else
        idx1 = 0;
    end
    LFP_tend(1, elec_iter) = ds_tstamp(positive_time(idx1+1));
    %%% Line for LFP amplitude
    [maxval,maxpos] = max(avg_LFP);
    [minval,minpos] = min(avg_LFP);
    
    %%% Line for LFP duration
    %%% Line for LFP amplitude
    [maxval,maxpos] = max(avg_LFP);
    [minval,minpos] = min(avg_LFP);
        
    LFP_duration(1, elec_iter) = LFP_tend(1, elec_iter) - LFP_tstart(1, elec_iter);
    LFP_amplitude(1, elec_iter) = (maxval - minval)*1e5;
    LFP_minamp(1, elec_iter) = minval*1e5;
    LFP_maxamp(1, elec_iter) = maxval*1e5;
    
    %% Plotting    
    if fig_bool & ~ isempty(temp_LFP_w_threshold)
        figure; hold on;
        plot(ds_tstamp, avg_LFP*1e5, 'k-', 'LineWidth', 2);
        line([-200 200], [uppth uppth]*1e5, 'Color', 'b');
        line([-200 200], [lowth lowth]*1e5, 'Color', 'b');
        plot(ds_tstamp(exceed_index), max(avg_LFP)*1e5*1.1, 'r.');
        xlim([-150 150]);
        ylim([min(avg_LFP)*1e5*1.3 max(avg_LFP)*1e5*1.3]);
        line([LFP_tstart(1, elec_iter) LFP_tstart(1, elec_iter)], [min(avg_LFP)*1e5*1.3 max(avg_LFP)*1e5*1.3], 'Color', 'k', 'LineStyle', '--');
        line([LFP_tend(1, elec_iter) LFP_tend(1, elec_iter)], [min(avg_LFP)*1e5*1.3 max(avg_LFP)*1e5*1.3], 'Color', 'k', 'LineStyle', '--');
        line([ds_tstamp(1) ds_tstamp(maxpos)], [maxval*1e5 maxval*1e5], 'Color', 'g', 'LineStyle', '--');
        line([ds_tstamp(1) ds_tstamp(minpos)], [minval*1e5 minval*1e5], 'Color', 'g', 'LineStyle', '--');
        xlabel('Time (ms)');
        ylabel('uV');
        title(['Average LFP event shape, electrode #', num2str(elec_iter)]);
    end
end

%% Last, LFP event frequency
LFP_rate = nan(1, size(vm,2)); %% LFP minamp
for elec_iter = 1:size(vm,2)
    LFP_rate(1,elec_iter) = numel(LFP_tstamp{1,elec_iter}) ./ size(vm,1) * SF;
end
%% Show in command window
LFP_rate
LFP_duration
LFP_amplitude
LFP_minamp
LFP_maxamp
