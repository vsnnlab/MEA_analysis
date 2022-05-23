clc;clear;close all;

load('MEA_sample_data.mat'); %% 30s sample data
%%% vm = [t x ch]
%%% tstamp = [t x 1]

SF = 1/(tstamp(2) - tstamp(1)) * 1000; %% Hz
intv = tstamp(2) - tstamp(1); %% ms

wo = 60/(SF/2);
bw = wo/35;
[b,a] = iirnotch(wo,bw, -25); %% parameters for Notch filter

fsig = 0.075/intv; %% ms
fsize = floor(fsig*4)*2+1;
xx = linspace(-fsize/2, fsize/2, fsize);
gf = exp(-xx.^2/(2*fsig^2));
gf = gf/sum(gf); %% gaussian filter

thres_candidate_fold = [4.5 5]; %% candidate for threshold
thN = numel(thres_candidate_fold);
burst_ISI = 10; %% threshold for burst detection

elNum = size(vm, 2);

spknum = nan(elNum, thN);%% pre-defined spk number in each el.
spktiming = cell(elNum, thN);%% pre-defined cell for spk timing
burst_num = nan(elNum, thN);%% pre-defined burst num in each el.
burstRatio = nan(elNum, thN); %% pre-defined burst ratio in each el.

% figure;
%%% if you want to see the trace
%%% you may uncomment figure and subplots in for loop below
for el_i = 1:elNum    
    ddata = vm(:,el_i); %% data in each electrode
    
    filtdata = filter(gf, 1, ddata); %% gaussian filter
    filtered_signal = filter(b,a,filtdata); %% notch filter
    mean_ch = mean(filtered_signal);
    std_ch = std(filtered_signal);
    
    thres_candidate = thres_candidate_fold * std_ch;
    %     quotient = floor(el_i/8)+1;
    %     rem = mod(el_i,8); if rem==0 rem = 8; quotient = quotient - 1; end
    %     posidx = (8-rem)*8 + quotient;
    %     subplot(8,8, posidx); box off; hold on; axis off;
    %     plot(ddata, 'k-');    
    %     ylim([mean_ch-30*std_ch mean_ch+30*std_ch]);
    
    for thres_i = 1:numel(thres_candidate_fold)
        detthr = mean_ch + thres_candidate(thres_i);
        hypertthr = mean_ch - thres_candidate(thres_i);
        interv = 1.5;
        
        tempspk = spike_det_avoid_close_with_hyper_v220418(filtered_signal, detthr, hypertthr, interv, SF);
        
        spknum(el_i, thres_i) = numel(find(tempspk));
        spktiming(el_i, thres_i) = {find(tempspk)};
        
        temp_spkt = find(tempspk)/SF * 1000;
        ISI = diff(temp_spkt); %% inter-spike-interval in each electrode
        burst_idx = find(ISI<burst_ISI);
        burst_idx_finalize = [burst_idx, burst_idx+1];
        temp_i = 1;
        while numel(burst_idx_finalize) ~= numel(unique(burst_idx_finalize))
            if burst_idx_finalize(temp_i,2) == burst_idx_finalize(temp_i+1, 1)
                burst_idx_finalize(temp_i,2) = burst_idx_finalize(temp_i+1,2);
                burst_idx_finalize(temp_i+1,:) = [];
                %             temp_i = temp_i;
            else
                temp_i = temp_i+1;
            end
        end                         %% End until inter-spike-interval
        
        bursting_spk_t = [];
        bursting_num = zeros(size(burst_idx_finalize,1), 1);
        for n_i = 1:size(burst_idx_finalize,1)
            bursting_spk_t = [bursting_spk_t; temp_spkt( (burst_idx_finalize(n_i,1)) : (burst_idx_finalize(n_i, 2)))];
            bursting_num(n_i, 1) = burst_idx_finalize(n_i,2) + 1 - burst_idx_finalize(n_i,1); %% number of burst spikes in a single burst
        end
        bursting_ratio = sum(bursting_num) / numel(temp_spkt); %% Burst ratio among all spikes
        
        burst_num(el_i, thres_i) = sum(bursting_num); 
        burstRatio(el_i, thres_i) = bursting_ratio;
    end
end
%% Firing rate Plotting
for thres_i = 1:numel(thres_candidate_fold)
    fr_in2d = reshape(spknum(:, thres_i), 8, 8) ./ (tstamp(end) / 1000);    
    %%% Firing rate Plotting
    figure;
    imagesc(fr_in2d); colormap jet; axis equal image tight; colorbar;
    title(['Firing rate (Hz) when thres = ', num2str(thres_candidate_fold(thres_i)), 'SD']);
    set(gca, 'XTick', [], 'YTick', []);
end

%% Burst parameters Plotting
for thres_i = 1:numel(thres_candidate_fold)
    burst_in2d = reshape(burst_num(:, thres_i), 8, 8) ./ (tstamp(end) / 1000);
    burstNum_in2d = reshape(burst_num(:, thres_i), 8, 8);
    burstRatio_in2d = reshape(burstRatio(:, thres_i), 8, 8);
    %%% Burst rate plotting
    figure;
    imagesc(burst_in2d); colormap jet; axis equal image tight; colorbar;
    title(['Burst rate (Hz) when thres = ', num2str(thres_candidate_fold(thres_i)), 'SD']);
    set(gca, 'XTick', [], 'YTick', []);
        
    %%% Number of burst spikes in a single burst
    figure;
    imagesc(burstNum_in2d); colormap jet; axis equal image tight; colorbar;
    title(['# of spikes in a single burst when thres = ', num2str(thres_candidate_fold(thres_i)), 'SD']);
    set(gca, 'XTick', [], 'YTick', []);
    
    %%% Burst ratio among all spikes
    figure;
    imagesc(burstRatio_in2d); colormap jet; axis equal image tight; colorbar;
    title(['Burst ratio among all spikes when thres = ', num2str(thres_candidate_fold(thres_i)), 'SD']);
    set(gca, 'XTick', [], 'YTick', []);
end