function [A] = spike_det_avoid_close_with_hyper_v220418(ddata, detthr, hypertthr, interv, SF)
% Phase 2
% Spike detector using point above +thr and below -thr in 1ms

% interv = 1;         % Interval between depolarization and hyperpolarization
iInt = interv/1000 * SF;      % Interval in index

temp = (imregionalmax(ddata).*(ddata>detthr)) + (imregionalmin(ddata).*(ddata<hypertthr));
ftemp = find(temp);

for ff = 1:sum(temp)
    temp(ftemp(ff)+1:(min(ftemp(ff)+iInt, length(ddata)))) = 0;
end

A = temp;
end

