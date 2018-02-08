% Calibration

[ rawData ] = getLabJackRawData( {'AIN12'}, 5000, 2, 5 );
 
plot(rawData)

% ylim([0 6])


try
    n = length(avg)+1;
catch
    n = 1;
end

mmH2O(n,:) = 0;
raw(n,:) = rawData(:);
avg(n,:) = mean(rawData);


save('calibration.mat')


scatter(mmH2O,avg,'filled')
grid minor


%%
inH2O = mmH2O/25.4;
scatter(avg,inH2O,'filled')
xlabel('Signal Voltage [V]')
ylabel('Pressure Difference [inH2O]')
grid minor