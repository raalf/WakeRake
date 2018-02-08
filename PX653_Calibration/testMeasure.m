N = 600;

avg = nan(N,1);
inH2O = nan(N,1);

for x = 1:N

[ rawData ] = getLabJackRawData( {'AIN12'}, 5000, 1, 5 );


avg(x) = mean(rawData);

inH2O(x) = 0.06928*avg(x) - 0.06843;

idx = (x-30:x);

try
    plot(idx, inH2O(idx), 'o-')
    title(sprintf('%.3f inH2O',inH2O(x)),'FontSize',36)
catch
    plot(inH2O, 'o-')
    title(sprintf('%.3f inH2O',inH2O(x)),'FontSize',36)
end
pause(0.001)


end
