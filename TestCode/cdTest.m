N = 600;

cd = nan(N,1);

for x = 1:N

[ rawData ] = getLabJackRawData( {'AIN6','AIN7','AIN12'}, 500, 1, 5 );
clc

meanData = mean(rawData,1);
AIN6 = meanData(1);
AIN7 = meanData(2);
AIN12= meanData(3);



T0 = AIN6*100;
TC = (T0-32)/1.8;

q0Pa = ((AIN7-0.136)*60.3547)-4.10261;

inHg = 29.79;
Pa = inHg * 3386.39;
rho = Pa/287.058/(TC+273);

rho0 = rho*0.00194032;
dptinH2O = 0.06928*AIN12 - 0.06843;


q0Psi = q0Pa * 0.000145038;
q0inH2O = q0Pa * 0.00401865;

dptPsi = dptinH2O*0.0360912;


Cd_measured = fcnWakeRake( dptinH2O, T0, q0inH2O, rho0 );
cd(x) = Cd_measured(end);


try
    plot(idx, cd(idx), 'o-')
    title(sprintf('%.3f inH2O, Cd %.3f ',dptinH2O, cd(x)),'FontSize',48)
catch
    plot(cd, 'o-')
    title(sprintf('%.3f inH2O, Cd %.3f ',dptinH2O ,cd(x)),'FontSize',48)
end

ylim([-1 1])
grid minor
pause(0.001)


end










