


tare = [0.00259392492305487];


[ rawData ] = getLabJackRawData( {'AIN7','AIN8','AIN12'}, 2000, 1, 5 );


AIN7 = mean(rawData(:,1));
AIN8 = mean(rawData(:,2));
AIN12 = mean(rawData(:,3));

Vinf =sqrt((((AIN7-tare)*58.3)+11.9)*2/(100270.97/(287.14*(((AIN8*100)+459.67)*5/9))));

inH2O = 0.06928*AIN12 - 0.06843;



