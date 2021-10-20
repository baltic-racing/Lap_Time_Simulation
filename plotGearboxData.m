function plotGearboxData(gearboxData, n_shift, t_shift, n_max)

clf

TIRparam = loadTIR('C19_CONTINENTAL_FORMULASTUDENT_205_470_R13_65kPa.tir');

R0 = TIRparam.UNLOADED_RADIUS;

resultData(1,1) = 0;
resultData(1,2) = 0;

for i = 2:size(gearboxData,1)+1
    speed(i) = (3.6 * n_max * pi * R0) / (30 * gearboxData(i-1,2));  % Max Speed km/h
    
    resultData(end+1,1) = speed(i);  
    resultData(end,2) = n_shift;                             % shift rpm
    
    if i <= size(gearboxData,1) 
        resultData(end+1,1) = speed(i); 
        resultData(end,2) = (speed(i) * 30 * gearboxData(i,2))/(3.6 * pi * R0);
    end

end

plot(resultData(:,1),resultData(:,2));

xlabel('car speed [km/h]')
ylabel('engine rpm [1/min]')

annotation('textbox', [0.2, 0.8, 0.1, 0.1],...
                        'String', {['Theoretical topspeed = ' num2str(resultData(end,1)) ' km/h']})
    
end