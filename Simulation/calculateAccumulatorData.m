function [Accumulator] = calculateAccumulatorData(Accumulator, P_el, P_M, P_Bh, t, i)
    
    %% Initialisation of all variable
    load('Emraxefficiencydata_interpolated.mat');
    load('CorrectedDischargeInterpolated.mat');
    load('RandomizedCellData.mat');
    load('CellparametersVoltageInterpolation.mat');

    Accumulator.V_i(i) = sum(Accumulator.Voltage_Cellpack(:,i));
    
    % Battery Currents (Akkustroeme)
    Accumulator.A_accu_cell(i) = P_el(i) / Accumulator.V_i(i) / setup.nZellen_Parallel;  
    
    Accumulator.Current_Cellpack_Pointer(i) = P_M(i) / Accumulator.V_i(i) * 10 ; %Strombelastung eines %er Parrallel Paketes in 0,1A parameter fuer die berechnung der korrigierten belastung mit hoehren verlusten durch hoehere zellstroeme
    if Accumulator.Current_Cellpack_Pointer(i) <= 1
        Accumulator.Current_Cellpack_Pointer(i)=1;
    end
    
    if Accumulator.Current_Cellpack_Pointer(i) >= 1500 %begrenzen des max Zellstromes auf 30A pro Zelle im 5er parralelverbund also 150A
        Accumulator.Current_Cellpack_Pointer(i)=1500;
    end
    
    Accumulator.VirtualCurrent_Cellpack(i) = CorrectedDischargeInterpolated(1,round(Accumulator.Current_Cellpack_Pointer(i))); %Berechnung der Virtuell hoeheren zellstr�me basierend auf den h�heren verlsuten durch h�here Str�me
    
    Accumulator.Energy_Cellpack(i) = (Accumulator.VirtualCurrent_Cellpack(i)*(t(i+1)-t(i))) - ((P_Bh(i) / Accumulator.V_i(i)) * (t(i+1)-t(i))) ; %Energieverbrauch in As f�r ein 5erpacket an akkuzellen -> Akkustrom zum zeitpunkt i mal Zeitdifferenz zwischen i und i+1
    Accumulator.Energy_Cellpack_Total(i+1) = Accumulator.Energy_Cellpack_Total(i) + Accumulator.Energy_Cellpack(i); % �ber Endurance Run Integrierte Energieverbrauch in As f�r ein 5erpacket an akkuzellen
    
    Accumulator.Capacity_Cellpack(1:131,i+1) =  Accumulator.Capacity_Cellpack(1:131,i) - Accumulator.Energy_Cellpack(i); 
    
    Accumulator.SOC_Cellpack(1:131,i+1) = Accumulator.Capacity_Cellpack(1:131,i) ./ Accumulator.Capacity_Cellpack(1:131,1); %Berechnung des SOC f�r den n�chsten tick basierend auf der aktuellen cellcapacity und der im n�chsten tick
    
    Accumulator.SOC_Pointer(1:131,i+1) = round(Accumulator.SOC_Cellpack(1:131,i+1)*1000);
    Accumulator.Current_Cellpack_Pointer_Voltage(1,i+1) = round(Accumulator.Current_Cellpack_Pointer(i)/5);
    
    if Accumulator.Current_Cellpack_Pointer_Voltage(i) <= 3
        Accumulator.Current_Cellpack_Pointer_Voltage(i)=3;
    end
    
    if size(Track,1) < Accumulator.SOC_Pointer(1:131,i+1)
        Accumulator.SOC_Pointer(1:131,i+1) = size(Track,1); 
    end
    
    if size(Track,1) < Accumulator.Current_Cellpack_Pointer_Voltage(1,i)
        Accumulator.Current_Cellpack_Pointer_Voltage(1,i) = size(Track,1);
    end
    
    try % Catch Exception when Accumulator is empty 
        Accumulator.Voltage_Cellpack(1:131,i+1) = Voltage_inter(Accumulator.Current_Cellpack_Pointer_Voltage(1,i),Accumulator.SOC_Pointer(1,i+1));
    catch
    
    end
end