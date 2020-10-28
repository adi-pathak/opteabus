function [P_EL_MIN, DATA] = Simulation_energy_consumption(dc, vehicle, architecture)

%% Calculation Resistance
v = dc.speed';              %determine velocity vector, v in m/s
a = [0 diff(v)];            %determine acceleration, a in m/s²
alpha = zeros(size(v));     %elevation angle, in this case 0° for all steps

[V]=Calc_Resistance(v,a,alpha, vehicle);

%% Calculation of Power requirement - different for each powertrain architecture

 Simulation=str2func(strcat('Simulation_', architecture));
 [P_EL_MIN, DATA]=Simulation(V, v, a, vehicle);


end

