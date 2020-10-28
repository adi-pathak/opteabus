function [Econs,range,vmax,amax,grade,throughput] = EnergyConsumption( Param, driving_cycle )
whichSim='E';
v_max_sim = 100;    %maximum simulated velocity in km/h. When vehicle reaches this velocity, the simulation ends.
t_sim = 100;        %simulated time in s. Simulation ends after this time.
visualize_ECS = 0;

%% Environment
vehicle.environment.roh_L   = 1.184;                                       % in kg/m^3
vehicle.environment.g       = 9.81;                                        % in N/kg

%% Vehicle parameters
vehicle.m                 = Param.Grossvehiclemass ;%Param.mass.adaptedmass;                             % in kg - Use Mass_calculation() regression functions
% if Param.chassis.numaxles == 3
%     vehicle.rotatingMass  = vehicle.m*0.11;                              % in kg - estimation of inertia from rotating parts
% else
%     vehicle.rotatingMass  = vehicle.m*0.07;
% end
vehicle.rotatingMass = 0;
vehicle.c_w             = 0.70;                                            % in - Use regression functions in Excel sheet
vehicle.A_front         = Param.Body.height/1000 *...
                          Param.Body.width/1000;                          % in m^2 -  Use: "0.85*Width*Height" to estimate it
vehicle.r_tyre          = 0.316;                                           % in m
vehicle.f_R             = 0.008;                                           % in - Estimation
vehicle.battery_cap     = Param.Battery.capacity;                          % in kWh - Use battery design description
vehicle.auxiliary       = Param.Passengercapacity * 1000 * 0.1575;       % in W - Estimation
vehicle.GEARBOX         = cell(1,4);                                       % initialize empty gearbox struct
vehicle.MOTOR           = cell(1,4);                                       % initialize empty motor struct



switch Param.Powertrain.topology
   case 1
        vehicle.MOTOR{1} = Param.Powertrain.RA1.motor.efficiencymap;
        vehicle.GEARBOX{1}.gear_ratio       =...
            [7 20]; %[Param.motor.FAGearRatio1 Param.motor.FAGearRatio2];
        vehicle.GEARBOX{1}.eff              = [1 1];
        architecture = 'GM_X';
    case 2
        vehicle.GEARBOX{1}.gear_ratio       =...
            [Param.motor.FAGearRatio1 Param.motor.FAGearRatio2];
        vehicle.GEARBOX{1}.eff              = [1 1];
        vehicle.GEARBOX{3}.gear_ratio       =...
            [Param.motor.RAGearRatio1 Param.motor.RAGearRatio2];
        vehicle.GEARBOX{3}.eff              = [1 1];
        vehicle.MOTOR{1} = motor;
        vehicle.MOTOR{1} = motor;
        vehicle.MOTOR{3} = motor;
        vehicle.MOTOR{3} = motor;
        architecture = 'GM_GM';
    case 3
        vehicle.GEARBOX{1}.gear_ratio       =...
            [Param.motor.FAGearRatio1 Param.motor.FAGearRatio2];
        vehicle.GEARBOX{1}.eff              = [1 1];
        vehicle.GEARBOX{2}.gear_ratio       =...
            [Param.motor.RAGearRatio1 Param.motor.RAGearRatio2];
        vehicle.GEARBOX{2}.eff              = [1 1];
        vehicle.MOTOR{1} = motor;
        vehicle.MOTOR{1} = motor;
        vehicle.MOTOR{2} = motor;
        vehicle.MOTOR{2} = motor;
        architecture = '2G2M_X';
    case 4
        vehicle.GEARBOX{1}.gear_ratio       =...
            [Param.motor.FAGearRatio1 Param.motor.FAGearRatio2];
        vehicle.GEARBOX{1}.eff              = [1 1];
        vehicle.GEARBOX{2}.gear_ratio       =...
            [Param.motor.FAGearRatio1 Param.motor.FAGearRatio2];
        vehicle.GEARBOX{2}.eff              = [1 1];
        vehicle.GEARBOX{3}.gear_ratio       =...
            [Param.motor.RAGearRatio1 Param.motor.RAGearRatio2];
        vehicle.GEARBOX{3}.eff              = [1 1];
        vehicle.MOTOR{1} = motor;
        vehicle.MOTOR{1} = motor;
        vehicle.MOTOR{2} = motor;
        vehicle.MOTOR{2} = motor;
        vehicle.MOTOR{3} = motor;
        vehicle.MOTOR{3} = motor;
        architecture = '2G2M_GM';
    case 5
        vehicle.GEARBOX{1}.gear_ratio       =...
            [Param.motor.FAGearRatio1 Param.motor.FAGearRatio2];
        vehicle.GEARBOX{1}.eff              = [1 1];
        vehicle.GEARBOX{2}.gear_ratio       =...
            [Param.motor.FAGearRatio1 Param.motor.FAGearRatio2];
        vehicle.GEARBOX{2}.eff              = [1 1];
        vehicle.GEARBOX{4}.gear_ratio       =...
            [Param.motor.RAGearRatio1 Param.motor.RAGearRatio2];
        vehicle.GEARBOX{4}.eff              = [1 1];
        vehicle.GEARBOX{3}.gear_ratio       =...
            [Param.motor.RAGearRatio1 Param.motor.RAGearRatio2];
        vehicle.GEARBOX{3}.eff              = [1 1];
        vehicle.MOTOR{1} = motor;
        vehicle.MOTOR{1} = motor;
        vehicle.MOTOR{4} = motor;
        vehicle.MOTOR{4} = motor;
        vehicle.MOTOR{3} = motor;
        vehicle.MOTOR{3} = motor;
        vehicle.MOTOR{2} = motor;
        vehicle.MOTOR{2} = motor;
        architecture = '2G2M_2G2M';
end
dc.speed=driving_cycle(:,2);
dc.time=driving_cycle(:,1);
[P_el_min, data_R] =...
    Simulation_energy_consumption(dc, vehicle, architecture);
 [vmax,amax,grade] =...
     Simulation_acceleration(Param.Grossvehiclemass,vehicle, v_max_sim, t_sim);
Results.E = cumsum(P_el_min);                                              % Required Energy in Ws
Results.S = cumsum(dc.speed);                                   % Distance in m
Econs = Results.E(end)/(1000*3600)*1/(Results.S(end)/1000); %kwh/km 
range = (Param.Battery.capacity/Econs);% Range of vehicle in km
% % Adaptation to Teichert
 Pregen(P_el_min<0) = P_el_min(P_el_min<0)/1000;
 Pdriving(P_el_min>=0) = P_el_min(P_el_min>=0)/1000;
 Echarge=-trapz(Pregen)/3600;
 Edriving=trapz(Pdriving)/3600;
 Etotal=(trapz(P_el_min/1000)/3600);
b_h_regen = (Echarge/driving_cycle(end,1))*3600;
b_h_discharge = (Edriving/driving_cycle(end,1))*3600;
b_h = (Etotal/driving_cycle(end,1))*3600;
b_h_regen = (Echarge/(Results.S(end)/1000));
b_h_discharge = (Edriving/(Results.S(end)/1000));
b_h = (Etotal/(Results.S(end)/1000));

throughput=b_h+b_h_regen+b_h_discharge;
% 
% % Energy consumption in kWh/100km;
% if isnan(Param.consumption.ECons)
%     Param.consumption.ECons=10000000;
% elseif Param.consumption.ECons > 1000000
%     Param.consumption.ECons=10000000;
% elseif Param.consumption.ECons < 0
%     Param.ECons=10000000;
% end
% 
% Param.consumption.range = (vehicle.battery_cap/Param.consumption.ECons)*100;% Range of vehicle in km
% Param.motor.architecture = architecture;
% dynamicplot_motor(vehicle, driving_cycle, data_R.Tn);
end


