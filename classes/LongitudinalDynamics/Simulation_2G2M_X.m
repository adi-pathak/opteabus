%% Simulation vehicle

function [P_EL_MIN, DATA_R]=Simulation_GM_X(V, v, a, vehicle)

%SIMULATION_GM_X calculates the required energy consumption of a vehicle
%with the powertrain GM_X (Gearbox and Motor at one axle, other axle
%powerless)

% Author:   Alexander Koch, FTM, TUM
% Date:     12/10/2018

%% Inputs
%   V:             Matrix with required Torque of vehicle in first line and revolutional speed of wheels in second line
%   vehicle:    struct with the vehicle parameters


%% Outputs
%   P_EL_MIN:      Vektor with minimal electrical power requirement for every time step
%   DATA_R:         Struct with Data about the motor torque, revolutional speed and Gear

 %V_Axle_Front = [(vehicle.MOTOR{1}.eff_T_axis(end)*max(vehicle.GEARBOX{1}.gear_ratio)+(vehicle.MOTOR{2}.eff_T_axis(end)*max(vehicle.GEARBOX{2}.gear_ratio)));(vehicle.MOTOR{1}.eff_T_axis(1)*min(vehicle.GEARBOX{1}.gear_ratio)+(vehicle.MOTOR{2}.eff_T_axis(1)*min(vehicle.GEARBOX{2}.gear_ratio)))];

%% Powertrain

 V_Axle_Front_Left_Right=V/2;
 
 [V_Gearbox_Front_Left]= Partition_Gearbox(vehicle.GEARBOX{1,1}, V_Axle_Front_Left_Right);
 
 [V_Gearbox_Front_Right]= Partition_Gearbox(vehicle.GEARBOX{1,2}, V_Axle_Front_Left_Right);

 [P_EL_Motor_Front_Left, DATA_Front_Left] = Determination_Power (V_Gearbox_Front_Left, vehicle.MOTOR{1,1});
 
 [P_EL_Motor_Front_Right, DATA_Front_Right] = Determination_Power (V_Gearbox_Front_Right, vehicle.MOTOR{1,2});
 
 [P_EL_Gearbox_Front_Left, DATA_Front_Left] = Determination_Gear (P_EL_Motor_Front_Left, vehicle.GEARBOX{1,1}, DATA_Front_Left);

 [P_EL_Gearbox_Front_Right, DATA_Front_Right] = Determination_Gear (P_EL_Motor_Front_Right, vehicle.GEARBOX{1,2}, DATA_Front_Right);

 P_EL_Gearbox_Front = P_EL_Gearbox_Front_Left + P_EL_Gearbox_Front_Right;
 
 DATA=[DATA_Front_Left;DATA_Front_Right];
 
%  P_EL_Gearbox_Front=P_EL_Gearbox_Front_Left(1,:)+P_EL_Gearbox_Front_Right(1,:);

%  [P_EL_MIN, DATA] = Determination_Distribution (P_EL_Gearbox_Front_Left, P_EL_Gearbox_Front_Right, [DATA_Front_Left;DATA_Front_Right]);
 
 P_EL_MIN=P_EL_Gearbox_Front+vehicle.auxiliary;

%% Battery (Simple Eff)

[P_EL_MIN]=BatteryEfficiency(P_EL_MIN);

%% Rearangement of DATA

% if v and a == 0 No Torque is required
DATA(1,(v==0&a==0))=0;
DATA(4,(v==0&a==0))=0;

DATA_R.Tn=zeros(4,length(DATA(1,:)));
 
 DATA_R.Tn(1,:)=DATA(1,:);
 DATA_R.Tn(2,:)=DATA(2,:);
 DATA_R.Tn(3,:)=DATA(4,:);
 DATA_R.Tn(4,:)=DATA(5,:);
 DATA_R.Gear_1=DATA(3,:);
 DATA_R.Gear_2=DATA(6,:);

end
 
 
