clc 
clear
% Generates efficiency maps for a specific power range
    
motor.type = 'PSM';         % PSM or ASM
motor.power=90;              % in kW
motor.n_n=1500;         	% in 1/min
motor.n_max=7500;          % in 1/min
motor.U_n=400;              % Motor voltage

motor = GENERATE_motor_efficiency_map(motor);
MData(90)=motor;            % Stores Data in MData
clear motor
%disp('All motors stored');
