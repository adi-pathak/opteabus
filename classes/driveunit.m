classdef driveunit%<motor & inverter & gearbox
    properties
       driveunitcost
       driveunitpower
       motor
       inverter
       transmission
%         powermap
%         M1T
%         M2T
%         Trc
%         n1
%         n2
    end
     
    methods
        function obj=driveunit(motorpower,numberofgears)
           obj.driveunitpower=motorpower;
           obj.motor=motor(motorpower);
           obj.inverter=inverter(motorpower);
           obj.transmission=gearbox(numberofgears);
           obj.driveunitcost=obj.motor.cost+ obj.inverter.cost+obj.transmission.cost;
        end

     end
    
end