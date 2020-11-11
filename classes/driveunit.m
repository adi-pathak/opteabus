classdef driveunit%<motor & inverter & gearbox
    properties
       driveunitcost
       driveunitpower
       motor
       inverter
       transmission
    end
     
    methods
        function obj=driveunit(motorpower,gearratio)
           obj.driveunitpower=motorpower;
           obj.motor=motor(motorpower);
           obj.inverter=inverter(motorpower);
           obj.transmission=gearbox(gearratio);
           obj.driveunitcost=obj.motor.cost+ obj.inverter.cost+obj.transmission.cost;
        end

     end
    
end