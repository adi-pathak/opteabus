classdef powertrain %< driveunit & battery & wheel
    properties
        
        topology
        mass
        cost
        totalpower
        powersplit
        FA1
        FA2
        RA1
        RA2
        emissions
    end
    properties (Constant)
        structure=[0 0 1 1 2
            1 2 1 2 2];
    end
    methods
        function obj=powertrain(topology,totalpower,powersplit,numberofgears,vehicle)
            obj.topology=topology;
            obj.totalpower=totalpower;
            obj.powersplit=powersplit;
            
            if  obj.structure(1,obj.topology)>0
                if  obj.structure(1,obj.topology)<2
                    motorpower=obj.totalpower*obj.powersplit;
                    obj.FA1=driveunit(motorpower,numberofgears);
                    obj.FA2=0;
                    FAcost=obj.FA1.driveunitcost;
                    FAmass=obj.FA1.motor.mass+obj.FA1.inverter.mass+obj.FA1.transmission.mass;
                else
                    motorpower=obj.totalpower*obj.powersplit/2;
                    obj.FA1=driveunit(motorpower,numberofgears);
                    obj.FA2=obj.FA1;
                    FAcost=2*(obj.FA1.driveunitcost);
                    FAmass=2*(obj.FA1.motor.mass+obj.FA1.inverter.mass+obj.FA1.transmission.mass);
                end
            else
                obj.FA1=0; obj.FA2=0; FAcost=0; FAmass=0;
            end
            if  obj.structure(2,obj.topology)<2
                motorpower=obj.totalpower*(1-obj.powersplit);
                obj.RA1=driveunit(motorpower,numberofgears);
                obj.RA2=0;
                RAcost=obj.RA1.driveunitcost;
                RAmass=obj.RA1.motor.mass+obj.RA1.inverter.mass+obj.RA1.transmission.mass;
            else
                motorpower=(obj.totalpower*(1-obj.powersplit))/2;
                obj.RA1=driveunit(motorpower,numberofgears);
                obj.RA2=obj.RA1;
                RAcost=2*(obj.RA1.driveunitcost);
                RAmass=2*(obj.RA1.motor.mass+obj.RA1.inverter.mass+obj.RA1.transmission.mass);
            end
            obj.cost=RAcost+FAcost;
            obj.mass=RAmass+FAmass;
            obj.emissions.CO2=obj.mass*18.664888;
            obj.emissions.DCB=obj.mass*84.018099;
            obj.emissions.PM10=obj.mass*0.066884;
            
        end
        function obj=updatepowertrain(obj)
            
            obj=updatebattery(obj);
            if  obj.structure(1,obj.topology)>0
                if  obj.structure(1,obj.topology)<2
                    obj.FA1.motorpower=obj.totalpower*obj.powersplit;
                    obj.FA2=0;
                else
                    obj.FA1.motorpower=obj.totalpower*obj.powersplit/2;
                    obj.FA2=obj.FA1;
                end
            else
                obj.FA1=0; obj.FA2=0;
            end
            if  obj.structure(2,obj.topology)<2
                obj.RA1.motorpower=obj.totalpower*(1-obj.powersplit);
                obj.RA2=0;
            else
                obj.RA1.motorpower=(obj.totalpower*(1-obj.powersplit))/2;
                obj.RA2=obj.RA1;
            end
            
            obj=updatedriveunit(obj);
            
        end
    end
end
