classdef inverter
    properties
      cost
      mass
      efficiency
    end
    methods
        function obj=inverter(motorpower)
           obj.mass= (0.0886*motorpower)+3.357; 
           obj.cost=(7.5*motorpower)+145;
        end
    end
    
end