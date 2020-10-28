classdef gearbox
    properties
        cost
        type
        numberofgears
        mass
        transmissionefficiency=.98;
        ratio=15;
    end
    
    methods
        
        
        function obj=gearbox(numberofgears)
            obj.cost=300;
            obj.mass=0;
        end
       
    end
    
end