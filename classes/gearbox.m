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
        
        
        function obj=gearbox(ratio)
            obj.cost=300;
            obj.mass=0;
            obj.ratio=ratio;
        end
       
    end
    
end