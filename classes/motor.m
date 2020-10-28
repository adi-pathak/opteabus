classdef motor
    properties
        type='PSM'
        power %nominal/peak?
        cost
        mass
        inertia
        efficiencymap
        diameter
        length
    end
    
    methods
        
        function obj=motor(motorpower)
            obj.power=motorpower;
            %obj.path='C:\Users\aditya.pathak\Documents\Packaging Code\Components\PSM_i3_75_4800_12000_360.mat';
            
            basemap=load(fullfile(pwd,'Inputs',(filesep),'basemap.mat'));
            obj.efficiencymap=obj.scalemap(motorpower,basemap);
            T=max(obj.efficiencymap.T_max);
            Tn= motorpower*1000/(obj.efficiencymap.n_n*2*pi/60);
            obj.cost=(9.2*motorpower)+250; % rated motor power
            obj.mass=((13.847*log(motorpower))-13.003+(10.979*log(Tn))-17.908)*0.5;
            obj.inertia=0.0002*Tn-0.0029;
            motorvolume=5.7+0.1*Tn*100^3; % felgenhauer
            maxmotorspeed=obj.efficiencymap.n_max;
            lengthtodiameterratio=0.26+7.02E-5*maxmotorspeed; %revspermin
            obj.diameter=round((motorvolume*4/(pi*lengthtodiameterratio))^(1/3));
            obj.length=round(obj.diameter*lengthtodiameterratio);
          obj.diameter=250;
        end
       
        function efficiencymap=scalemap(obj,power,basemap)
            %map - base motor
            map=basemap.motor;
            scalefactor=power/map.power; %basemotor -100kW
            efficiencymap=map;
            efficiencymap.eff_n_axis= map.eff_n_axis;
            efficiencymap.eff_T_axis=map.eff_T_axis.*scalefactor;
            efficiencymap.T_max=map.T_max.*scalefactor;
            efficiencymap.power=power;
        end
         function plotmotor(obj,position,handle)
             diameter=obj.diameter;
             length=obj.length;
                  vo1=[0.5*diameter*cos(linspace(0,pi,40));-0.5*length+0*(linspace(0,pi,40));0.5*diameter*sin(linspace(0,pi,40))];
        vo2=[0.5*diameter*cos(linspace(0,pi,40));0.5*length+0*(linspace(0,pi,40));0.5*diameter*sin(linspace(0,pi,40))];
        vo3=[0.5*diameter*cos(linspace(pi,2*pi,40));-0.5*length+0*(linspace(0,pi,40));0.5*diameter*sin(linspace(pi,2*pi,40))];
        vo4=[0.5*diameter*cos(linspace(pi,2*pi,40));0.5*length+0*(linspace(0,pi,40));0.5*diameter*sin(linspace(pi,2*pi,40))];
        V1=obj.translate([vo1,flip(vo2,2)],position);
        V2=obj.translate([vo3,flip(vo4,2)],position);
        V3=obj.translate([vo1,vo3],position);
        V4=obj.translate([vo2,vo4],position);
       
        colr=[0 0 1];
        patch(handle,'Faces',[1:80],'Vertices',V1','FaceColor', colr)
        patch(handle,'Faces',[1:80],'Vertices',V2','FaceColor', colr)
        patch(handle,'Faces',[1:80],'Vertices',V3','FaceColor', colr)
        patch(handle,'Faces',[1:80],'Vertices',V4','FaceColor', colr)
        
         end
         function vertices=translate(obj,vertices,position)
            T= [1 0 0  position(1);0 1 0  position(2);0 0 1  position(3);0 0 0 1] ;
            vertices(4,:)=1;
            vertices=(T*vertices);
            vertices=vertices(1:3,:);
        end
       
         
    end
    
    
end

