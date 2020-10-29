classdef interior
    
    properties
        interiorlayout; % urban 1, coach, urban 2
        standingspace;
        seatwidth=1;
        seatpitch=1;
        gangwaywidth=0.5;
        wheelchairzones=0;
        seatingratio
        costs
        mass
        standingpassengers
        passengercapacity
        numberseats
        loadedmass
        numberofdoors
        floorheight=370;
        interiorlength
        interiorwidth
        frontcrashlength=450;
        rearcrashlength=450;
    end
    properties (Access=private)
        passengermass=65;
        w_door = 1200;
        d_wheelhouse = 767;
        seatlength = 450;
        seatheight=450;
        backrestheight=500;
        l_seat = 500;
        t_seat = 40;
        w_seat_gap = 30;
        p_shoulder_width = 530;
        p_body_depth = 450;
        w_tire = 215;
        d_tire = 767;
        d_spring = 300;
        t_body = 5;
        curvature=50;
        doorwidth=1250;
        wheelchairzonelength=1250+250;
    end
    methods
        
        function obj=interior(layout,standingspace,seatwidth,seatpitch,wheelchairzones,Vehicle)
            if nargin==6
                obj.interiorlayout=layout;
                obj.standingspace=standingspace;
                obj.seatwidth=seatwidth;
                obj.seatpitch=seatpitch;
                obj.wheelchairzones=1; %wheelchairzones;
                if Vehicle.Body.length>8500
                    obj.numberofdoors=2;
                else
                    obj.numberofdoors=1;
                end
                         
            end
            [obj]=updateinterior(obj,Vehicle);
            obj.loadedmass=obj.passengermass*obj.passengercapacity;
        end
        function [obj]=updateinterior(obj,Vehicle)
%                      
            if obj.interiorlayout==1 % coach layout
                vehiclelength=Vehicle.Body.length;
                interiorlength=vehiclelength-300-300; % reduce interior length for front and rear zones
                obj.interiorlength=interiorlength;
                seatpitch=obj.seatpitch;
                seatwidth=obj.seatwidth;
                seatlength=obj.seatlength;
                vehiclewidth=Vehicle.Body.width;
                interiorwidth=vehiclewidth-200;
                obj.interiorwidth=interiorwidth;
                seatgap=obj.w_seat_gap;
                tgap=obj.t_body;
                doorwidth=obj.doorwidth;
                wheelchairs=0;
                doors=obj.numberofdoors;
                passengerdensity=1/4;% 4p/m^2
                if doors==1
                    seatarea1_length=(interiorwidth-tgap); % calculate numseat at the front, back
                    numseats1=floor(seatarea1_length/(seatwidth+2*seatgap));
                    %numseats2=numseats1;
                    legroom=seatpitch-seatlength;
                   
                   
                    if wheelchairs==0
                          seatarea2_length=(interiorlength-2*(tgap))-(seatpitch); % calculate numseats at side
                          numseats2=2*floor(seatarea2_length/(seatpitch));
                   
                    seatarea3_length=(interiorlength-2*(tgap)-doorwidth)/2-(seatpitch); % calculate numseats behind door
                    numseats3=2*floor(seatarea3_length/(seatpitch));
                    seatarea4_length=(interiorlength-2*(tgap)-doorwidth)/2;%
                    numseats4=2*floor(seatarea4_length/(seatpitch));
                    numseats=numseats1+numseats2+numseats3+numseats4;
                    % calculate standing space
                    aislewidth=interiorwidth-4*(seatwidth+seatgap);
                    standingarea= ((seatarea2_length-150) *aislewidth)/(1000^2);
                    standingpassengers=floor(standingarea/passengerdensity);
                    standingarea2=doorwidth*((interiorwidth-aislewidth-2*(seatwidth+seatgap)))/(1000^2);
                    standingpassengers2=floor(standingarea2/passengerdensity);
                    standingpassengers=standingpassengers+standingpassengers2;
                    obj.numberseats=numseats;
                    obj.passengercapacity= standingpassengers+numseats;
                    obj.standingpassengers= standingpassengers;
                    obj.seatingratio=obj.numberseats./obj.passengercapacity;
                    obj.costs=obj.numberseats*265;
                    obj.mass=obj.numberseats*15;
                     else
                        %
                    end
                else
                       wheeldiameter=900; % initial estimate to calculate overhangs
                       wheelbase=Vehicle.Body.wheelbase;
                  frontoverhang=400+obj.doorwidth+1.5*wheeldiameter/2;
                  rearoverhang=vehiclelength-frontoverhang-wheelbase;
                    seatarea1_length=(interiorwidth-tgap); % calculate numseat at the front, back
                    numseats1=floor(seatarea1_length/(seatwidth+2*seatgap));
                    seatarea2_length=wheelbase/2+rearoverhang-obj.doorwidth/2;
                     seatarea2_length= seatarea2_length-(seatpitch-seatlength); % number of seats in the rear left
                     numseats2=floor(seatarea2_length/(seatpitch+seatlength/2))*2; % 2 seats/row
                     
                     numseats3=numseats2; % on thre r side
                     seatarea4_length=wheelbase/2+frontoverhang-obj.doorwidth/2-...
                         (obj.wheelchairzonelength-obj.doorwidth);
                     numseats4=floor(seatarea4_length/(seatpitch+seatlength/2))*2; % front seats
                     if obj.wheelchairzones<2
                    seatarea5_length=wheelbase/2+wheeldiameter/2;
                    numseats5=floor(seatarea5_length/(seatpitch+seatlength/2))*2; % front seats
                    numseats=numseats1+numseats2+numseats3+numseats4+numseats5;
                     end
                     % calculate standing passengers
                   passengerdensity=1/4;% 4p/m^2
                   aislelength=interiorlength-(seatpitch+seatlength);
                   aislewidth=(interiorwidth-4*(seatwidth+seatgap));
                    standingarea= aislelength *aislewidth/(1000^2);
                    standingpassengers=floor(standingarea/passengerdensity); %aisle
                    standingarea2=doorwidth*(interiorwidth-aislewidth)/(1000^2); % rear door standing
                    standingpassengers2=floor(standingarea2/passengerdensity);
                    standingarea3=doorwidth*(interiorwidth-aislewidth-2*(seatwidth+seatgap))/(1000^2); % rear door standing
                    standingpassengers3=floor(standingarea3/passengerdensity);
                    standingpassengers=standingpassengers+standingpassengers2+standingpassengers3;
                    obj.passengercapacity= standingpassengers+numseats;
                    obj.numberseats=numseats;
                      obj.standingpassengers= standingpassengers;
                    obj.seatingratio=obj.numberseats./obj.passengercapacity;
                    obj.costs=obj.numberseats*264;
                    obj.mass=obj.numberseats*15;
                end    
              
               
           obj.seatingratio=obj.numberseats./obj.passengercapacity;
            obj.costs=obj.numberseats*264; %derive better interior cost model
            obj.mass=obj.numberseats*15; % weight of seats - derive better weight model
           
                
                
            elseif obj.interiorlayout==2 % Urban Layout 1
                
                vehiclelength=Vehicle.Body.length;
                interiorlength=vehiclelength-obj.frontcrashlength-obj.rearcrashlength; % reduce interior length for front and rear zones
                obj.interiorlength=interiorlength;
             
                seatpitch=obj.seatpitch;
                seatwidth=obj.seatwidth;
                seatlength=obj.seatlength;
                vehiclewidth=Vehicle.Body.width;
                interiorwidth=vehiclewidth-200;
                obj.interiorwidth=interiorwidth;
                seatgap=obj.w_seat_gap;
                tgap=obj.t_body;
                doorwidth=obj.doorwidth;
                wheelchairs=0;
                doors=1;
                passengerdensity=1/4;% 4p/m^2
                if doors==1
                    seatarea1_length=(interiorwidth-tgap); % calculate numseat at the front, back
                    numseats1=floor(seatarea1_length/(seatwidth+2*seatgap));
                    numseats2=numseats1;
                    legroom=250+150;
                    seatarea3_length=(interiorlength-2*(tgap))-2*(seatlength+legroom); % calculate numseats at side
                    numseats3=floor(seatarea3_length/(seatwidth+2*seatgap));
                    seatarea4_length=(seatarea3_length-(doorwidth+200))/2;
                    numseats4=floor(seatarea4_length/(seatwidth+2*seatgap))*(seatarea4_length>0);
                    if wheelchairs==0
                        numseats5=numseats4;
                    else
                        %
                    end
                    numseats=numseats1+numseats2+numseats3+numseats4+numseats5;
                    % calculate standing space
                    standingarea= (seatarea3_length) *(seatarea1_length-2*(seatlength+legroom))/(1000^2);
                    standingpassengers=floor(standingarea/passengerdensity);
                    standingarea2=doorwidth*(seatlength)/(1000^2);
                    standingpassengers2=floor(standingarea2/passengerdensity);
                    standingpassengers=standingpassengers+standingpassengers2;
                    obj.numberseats=numseats;
                    obj.passengercapacity= standingpassengers+numseats;
                    obj.standingpassengers= standingpassengers;
                    obj.seatingratio=obj.numberseats./obj.passengercapacity;
                    obj.costs=obj.numberseats*264;
                    obj.mass=obj.numberseats*15;
                else
                    
                end    
              
                end
           obj.seatingratio=obj.numberseats./obj.passengercapacity;
            obj.costs=obj.numberseats*264; %derive better interior cost model
            obj.mass=obj.numberseats*15; % weight of seats - derive better weight model
        end
      
        function plotinterior(obj,Vehicle,handle)
             if obj.interiorlayout==1 % coach layout
                 plotinterior1(obj,Vehicle,handle);
             elseif obj.interiorlayout==2
                 plotinterior2(obj,Vehicle,handle);
                 
             end
        end
        function plotinterior1(obj,Vehicle,handle) %coach layout
                vehiclelength=Vehicle.Body.length;
                 interiorlength=vehiclelength-obj.frontcrashlength-obj.rearcrashlength; % reduce interior length for front and rear zones
                obj.interiorlength=interiorlength;
                seatpitch=obj.seatpitch;
                seatwidth=obj.seatwidth;
                seatlength=obj.seatlength;
                vehiclewidth=Vehicle.Body.width;
                interiorwidth=vehiclewidth-200;
                obj.interiorwidth=interiorwidth;
                seatgap=obj.w_seat_gap;
                tgap=obj.t_body;
                doorwidth=obj.doorwidth;
                wheelchairs=0;
                doors=obj.numberofdoors;
                passengerdensity=1/4;% 4p/m^2
                if doors==1
                    seatarea1_length=(interiorwidth-tgap); % calculate numseat at the front, back
                    numseats1=floor(seatarea1_length/(seatwidth+2*seatgap));
                    legroom=seatpitch-seatlength;
                    xback=-0.5*interiorlength+50+seatlength/2;
                    y=linspace(-seatarea1_length/2+50+seatwidth/2,seatarea1_length/2-50-seatwidth/2,numseats1);
                    z=obj.seatheight+obj.floorheight;
                    orientback=[0 0 0];
                    for i=1:numseats1
                        position=[xback y(i) z];
                        obj.plotseats(position,orientback,obj.seatwidth,obj.seatlength,obj.backrestheight,handle)
                    end
                    
                    if wheelchairs==0
                        seatarea2_length=(interiorlength-2*(tgap))-(seatpitch); % calculate numseats at side
                        numseats2=2*(floor(seatarea2_length/(seatpitch)));
                         x=ones(1,numseats2/2)*seatpitch; % seat centre pattern
                         x(1)=xback+seatpitch;
                         x=cumsum(x); % linear pattern
                         y=-(interiorwidth/2-seatwidth/2);
                         y1=y+seatwidth+seatgap;
                         z=obj.seatheight+obj.floorheight;
                         for i=1:numseats2/2
                        position=[x(i) y z];
                        obj.plotseats(position,orientback,obj.seatwidth,obj.seatlength,obj.backrestheight,handle)
                        position=[x(i) y1 z];
                        obj.plotseats(position,orientback,obj.seatwidth,obj.seatlength,obj.backrestheight,handle);
                        end
                        seatarea3_length=(interiorlength-2*(tgap)-doorwidth)/2-(seatpitch); % calculate numseats behind door
                        numseats3=2*floor(seatarea3_length/(seatpitch));
                        
                         x1=ones(1,numseats3/2)*seatpitch;
                         x1(1)=x(1);
                         x1=cumsum(x1); % linear pattern
                         y=(interiorwidth/2-seatwidth/2);
                         y1=y-seatwidth-seatgap;
                         z=obj.seatheight+obj.floorheight;
                        for i=1:numseats3/2
                        position=[x1(i) y z];
                        obj.plotseats(position,orientback,obj.seatwidth,obj.seatlength,obj.backrestheight,handle)
                        position=[x1(i) y1 z];
                        obj.plotseats(position,orientback,obj.seatwidth,obj.seatlength,obj.backrestheight,handle);
                        end
                        
                        seatarea4_length=(interiorlength-2*(tgap)-doorwidth)/2;%
                        numseats4=2*floor(seatarea4_length/(seatpitch));
                         x=ones(1,numseats4/2)*seatpitch;
                         x(1)=doorwidth/2+seatlength;
                         x=cumsum(x); % linear pattern
                         y=(interiorwidth/2-seatwidth/2);
                         y1=y-seatwidth-seatgap;
                         z=obj.seatheight+obj.floorheight;
                        for i=1:numseats4/2
                        position=[x(i) y z];
                        obj.plotseats(position,orientback,obj.seatwidth,obj.seatlength,obj.backrestheight,handle)
                        position=[x(i) y1 z];
                        obj.plotseats(position,orientback,obj.seatwidth,obj.seatlength,obj.backrestheight,handle);
                        end
                        numseats=numseats1+numseats2+numseats3+numseats4;
                        % calculate standing space
                        aislewidth=interiorwidth-4*(seatwidth+seatgap);
                        standingarea= ((seatarea2_length-150) *aislewidth)/(1000^2);
                        standingpassengers=floor(standingarea/passengerdensity);
                        standingarea2=doorwidth*((interiorwidth-aislewidth)/2)/(1000^2);
                        standingpassengers2=floor(standingarea2/passengerdensity);
                        standingpassengers=standingpassengers+standingpassengers2;
                        obj.numberseats=numseats;
                        obj.passengercapacity= standingpassengers+numseats;
                        obj.standingpassengers= standingpassengers;
                        obj.seatingratio=obj.numberseats./obj.passengercapacity;
                        obj.costs=obj.numberseats*264;
                        obj.mass=obj.numberseats*15;
                    else
                        %
                    end
                else
                end     
            
        end
        function plotinterior2(obj,Vehicle,handle)
            
            %% plot urban layout1
            
            vehiclelength=Vehicle.Body.length;
            seatpitch=obj.seatpitch;
            seatwidth=obj.seatwidth;
            seatlength=obj.seatlength;
            interiorlength=obj.interiorlength;
            interiorwidth=obj.interiorwidth;
            vehiclewidth=Vehicle.Body.width;
            seatgap=obj.w_seat_gap;
            tgap=obj.t_body;
            doorwidth=obj.doorwidth;
            wheelchairs=0;
            doors=1;
            passengerdensity=1/4;% 4p/m^2
            floorheight=obj.floorheight;
            legroom=350;%obj.legroom;
            if doors==1
                seatarea1_length=(interiorwidth-tgap); % calculate numseat at the front, back
                numseats1=floor(seatarea1_length/(seatwidth+2*seatgap));
                xback=-0.5*interiorlength+50+seatlength/2;
                xfront=0.5*interiorlength-50-seatlength/2;
                y=linspace(-seatarea1_length/2+50+seatwidth/2,seatarea1_length/2-50-seatwidth/2,numseats1);
                z=obj.seatheight+floorheight;
                orientback=[0 0 0];
                orientfront=[pi 0 0];
                for i=1:numseats1
                    position=[xback y(i) z];
                    obj.plotseats(position,orientback,obj.seatwidth,obj.seatlength,obj.backrestheight,handle)
                    position=[xfront y(i) z];
                    obj.plotseats(position,orientfront,obj.seatwidth,obj.seatlength,obj.backrestheight,handle)
                end
                
                seatarea3_length=(interiorlength-2*(tgap))-2*(seatlength+legroom); % calculate numseats at side
                numseats3=floor(seatarea3_length/(seatwidth+2*seatgap));
                
                x=linspace(-seatarea3_length/2+seatwidth/2,seatarea3_length/2-seatwidth/2,numseats3);
                y=-0.5*vehiclewidth+50+seatlength/2;
                z=obj.seatheight+floorheight;
                orient=[-pi/2 0 0];
                for i=1:numseats3
                    position=[x(i) y z];
                    obj.plotseats(position,orient,obj.seatwidth,obj.seatlength,obj.backrestheight,handle)
                end
                seatarea4_length=(seatarea3_length-doorwidth-200)/2;
                numseats4=floor(seatarea4_length/(seatwidth+2*seatgap));
                if numseats4>0
                    x=linspace(seatarea4_length-seatwidth/2+doorwidth/2,doorwidth/2+200+seatwidth/2,numseats4);
                    y=0.5*vehiclewidth-50-seatlength/2;
                    z=obj.seatheight+floorheight;
                    orient=[pi/2 0 0];
                    for i=1:numseats4
                        position=[x(i) y z];
                        obj.plotseats(position,orient,obj.seatwidth,obj.seatlength,obj.backrestheight,handle)
                    end
                    if wheelchairs==0
                        numseats5=numseats4;
                        x=linspace(-doorwidth/2-seatwidth/2-200,-seatarea4_length+seatwidth/2-doorwidth/2,numseats4);
                        y=0.5*vehiclewidth-50-seatlength/2;
                        z=obj.seatheight+floorheight;
                        orient=[pi/2 0 0];
                        for i=1:numseats4
                            position=[x(i) y z];
                            obj.plotseats(position,orient,obj.seatwidth,obj.seatlength,obj.backrestheight,handle)
                        end
                    else
                    end
                    
                end
            end
            
            
           
        end
        
        function plotseats(obj,position,orient,seatwidth,seatlength,backrestheight,handle)
            
            seatthickness=50;
            
            
            colr = [0.1 0.1 0.1];
            alph = 1;
            n=10;
            r=0.1*backrestheight;
            v1=[0.5*seatlength;0.5*seatwidth;0.5*seatthickness];
            v2=[0.5*seatlength;0.5*seatwidth;-0.5*seatthickness];
            x= -0.5* seatlength+r*(1-(cos(linspace(0,pi/2,n))));
            y= -0.5*seatthickness+r*(1-(sin(linspace(0,pi/2,n))));
            fi11= [x;ones(1,n)*0.5*seatwidth;y];
            fi21=[x;-ones(1,n)*0.5*seatwidth;y];
            v3=[-0.5*seatlength;0.5*seatwidth;backrestheight-0.5*seatthickness];
            v4=[-0.5*seatlength+seatthickness;0.5*seatwidth;backrestheight-0.5*seatthickness];
            x= -0.5* seatlength+seatthickness+r*(1-(cos(linspace(0,pi/2,n))));
            y= 0.5*seatthickness+r*(1-(sin(linspace(0,pi/2,n))));
            fi12=[x;ones(1,n)*0.5*seatwidth;y];
            fi22=[x;-ones(1,n)*0.5*seatwidth;y];
            v6=[-0.5*seatlength;-0.5*seatwidth;backrestheight-0.5*seatthickness];
            v5=[-0.5*seatlength+seatthickness;-0.5*seatwidth;backrestheight-0.5*seatthickness];
            v8=[0.5*seatlength;-0.5*seatwidth;0.5*seatthickness];
            v7=[0.5*seatlength;-0.5*seatwidth;-0.5*seatthickness];
            V1=[v2,flip(fi11,2),v3,v6,fi21,v7];
            V2=[v1,flip(fi12,2),v4,v5,fi22,v8];
            V3=[v1,v2,flip(fi11,2),v3,v4,fi12,v1];
            V4=[v8,v7,flip(fi21,2),v6,v5,fi22,v8];
            V5=[v1,v8,v7,v2];
            V6=[v3,v4,v5,v6];
            
            V1=obj.rotate(V1',orient)';
            V2=obj.rotate(V2',orient)';
            V3=obj.rotate(V3',orient)';
            V4=obj.rotate(V4',orient)';
            V5=obj.rotate(V5',orient)';
            V6=obj.rotate(V6',orient)';
            
            V1=obj.translate(V1,position);
            V2=obj.translate(V2,position);
            V3=obj.translate(V3,position);
            V4=obj.translate(V4,position);
            V5=obj.translate(V5,position);
            V6=obj.translate(V6,position);
            
            colr=[0.9 0.9 0.9];
            patch(handle,'Faces',[1:25],'Vertices',V3','FaceColor', colr);
            patch(handle,'Faces',[1:25],'Vertices',V4','FaceColor', colr);
            patch(handle,'Faces',[1:24],'Vertices',V1','FaceColor', colr);
            patch(handle,'Faces',[1:24],'Vertices',V2','FaceColor', colr);
            patch(handle,'Faces',[1:4],'Vertices',V5','FaceColor', colr);
            patch(handle,'Faces',[1:4],'Vertices',V6','FaceColor', colr);
            
            
            
        end
        function vertices=rotate(obj,vertices,orientation)
            %% Form ZYX Rotation Matrix
            % [From Wikipedia Oct-11-2009 1:30 AM]
            % http://en.wikipedia.org/wiki/Euler_Angles
            % yaw pitch roll
            % Calculate Sines and Cosines
            c1 = cos(orientation(1));	s1 = sin(orientation(1));
            c2 = cos(orientation(2));	s2 = sin(orientation(2));
            c3 = cos(orientation(3));	s3 = sin(orientation(3));
            
            % Calculate rotation Matrix
            R = [c1*c2           -c2*s1          s2
                c3*s1+c1*s2*s3  c1*c3-s1*s2*s3  -c2*s3
                s1*s3-c1*c3*s2  c3*s1*s2+c1*s3  c2*c3]';
            vertices=(R*vertices')';
        end
        
        function vertices=translate(obj,vertices,position)
            T= [1 0 0  position(1);0 1 0  position(2);0 0 1  position(3);0 0 0 1] ;
            vertices(4,:)=1;
            vertices=(T*vertices);
            vertices=vertices(1:3,:);
        end
    end
    
    
    
    
end
