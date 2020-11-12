classdef body
    properties
        length;
        width;
        height;
        wheelbase;
        numberdecks
        groundclearance=150;
        doors
        cost
        mass
        wheelhousingheight
        wheelhousingwidth
        frontoverhang
        rearoverhang
        doorwidth=1250;
        doorheight=1900;
        superstructuremass
        superstructurecost
        sectionheight=60; %100 mm square section
        sectionwidth=120; %square section for ladder frame;
        sectionthickness=6; % 4mm thickness %c section - 210x76x6
    end
    properties (Access=private)
        aluminiumrawprice=3.5; %material costs
        structure_materialutilisation=0.52; %production cost Fuchs
        aluminiumdensity=2710; %kg/m3
        steeldensity=7800; %kg/m3
        steelrawprice=6.67; %EUR/kg
    end
    methods
        
        function obj=body(length,width,height,wheelbase,numberdecks)
            obj.length=length;
            obj.height=height;
            obj.width=width;
            obj.wheelbase=wheelbase;
            obj.numberdecks=numberdecks;
            obj.doors=((obj.length>8500)*2)+((obj.length<8500)*1);
            if obj.doors==1
                obj.frontoverhang=obj.length/2-obj.wheelbase/2;
                obj.rearoverhang=obj.frontoverhang;
            else
                wheeldiameter=950; % initial estimate to calculate FOH
                obj.frontoverhang=0.2*(obj.height-obj.groundclearance)+obj.doorwidth+1.5*wheeldiameter/2;
                obj.rearoverhang=obj.length-obj.frontoverhang-obj.wheelbase;
                
            end
            obj=updatebodyexterior(obj);
            [obj.superstructuremass,obj.superstructurecost]=superstructureestimation(obj);
            
            
        end
        function obj=update_body(obj,wheeldiameter,wheelwidth,airspringdiameter)
            % this function updates the wheelhouse dimensions and overhangs based on the
            % axle/tyre selection
            obj.doors=((obj.length>8500)*2)+((obj.length<8500)*1);
            if obj.doors==1
                obj.frontoverhang=obj.length/2-obj.wheelbase/2;
                obj.rearoverhang=obj.frontoverhang;
            else
                obj.frontoverhang=0.2*(obj.height-obj.groundclearance)+obj.doorwidth+1.5*wheeldiameter/2;
                obj.rearoverhang=obj.length-obj.frontoverhang-obj.wheelbase;
                
            end
            obj.wheelhousingheight=wheeldiameter+40;
            obj.wheelhousingwidth=wheelwidth+50+airspringdiameter;
        end
        function obj=updatebodyexterior(obj)
            
            L = obj.length/1000;
            W = obj.width/1000;
            H =obj.height/1000;
            panelthickness=5/1000;% 5mm thickness
            Volume = 2 * panelthickness * ((L * W) + (L * H) + (W * H)); % assuming 5mm thickness
            aluminiumdensity=2710; %kg/m3
            obj.mass=Volume*aluminiumdensity;
            obj.cost=obj.mass*...
                (obj.aluminiumrawprice/...
                obj.structure_materialutilisation)*1.53; % cost of body surfaces
        end
        function obj=updatebody(length,width,height,wheelbase,numberdecks)
            obj.length=length;
            obj.height=height;
            obj.width=width;
            obj.wheelbase=wheelbase;
            obj.numberdecks=obj.numberdecks;
            
        end
        function plotbody(obj,handle,vehicle)
            length=obj.length;
            width=obj.width;
            height=obj.height;
            groundclearance=obj.groundclearance;
            floorheight=vehicle.Interior.floorheight-groundclearance;
            
            
            if obj.doors==1
                %% front section
                position=[length/2-0.2*(height-groundclearance)
                    0
                    (height-groundclearance)/2];
                obj.section_end(length,width,height,groundclearance,position,floorheight,handle);
                %% door
                doorwidth=obj.doorwidth;
                doorheight=obj.doorheight;
                position=[0
                    0
                    (height-groundclearance)/2];
                
                obj.section_door(length,width,height,doorwidth,doorheight,...
                    floorheight,groundclearance,position,handle);
                
                %% wheel housings
                % if numaxles=2
                wheelbase=obj.wheelbase;
                diameter=vehicle.Chassis.tyrediameter;
                wheelwidth=vehicle.Chassis.tyrewidth;
                wheelhousingwidth=wheelwidth+50+vehicle.Chassis.airspringdiameter;
                position=[wheelbase/2 0 (height-groundclearance)/2];
                obj.section_housing(width,diameter,wheelwidth,groundclearance,height,floorheight,position,wheelhousingwidth,handle);
                position=[-wheelbase/2 0 (height-groundclearance)/2];
                obj.section_housing(width,diameter,wheelwidth,groundclearance,height,floorheight,position,wheelhousingwidth,handle);
                %% back section
                position=[-length/2+0.2*(height-groundclearance)
                    0
                    (height-groundclearance)/2];
                obj.section_end(length,width,height,groundclearance,position,floorheight,handle)
                %% sections
                % section 1
                length=(length/2-0.2*(height-groundclearance))-(wheelbase/2+1.5*diameter/2);
                if length>0
                    position=[(wheelbase/2+1.5*diameter/2+length/2),0,(height-groundclearance)/2];
                    obj.section_window(length,width,height,groundclearance,floorheight,position,handle);
                    % section 2
                    position=[-(wheelbase/2+1.5*diameter/2+length/2),0,(height-groundclearance)/2];
                    obj.section_window(length,width,height,groundclearance,floorheight,position,handle);
                end
                % section 3
                length=(wheelbase/2-1.5*diameter/2)-doorwidth/2;
                if length>0
                    position=[(doorwidth/2+length/2),0,(height-groundclearance)/2];
                    obj.section_window(length,width,height,groundclearance,floorheight,position,handle);
                    % section 4
                    length=(wheelbase/2-1.5*diameter/2)-doorwidth/2;
                    position=[-(doorwidth/2+length/2),0,(height-groundclearance)/2];
                    obj.section_window(length,width,height,groundclearance,floorheight,position,handle);
                end
                
                length=obj.length-0.4*(height-groundclearance);
                position=[0 width/2-wheelhousingwidth 0];
                plotstucture_frame(obj,handle,position,length)
                position=[0 -width/2+wheelhousingwidth 0];
                plotstucture_frame(obj,handle,position,length)
                %             axis off
                %              view(handle,[0 00])
            else
                %% section front
                position=[obj.wheelbase/2+obj.frontoverhang-0.2*(height-groundclearance)
                    0
                    (height-groundclearance)/2];
                obj.section_end(length,width,height,groundclearance,position,floorheight,handle);
                %% door 1 is at centre
                doorwidth=obj.doorwidth;
                doorheight=obj.doorheight;
                position=[0
                    0
                    (height-groundclearance)/2];
                
                obj.section_door(length,width,height,doorwidth,doorheight,...
                    floorheight,groundclearance,position,handle)
                %% door 2 is at front
                doorwidth=obj.doorwidth;
                doorheight=obj.doorheight;
                position=[obj.wheelbase/2+obj.frontoverhang-0.2*(height-groundclearance)-doorwidth/2
                    0
                    (height-groundclearance)/2];
                
                obj.section_door(length,width,height,doorwidth,doorheight,...
                    floorheight,groundclearance,position,handle)
                %% wheel housings
                % if numaxles=2
                wheelbase=obj.wheelbase;
                diameter=vehicle.Chassis.tyrediameter;
                wheelwidth=vehicle.Chassis.tyrewidth;
                wheelhousingwidth=wheelwidth+50+vehicle.Chassis.airspringdiameter;
                position=[wheelbase/2 0 (height-groundclearance)/2];
                obj.section_housing(width,diameter,wheelwidth,groundclearance,height,floorheight,position,wheelhousingwidth,handle);
                position=[-wheelbase/2 0 (height-groundclearance)/2];
                obj.section_housing(width,diameter,wheelwidth,groundclearance,height,floorheight,position,wheelhousingwidth,handle);
                %% back section
                position=[-obj.wheelbase/2-obj.rearoverhang+0.2*(height-groundclearance)
                    0
                    (height-groundclearance)/2];
                obj.section_end(length,width,height,groundclearance,position,floorheight,handle)
                %% sections
                % section 1
                length=(wheelbase/2-1.5*diameter/2)-doorwidth/2;
                if length>0
                    position=[(wheelbase/2-1.5*diameter/2-length/2),0,(height-groundclearance)/2];
                    obj.section_window(length,width,height,groundclearance,floorheight,position,handle);
                    % section 2
                    position=[-(wheelbase/2-1.5*diameter/2-length/2),0,(height-groundclearance)/2];
                    obj.section_window(length,width,height,groundclearance,floorheight,position,handle);
                end
                % section 3
                length= obj.rearoverhang-0.2*(height-groundclearance)-1.5*diameter/2;
                if length>0
                    position=[(-wheelbase/2-1.5*diameter/2-length/2),0,(height-groundclearance)/2];
                    obj.section_window(length,width,height,groundclearance,floorheight,position,handle);
                    
                end
                
                length=obj.length-0.4*(height-groundclearance);
                position=[ -obj.length/2+(wheelbase/2+obj.frontoverhang)...
                    width/2-wheelhousingwidth 0];
                plotstucture_frame(obj,handle,position,length)
                position=[ -obj.length/2+(wheelbase/2+obj.frontoverhang) ...
                    -width/2+wheelhousingwidth 0];
                plotstucture_frame(obj,handle,position,length)
                
                
            end
        end
        function [mass,cost]=superstructureestimation(obj)
            
            % This function estimates the mass of the bus superstructure
            % using the density of the material and the volume of the rail
            % sections. The function returns the total mass of the
            % structure
            
            %%
            sectionheight=60/1000; %in mm
            sectionwidth=50/1000;%100/1000; %in mm
            sectionthickness=4/1000; %in mm
            length=obj.length/1000; %vehicle length
            width=obj.width/1000; %vehicle width
            height=obj.height/1000; %vehicle height
            numwaistrails=2; % 2 waist rails
            numcantrails=2; % 2 cant rails
            numseatrails=2;
            numskirtrails=2;
            density=obj.aluminiumdensity;
            density=7800;
            %
            sectionarea=((sectionheight*sectionwidth) - ((sectionheight-sectionthickness*2)*...
                (sectionwidth-sectionthickness*2)));
            %
            % cant rail - runs longitudinally on the roof edges
            cantrailvolume=sectionarea*length;
            cantrailmass=cantrailvolume*density*numcantrails;
            
            % waist rail - longitudinally besides the seats
            waistrailvolume=sectionarea*length;
            waistrailmass=waistrailvolume*density*numwaistrails;
            % seat rail
            seatrailvolume=sectionarea*length;
            seatrailmass=seatrailvolume*density*numseatrails;
            
            % skirt rail
            skirtrailvolume=sectionarea*length;
            skirtrailmass=skirtrailvolume*density*numskirtrails;
            
            %%
            % A Pillar
            apillarvolume=sectionarea*height;
            apillarmass=apillarvolume*density*2;
            
            % Window  vertical pillars
            intercrossmemberdistance=900/1000; %mm
            numcrossmembers=floor(length/intercrossmemberdistance);
            vertpillarvolume=sectionarea*height;
            vertpillarmass=vertpillarvolume*density*2*numcrossmembers;
            
            
            % Roof arches
            roofarchvolume=sectionarea*width;
            roofarchmass=roofarchvolume*density*numcrossmembers;
            
            %roof rail
            roofarchvolume=sectionarea*length;
            lateralcrossmemberdistance=600/1000;
            numlatcrossmembers=floor(width/lateralcrossmemberdistance);
            roofrailmass=roofarchvolume*density*numlatcrossmembers;
            
            % floor rail
            intercrossmemberdistance=600/1000; %mm
            numcrossmembers=floor(length/intercrossmemberdistance);
            floorrailvolume=sectionarea*width;
            floorrailmass=floorrailvolume*density*numcrossmembers;
            
            %floor crossmember
            floorlatrailvolume=sectionarea*length;
            lateralcrossmemberdistance=600/1000;
            numlatcrossmembers=floor(width/lateralcrossmemberdistance);
            floorlatrailmass=floorlatrailvolume*density*numlatcrossmembers;
            
            
            %
            
            mass=roofarchmass+vertpillarmass+apillarmass+skirtrailmass+waistrailmass...
                +seatrailmass+cantrailmass+roofrailmass+floorrailmass+floorlatrailmass;
            cost= mass*(obj.steelrawprice/...
                obj.structure_materialutilisation)*1.53; % cost of superstructure
        end
        function plotstucture_frame(obj,handle,position,vehiclelength)
            %             position =[0 0 0];
            length=vehiclelength;
            width= obj.sectionwidth;
            height = obj.sectionheight;
            position(3)=position(3)+height/2+5;
            position(2)=((position(2)>0)* (position(2)-width/2))+...
                ((position(2)<0)* (position(2)+width/2));
            thickness=obj.sectionthickness;
            orient=[0 0 0];
            colr = [.8 .8 .8];
            alph = 1;
            %% Calculate vertices at origin
            v1=[-0.5* length;0.5* width;-0.5* height];
            v2=[0.5* length;0.5* width;-0.5* height];
            v3=[0.5* length;0.5* width;0.5* height];
            v4=[-0.5* length;0.5* width;0.5* height];
            v5=[-0.5* length;-0.5* width;0.5* height];
            v6=[0.5* length;-0.5* width;0.5* height];
            v7=[0.5* length;-0.5* width;-0.5* height];
            v8=[-0.5* length;-0.5* width;-0.5* height];
            
            v9=[-0.5* length;0.5* width-thickness;-0.5* height+thickness];
            v10=[0.5* length;0.5* width-thickness;-0.5* height+thickness];
            v11=[0.5* length;0.5* width-thickness;0.5* height-thickness];
            v12=[-0.5* length;0.5* width-thickness;0.5* height-thickness];
            v13=[-0.5* length;-0.5* width+thickness;0.5* height-thickness];
            v14=[0.5* length;-0.5* width+thickness;0.5* height-thickness];
            v15=[0.5* length;-0.5* width+thickness;-0.5* height+thickness];
            v16=[-0.5* length;-0.5* width+thickness;-0.5* height+thickness];
            
            
            Vertices=[v1,v2,v3,v4,v5,v6,v7,v8];
            Vertices=[Vertices,v9,v10,v11,v12,v13,v14,v15,v16];
            Vertices=obj.translate(Vertices,position);
            Faces=[1 2 3 4
                5 6 7 8
                4 3 6 5
                %3 2 7 6
                2 1 8 7
                % 1 4 5 8
                ];
            Faces=[Faces;8+Faces;[1 9 12 4;12 4 5 13;...
                13 5 8 16; 16 8 1 9];[2 10 11 3;11 3 6 14;...
                14 6 7 15; 15 7 2 10]];
            patch(handle,'Faces',Faces,'Vertices',Vertices','FaceColor', colr,'FaceAlpha',alph);
            
        end
        function section_housing(obj,width,wheeldiameter,wheelwidth,groundclearance,height,floorheight,position,wheelhousingwidth,handle)
            
            height=height-groundclearance;
            orient=[0 0 0];
            colr = [1 1 1];
            alph = 0.5;
            length = 1.5*wheeldiameter;
            bottomsidewallheight=950; %
            windowheight=0.75*(height-bottomsidewallheight);
            sidewallheight=height-windowheight-bottomsidewallheight;
             hw=windowheight;
            %% roof
            % side right surface
            n=10; % points in the fillet radi
            r=0.05*height;
            h=sidewallheight; %height/3;
            v1=[0.5* length;0.5* width; -0.5* h];
            v2=[-0.5* length;0.5* width;-0.5* h];
            v3=[-0.5* length;0.5* width;(0.5)* h-r];
            v4=[0.5* length;0.5* width;(0.5)* h-r];
            x= (0.5)* width-r*(1-(cos(linspace(0,pi/2,n))));
            y= (0.5)* h-r*(1-(sin(linspace(0,pi/2,n))));
            fi11= [ones(1,n)*0.5*length;x;y];
            fi12=[-ones(1,n)*0.5*length;x;y];
            fil1=sortrows([fi11  fi12]',3,'ascend')';
            % top surface
            v5=[-0.5*length; (0.5)* width-r;0.5* h];
            v6=[0.5* length;(0.5)* width-r;0.5* h];
            v7=[0.5* length;-(0.5)* width+r;0.5* h];
            v8=[-0.5* length;-(0.5)* width+r;0.5* h];
            % left wall
            v9=[0.5* length;-0.5* width; -0.5* h];
            v10=[-0.5* length;-0.5* width;-0.5* h];
            v11=[-0.5* length;-0.5* width;(0.5)* h-r];
            v12=[0.5* length;-0.5* width;(0.5)* h-r];
            % right fillet
            x= -(0.5)* width+r*(1-(cos(linspace(0,pi/2,n))));
            y= (0.5)* h-r*(1-(sin(linspace(0,pi/2,n))));
            fi21= [ones(1,n)*0.5*length;x;y];
            fi22=[-ones(1,n)*0.5*length;x;y];
            fil2=sortrows([fi21  fi22]',3,'descend')';
            
            roofVertices=[v1,v2,v3,v4,fil1,v5,v6,v7,v8,fil2,v9,v10,v11,v12];
            rf1Vertices=[v1,v2,v3,v4];
            rf2Vertices=[fi12,flip(fi11,2)];
            rf3Vertices=[v9,v10,v11,v12];
            rf4Vertices=[fi22,flip(fi21,2),v12];
            rf5Vertices=[v5,v6,v7,v8];
            
            rf1Vertices=obj.translate(rf1Vertices,[0 0 0.5*hw+h/2]); % translate to origin 0 0 0
            rf2Vertices=obj.translate(rf2Vertices,[0 0 0.5*hw+h/2]);
            rf3Vertices=obj.translate(rf3Vertices,[0 0 0.5*hw+h/2]);
            rf4Vertices=obj.translate(rf4Vertices,[0 0 0.5*hw+h/2]);
            rf5Vertices=obj.translate(rf5Vertices,[0 0 0.5*hw+h/2]);
            
            
            %% create bottom surfaces
            % side right surface
            n=10; % points in the fillet radi
            r=0.05*height;
           
            h=bottomsidewallheight; %height/3; % height of right surface
            l1=0.5*length-0.2*wheeldiameter;
            v1=[0.5* length;0.5* width;(0.5)* h];
            v2=[0.5* length;0.5* width;-(0.5)* h+r];
            v3=[l1;0.5* width;-(0.5)* h+r];
            v4=[l1;0.5* width;-(0.5)*h+0.25*wheeldiameter];
            
            x= (0.5)* width-r*(1-(cos(linspace(0,pi/2,n))));
            y= -(0.5)* h+r*(1-(sin(linspace(0,pi/2,n))));
            
            fi11= [ones(1,n)*0.5*length;x;y];
            fi13=[-ones(1,n)*0.5*length;x;y];
            fi12= [ones(1,n)*l1;x;y];
            fi14=[-ones(1,n)*l1;x;y];
            fil1=sortrows([fi11  fi12 fi13 fi14]',3,'descend')';
            
            n=20;
            
            x= 1.1*wheeldiameter/2*(cos(linspace(0,pi,n)));
            y= -(0.5)*h+floorheight+0.9*wheeldiameter/2*(sin(linspace(0,pi,n))); %% would need to modify for large wheel sizes
            hsg1=[x;0.5* width+0*y;y];
            v5=[-l1;0.5* width;-(0.5)*h+0.25*wheeldiameter];
            v6=[-l1;0.5* width;-(0.5)* h+r];
            v8=[-0.5* length;0.5* width;(0.5)* h];
            v7=[-0.5* length;0.5* width;-(0.5)* h+r];
            rwallVertices=[v1,v2,v3,v4,hsg1,v5,v6,v7,v8];
            rwallVertices=obj.translate(rwallVertices,[0 0 -0.5*hw-h/2]);
            floorheight=-height/2+floorheight;
            fl1=[0.5*length;0.5* width;floorheight];
            fl2=[0.5*length;-0.5* width;floorheight];
            fl3=[l1;-0.5* width;floorheight];
            fl4=[l1;0.5* width;floorheight];
            fl5=[-0.5*length;0.5* width;floorheight];
            fl6=[-0.5*length;-0.5* width;floorheight];
            fl7=[-l1;-0.5* width;floorheight];
            fl8=[-l1;0.5* width;floorheight];
            fl9=[l1;-0.5* width+wheelhousingwidth;floorheight];
            fl10=[l1;0.5* width-wheelhousingwidth;floorheight];
            fl12=[-l1;-0.5* width+wheelhousingwidth;floorheight];
            fl11=[-l1;0.5* width-wheelhousingwidth;floorheight];
            fl13=[0.5*length;-0.5* width+wheelhousingwidth;-0.5*height];
            fl14=[0.5*length;0.5* width-wheelhousingwidth;-0.5*height];
            fl16=[-0.5*length;-0.5* width+wheelhousingwidth;-0.5*height];
            fl15=[-0.5*length;0.5* width-wheelhousingwidth;-0.5*height];
            
            F1=[fl1,fl2,fl3,fl4];
            F2=[fl5,fl6,fl7,fl8];
            F3=[fl9,fl10,fl11,fl12];
            F4=[fl13,fl14,fl15,fl16];
            
            
            % left surface
            leftvertices=[rwallVertices(1,:);-rwallVertices(2,:);rwallVertices(3,:)];
            
            
            % first fillet
            fi21=fi11;
            fi21(2,:)=-fi21(2,:);
            fi22=fi12;
            fi22(2,:)=-fi22(2,:);
            strip1Vertices=[fi11,[fi11(1,end);0.5*width-wheelhousingwidth;fi11(3,end)],...
                [fi12(1,end);0.5*width-wheelhousingwidth;fi12(3,end)],flip(fi12,2)];
            strip2Vertices=[fi22,[fi22(1,end);-0.5*width+wheelhousingwidth;fi22(3,end)],...
                [fi21(1,end);-0.5*width+wheelhousingwidth;fi21(3,end)],flip(fi21,2)];
            
            strip1Vertices=obj.translate(strip1Vertices,[0 0 -0.5*hw-h/2]);
            strip2Vertices=obj.translate(strip2Vertices,[0 0 -0.5*hw-h/2]);
            
            % 2nd fillet
            fi21=fi13;
            fi21(2,:)=-fi21(2,:);
            fi22=fi14;
            fi22(2,:)=-fi22(2,:);
            strip3Vertices=[fi13,[fi13(1,end);0.5*width-wheelhousingwidth;fi13(3,end)],...
                [fi14(1,end);0.5*width-wheelhousingwidth;fi14(3,end)],flip(fi14,2)];
            strip4Vertices=[fi22,[fi22(1,end);-0.5*width+wheelhousingwidth;fi22(3,end)],...
                [fi21(1,end);-0.5*width+wheelhousingwidth;fi21(3,end)],flip(fi21,2)];
            strip3Vertices=obj.translate(strip3Vertices,[0 0 -0.5*hw-h/2]);
            strip4Vertices=obj.translate(strip4Vertices,[0 0 -0.5*hw-h/2]);
            % wheel housing surface
            hsgv2=hsg1;
            hsgv2(2,:)=hsgv2(2,:)-wheelhousingwidth;
            housing1Vertices=[hsg1,flip(hsgv2,2)];
            housing1Vertices=obj.translate(housing1Vertices,[0 0 -0.5*hw-h/2]);
            V=[hsg1(:,1),fi12,[hsgv2(1:2,1);fi14(3,end)],[hsgv2(1:2,1);hsg1(3,1)]];
            V1=[-V(1,:);V(2:3,:)];
            V2=[V(1,:);-V(2,:);V(3,:)];
            V3=[V1(1,:);-V1(2,:);V1(3,:)];
            V4=[hsgv2];
            V5=[V4(1,:);-V4(2,:);V4(3,:)];
            
            hsgv2=hsg1;
            hsg1(2,:)=-hsg1(2,:);
            hsgv2(2,:)=-hsgv2(2,:)+wheelhousingwidth;
            housing2Vertices=[hsg1,flip(hsgv2,2)];
            housing2Vertices=obj.translate(housing2Vertices,[0 0 -0.5*hw-h/2]);
            
            V=obj.translate(V,[0 0 -0.5*hw-h/2]);
            V1=obj.translate(V1,[0 0 -0.5*hw-h/2]);
            V2=obj.translate(V2,[0 0 -0.5*hw-h/2]);
            V3=obj.translate(V3,[0 0 -0.5*hw-h/2]);
            V4=obj.translate(V4,[0 0 -0.5*hw-h/2]);
            V5=obj.translate(V5,[0 0 -0.5*hw-h/2]);
            
            
            %%
            %% left window
            h=windowheight; %height/3;
            v1=[0.5* length;-0.5* width; 0.5* h];
            v2=[-0.5* length;-0.5* width;0.5* h];
            v3=[-0.5* length;-0.5* width;-0.5* h];
            v4=[0.5* length;-0.5* width;-0.5* h];
            lwindowVertices=[v1,v2,v3,v4];
            lwindowFaces=[1 2 3 4];
            
            %'LineStyle','none'
            %% right window
            
            v1=[0.5* length;0.5* width; 0.5* h];
            v2=[-0.5* length;0.5* width;0.5* h];
            v3=[-0.5* length;0.5* width;-0.5* h];
            v4=[0.5* length;0.5* width;-0.5* h];
            rwindowVertices=[v1,v2,v3,v4];
            rwindowFaces=[1 2 3 4];
            
            
            
            %% plot
             F1=obj.translate(F1,position); %floor surfaces
            F2=obj.translate(F2,position);
            F3=obj.translate(F3,position);
            F4=obj.translate(F4,position);
            
            position(3)=position(3)+hw/2+ bottomsidewallheight- height/2;
            roofVertices=obj.translate(roofVertices,position);
            rwallVertices=obj.translate(rwallVertices,position);
            leftvertices=obj.translate(leftvertices,position);
            strip1Vertices=obj.translate(strip1Vertices,position);
            strip2Vertices=obj.translate(strip2Vertices,position);
            strip3Vertices=obj.translate(strip3Vertices,position);
            strip4Vertices=obj.translate(strip4Vertices,position);
            housing1Vertices=obj.translate(housing1Vertices,position);
            housing2Vertices=obj.translate(housing2Vertices,position);
            V=obj.translate(V,position);
            V1=obj.translate(V1,position);
            V2=obj.translate(V2,position);
            V3=obj.translate(V3,position);
            V4=obj.translate(V4,position);
            V5=obj.translate(V5,position);
            lwindowVertices=obj.translate(lwindowVertices,position);
            rwindowVertices=obj.translate(rwindowVertices,position);
            rf1Vertices=obj.translate(rf1Vertices,position); %roof surfaces
            rf2Vertices=obj.translate(rf2Vertices,position);
            rf3Vertices=obj.translate(rf3Vertices,position);
            rf4Vertices=obj.translate(rf4Vertices,position);
            rf5Vertices=obj.translate(rf5Vertices,position);
           
            alph1=0.8;
            rightwall=patch(handle,'Faces', [1:28],'Vertices', rwallVertices','FaceColor', colr,'FaceAlpha',alph1);
            leftwall=patch(handle,'Faces', [1:28],'Vertices', leftvertices','FaceColor', colr,'FaceAlpha',alph1);
            strip1=patch(handle,'Faces', [1:22],'Vertices',strip1Vertices','FaceColor', colr);
            strip2=patch(handle,'Faces', [1:22],'Vertices',strip2Vertices','FaceColor', colr);
            strip3=patch(handle,'Faces', [1:22],'Vertices',strip3Vertices','FaceColor', colr);
            strip4=patch(handle,'Faces', [1:22],'Vertices',strip4Vertices','FaceColor', colr);
            
            housing1=patch(handle,'Faces', [1:40],'Vertices',housing1Vertices','FaceColor', colr);
            patch(handle,'Faces', [1:40],'Vertices',housing2Vertices','FaceColor', colr);
            patch(handle,'Faces', [1:13],'Vertices',V','FaceColor', colr);
            patch(handle,'Faces', [1:13],'Vertices',V1','FaceColor', colr);
            patch(handle,'Faces', [1:13],'Vertices',V2','FaceColor', colr);
            patch(handle,'Faces', [1:13],'Vertices',V3','FaceColor', colr);
            
            
            patch(handle,'Faces', [1:20],'Vertices',V4','FaceColor', colr);
            patch(handle,'Faces', [1:20],'Vertices',V5','FaceColor', colr);
            
            
            patch(handle,'Faces', [1:4],'Vertices', rf1Vertices','FaceColor', colr);
            patch(handle,'Faces', [1:20],'Vertices', rf2Vertices','FaceColor', colr);
            patch(handle,'Faces', [1:4],'Vertices', rf3Vertices','FaceColor', colr);
            patch(handle,'Faces', [1:20],'Vertices', rf4Vertices','FaceColor', colr);
            patch(handle,'Faces', [1:4],'Vertices', rf5Vertices','FaceColor', colr);
            
            patch(handle,'Faces', [1:4],'Vertices', F1','FaceColor', colr);
            patch(handle,'Faces', [1:4],'Vertices', F2','FaceColor', colr);
            patch(handle,'Faces', [1:4],'Vertices', F3','FaceColor', colr);
            patch(handle,'Faces', [1:4],'Vertices', F4','FaceColor', colr,'FaceAlpha',0.5);
            
            colr=[0 0 1];
            alph=0.2;
            lwindow=patch(handle,'Faces', lwindowFaces,'Vertices', lwindowVertices','FaceColor', colr,'FaceAlpha',alph);
            rwindow=patch(handle,'Faces', rwindowFaces,'Vertices', rwindowVertices','FaceColor', colr,'FaceAlpha',alph);
            
            
            
            
            
        end
        function section_end(obj,length,width,height,groundclearance,position,floorheight,handle)
            
            height=height-groundclearance;
            orient=[0 0 0];
            colr = [1 1 1];
            alph = 0.5;
            %% define profile
            %% create profile to sweep
            n=10; % points in the fillet radi
            r=0.05*height;
            v1=[0; (0.5)* width-r;0];
            v2=[0;(0.5)* width-r;0];
            v3=[0;-(0.5)* width+r;0];
            v4=[0;-(0.5)* width+r;0];
            x= (0.5)* width-r*(1-(cos(linspace(0,pi/2,n))));
            y= -r*(1-(sin(linspace(0,pi/2,n))));
            fi11= [ones(1,n)*0;x;y];
            fil1=sortrows(fi11',3,'ascend')';
            x= -(0.5)* width+r*(1-(cos(linspace(0,pi/2,n))));
            y= 0-r*(1-(sin(linspace(0,pi/2,n))));
            fi21= [ones(1,n)*0;x;y];
            profile=[fil1 flip(fi21,2)]';
            
            %%  Define path
            r=0.10*height;
            n=10;
            x= 0.10*height+r*((cos(linspace(pi/2,0,n))));
            z= -0.5*height+r*(1-(sin(linspace(pi/2,0,n))));
            fi11= [x;ones(1,n)*0;z];
            fi12=[x;-ones(1,n)*0.5*length;z];
            fil1=sortrows([fi11]',3,'ascend')';
            fil1=[fil1,[x(end);0; -0.5*height+floorheight]];
            
            
            v5=[fil1(1,end);0.5*width;fil1(3,end)];
            v6=[fil1(1,end);-0.5*width;fil1(3,end)];
            
            v7=[0.15*height;0.5*width;0.5*height-0.15*height];
            v8=[0.15*height;-0.5*width;0.5*height-0.15*height];
            %
            Vertarcs=arc2points(0,[v6(1) v6(3)],[v8(1) v8(3)]); %returns tangent curver to 2 points
            Vertarcs(:,3)=Vertarcs(:,2);
            Vertarcs(:,2)=0;
            
            
            Vn=[0 ; 0; 0.5*height];
            Vn_1=[0.15*height/2;0;0.5*height];
            P=[0.15*height/2;0;0.5*height];
            dA=([0.15*height/2-0 ; 0 ; 0.5*height-0.5*height])./(0.15*height/2); % unit vector
            dB=([(0.15*height/2-Vertarcs(end,1));0;(0.5*height-Vertarcs(end,3))])./sqrt(...
                (0.15*height/2-Vertarcs(end,1))^2 +(0.5*height-Vertarcs(end,3))^2); % unit vector
            R=210;%/1000;
            a0=(2*pi-acos(dot(dA,dB)/(norm(dA)*norm(dB)))-pi); % angle
            Cx=P(1)-R*tan(a0/2);
            Cy=P(3)-R;
            n=10;
            x= Cx+R*((cos(linspace(pi/2-a0,pi/2,n))));
            z= Cy+R*((sin(linspace(pi/2-a0,pi/2,n))));
            path=[Vn(1) Vn(2) Vn(3);flip(x)' (x)'*0 flip(z)';flip(Vertarcs);flip(fil1(:,1:end-1),2)' ;0 0 -height*0.5];
            [facel,facer,roof,ws,bot,vertices]=sweep(path,profile,-0.5*height+floorheight);
            
            
            % define faces for windshield
            facews=[];
            facerf=[];
            facebt=[];
            n=1;
            n1=1;
            for i=1:size(roof,1)/20-1
                facerf=[facerf; [linspace(n,n+20-1,20) linspace(n+40-1,n+20,20)]];
                if i<=size(bot,1)/20-1
                    facebt=[facebt; [linspace(n,n+20-1,20) linspace(n+40-1,n+20,20)]];
                    if i<=size(ws,1)/20-1
                        facews=[facews; [linspace(n,n+20-1,20) linspace(n+40-1,n+20,20)]];
                    end
                end
                n=n+20;
            end
            floor=vertices(vertices(:,3)==-0.5*height+floorheight,:);
            if position(1)<0 % flip if less than 0
                facel(:,1)=-facel(:,1);
                facer(:,1)=-facer(:,1);
                roof(:,1)=-roof(:,1);
                bot(:,1)=-bot(:,1);
                ws(:,1)=-ws(:,1);
                floor(:,1)=-floor(:,1);
            end
            
            facel=obj.translate(facel',position);
            facer=obj.translate(facer',position);
            roof=obj.translate(roof',position);
            bot=obj.translate(bot',position);
            ws=obj.translate(ws',position);
            floor=obj.translate(floor',position);
            floor=[[facel(1:2,1);floor(3,1)],floor,[facer(1:2,1);floor(3,1)]];
            leftwall=patch(handle,'Faces', [1:32],'Vertices', facel','FaceColor', colr);
            rightwall=patch(handle,'Faces', [1:32],'Vertices', facer','FaceColor', colr);
            roofpatch=patch(handle,'Faces', facerf,'Vertices', roof','FaceColor', colr,'LineStyle','none');
            bottom=patch(handle,'Faces', facebt,'Vertices', bot','FaceColor', colr,'LineStyle','none');
            floor=patch(handle,'Faces', [1:22],'Vertices', floor','FaceColor', colr);
            colr=[0 0 1];
            alph=0.2;
            windShield=patch(handle,'Faces', facews,'Vertices', ws','FaceColor', colr,'FaceAlpha',alph,'LineStyle','none');
            
            
            
            function [facel,facer,roof,ws,bot,vertices]=sweep(path,profile,floorheight)
                path=round(path,4);
                profile=round(profile,4);
                
                % cog=[sum(1 .* profile(:,1)) / sum(profile(:,1)*0+1) ...
                %     sum(1 .* profile(:,2)) / sum(profile(:,1)*0+1)...
                %     sum(1 .* profile(:,3)) / sum(profile(:,1)*0+1)]; Centre of gravity
                vin=obj.translate(profile',path(1,:));
                vertices=vin';
                % figure
                %
                % plot3(path(:,1),path(:,2),path(:,3))
                % hold on
                % plot3(vertices(:,1),vertices(:,2),vertices(:,3))
                % axis equal
                facel=vertices(1,:);
                facer=vertices(end,:);
                for i=2:size(path,1)-1
                    orientation=[((path(i,2)-path(i-1,2))/(path(i,1)-path(i-1,1))) ... %yaw
                        ((path(i,3)-path(i-1,3))/(path(i,1)-path(i-1,1))) ... % pitch
                        (((path(i,3)-path(i-1,3))/(path(i,2)-path(i-1,2)))) ... %roll
                        ];
                    orientation( isinf(orientation) | isnan(orientation))=0;
                    orientation=atan(orientation);
                    if orientation(2)>0
                        orientation(2)=orientation(2)-pi;
                        
                    end
                    if orientation(1)~=0 | orientation(2)~=0 | orientation(3)~=0
                        vin=obj.rotate(profile,orientation);
                        vin=obj.translate(vin',path(i,:));
                    else
                        vin=vin';
                        vin=obj.translate(vin',path(i,:)-path(i-1,:));
                    end
                    
                    %     plot3(vin(1,:)',vin(2,:)',vin(3,:)')
                    facel=[facel;vin(:,1)'];
                    facer=[facer;vin(:,end)'];
                    vertices=[vertices;vin'];
                    % find orientation
                    %  roll - dzdy
                    %  pitch - dxdz
                    %  yaw  - dydx
                    
                end
                vin=obj.rotate(profile,[0 pi 0]);
                vin=obj.translate(vin',path(end,:));
                facel=[facel;vin(:,1)'];
                facer=[facer;vin(:,end)'];
                vertices(401:420,3)=floorheight;
                vertices=[vertices;vin'];
                roof=vertices(1:240,:);
                ws=vertices(221:420,:);
                bot=vertices(401:end,:);
            end
            
            function vertices=arc2points(a0,c1,c2)
                % a0 is the initial starting angle in radians % angle is 2 * tan opp/ adj
                % c1 1st coordinate
                % c2 2nd coordinate
                
                x=[];
                r=((c1(1)-c2(1))^2 + (c1(2)-c2(2))^2)/(2*(abs(c1(1)-c2(1))*cos(a0)) + 2*(abs(c1(2)-c2(2))*sin(a0)));
                
                % find centre
                coeff=-(c1*2)-(-c2*2);
                coeffc=sum(c1.^2-c2.^2);
                
                a= 1+(coeff(1)/coeff(2))^2;
                b= 2*-c2(1) +(-(coeffc(1)/coeff(2))-c2(2))*2*-(coeff(1)/coeff(2));
                c= c2(1)^2 + (-(coeffc(1)/coeff(2))-c2(2))^2 - (r^2);
                x=roots([a b c]);
                y=-((coeff(1).*x)+coeffc)/coeff(2);
                a1=atan((((c1(1)-c2(1))^2-(c1(2)-c2(2))^2)*sin(a0)-2*(((c1(1)-c2(1))*(c1(2)-c2(2)))*cos(a0)))/ ...
                    (((c1(2)-c2(2))^2-(c1(1)-c2(1))^2)*cos(a0)-2*(((c1(1)-c2(1))*(c1(2)-c2(2)))*sin(a0)))); % end angle
                c=[x y];
                
                
                %  hold on
                x=r*cos(linspace(0,a1,10));
                y=r*sin(linspace(0,a1,10));
                x=x+c(2,1);
                y=y+c(2,2);
                % plot3(x,x*0,y);
                vertices=[x' y'];
            end
            
        end
        function section_door(obj,length,width,height,doorwidth,doorheight,...
                floorheight,groundclearance,position,handle)
             orient=[0 0 0];
            length=doorwidth;
            colr = [1 1 1];
            alph = 0.2;
            height=height-groundclearance;
            bottomsidewallheight=950; %
            windowheight=0.75*(height-bottomsidewallheight);
            sidewallheight=height-windowheight-bottomsidewallheight;
            hw=windowheight;
            % position =[0 0 height/2];
            %% entrance height
            % bottom surface
            
            v1=[-0.5*length; (0.5)* width;-height/2+floorheight];
            v2=[0.5* length;(0.5)* width;-height/2+floorheight];
            v3=[0.5* length;-(0.5)* width;-height/2+floorheight];
            v4=[-0.5* length;-(0.5)* width;-height/2+floorheight];
            floorVertices=[v1,v2,v3,v4];
            floorfaces=[1 2 3 4];
            
            %% create bottom bucket
            h1=(height-doorheight)/2;
            h2=height/3;
            % side right surface
            n=10; % points in the fillet radi
            r=0.05*height;
            v1=[0.5* length;0.5* width;  -(0.5)* h2+ floorheight];
            v2=[-0.5* length;0.5* width; -(0.5)* h2+ floorheight];
            v3=[-0.5* length;0.5* width;-(0.5)* h2+r];
            v4=[0.5* length;0.5* width;-(0.5)* h2+r];
            x= (0.5)* width-r*(1-(cos(linspace(0,pi/2,n))));
            y= -(0.5)* h2+r*(1-(sin(linspace(0,pi/2,n))));
            fi11= [ones(1,n)*0.5*length;x;y];
            fi12=[-ones(1,n)*0.5*length;x;y];
            fil1=sortrows([fi11  fi12]',3,'descend')';
            % bottom surface
            v5=[-0.5*length; (0.5)* width-r;-0.5* h2];
            v6=[0.5* length;(0.5)* width-r;-0.5* h2];
            v7=[0.5* length;-(0.5)* width+r;-0.5* h2];
            v8=[-0.5* length;-(0.5)* width+r;-0.5* h2];
            % left wall
            v9=[0.5* length;-0.5* width; 0.5* h2];
            v10=[-0.5* length;-0.5* width;0.5* h2];
            v11=[-0.5* length;-0.5* width;-(0.5)* h2+r];
            v12=[0.5* length;-0.5* width;-(0.5)* h2+r];
            % right fillet
            x= -(0.5)* width+r*(1-(cos(linspace(0,pi/2,n))));
            y= -(0.5)* h2+r*(1-(sin(linspace(0,pi/2,n))));
            fi21= [ones(1,n)*0.5*length;x;y];
            fi22=[-ones(1,n)*0.5*length;x;y];
            fil2=sortrows([fi21  fi22]',3)';
            bottomVertices=[v1,v2,v3,v4,fil1,v5,v6,v7,v8,fil2,v9,v10,v11,v12];
            bot1Vertices=[v1,v2,v3,v4];
            bot2Vertices=[fi12,flip(fi11,2)];
            bot3Vertices=[v9,v10,v11,v12];
            bot4Vertices=[fi22,flip(fi21,2)];
            bot5Vertices=[v5,v6,v7,v8];
            
            bottomVertices=obj.translate(bottomVertices,[0 0 -h2]);
            bot1Vertices=obj.translate(bot1Vertices,[0 0 -h2]);
            bot2Vertices=obj.translate(bot2Vertices,[0 0 -h2]);
            bot3Vertices=obj.translate(bot3Vertices,[0 0 -h2]);
            bot4Vertices=obj.translate(bot4Vertices,[0 0 -h2]);
            bot5Vertices=obj.translate(bot5Vertices,[0 0 -h2]);
            
            %% roof
            % side right surface
            n=10; % points in the fillet radi
            h3=height/3;
            h4=775-h3/2;%height-h3-doorheight-floorheight; % remaining panel height
            h4= (doorheight/2)-((height/2-floorheight)-doorheight/2);
            r=0.05*height;
            v1=[0.5* length;0.5* width;  h4];
            v2=[-0.5* length;0.5* width;  h4];
            v3=[-0.5* length;0.5* width; height/2-r];
            v4=[0.5* length;0.5* width;height/2-r];
            x= (0.5)* width-r*(1-(cos(linspace(0,pi/2,n))));
            y=  height/2-r*(1-(sin(linspace(0,pi/2,n))));
            fi11= [ones(1,n)*0.5*length;x;y];
            fi12=[-ones(1,n)*0.5*length;x;y];
            fil1=sortrows([fi11  fi12]',3,'ascend')';
            % top surface
            v5=[-0.5*length; (0.5)* width-r;height/2];
            v6=[0.5* length;(0.5)* width-r;height/2];
            v7=[0.5* length;-(0.5)* width+r;height/2];
            v8=[-0.5* length;-(0.5)* width+r;height/2];
            % left wall
            v9=[0.5* length;-0.5* width; height/2-h3];
            v10=[-0.5* length;-0.5* width;height/2-h3];
            v11=[-0.5* length;-0.5* width;height/2-r];
            v12=[0.5* length;-0.5* width; height/2-r];
            % right fillet
            x= -(0.5)* width+r*(1-(cos(linspace(0,pi/2,n))));
            y= height/2-r*(1-(sin(linspace(0,pi/2,n))));
            fi21= [ones(1,n)*0.5*length;x;y];
            fi22=[-ones(1,n)*0.5*length;x;y];
            fil2=sortrows([fi21  fi22]',3,'descend')';
            roofVertices=[v1,v2,v3,v4,fil1,v5,v6,v7,v8,fil2,v9,v10,v11,v12];
            rf1Vertices=[v1,v2,v3,v4];
            rf2Vertices=[fi12,flip(fi11,2)];
            rf3Vertices=[v9,v10,v11,v12];
            rf4Vertices=[fi22,flip(fi21,2),v12];
            rf5Vertices=[v5,v6,v7,v8];
            
            
            %roofVertices=translate(roofVertices,[0 0 h3]);
            %Vertices=[Vertices fil1' fil2']
            
            
            %roofVertices=translate(roofVertices,[0 0 2*h3])
            
            
            %% left window
            h4=height/3;
            v1=[0.5* length;-0.5* width; 0.5* h4];
            v2=[-0.5* length;-0.5* width;0.5* h4];
            v3=[-0.5* length;-0.5* width;-0.5* h4];
            v4=[0.5* length;-0.5* width;-0.5* h4];
            lwindowVertices=[v1,v2,v3,v4];
            lwindowFaces=[1 2 3 4];
            
            
            
            %% right door
            
            v1=[0.5* doorwidth;0.5* width; 0.5* doorheight];
            v2=[-0.5* doorwidth;0.5* width;0.5* doorheight];
            v3=[-0.5* doorwidth;0.5* width;-0.5* doorheight];
            v4=[0.5* doorwidth;0.5* width;-0.5* doorheight];
            doorVertices=[v1,v2,v3,v4];
            doorFaces=[1 2 3 4];
            
            doorVertices=obj.translate(doorVertices,[0 0 -((height/2-floorheight)-doorheight/2)]);
            
            
            
            floorVertices=obj.translate(floorVertices,position);
            doorVertices=obj.translate(doorVertices,position);
            roofVertices=obj.translate(roofVertices,position);
            bottomVertices=obj.translate(bottomVertices,position);
            lwindowVertices=obj.translate(lwindowVertices,position);
            rf1Vertices=obj.translate(rf1Vertices,position);
            rf2Vertices=obj.translate(rf2Vertices,position);
            rf3Vertices=obj.translate(rf3Vertices,position);
            rf4Vertices=obj.translate(rf4Vertices,position);
            rf5Vertices=obj.translate(rf5Vertices,position);
            bot1Vertices=obj.translate(bot1Vertices,position);
            bot2Vertices=obj.translate(bot2Vertices,position);
            bot3Vertices=obj.translate(bot3Vertices,position);
            bot4Vertices=obj.translate(bot4Vertices,position);
            bot5Vertices=obj.translate(bot5Vertices,position);
            
            floor=patch(handle,'Faces', floorfaces,'Vertices', floorVertices','FaceColor', colr); %,'FaceAlpha',alph)
            %bucket=patch('Faces', bottomfaces,'Vertices', bottomVertices','FaceColor', colr) %,'FaceAlpha',alph)
            %roof=patch('Faces', rooffaces,'Vertices', roofVertices','FaceColor', colr) %,'FaceAlpha',alph)
            patch(handle,'Faces', [1:4],'Vertices', rf1Vertices','FaceColor', colr);
            patch(handle,'Faces', [1:20],'Vertices', rf2Vertices','FaceColor', colr);
            patch(handle,'Faces', [1:4],'Vertices', rf3Vertices','FaceColor', colr);
            patch(handle,'Faces', [1:20],'Vertices', rf4Vertices','FaceColor', colr);
            patch(handle,'Faces', [1:4],'Vertices', rf5Vertices','FaceColor', colr);
            alph1=0.8;
            patch(handle,'Faces', [1:4],'Vertices', bot1Vertices','FaceColor', colr,'FaceAlpha',alph1);
            patch(handle,'Faces', [1:20],'Vertices', bot2Vertices','FaceColor', colr,'FaceAlpha',alph1);
            patch(handle,'Faces', [1:4],'Vertices', bot3Vertices','FaceColor', colr,'FaceAlpha',alph1);
            patch(handle,'Faces', [1:20],'Vertices', bot4Vertices','FaceColor', colr,'FaceAlpha',alph1);
            patch(handle,'Faces', [1:4],'Vertices', bot5Vertices','FaceColor', colr,'FaceAlpha',0.5);
            
            
            colr=[0 0 1];
            alph=0.2;
            lwindow=patch(handle,'Faces', lwindowFaces,'Vertices', lwindowVertices','FaceColor', colr,'FaceAlpha',alph);
            
            colr=[0 0 0];
             alph=0.52;
            door=patch(handle,'Faces',doorFaces,'Vertices', doorVertices','FaceColor', colr,'FaceAlpha',alph);
            
            
        end
        function section_window(obj,length,width,height,groundclearance,floorheight,position,handle)
            
            height=height-groundclearance;
            orient=[0 0 0];
            colr = [1 1 1];
            alph = 0.5;
            %position =[2000 0 height/2];
            bottomsidewallheight=950; %
            windowheight=0.75*(height-bottomsidewallheight);
            sidewallheight=height-windowheight-bottomsidewallheight;
            
            %% floor
            % bottom surface
            v1=[-0.5*length; (0.5)* width;-0.5*height+floorheight];
            v2=[0.5* length;(0.5)* width;-0.5*height+floorheight];
            v3=[0.5* length;-(0.5)* width;-0.5*height+floorheight];
            v4=[-0.5* length;-(0.5)* width;-0.5*height+floorheight];
            floorVertices=[v1,v2,v3,v4];
            floorfaces=[1 2 3 4];
            
             
            
            %% create bottom bucket
            % side right surface
            n=10; % points in the fillet radi
            r=0.05*height;
            h=bottomsidewallheight;%height/3; % height of right surface
            v1=[0.5* length;0.5* width; 0.5* h];
            v2=[-0.5* length;0.5* width;0.5* h];
            v3=[-0.5* length;0.5* width;-(0.5)* h+r];
            v4=[0.5* length;0.5* width;-(0.5)* h+r];
            x= (0.5)* width-r*(1-(cos(linspace(0,pi/2,n))));
            y= -(0.5)* h+r*(1-(sin(linspace(0,pi/2,n))));
            fi11= [ones(1,n)*0.5*length;x;y];
            fi12=[-ones(1,n)*0.5*length;x;y];
            fil1=sortrows([fi11  fi12]',3,'descend')';
            % bottom surface
            v5=[-0.5*length; (0.5)* width-r;-0.5* h];
            v6=[0.5* length;(0.5)* width-r;-0.5* h];
            v7=[0.5* length;-(0.5)* width+r;-0.5* h];
            v8=[-0.5* length;-(0.5)* width+r;-0.5* h];
            % left wall
            v9=[0.5* length;-0.5* width; 0.5* h];
            v10=[-0.5* length;-0.5* width;0.5* h];
            v11=[-0.5* length;-0.5* width;-(0.5)* h+r];
            v12=[0.5* length;-0.5* width;-(0.5)* h+r];
            % right fillet
            x= -(0.5)* width+r*(1-(cos(linspace(0,pi/2,n))));
            y= -(0.5)* h+r*(1-(sin(linspace(0,pi/2,n))));
            fi21= [ones(1,n)*0.5*length;x;y];
            fi22=[-ones(1,n)*0.5*length;x;y];
            fil2=sortrows([fi21  fi22]',3)';
            bot1Vertices=[v1,v2,v3,v4];
            bot2Vertices=[fi12,flip(fi11,2)];
            bot3Vertices=[v9,v10,v11,v12];
            bot4Vertices=[fi22,flip(fi21,2)];
            bot5Vertices=[v5,v6,v7,v8];
            bottomVertices=[v1,v2,v3,v4,fil1,v5,v6,v7,v8,fil2,v9,v10,v11,v12];
            %Vertices=[Vertices fil1' fil2']
            bot1Vertices=obj.translate(bot1Vertices,[0 0 -0.5*height+h/2]);
            bot2Vertices=obj.translate(bot2Vertices,[0 0 -0.5*height+h/2]);
            bot3Vertices=obj.translate(bot3Vertices,[0 0 -0.5*height+h/2]);
            bot4Vertices=obj.translate(bot4Vertices,[0 0 -0.5*height+h/2]);
            bot5Vertices=obj.translate(bot5Vertices,[0 0 -0.5*height+h/2]);
            %%
            bottomVertices=obj.translate(bottomVertices,[0 0 -0.5*height+h/2]);
            
            %% roof
            % side right surface
            n=10; % points in the fillet radi
            r=0.05*height;
            h= sidewallheight;%height/3;
            v1=[0.5* length;0.5* width; -0.5* h];
            v2=[-0.5* length;0.5* width;-0.5* h];
            v3=[-0.5* length;0.5* width;(0.5)* h-r];
            v4=[0.5* length;0.5* width;(0.5)* h-r];
            x= (0.5)* width-r*(1-(cos(linspace(0,pi/2,n))));
            y= (0.5)* h-r*(1-(sin(linspace(0,pi/2,n))));
            fi11= [ones(1,n)*0.5*length;x;y];
            fi12=[-ones(1,n)*0.5*length;x;y];
            fil1=sortrows([fi11  fi12]',3,'ascend')';
            % top surface
            v5=[-0.5*length; (0.5)* width-r;0.5* h];
            v6=[0.5* length;(0.5)* width-r;0.5* h];
            v7=[0.5* length;-(0.5)* width+r;0.5* h];
            v8=[-0.5* length;-(0.5)* width+r;0.5* h];
            % left wall
            v9=[0.5* length;-0.5* width; -0.5* h];
            v10=[-0.5* length;-0.5* width;-0.5* h];
            v11=[-0.5* length;-0.5* width;(0.5)* h-r];
            v12=[0.5* length;-0.5* width;(0.5)* h-r];
            % right fillet
            x= -(0.5)* width+r*(1-(cos(linspace(0,pi/2,n))));
            y= (0.5)* h-r*(1-(sin(linspace(0,pi/2,n))));
            fi21= [ones(1,n)*0.5*length;x;y];
            fi22=[-ones(1,n)*0.5*length;x;y];
            fil2=sortrows([fi21  fi22]',3,'descend')';
            rf1Vertices=[v1,v2,v3,v4];
            rf2Vertices=[fi12,flip(fi11,2)];
            rf3Vertices=[v9,v10,v11,v12];
            rf4Vertices=[fi22,flip(fi21,2),v12];
            rf5Vertices=[v5,v6,v7,v8];
            roofVertices=[v1,v2,v3,v4,fil1,v5,v6,v7,v8,fil2,v9,v10,v11,v12];
            rf1Vertices=obj.translate(rf1Vertices,[0 0 0.5*height-h/2]);
            rf2Vertices=obj.translate(rf2Vertices,[0 0 0.5*height-h/2]);
            rf3Vertices=obj.translate(rf3Vertices,[0 0 0.5*height-h/2]);
            rf4Vertices=obj.translate(rf4Vertices,[0 0 0.5*height-h/2]);
            rf5Vertices=obj.translate(rf5Vertices,[0 0 0.5*height-h/2]);
            
            %% left window
            h=windowheight;
            v1=[0.5* length;-0.5* width; 0.5* h];
            v2=[-0.5* length;-0.5* width;0.5* h];
            v3=[-0.5* length;-0.5* width;-0.5* h];
            v4=[0.5* length;-0.5* width;-0.5* h];
            lwindowVertices=[v1,v2,v3,v4];
            lwindowFaces=[1 2 3 4];
            
            
            %'LineStyle','none'
            %% right window
            
            v1=[0.5* length;0.5* width; 0.5* h];
            v2=[-0.5* length;0.5* width;0.5* h];
            v3=[-0.5* length;0.5* width;-0.5* h];
            v4=[0.5* length;0.5* width;-0.5* h];
            rwindowVertices=[v1,v2,v3,v4];
            rwindowFaces=[1 2 3 4];
            
            lwindowVertices=obj.translate(lwindowVertices,[0 0 ...
                0.5*height-(sidewallheight+h/2)]);
            rwindowVertices=obj.translate(rwindowVertices,[0 0 ... 
                0.5*height-(sidewallheight+h/2)]);
            
            floorVertices=obj.translate(floorVertices,position);
            lwindowVertices=obj.translate(lwindowVertices,position);
            rwindowVertices=obj.translate(rwindowVertices,position);
            rf1Vertices=obj.translate(rf1Vertices,position);
            rf2Vertices=obj.translate(rf2Vertices,position);
            rf3Vertices=obj.translate(rf3Vertices,position);
            rf4Vertices=obj.translate(rf4Vertices,position);
            rf5Vertices=obj.translate(rf5Vertices,position);
            bot1Vertices=obj.translate(bot1Vertices,position);
            bot2Vertices=obj.translate(bot2Vertices,position);
            bot3Vertices=obj.translate(bot3Vertices,position);
            bot4Vertices=obj.translate(bot4Vertices,position);
            bot5Vertices=obj.translate(bot5Vertices,position);
            
            
            floor=patch(handle,'Faces', floorfaces,'Vertices', floorVertices','FaceColor', colr); %,'FaceAlpha',alph)
            patch(handle,'Faces', [1:4],'Vertices', rf1Vertices','FaceColor', colr);
            patch(handle,'Faces', [1:20],'Vertices', rf2Vertices','FaceColor', colr);
            patch(handle,'Faces', [1:4],'Vertices', rf3Vertices','FaceColor', colr);
            patch(handle,'Faces', [1:20],'Vertices', rf4Vertices','FaceColor', colr);
            patch(handle,'Faces', [1:4],'Vertices', rf5Vertices','FaceColor', colr);
            alph1=0.8;
            patch(handle,'Faces', [1:4],'Vertices', bot1Vertices','FaceColor', colr,'FaceAlpha',alph1);
            patch(handle,'Faces', [1:20],'Vertices', bot2Vertices','FaceColor', colr,'FaceAlpha',alph1);
            patch(handle,'Faces', [1:4],'Vertices', bot3Vertices','FaceColor', colr,'FaceAlpha',alph1);
            patch(handle,'Faces', [1:20],'Vertices', bot4Vertices','FaceColor', colr,'FaceAlpha',alph1);
            patch(handle,'Faces', [1:4],'Vertices', bot5Vertices','FaceColor', colr,'FaceAlpha',0.5);
            colr=[0 0 1];
            alph=0.2;
            lwindow=patch(handle,'Faces', lwindowFaces,'Vertices', lwindowVertices','FaceColor', colr,'FaceAlpha',alph);
            rwindow=patch(handle,'Faces', rwindowFaces,'Vertices', rwindowVertices','FaceColor', colr,'FaceAlpha',alph);
            
            
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
