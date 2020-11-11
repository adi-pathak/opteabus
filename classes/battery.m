classdef battery
    properties
        parallelcells
        seriescells
        mass %in kg
        cost %battery acquisition cost
        capacity %in kwh
        emissions
        replacements
        cycles=3000; % cycle life
        nominalcellvoltage=4.2; % 3V for calculation
        motorvoltage=450; %voltage of the motor
        life
        maxcapacity %maximum energy capacity that can be stored in the available space
        numberofpacks
        width
        height
        length
        position=1;
    end
    properties (Constant)
        C_ratemax=3; %Max discharge rate
        C_ratemaxc=2; %max charge rate
        celltemperature=298.15; % celltemperature in K
        cellcapacity=2.6; %2.05; %Ah
        
        specificbatteryweight=151; % Wh/kg https://www.researchgate.net/publication/323486141_Status-Elektromobilitaet-2018-HL
        volumetricenergydensity=300;%200;%170; %230-300 Wh/L - prottera,akasol
        specificbatterycost=148; % EUR/kwh - Fries cost paper ftm
         end
    
    
    methods
        function obj=battery(batterycapacity,body)
            
            obj.capacity=batterycapacity;
            obj.seriescells=ceil(obj.motorvoltage/obj.nominalcellvoltage);  %number of series 4.2V cells to achieve 420V for motor
            obj.parallelcells=ceil(batterycapacity*1000/(obj.nominalcellvoltage*obj.cellcapacity*obj.seriescells));
            obj.mass=(batterycapacity/(obj.specificbatteryweight))*1000;
            obj.cost=obj.specificbatterycost*obj.capacity;
            
        end
        
        function obj=maxstoragecapacity(obj,wheelbase,bodywidth,wheelhousingwidth,...
               rearoverhang,tyrediameter,floorheight,sectionwidth,sectionheight,groundclearance)
            %maxheight=175;
            if obj.position==1
            maxlength=wheelbase-1.5*tyrediameter;
            maxwidth= bodywidth-2*wheelhousingwidth-2*sectionwidth-5; % section height and width are equal
            maxheight=floorheight-groundclearance-30; % 30 mm for ramp
            elseif obj.position==2
            maxlength=wheelbase-1.5*tyrediameter;
            maxwidth= bodywidth-wheelhousingwidth;
            maxheight=floorheight-sectionheight-groundclearance-10;
            end
            maxvolume=maxlength/1000*maxwidth/1000*maxheight/1000 *1000;
            obj.maxcapacity=maxvolume*obj.volumetricenergydensity/1000;
            obj.width=maxwidth;
            obj.height=maxheight;
            
            % if battery capacity exceeds and rear overhang has space for
            % 2nd pack
            if obj.capacity>obj.maxcapacity & rearoverhang>1.5*tyrediameter
                obj.numberofpacks=2;
                %check maximum space for 2nd pack
                maxlength=rearoverhang-1.5*tyrediameter-200;
                maxwidth= bodywidth-2*wheelhousingwidth;
                maxheight=floorheight-sectionheight-groundclearance-5;
                maxvolume=maxlength/1000*maxwidth/1000*maxheight/1000 *1000;
                obj.maxcapacity=maxvolume*obj.volumetricenergydensity/1000;
            
            end
                
        end
        function emissions=batteryemissions(obj)
            emissions.CO2=(obj.mass*obj.replacements*.63*27.4137)+(obj.mass*obj.replacements*.37*24.5588); % CO2 emissions of cell+housing %24.79118716 %21.937479
            emissions.DCB=(obj.mass*obj.replacements*.63*29.11545727)+(obj.mass*obj.batterylife*.37*26.27587); % DCB emissions of cell+housing
            emissions.PM10=(obj.mass*obj.replacements*.63*0.1269670696)+(obj.mass*obj.batterylife*.37*0.050547); % PM-10 
        end
        function plotbattery(obj,position,handle)
          %  calculate the battery pack dimensions
            if ~isempty(obj.width)
                if obj.maxcapacity<obj.capacity
                    disp('Battery capacity exceeds available space')
                  return  
                end
                if obj.position==1
                      packwidth=obj.width;
                packheight=obj.height;
                packlength=(obj.capacity*1000/obj.volumetricenergydensity)/(packwidth/1000*packheight/1000);
              position=[0 0 packheight/2+2]; % 2mm tolerance
                     
                else
              
                packwidth=obj.width;
                packheight=obj.height;
                packlength=(obj.capacity*1000/obj.volumetricenergydensity)/(packwidth/1000*packheight/1000);
                  position=[0 0 packheight/2+60]; % 60 is frameheight (battery on top of frame)
                end
            else
            cellheight=65; %18650 cell
            cellwidth=18; %diameter 18650
            spacebetweencells=2; %mm
            bmsthickness=10; %mm
            
            packheight=cellheight+bmsthickness;
            packwidth=obj.parallelcells*(cellwidth+spacebetweencells);
            packlength=obj.seriescells*(cellwidth+spacebetweencells);
            if packlength<packwidth %rotate pack 
               l=packlength;
               packlength=packwidth;
               packwidth=l;
            end
            end
            % define vertices
            v1=[-0.5*packlength;0.5*packwidth;-0.5*packheight];
            v2=[0.5*packlength;0.5*packwidth;-0.5*packheight];
            v3=[0.5*packlength;0.5*packwidth;0.5*packheight];
            v4=[-0.5*packlength;0.5*packwidth;0.5*packheight];
            v5=[-0.5*packlength;-0.5*packwidth;0.5*packheight];
            v6=[0.5*packlength;-0.5*packwidth;0.5*packheight];
            v7=[0.5*packlength;-0.5*packwidth;-0.5*packheight];
            v8=[-0.5*packlength;-0.5*packwidth;-0.5*packheight];
            Vertices=[v1,v2,v3,v4,v5,v6,v7,v8];
           
            Vertices=obj.translate(Vertices,position);
            colr=[1 0 0];
            patch(handle,'Faces',[1 2 3 4
                5 6 7 8
                4 3 6 5
                3 2 7 6
                2 1 8 7
                1 4 5 8],'Vertices',Vertices','FaceColor',colr);
%             text(handle,position(1),position(2),position(3),'Battery');
        end
        function vertices=translate(obj,vertices,position)
            T= [1 0 0  position(1);0 1 0  position(2);0 0 1  position(3);0 0 0 1] ;
            vertices(4,:)=1;
            vertices=(T*vertices);
            vertices=vertices(1:3,:);
        end
    end
    
end