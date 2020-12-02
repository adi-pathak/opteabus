classdef vehicleconcept
    properties
        Acquisitioncost %in SGD
        BOMcost % in SGD
        Cd=0.5; %drag coefficient
        Af %Frontal Area
        Cr=0.008; %rolling resistance
        rho=1.225;
        Unladenmass
        Glidermass
        Grossvehiclemass
        Body
        Interior
        Powertrain
        Battery
        Passengercapacity
        %Cooling
        Chassis
        Productionemissions
        Distributionemissions
        EOLemissions
        WTTemissions
        LCAemissions
        TTWemissions
        Energyconsumption
        Emptyenergyconsumption
        Energythroughput
        Lifetime=17; % 17 years
        TCO
        TCO_km
        TCO_passengerkm
        av=1;
        diesel=0;
        Properties
        Groundclearance=150;
        constraints=0;
        Range
        Topspeed
        Acceleration
        Gradeability
    end
    properties (Access = private)
        exchangerate=1.56; % EUR to SGD
        discountfactor=0.03; %3 % Ongel 2018
        hourlywages_assembly=25; % 25 SGD/hr
        GST=0.19; % GST in Germany
        GST_SG=0.07; % GST in SG
        markup=2;
        registration_fee=140; % SGD; Ongel 2018
        additional_registration_fee=0.05; % Ongel 2018
        roadtax=524; % SGD
        daysofoperation=329; %347;
        cleaninginterval=1;
        averagecleaningtimeperbus=0.25; % hours
        labourcosts=6.9; % SGD per hour
        electricitycost=0.2065; %SGD/kWh
    end
    methods
        function Vehicle=vehicleconcept(vehicle_parameters)
            
            numberdecks=vehicle_parameters(1);
            interiorlayout=vehicle_parameters(2);
            length=vehicle_parameters(3);
            width=vehicle_parameters(4);
            height=vehicle_parameters(5);
            wheelbase=vehicle_parameters(6);
            %gangwaywidth=vehicle_parameters(7);
            seatpitch=vehicle_parameters(7);
            seatwidth=vehicle_parameters(8);
            standingspace=vehicle_parameters(9);
            powertraintopology=vehicle_parameters(10);
            totalpower=vehicle_parameters(11);
            powersplit=vehicle_parameters(12);
            gearratio=vehicle_parameters(13);
            batterycapacity=vehicle_parameters(14);
            
            
            Vehicle.Body=body(length,width,height,wheelbase,numberdecks);
            Vehicle.Interior=interior(interiorlayout,standingspace,seatwidth,seatpitch,0,Vehicle);%,wheelchairzones)
            Vehicle.Passengercapacity=Vehicle.Interior.passengercapacity;
            Vehicle.Battery=battery(batterycapacity,Vehicle.Body);
            Vehicle.Powertrain=powertrain(powertraintopology,totalpower,powersplit,gearratio);
            [Chassis,Vehicle,Vehicle.Grossvehiclemass,Vehicle.Unladenmass...
                ,Vehicle.Glidermass]=chassis(Vehicle);
            Vehicle.Chassis=Chassis;
            [Vehicle.Acquisitioncost,Vehicle.BOMcost]=BOM(Vehicle);
            Vehicle.constraints=Vehicle.constraintscheck();
        end
        
        function [acquisitioncost,bomcost]=BOM(Vehicle)
            BOMcost=(Vehicle.Body.cost+Vehicle.Body.superstructurecost+...
                +Vehicle.Chassis.axlecost+Vehicle.Chassis.ladderframecost+Vehicle.Interior.costs+Vehicle.Battery.cost+Vehicle.Powertrain.cost)*Vehicle.exchangerate;
            if Vehicle.Body.length==120003
                glidercost=210000;
                BOMcost=(glidercost+Vehicle.Battery.cost+Vehicle.Powertrain.cost)*Vehicle.exchangerate;
            end
            BOMcost=BOMcost+35000*(Vehicle.av);
            % assembly costs
            man_hours=(-0.0776*(Vehicle.Interior.numberseats^2))+7.4516*Vehicle.Interior.numberseats+21.435; % Ongel 2018 % man hours for assembly
            assemblycost=man_hours*Vehicle.hourlywages_assembly;
            acquisitioncost=((BOMcost+assemblycost)*Vehicle.markup)/(1+Vehicle.GST); % Ongel 2018
            bomcost.powertrain=(Vehicle.Battery.cost+Vehicle.Powertrain.cost)*Vehicle.exchangerate;
            bomcost.chassis=(Vehicle.Chassis.axlecost+Vehicle.Chassis.ladderframecost)*Vehicle.exchangerate;
            bomcost.body=(Vehicle.Body.cost+Vehicle.Body.superstructurecost)*Vehicle.exchangerate;
            bomcost.av=34628;
            bomcost.total=BOMcost;
            bomcost.interior=Vehicle.Interior.costs;
            
            %plot bomcost pie chart 
         %   pie([ bomcost.powertrain, bomcost.chassis, bomcost.body,bomcost.interior,bomcost.av])
          %  legend({'powertrain','chassis','body','interior','av'})
        end
        function [batteryreplacements,TCO]=LCC(Vehicle,dailykm,blocks,obj)
            year0_cost=Vehicle.Acquisitioncost+(Vehicle.Acquisitioncost*Vehicle.GST_SG)...
                +Vehicle.registration_fee+(Vehicle.Acquisitioncost*Vehicle.additional_registration_fee);
            distance_to_failure=Vehicle.Battery.cycles*Vehicle.Battery.capacity/Vehicle.Energythroughput;
            batterylife=distance_to_failure./(mean(dailykm)*Vehicle.daysofoperation); % years
            %  batterylife=mean(batterylife);
            batteryreplacements=(Vehicle.Lifetime/batterylife);
            year=[];
            batteryEOLcost=zeros(1,Vehicle.Lifetime);
            batteryreplacementcost=zeros(1,Vehicle.Lifetime);
            for i=1:floor(batteryreplacements)
                year=[year floor(batterylife*i)];
            end
            personnelcost=(43125*1.8)*ones(1,Vehicle.Lifetime)*(1-Vehicle.av);
            electricitycost=(mean(dailykm)*Vehicle.daysofoperation*Vehicle.Energyconsumption/.95*Vehicle.electricitycost)*ones(1,Vehicle.Lifetime);
            roadtax=Vehicle.roadtax*ones(1,Vehicle.Lifetime);
            maintenancecost=(0.55*mean(dailykm)*Vehicle.daysofoperation*(Vehicle.Body.length/1000)/12*(0.5*(Vehicle.diesel+1)))*ones(1,Vehicle.Lifetime);
            insurancecost=((4395*Vehicle.Acquisitioncost/403340)/(1+Vehicle.av))*ones(1,Vehicle.Lifetime);
            cleaningcost=(Vehicle.daysofoperation*Vehicle.cleaninginterval...
                *Vehicle.averagecleaningtimeperbus*...
                (Vehicle.labourcosts*1)/2)*ones(1,Vehicle.Lifetime);
            batteryEOLcost(year)=-(Vehicle.Battery.cost/1.05)*.50; % assumes 50% resale value of battery pack
            batteryreplacementcost(year)=Vehicle.Battery.cost+1000; % assume cost of battery replacement assembly 1000 SGD
            if exist('nchargers')==0
                nchargers=obj.fleetsize;
                chargingpower=50;
            end
            chargerfixedcost=20000;
            variablechargercosts=455;
            
            chargercosts=zeros(1,Vehicle.Lifetime);
            chargerlife=10; %10years
            chargercosts(1)=nchargers*(chargerfixedcost+(variablechargercosts*chargingpower));
            chargercosts(chargerlife)=nchargers*(chargerfixedcost+(variablechargercosts*chargingpower));
            TCO.total=sum([year0_cost+chargercosts(1) electricitycost+roadtax+maintenancecost+personnelcost+insurancecost+cleaningcost+batteryEOLcost+batteryreplacementcost]./((1+Vehicle.discountfactor).^(0:(Vehicle.Lifetime))));
            
            TCO.yearly.total=([year0_cost+chargercosts(1) electricitycost+roadtax+maintenancecost+personnelcost+insurancecost+cleaningcost+batteryEOLcost+batteryreplacementcost]./((1+Vehicle.discountfactor).^(0:(Vehicle.Lifetime))));
            TCO.year0=year0_cost;
            TCO.yearly.cleaning=cleaningcost./((1+Vehicle.discountfactor).^(1:(Vehicle.Lifetime)));
            TCO.yearly.maintenance=maintenancecost./((1+Vehicle.discountfactor).^(1:(Vehicle.Lifetime)));
            TCO.yearly.EOL=batteryEOLcost./((1+Vehicle.discountfactor).^(1:(Vehicle.Lifetime)));
            TCO.yearly.energy=electricitycost./((1+Vehicle.discountfactor).^(1:(Vehicle.Lifetime)));
            TCO.yearly.personnel=personnelcost./((1+Vehicle.discountfactor).^(1:(Vehicle.Lifetime)));
            TCO.yearly.batteryreplacement=batteryreplacementcost./((1+Vehicle.discountfactor).^(1:(Vehicle.Lifetime)));
            TCO.yearly.insurance=insurancecost./((1+Vehicle.discountfactor).^(1:(Vehicle.Lifetime)));
            TCO.yearly.roadtax=roadtax./((1+Vehicle.discountfactor).^(1:(Vehicle.Lifetime)));
            TCO.cleaning=sum(TCO.yearly.cleaning);
            TCO.maintenance=sum(TCO.yearly.maintenance);
            TCO.EOL=sum(TCO.yearly.EOL);
            TCO.energy=sum(TCO.yearly.energy);
            TCO.personnel=sum(TCO.yearly.personnel);
            TCO.batteryreplacement=sum(TCO.yearly.batteryreplacement);
            TCO.insurance=sum(TCO.yearly.insurance);
            TCO.roadtax=sum(TCO.yearly.roadtax);
            TCO.overallcost=TCO.total*obj.fleetsize;
            TCO.TCO_km=TCO.overallcost/(sum(dailykm)*Vehicle.daysofoperation*Vehicle.Lifetime);
            TCO.passengerkm= (TCO.total*obj.fleetsize)/((obj.dailypassengerkm)*Vehicle.daysofoperation*Vehicle.Lifetime); %SGD per km
        end
        function [LCAemissions,EOLemissions,WTTemissions,distributionemissions...
                productionemissions,gCO2_km,gCO2_passengerkm]=LCemissions(Vehicle,distance,passengerkm,obj);
            Vehicle.Battery.emissions=batteryemissions(Vehicle.Battery);
            glideremissions=(Vehicle.Glidermass*5.29473200045908)+(-0.170205385056812*Vehicle.Glidermass);
            productionemissions.CO2=glideremissions+Vehicle.Powertrain.emissions.CO2+Vehicle.Battery.emissions.CO2;
            distributionemissions.CO2=(0.32854*(Vehicle.Glidermass+Vehicle.Powertrain.mass+(Vehicle.Battery.replacements*Vehicle.Battery.mass)));
            WTTemissions=0.4845515892*distance*Vehicle.Lifetime*Vehicle.daysofoperation*Vehicle.Energyconsumption/.95; %well to tank emissions
            EOLemissions.CO2=(Vehicle.Glidermass*0.90600252)+(Vehicle.Battery.replacements*Vehicle.Battery.mass*1.410834)+...
                (Vehicle.Powertrain.mass*0.052409)+(Vehicle.Glidermass*-2.34086996)+(Vehicle.Battery.replacements*Vehicle.Battery.mass*-7.664564);
            LCAemissions.CO2=WTTemissions+EOLemissions.CO2+distributionemissions.CO2+productionemissions.CO2;
            gCO2_km=(sum(LCAemissions.CO2))/sum(distance*Vehicle.Lifetime*Vehicle.daysofoperation)*1000;
            gCO2_passengerkm=(sum(LCAemissions.CO2)*obj.fleetsize)/(passengerkm*Vehicle.Lifetime*Vehicle.daysofoperation*obj.fleetsize)*1000;
        end
        function package(Vehicle,handle)
            if exist('handle','var')==0
            handle=axes;
            end
            cla(handle)
            axis(handle,'equal')
            Vehicle.Body.plotbody(handle,Vehicle);
            vehiclewidth=Vehicle.Body.width;
            wheelbase=Vehicle.Body.wheelbase;
            groundclearance=Vehicle.Groundclearance;
           
            Vehicle.Chassis.plotchassis(handle,wheelbase,vehiclewidth,groundclearance);
            Vehicle.Battery.plotbattery([],handle);
            Vehicle.Interior.plotinterior(Vehicle,handle);
            Vehicle.Powertrain.RA1.motor.plotmotor([-wheelbase/2 0 (Vehicle.Chassis.tyrediameter/2-groundclearance)-55],handle);
            view(handle,[180 0]);
        end
        function flag=constraintscheck(Vehicle)
            % constraint 1 - battery capacity should be less than max
            % battery
            % groundclearance <250 & >150
            %floorheight <40 >30
            %wheelbase/2+foverhang<?
            if Vehicle.Battery.maxcapacity<=Vehicle.Battery.capacity
                flag=-1;
                disp('battery too big')
                return
            
            % constraint2 wheelbase should be greater than 0.5 *length
            elseif Vehicle.Body.wheelbase<0.4*Vehicle.Body.length
                flag=-1;
                disp('wheelbase too small')
                return
            
            % constraint 3 wheelbase +tire diameter<vehicle length
            elseif (5+Vehicle.Body.wheelbase+(1.05*Vehicle.Chassis.tyrediameter))>Vehicle.Body.length
            flag=-1;
            disp('wheelbase too big')
            return
            elseif Vehicle.Interior.seatpitch<650 % seat pitch cannot be less than 650
                flag=-1;
                disp('legroom too small')
                return
            elseif Vehicle.Interior.aislewidth<350 % aisle width cannot be less than 400
                flag=-1;
                disp('Aisle width too small')
                return
            else
                flag=0;
            end
            
            
            
        end
      
        
        function Vehicle=PropertyEvaluation(Vehicle,OperatingPerformance)
            % Evaluate the properties
            
            %% Usability
            % Road Damage
            AL = [300 890 4450 6230 8000 8900 13340];
            RD = [0 0.0003 0.118 0.399 1 1.4 7.9];
            ESAL = interp1(AL, RD, Vehicle.Chassis.axleload);
            Phy.rd = [9e9 10 7.5 0];
            Not.rd = [0 0 5 10];
            Property.RoadDamage.Score = interp1(Phy.rd, Not.rd, ESAL);
            Property.RoadDamage.Value=ESAL;
            %Target(1) = MinReq(1);
            % Space Usage
            Phy.rsu = [9e9 6 5.2 2.6 0];
            Not.rsu = [0 0 5 10 10];
            rsu = (Vehicle.Body.width/1000 *...
                Vehicle.Body.length/1000)/ (Vehicle.Passengercapacity * 0.17);
            Property.SpaceUsage.Score = interp1(Phy.rsu, Not.rsu, rsu);
            Property.SpaceUsage.Value=rsu;
            % Target(3)=MinReq(1);
            % Turning Radius
            if (Vehicle.Chassis.numberofaxles == 2) && (Vehicle.Chassis.RA.RearWheelSteering == 1) % 4-WheelSteering
                tr = (Vehicle.Body.wheelbase/2000)^2 + ((Vehicle.Body.wheelbase/1000)*...
                    ((cotd(46)+cotd(35))/2))^2;
            else
                tr = (Vehicle.Body.wheelbase/1000)/sind(22);
            end
            Phy.tr = [9e9 18.99 19 13 0];
            Not.tr =[0 0 5 10 10];
            Property.TurningRadius.Score = interp1(Phy.tr, Not.tr, tr);
            Property.TurningRadius.Value=tr;
            %Target(2)=MinReq(1);
            %% Design
            % Aspect Ratio
            Phy.ar = [0 0.57 0.75 9e9];
            Not.ar = [0 5 10 10];
            ar =(Vehicle.Body.width/1000 *Vehicle.Body.height/1000);
            Property.AspectRatio.Score = interp1(Phy.ar, Not.ar, ar);
            Property.AspectRatio.Value=ar;
            %Target(4)=MinReq(2);
            
            % Wheelbase Ratio
            Phy.wr = [0 0.4 0.6 9e9];
            Not.wr = [0 5 10 10];
            wr = Vehicle.Body.wheelbase/Vehicle.Body.length;
            Property.WheelBaseRatio.Score = interp1(Phy.wr, Not.wr, wr);
            Property.WheelBaseRatio.Value=wr;
            %  Target(5)=MinReq(2);
            
            
            %% Longitudinal Dynamics
            % Top Speed
            Phy.vmax = [0 59 60 90 9e9];
            Not.vmax = [0 0 5 10 10];
            Property.TopSpeed.Score = interp1(Phy.vmax, Not.vmax, Vehicle.Properties.vmax);
            Property.TopSpeed.Value=Vehicle.Properties.vmax;
            %    Target(6)=MinReq(3);
            
            % Acceleration
            Phy.amax = [0 0.6 1.5 1.75 9e9];
            Not.amax = [0 5 10 10 10];
            Property.Acceleration.Score = interp1(Phy.amax, Not.amax, Vehicle.Properties.amax);
            Property.Acceleration.Value=Vehicle.Properties.amax;
            % Target(7)=MinReq(3);
            
            % Gradeability
            Phy.grade = [0 9.9 10 20 9e9];
            Not.grade = [0 0 5 10 10];
            Property.Gradeability.Score = interp1(Phy.grade, Not.grade, Vehicle.Properties.gradeability);
            Property.Gradeability.Value=Vehicle.Properties.gradeability;
            %   Target(8)=MinReq(3);
            %% Comfort
            % Crowding
            %   Prop(9) = Basisfahrzeug.interior.CrowdingComfort;
            %  Target(9)=MinReq(4);
            
            % Seat comfort
            Phy.sw = [0 0.39 0.4 0.45 9e9];
            Not.sw = [0 0 5 10 10];
            sw = Vehicle.Interior.seatwidth/1000;
            Property.SeatWidth.Score = interp1(Phy.sw, Not.sw, sw);
            Property.SeatWidth.Value=sw;
            % Target(11)=MinReq(4);
            % Leg room
            Phy.sp = [0 0.64 0.65 0.71 0.76 0.81 0.91 9e9];
            Not.sp = [0 0 5 5.8 7.2 8.8 10 10];
            sp = Vehicle.Interior.seatpitch;
            Property.Legroom.Score = interp1(Phy.sp, Not.sp, sp);
            Property.Legroom.Value=sp;
            % Target(12)=MinReq(4);
            % Standing space
            
            
            % climate comfort
            Property.Climatecomfort.Value=010;
            Property.Ridecomfort.Value=10;
            %% Availability
            % Passenger Capacity
            
            % Waiting Time
            
            % Travel Time
            %            Phy.dwt = [9e9 61 60 20 0];
            %            Not.dwt = [0 0 5 10 10];
            %            Prop(15) = interp1(Phy.dwt, Not.dwt, Basisfahrzeug.interior.dwellingtime);
            %            Target(15)=MinReq(6);
            % Frequency
            Phy.phf = [0 3 6 9e9];
            Not.phf = [0 5 10 10];
            phf = mean(mean(60./[OperatingPerformance.lines(:).headway]));
            Property.Frequency.Score = interp1(Phy.phf, Not.phf, phf);
            Property.Frequency.Value=phf;
            %Target(13)=MinReq(5);
            % Seat availability
            
            
            %% Accessibility
            % Low entrance
            
            % Number of Doors
            
            % Number of Decks
            
            % Wheelchair zones
            Phy.wcu = [0 1 2 9e9];
            Not.wcu = [0 5 10 10];
            
            Property.WheelChairZones.Score = interp1(Phy.wcu, Not.wcu, Vehicle.Interior.wheelchairzones);
            Property.WheelChairZones.Value=Vehicle.Interior.wheelchairzones;
            %Target(16)=MinReq(6);
            Vehicle.Properties=Property;
        end
        
    end
end