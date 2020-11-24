classdef chassis
    properties
        axleload
        numberofaxles
        axlecost
        ladderframemass
        ladderframecost
        mass
        FA
        RA
        TA
        airspring
        brakes
        steering
        ESP
        tyrewidth
        tyrediameter
        tyreinnerdiameter
        rimsize
        airspringdiameter=300;
        airspringheight=285;
        airtankdiameter=210;
        airtanklength=450;
    end
    properties (Access=private)
        passengermass=70; %kg
        sectionheight=60; %100 mm square section
        sectionwidth=120; %square section for ladder frame; % source at.govt.nzmedia/1976524/enviro200ev-brochure.pdf
        sectionthickness=6; % 6mm thickness
        steeldensity=7750; %kg/m3 density of steel used in ladder frame
        steelrawprice=6.67; %EUR/kg
        structure_materialutilisation=0.52; %production cost Fuchs
        
    end
    methods
        function [obj,Vehicle,gvm,unladenmass,glidermass]=chassis(Vehicle)
            battery=Vehicle.Battery;
            powertrain=Vehicle.Powertrain;
            body=Vehicle.Body;
            interior=Vehicle.Interior;
            passengermass=interior.passengercapacity*obj.passengermass;
            [framemass,framecost]=ladderframe(obj,body.length,body.width);
            obj.ladderframemass=framemass;
            obj.ladderframecost=framecost;
            sprungmass = (battery.mass + powertrain.mass +...
                body.mass + body.superstructuremass+interior.mass +...
                passengermass+framemass);                              % Total Axle Load
            
            
            if  (sprungmass>23500)             % Add a third axle if load exceeds
                obj.numberofaxles=3;
            else
                obj.numberofaxles=2;
            end
            
            AxleLoad = sprungmass/obj.numberofaxles; % Axle Load for chassis decision assuming 50% distribution
            obj.axleload=AxleLoad;
            %% Axle correlation (Based on the ZF Database)
            % Front axle selection
            obj=obj.ZFchassis(string('FA'),AxleLoad);
            
            % Rear axle selection
            obj=obj.ZFchassis('RA',AxleLoad);
            
            % Third axle selection (if applicable)
            if obj.numberofaxles==3
                obj=obj.ZFchassis('TA',AxleLoad);
                
                obj.mass= obj.TA.mass + obj.FA.mass + obj.RA.mass;
            else
                obj.mass= obj.FA.mass + obj.RA.mass;
            end
            obj = obj.axlecosts(AxleLoad);% Chassis cost
            glidermass = (body.mass + body.superstructuremass+ interior.mass +...
                obj.mass+ powertrain.mass+framemass);                              % glider mass
            unladenmass=(battery.mass + powertrain.mass +...
                obj.mass+  body.mass + body.superstructuremass+ interior.mass+framemass);
            gvm=obj.mass+sprungmass;
            obj=tyresizing(obj,unladenmass,gvm);
           % obj=airspringsizing(obj,unladenmass,gvm);
           % add tire mass to gvm
           %% update body dimensions - tire size,front and rear overhang
           Vehicle.Body=Vehicle.Body.update_body(obj.tyrediameter,obj.tyrewidth...
               ,obj.airspringdiameter);
           Vehicle.Battery=Vehicle.Battery.maxstoragecapacity(body.wheelbase,body.width,Vehicle.Body.wheelhousingwidth,...
               Vehicle.Body.rearoverhang,obj.tyrediameter,interior.floorheight,body.sectionwidth,...
               body.sectionheight, body.groundclearance);
        
        end
        function cog(body,chassis,battery,powertrain)
            
            
        end
        function obj=ZFchassis(obj,type,load)
            if strcmp(type,'FA')
                
                if load < 5800
                    obj.FA.name = 'RL 55 EC';
                    obj.FA.maxLoad = 5800;
                    obj.FA.maxSteering = 55;
                    obj.FA.mass = 425;
                    obj.FA.RearWheelSteering = 1;
                    obj.FA.wheels=2; %number of wheels
                elseif load < 7500
                    obj.FA.name = 'RL 75 A';
                    obj.FA.maxLoad = 7500;
                    obj.FA.maxSteering = 60;
                    obj.FA.mass = 527;
                    obj.FA.wheels=2; %number of wheels
                elseif load < 8500
                    obj.RA.name = 'RL 82 EC';
                    obj.RA.maxLoad = 8200;
                    obj.RA.maxSteering = 56;
                    obj.RA.mass = 482;
                    obj.RA.RearWheelSteering = 1; 
                    obj.FA.wheels=2; %number of wheels
                else
                    obj.FA.name = 'RL 82 A';
                    obj.FA.maxLoad = 8200;
                    obj.FA.maxSteering = 55;
                    obj.FA.mass = 527;
                    obj.FA.wheels=2; %number of wheels
                end
            elseif strcmp(type, 'TA')
                if load < 7500
                    obj.TA.name = 'RL 75 A';
                    obj.TA.maxLoad = 7500;
                    obj.TA.maxSteering = 60;
                    obj.TA.mass = 527;
                else
                    obj.TA.name = 'RL 82 A';
                    obj.TA.maxLoad = 8200;
                    obj.TA.maxSteering = 55;
                    obj.TA.mass = 527;
                end
            else
                if load < 5800
                    obj.RA.name = 'RL 55 EC';
                    obj.RA.maxLoad = 5800;
                    obj.RA.maxSteering = 55;
                    obj.RA.mass = 425;
                    obj.RA.RearWheelSteering = 1;
                    obj.RA.wheels=2; %number of wheels
                elseif load < 7500
                    obj.RA.name = 'RL 75 E';
                    obj.RA.maxLoad = 7500;
                    obj.RA.maxSteering = 60;
                    obj.RA.mass = 496;
                    obj.RA.RearWheelSteering = 1;
                    obj.RA.wheels=4; %number of wheels
                elseif load < 13000
                    obj.RA.name = 'RL 82 EC';
                    obj.RA.maxLoad = 8200;
                    obj.RA.maxSteering = 56;
                    obj.RA.mass = 482;
                    obj.RA.RearWheelSteering = 1;
                    obj.RA.wheels=4; %number of wheels
                elseif load < 8200
                end
            end
        end
        function obj=axlecosts(obj,axleload)
            if axleload>8200
                obj.RA.cost=(obj.numberofaxles-1)*obj.RA.mass*((0.9021*...
                    6.67) +... %Materialkosten.Stahl.Schmide_Fahrwerk.Preis_Verarbeitet
                    (0.0818 * 3.5) + (0.0038 *...%Materialkosten.Aluminium.Rohpreis
                    2) + (0.0105 *... %Materialkosten.Kunststoff.Elastomere.Rohpreis
                    4)); %Materialkosten.Kunststoff.PA.Rohpreis
            else
                if obj.numberofaxles==3
                    obj.TA.cost = obj.TA.mass *...
                        ((0.9021*...
                        6.67) +... %Materialkosten.Stahl.Schmide_Fahrwerk.Preis_Verarbeitet
                        (0.0818 * 3.5) + (0.0038 *...%Materialkosten.Aluminium.Rohpreis
                        2) + (0.0105 *... %Materialkosten.Kunststoff.Elastomere.Rohpreis
                        4)); %Materialkosten.Kunststoff.PA.Rohpreis
                    obj.RA.cost = (obj.RA.mass *...
                        ((0.8331 * 6.67) +... %Materialkosten.Stahl.Schmide_Fahrwerk.Preis_Verarbeitet
                        (0.1425 * 3.5) +... %Materialkosten.Aluminium.Rohpreis
                        (0.019 * 2))); %Materialkosten.Kunststoff.Elastomere.Rohpreis
                else
                    obj.RA.cost= (obj.RA.mass *...
                        ((0.8331 * 6.67)+  ... %Materialkosten.Stahl.Schmide_Fahrwerk.Preis_Verarbeitet
                        (0.1425 *3.5) ...  % Materialkosten.Aluminium.Rohpreis
                        + (0.019 *2))); % Materialkosten.Kunststoff.Elastomere.Rohpreis
                end
            end
            obj.FA.cost = obj.FA.mass *...
                ((0.8331 * 6.67) +...
                (0.1425 * 3.5) +...
                (0.019 * 2));
            
            extracost = obj.numberofaxles * 2 * 2700;
            
            if obj.numberofaxles>2
                obj.axlecost=obj.FA.cost + obj.RA.cost+obj.TA.cost+ extracost;
            else
                obj.axlecost=obj.FA.cost + obj.RA.cost+ extracost;
            end
        end
        function obj=tyresizing(obj,unladenmass,grossvehiclemass)
            %% Mass configurations vehicle (unladen, unladen +25% payload, unladen +50% payload, Fully laden, Laden + 25% Payload)
            % FAWL: Front Axle Weight Laden
            % FAWUL: Front Axle Weight Unladen
            % RAWL: Rear Axle Weight Laden
            % RAWUL: Rear Axle Weight Unladen
            % in kg
            
            FAWL = grossvehiclemass/2; % assuming 50-50 weight distribution
            FAWUL = unladenmass/2;
            
            RAWL = grossvehiclemass/2;
            RAWUL = unladenmass/2;
            
            
            % Pay Load (Front: Fpl || Rear: Rpl) starts with Unladen to Fully Laden*n in j steps
            j = 0.25;
            n = 1.25;
            
            Fpl = struct('Loads_Front', FAWL - FAWUL);
            Rpl = struct('Loads_Rear', RAWL - FAWUL);
            i=2;
            while j<=n
                Fpl(i).Loads_Front = FAWUL + j*Fpl(1).Loads_Front;
                Rpl(i).Loads_Rear = RAWUL + j*Rpl(1).Loads_Rear;
                j=0.25;
                j=j*i;
                i=i+1;
            end
            
            % Unsprung mass FVSM Front and RVSM Rear
            FUSM = 100; %!!!Hardcoded!!!
            RUSM = 100; %!!!Hardcoded!!!
            
            %Sprung mass FSM Front / RSM Rear per wheel
            FSM = struct('Sprung_Loads_Front', (Fpl(1).Loads_Front - FUSM)/2);
            RSM = struct('Sprung_Loads_Rear', (Rpl(1).Loads_Rear - RUSM)/2);
            i=2;
            while i<= length(Fpl(1,:))
                FSM(i).Sprung_Loads_Front = ((Fpl(i).Loads_Front - FUSM)/2); % in kg
                RSM(i).Sprung_Loads_Rear = ((Rpl(i).Loads_Rear - RUSM)/2); % in kg
                i=i+1;
            end
            
            %Index for fully laden +25%
            p = length(FSM(1,:));
            
            %% Tyre Selection
            % Load Data from Database
           load(fullfile(pwd,filesep,'classes',(filesep),'Tyredata.mat'));
            % 50:50 axle load!
            
            i =1;
            while i<=size(Tyre,1)
                
                if (Tyre.load_min_P(i) > FSM(p).Sprung_Loads_Front)
                    load_selection(i) = table2struct(Tyre(i,:));
                end
                
                i = i+1;
            end
            
            % Delete empty rows from load_selection
            i = 1;
            
            while i<=length(load_selection)
                if isempty(load_selection(i).tyre_id)
                    load_selection(i) = [];
                    i = 1;
                    
                else i = i+1;
                end
                
            end
            
            % Smallest Rim Size selection
            
            Afields = fieldnames(load_selection);
            Acell = struct2cell(load_selection);
            sz = size(Acell);
            
            
            % Convert to a matrix
            Acell = reshape(Acell, sz(1), []);      % Px(MxN)
            
            % Make each field a column
            Acell = Acell';                         % (MxN)xP
            
            % Sort by field 5 "rim_size"
            Acell = sortrows(Acell, 5);
            
            % Put back into original cell array format
            Acell = reshape(Acell', sz);
            
            % Convert to Struct
            load_selection = cell2struct(Acell, Afields, 1);
            
            % Select all Tyres with smallest Rim Size
            
            % Delete rows with bigger Tyre Sizes from load_selection
            min_rim_size = load_selection(1).rim_size;
            i = 2;
            
            while i<=length(load_selection)
                if load_selection(i).rim_size > min_rim_size
                    load_selection(i) = [];
                    i = 1;
                end
                i = i+1;
            end
            
            % Smallest Tyre Width Selection
            Afields = fieldnames(load_selection);
            Acell = struct2cell(load_selection);
            sz = size(Acell);
            
            
            % Convert to a matrix
            Acell = reshape(Acell, sz(1), []);      % Px(MxN)
            
            % Make each field a column
            Acell = Acell';                         % (MxN)xP
            
            % Sort by field 5 "width"
            Acell = sortrows(Acell, 2);
            
            % Put back into original cell array format
            Acell = reshape(Acell', sz);
            
            % Convert to Struct
            load_selection = cell2struct(Acell, Afields, 1);
            
            Tyre_selection = load_selection(1);
            obj.tyrewidth=Tyre_selection.width;
            obj.tyrediameter=Tyre_selection.outside_diameter;
            obj.tyreinnerdiameter=obj.tyrediameter-(obj.tyrewidth*Tyre_selection.aspect_ratio/100)*2;
            obj.rimsize=Tyre_selection.rim_size*25.4;
            
        end
        function obj=airspringsizing(obj,unladenmass,grossvehiclemass)
            %% Mass configurations vehicle (unladen, unladen +25% payload, unladen +50% payload, Fully laden, Laden + 25% Payload)
            % FAWL: Front Axle Weight Laden
            % FAWUL: Front Axle Weight Unladen
            % RAWL: Rear Axle Weight Laden
            % RAWUL: Rear Axle Weight Unladen
            
            FAWL = grossvehiclemass/2;
            FAWUL = unladenmass/2;
            
            RAWL =  grossvehiclemass/2;
            RAWUL =unladenmass/2;
          
            
            
            % Pay Load (Front: Fpl || Rear: Rpl) starts with Unladen to Fully Laden*n in j steps
            j = 0.25;
            n = 1.25;
            
            Fpl = struct('Loads_Front', FAWL - FAWUL);
            Rpl = struct('Loads_Rear', RAWL - FAWUL);
            i=2;
            while j<=n
                Fpl(i).Loads_Front = FAWUL + j*Fpl(1).Loads_Front;
                Rpl(i).Loads_Rear = RAWUL + j*Rpl(1).Loads_Rear;
                j=0.25;
                j=j*i;
                i=i+1;
            end
            
            % Unsprung mass FVSM Front and RVSM Rear
            FUSM = 165; %!!!Hardcoded!!!
            RUSM = 165; %!!!Hardcoded!!!
            
            %Sprung mass FSM Front / RSM Rear per wheel
            FSM = struct('Sprung_Loads_Front', (Fpl(1).Loads_Front - FUSM)/2);
            RSM = struct('Sprung_Loads_Rear', (Rpl(1).Loads_Rear - RUSM)/2);
            i=2;
            while i<= length(Fpl(1,:))
                FSM(i).Sprung_Loads_Front = ((Fpl(i).Loads_Front - FUSM)/2)*2.20462; % in lbs
                RSM(i).Sprung_Loads_Rear = ((Rpl(i).Loads_Rear - RUSM)/2)*2.20462; % in lbs
                i=i+1;
            end
            
            %Index for fully laden +25%
            p = length(FSM(1,:));
            
            % Load Data from Database
            %DBConnection; % calls script for the data base connection
            load(strcat(pwd,'\Classes\Airspringdata.mat'));
            AirSpring=table2struct(AirSpring);
            %###Suspension Selection###%
            % Selection criteria:
            % 1. Load
            % 2. Suspension travel min. 180mm (7.07in)
            % 3. Isolation effectivness
            % 4. Best natural frequency (close to 1Hz)
            
            % Front = Rear; #calculation for fully laden + 25%#
            
            i = 1; % 1. and 2. Loop Load selection @ 40 PSIG, @ 60 PSIG , @ 80 PSIG
            
            
            while i<=length(AirSpring)
                
                if (AirSpring(i).load_40psig > FSM(p).Sprung_Loads_Front) && ((AirSpring(i).design_height - AirSpring(i).min_height)>= 7.07)
                    load_selection(i) = AirSpring(i);
                    
                elseif (AirSpring(i).load_60psig > FSM(p).Sprung_Loads_Front) && ((AirSpring(i).design_height - AirSpring(i).min_height)>= 7.07)
                    load_selection(i) = AirSpring(i);
                    
                elseif (AirSpring(i).load_80psig > FSM(p).Sprung_Loads_Front) && ((AirSpring(i).design_height - AirSpring(i).min_height)>= 7.07)
                    load_selection(i) = AirSpring(i);
                end
                
                i = i+1;
            end
            
            % Delete empty rows from load_selection
            i = 1;
            
            while i<=length(load_selection)
                if isempty(load_selection(i).style_number)
                    load_selection(i) = [];
                    i = 1;
                    
                else i = i+1;
                end
                
            end
            
            
            % Overall Load smallest spring selection
            
            Afields = fieldnames(load_selection);
            Acell = struct2cell(load_selection);
            sz = size(Acell);
            
            
            % Convert to a matrix
            Acell = reshape(Acell, sz(1), []);      % Px(MxN)
            
            % Make each field a column
            Acell = Acell';                         % (MxN)xP
            
            % Sort by field 14 "Max_Diameter"
            Acell = sortrows(Acell, 14);
            
            % Put back into original cell array format
            Acell = reshape(Acell', sz);
            
            % Convert to Struct
            load_selection = cell2struct(Acell, Afields, 1);
            
            % Select first 3 elements
            if length(load_selection(1,:))<3
                
            else
                load_selection = [load_selection(1),load_selection(2),load_selection(3)];
            end
            
            % 3. Isolation effectivness
            Iso_eff_selection = load_selection;
            i = 1;
            while i<=length(load_selection)
                
                Iso_eff_400 = 100 - load_selection(i).isolation_400cpm; % Isolation effectivness @ 400 cpm
                Iso_eff_800 = 100 - load_selection(i).isolation_800cpm; % Isolation effectivness @ 800 cpm
                Iso_eff_1500 = 100 - load_selection(i).isolation_1500cpm; % Isolation effectivness @ 1500 cpm
                
                % Average percentage of vibriation transmission - small value is
                % better -> less vibration transmission. Vehicle related vibrations < 30Hz
                % (1680 cpm)
                
                Iso_eff_average = (Iso_eff_400 + Iso_eff_800 + Iso_eff_1500)/3; % Average Isolation effectivness from 400 cpm to 1500 cpm
                
                Iso_eff_selection(i).Isolation_effectivness_average = Iso_eff_average;
                i = i+1;
            end
            
            % Isolation effectivness smallest value selection
            if length(Iso_eff_selection)>1 % if length is greater than 1 select smallest 2
            Afields = fieldnames(Iso_eff_selection);
            Acell = struct2cell(Iso_eff_selection);
            sz = size(Acell);
            
            
            % Convert to a matrix
            Acell = reshape(Acell, sz(1), []);      % Px(MxN)
            
            % Make each field a column
            Acell = Acell';                         % (MxN)xP
            
            % Sort by field 15 "Isolation_effectivness_average"
            Acell = sortrows(Acell, 15);
            
            % Put back into original cell array format
            Acell = reshape(Acell', sz);
            
            % Convert to Struct
            Iso_eff_selection = cell2struct(Acell, Afields, 1);
            
            % Select first 2 elements
            if length(Iso_eff_selection(:))<2
                
            else
                Iso_eff_selection = [Iso_eff_selection(1),Iso_eff_selection(2)];
            end
            
            % Preselection of lowest possible pressure level
            default_load = [];
            default_pressure = [];
            
            i = 1;
            while i<=length(Iso_eff_selection)
                
                if Iso_eff_selection(i).load_40psig > FSM(p).Sprung_Loads_Front
                    
                    default_load(i) = Iso_eff_selection(i).load_40psig;
                    default_pressure(i) = 40; % 40 PSIG
                    
                elseif Iso_eff_selection(i).load_60psig > FSM(p).Sprung_Loads_Front
                    
                    default_load(i) = Iso_eff_selection(i).load_60psig;
                    default_pressure(i) = 60; % 60 PSIG
                    
                else
                    
                    default_load(i) = Iso_eff_selection(i).load_80psig;
                    default_pressure(i) = 80; % 80 PSIG
                    
                end
                
                i = i+1;
            end
            
            % Calculation of the effective area and required pressure
            
            i = 1;
            while i<=length(default_load)
                
                Iso_eff_selection(i).effective_area_design_height = default_load(i)/default_pressure(i); % calculation of the effective area at design hight with default values from table
                Iso_eff_selection(i).required_pressure = FSM(p).Sprung_Loads_Front/Iso_eff_selection(i).effective_area_design_height; % calculation of the required pressure
                
                i = i+1;
            end
            
            % Calculation of the dynamic spring rate
            p_force_table = []; % Vector to find the pressure value to interpolate with
            
            
            i = 1; % Loop number of springs
            j = 1; % Loop force table
            K_selection = Iso_eff_selection; % Struct for vertical spring rate selection
            
            while i<=length(Iso_eff_selection)
                
                p_force_table = [abs(20 - Iso_eff_selection(i).required_pressure); abs(40 - Iso_eff_selection(i).required_pressure); abs(60 - Iso_eff_selection(i).required_pressure); abs(80 - Iso_eff_selection(i).required_pressure); abs(100 - Iso_eff_selection(i).required_pressure)];
                [p_force_table,I] = min(p_force_table); % select pressure level in force table
                
                height_Ac = Iso_eff_selection(i).design_height - 0.5; % assembly height 1/2 inch below design height
                height_Ae = Iso_eff_selection(i).design_height + 0.5; % assembly height 1/2 inch above design height
                
                trigger_Ac = 0;
                trigger_Ae = 0;
                
                switch  I % Case 1: 20 PSIG Case 2: 40 PSIG Case 3: 60 PSIG Case 4: 80 PSIG Case 5: 100 PSIG
                    
                    
                    case 1 % pressure level @ 20 PSIG
                        
                        while j <= length(Iso_eff_selection(i).ForceTable)% check if height_Ac = assembly height
                            
                            % Height and Volume
                            if height_Ac == Iso_eff_selection(i).ForceTable(j).Assembly_Height % Ac loop
                                
                                Ac = Iso_eff_selection(i).ForceTable(j).Pounds_Force_20_PSIG/20; % Effective Area at 1/2 inch below design height
                                trigger_Ac = 1;
                                
                                Vc = Iso_eff_selection(i).ForceTable(j).Volume_100_PSIG; % Vc at 1/2 inch below design height
                            end
                            
                            if height_Ae == Iso_eff_selection(i).ForceTable(j).Assembly_Height % Ae loop
                                
                                Ae = Iso_eff_selection(i).ForceTable(j).Pounds_Force_20_PSIG/20; % Effective Area at 1/2 inch above design height
                                trigger_Ae = 1;
                                
                                Ve = Iso_eff_selection(i).ForceTable(j).Volume_100_PSIG; % Vc at 1/2 inch above design height
                            end
                            
                            a_Ace(j) = Iso_eff_selection(i).ForceTable(j).Assembly_Height; % Assembly Height values
                            
                            j = j+1;
                        end
                        
                        if trigger_Ac == 0 % if Ac/Ae =~ assembly height --> interpolation
                            
                            % find limits for height and volume interpolation
                            
                            b_Ac = height_Ac;
                            d_Ac = sort(abs(b_Ac - a_Ace));
                            lowest_Ac = find(abs(b_Ac - a_Ace) == d_Ac(1));
                            sec_lowest_Ac = find(abs(b_Ac - a_Ace) == d_Ac(2));
                            
                            % Ac height interpolation
                            Ac = (Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Pounds_Force_20_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Pounds_Force_20_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Pounds_Force_20_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height))*(height_Ac - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height))/20;
                            
                            % Vc volume interpolation
                            Vc = (Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Volume_100_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Volume_100_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Volume_100_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height))*(height_Ac - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height));
                            
                        end
                        
                        if trigger_Ae == 0 % if Ac/Ae =~ assembly height --> interpolation
                            
                            % find limits for height interpolation
                            
                            b_Ae = height_Ae;
                            d_Ae = sort(abs(b_Ae - a_Ace));
                            lowest_Ae = find(abs(b_Ae - a_Ace) == d_Ae(1));
                            sec_lowest_Ae = find(abs(b_Ae - a_Ace) == d_Ae(2));
                            
                            % Ae height interpolation
                            Ae = (Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Pounds_Force_20_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Pounds_Force_20_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Pounds_Force_20_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height))*(height_Ae - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height))/20;
                            
                            % Ve volume interpolation
                            Ve = (Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Volume_100_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Volume_100_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Volume_100_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height))*(height_Ae - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height));
                            
                        end
                        
                        
                    case 2 % pressure level @ 40 PSIG
                        
                        while j <= length(Iso_eff_selection(i).ForceTable)% check if height_Ac = assembly height
                            
                            % Height and Volume
                            if height_Ac == Iso_eff_selection(i).ForceTable(j).Assembly_Height % Ac loop
                                
                                Ac = Iso_eff_selection(i).ForceTable(j).Pounds_Force_40_PSIG/40; % Effective Area at 1/2 inch below design height
                                trigger_Ac = 1;
                                
                                Vc = Iso_eff_selection(i).ForceTable(j).Volume_100_PSIG; % Vc at 1/2 inch below design height
                            end
                            
                            if height_Ae == Iso_eff_selection(i).ForceTable(j).Assembly_Height % Ae loop
                                
                                Ae = Iso_eff_selection(i).ForceTable(j).Pounds_Force_40_PSIG/40; % Effective Area at 1/2 inch above design height
                                trigger_Ae = 1;
                                
                                Ve = Iso_eff_selection(i).ForceTable(j).Volume_100_PSIG; % Vc at 1/2 inch above design height
                            end
                            
                            a_Ace(j) = Iso_eff_selection(i).ForceTable(j).Assembly_Height; % Assembly Height values
                            
                            j = j+1;
                        end
                        
                        if trigger_Ac == 0 % if Ac/Ae =~ assembly height --> interpolation
                            
                            % find limits for height and volume interpolation
                            
                            b_Ac = height_Ac;
                            d_Ac = sort(abs(b_Ac - a_Ace));
                            lowest_Ac = find(abs(b_Ac - a_Ace) == d_Ac(1));
                            sec_lowest_Ac = find(abs(b_Ac - a_Ace) == d_Ac(2));
                            
                            % Ac height interpolation
                            Ac = (Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Pounds_Force_40_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Pounds_Force_40_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Pounds_Force_40_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height))*(height_Ac - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height))/40;
                            
                            % Vc volume interpolation
                            Vc = (Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Volume_100_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Volume_100_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Volume_100_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height))*(height_Ac - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height));
                            
                        end
                        
                        if trigger_Ae == 0 % if Ac/Ae =~ assembly height --> interpolation
                            
                            % find limits for height interpolation
                            
                            b_Ae = height_Ae;
                            d_Ae = sort(abs(b_Ae - a_Ace));
                            lowest_Ae = find(abs(b_Ae - a_Ace) == d_Ae(1));
                            sec_lowest_Ae = find(abs(b_Ae - a_Ace) == d_Ae(2));
                            
                            % Ae height interpolation
                            Ae = (Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Pounds_Force_40_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Pounds_Force_40_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Pounds_Force_40_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height))*(height_Ae - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height))/40;
                            
                            % Ve volume interpolation
                            Ve = (Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Volume_100_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Volume_100_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Volume_100_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height))*(height_Ae - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height));
                            
                        end
                        
                    case 3 % pressure level @ 60 PSIG
                        
                        while j <= length(Iso_eff_selection(i).ForceTable)% check if height_Ac = assembly height
                            
                            % Height and Volume
                            if height_Ac == Iso_eff_selection(i).ForceTable(j).Assembly_Height % Ac loop
                                
                                Ac = Iso_eff_selection(i).ForceTable(j).Pounds_Force_60_PSIG/60; % Effective Area at 1/2 inch below design height
                                trigger_Ac = 1;
                                
                                Vc = Iso_eff_selection(i).ForceTable(j).Volume_100_PSIG; % Vc at 1/2 inch below design height
                            end
                            
                            if height_Ae == Iso_eff_selection(i).ForceTable(j).Assembly_Height % Ae loop
                                
                                Ae = Iso_eff_selection(i).ForceTable(j).Pounds_Force_60_PSIG/60; % Effective Area at 1/2 inch above design height
                                trigger_Ae = 1;
                                
                                Ve = Iso_eff_selection(i).ForceTable(j).Volume_100_PSIG; % Vc at 1/2 inch above design height
                            end
                            
                            a_Ace(j) = Iso_eff_selection(i).ForceTable(j).Assembly_Height; % Assembly Height values
                            
                            j = j+1;
                        end
                        
                        if trigger_Ac == 0 % if Ac/Ae =~ assembly height --> interpolation
                            
                            % find limits for height and volume interpolation
                            
                            b_Ac = height_Ac;
                            d_Ac = sort(abs(b_Ac - a_Ace));
                            lowest_Ac = find(abs(b_Ac - a_Ace) == d_Ac(1));
                            sec_lowest_Ac = find(abs(b_Ac - a_Ace) == d_Ac(2));
                            
                            % Ac height interpolation
                            Ac = (Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Pounds_Force_60_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Pounds_Force_60_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Pounds_Force_60_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height))*(height_Ac - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height))/60;
                            
                            % Vc volume interpolation
                            Vc = (Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Volume_100_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Volume_100_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Volume_100_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height))*(height_Ac - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height));
                            
                        end
                        
                        if trigger_Ae == 0 % if Ac/Ae =~ assembly height --> interpolation
                            
                            % find limits for height interpolation
                            
                            b_Ae = height_Ae;
                            d_Ae = sort(abs(b_Ae - a_Ace));
                            lowest_Ae = find(abs(b_Ae - a_Ace) == d_Ae(1));
                            sec_lowest_Ae = find(abs(b_Ae - a_Ace) == d_Ae(2));
                            
                            % Ae height interpolation
                            Ae = (Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Pounds_Force_60_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Pounds_Force_60_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Pounds_Force_60_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height))*(height_Ae - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height))/60;
                            
                            % Ve volume interpolation
                            Ve = (Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Volume_100_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Volume_100_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Volume_100_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height))*(height_Ae - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height));
                            
                        end
                        
                    case 4 % pressure level @ 80 PSIG
                        
                        while j <= length(Iso_eff_selection(i).ForceTable)% check if height_Ac = assembly height
                            
                            % Height and Volume
                            if height_Ac == Iso_eff_selection(i).ForceTable(j).Assembly_Height % Ac loop
                                
                                Ac = Iso_eff_selection(i).ForceTable(j).Pounds_Force_80_PSIG/80; % Effective Area at 1/2 inch below design height
                                trigger_Ac = 1;
                                
                                Vc = Iso_eff_selection(i).ForceTable(j).Volume_100_PSIG; % Vc at 1/2 inch below design height
                            end
                            
                            if height_Ae == Iso_eff_selection(i).ForceTable(j).Assembly_Height % Ae loop
                                
                                Ae = Iso_eff_selection(i).ForceTable(j).Pounds_Force_80_PSIG/80; % Effective Area at 1/2 inch above design height
                                trigger_Ae = 1;
                                
                                Ve = Iso_eff_selection(i).ForceTable(j).Volume_100_PSIG; % Vc at 1/2 inch above design height
                            end
                            
                            a_Ace(j) = Iso_eff_selection(i).ForceTable(j).Assembly_Height; % Assembly Height values
                            
                            j = j+1;
                        end
                        
                        if trigger_Ac == 0 % if Ac/Ae =~ assembly height --> interpolation
                            
                            % find limits for height and volume interpolation
                            
                            b_Ac = height_Ac;
                            d_Ac = sort(abs(b_Ac - a_Ace));
                            lowest_Ac = find(abs(b_Ac - a_Ace) == d_Ac(1));
                            sec_lowest_Ac = find(abs(b_Ac - a_Ace) == d_Ac(2));
                            
                            % Ac height interpolation
                            Ac = (Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Pounds_Force_80_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Pounds_Force_80_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Pounds_Force_80_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height))*(height_Ac - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height))/80;
                            
                            % Vc volume interpolation
                            Vc = (Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Volume_100_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Volume_100_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Volume_100_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height))*(height_Ac - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height));
                            
                        end
                        
                        if trigger_Ae == 0 % if Ac/Ae =~ assembly height --> interpolation
                            
                            % find limits for height interpolation
                            
                            b_Ae = height_Ae;
                            d_Ae = sort(abs(b_Ae - a_Ace));
                            lowest_Ae = find(abs(b_Ae - a_Ace) == d_Ae(1));
                            sec_lowest_Ae = find(abs(b_Ae - a_Ace) == d_Ae(2));
                            
                            % Ae height interpolation
                            Ae = (Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Pounds_Force_80_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Pounds_Force_80_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Pounds_Force_80_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height))*(height_Ae - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height))/80;
                            
                            % Ve volume interpolation
                            Ve = (Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Volume_100_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Volume_100_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Volume_100_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height))*(height_Ae - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height));
                            
                        end
                        
                    case 5 % pressure level @ 100 PSIG
                        
                        while j <= length(Iso_eff_selection(i).ForceTable)% check if height_Ac = assembly height
                            
                            % Height and Volume
                            if height_Ac == Iso_eff_selection(i).ForceTable(j).Assembly_Height % Ac loop
                                
                                Ac = Iso_eff_selection(i).ForceTable(j).Pounds_Force_100_PSIG/100; % Effective Area at 1/2 inch below design height
                                trigger_Ac = 1;
                                
                                Vc = Iso_eff_selection(i).ForceTable(j).Volume_100_PSIG; % Vc at 1/2 inch below design height
                            end
                            
                            if height_Ae == Iso_eff_selection(i).ForceTable(j).Assembly_Height % Ae loop
                                
                                Ae = Iso_eff_selection(i).ForceTable(j).Pounds_Force_100_PSIG/100; % Effective Area at 1/2 inch above design height
                                trigger_Ae = 1;
                                
                                Ve = Iso_eff_selection(i).ForceTable(j).Volume_100_PSIG; % Vc at 1/2 inch above design height
                            end
                            
                            a_Ace(j) = Iso_eff_selection(i).ForceTable(j).Assembly_Height; % Assembly Height values
                            
                            j = j+1;
                        end
                        
                        if trigger_Ac == 0 % if Ac/Ae =~ assembly height --> interpolation
                            
                            % find limits for height and volume interpolation
                            
                            b_Ac = height_Ac;
                            d_Ac = sort(abs(b_Ac - a_Ace));
                            lowest_Ac = find(abs(b_Ac - a_Ace) == d_Ac(1));
                            sec_lowest_Ac = find(abs(b_Ac - a_Ace) == d_Ac(2));
                            
                            % Ac height interpolation
                            Ac = (Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Pounds_Force_100_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Pounds_Force_100_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Pounds_Force_100_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height))*(height_Ac - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height))/100;
                            
                            % Vc volume interpolation
                            Vc = (Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Volume_100_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Volume_100_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Volume_100_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ac(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height))*(height_Ac - Iso_eff_selection(i).ForceTable(sec_lowest_Ac(1)).Assembly_Height));
                            
                        end
                        
                        if trigger_Ae == 0 % if Ac/Ae =~ assembly height --> interpolation
                            
                            % find limits for height interpolation
                            
                            b_Ae = height_Ae;
                            d_Ae = sort(abs(b_Ae - a_Ace));
                            lowest_Ae = find(abs(b_Ae - a_Ace) == d_Ae(1));
                            sec_lowest_Ae = find(abs(b_Ae - a_Ace) == d_Ae(2));
                            
                            % Ae height interpolation
                            Ae = (Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Pounds_Force_100_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Pounds_Force_100_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Pounds_Force_100_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height))*(height_Ae - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height))/100;
                            
                            % Ve volume interpolation
                            Ve = (Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Volume_100_PSIG + ((Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Volume_100_PSIG - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Volume_100_PSIG)/(Iso_eff_selection(i).ForceTable(lowest_Ae(1)).Assembly_Height - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height))*(height_Ae - Iso_eff_selection(i).ForceTable(sec_lowest_Ae(1)).Assembly_Height));
                            
                        end
                        
                end
                
                % Calculation of K: Vertical Spring rate in lbs/inch
                V1 = Iso_eff_selection(i).volume_100psig; % Internal Volume at design height in in^3
                
                
                K_selection(i).Vertical_Spring_Rate = ((K_selection(i).required_pressure + 14.7)*(Ac*((V1/Vc)^1.38) - Ae*((V1/Ve)^1.38)) - 14.7*(Ac - Ae)); % in lbs/inch
                
                K_selection(i).Natural_frequency = (188*sqrt(K_selection(i).Vertical_Spring_Rate/FSM(p).Sprung_Loads_Front))*0.01666666666; % in Hz
                K_selection(i).Diff_to_1Hz = abs(1 - K_selection(i).Natural_frequency);
                
                
                j = 1;
                i = i+1;
            end
            
            AirSuspension_selection = K_selection;
            
            
            % Vertical Spring Rate smallest value selection
            
            Afields = fieldnames(AirSuspension_selection);
            Acell = struct2cell(AirSuspension_selection);
            sz = size(Acell);
            
            
            % Convert to a matrix
            Acell = reshape(Acell, sz(1), []);      % Px(MxN)
            
            % Make each field a column
            Acell = Acell';                         % (MxN)xP
            
            % Sort by field 15 "Differnce to 1Hz"
            Acell = sortrows(Acell, 20);
            
            % Put back into original cell array format
            Acell = reshape(Acell', sz);
            
            % Convert to Struct
            Iso_eff_selection = cell2struct(Acell, Afields, 1);
            
            % Select closest frequency to 1 Hz
            AirSuspension_selection = [AirSuspension_selection(1)];
            %%
            %%
            end
            
            
        end
        function [mass,cost]=ladderframe(obj,length,width)
           
            thickness=obj.sectionthickness;
            intercrossmemberdistance=700; %mm  
            numcrossmembers=floor(length/intercrossmemberdistance);
            raillength=length/1000;
            crossraillength=width/1000; %vehicle width
            height=obj.sectionheight/1000; % in m
            width=obj.sectionwidth/1000; % in m
            numrails=2; % 2 longitudinal rails
            railvolume=((height*width) - ((height-thickness*2/1000)*...
                (width-thickness*2/1000)))*raillength;
            crossrailvolume=((height*width) - ((height-thickness*2/1000)*...
                (width-thickness*2/1000)))*crossraillength;
            railmass=railvolume*obj.steeldensity;
            crossrailmass=crossrailvolume*obj.steeldensity;
            mass=numrails*railmass + numcrossmembers*crossrailmass;
            cost= mass*(obj.steelrawprice/...
                obj.structure_materialutilisation)*1.53; % cost of frame
        end
        function plotchassis(obj,handle,wheelbase,vehiclewidth,groundclearance)
            width=obj.tyrewidth;%274.32; % tire width
            outerdiameter=obj.tyrediameter;%817.88;
            innerdiameter=obj.tyreinnerdiameter;%431.8;
            springpositionz=(outerdiameter-groundclearance+50)-obj.airspringheight/2;
            springpositiony1=vehiclewidth/2-(width+50+obj.airspringdiameter/2-5) ;
            springpositionyr=-vehiclewidth/2+(width+50+obj.airspringdiameter/2-5);
            airspringpositionyleft=vehiclewidth/2;
            airspringpositionyright=-vehiclewidth/2;
            airspringpositionyfront=wheelbase/2-1.5*outerdiameter/2;
            airspringpositionyback=-airspringpositionyfront;
            
            % plot tire and spring positions
            position=[wheelbase/2 (vehiclewidth-width)/2 (outerdiameter/2-groundclearance)];
            obj.tireplot(position,width,outerdiameter,innerdiameter,handle);
            springpos=[wheelbase/2  springpositiony1  springpositionz];
            airspringplot(obj,springpos,handle);
            airtankplot(obj,[airspringpositionyfront airspringpositionyleft 0],handle);
            position=[-wheelbase/2  (vehiclewidth-width)/2  (outerdiameter/2-groundclearance)];
            obj.tireplot(position,width,outerdiameter,innerdiameter,handle);
            springpos=[-wheelbase/2 springpositiony1 springpositionz];
            airspringplot(obj,springpos,handle);
            airtankplot(obj,[airspringpositionyback airspringpositionyleft 0],handle);
            position=[wheelbase/2 -(vehiclewidth-width)/2 (outerdiameter/2-groundclearance)];
            obj.tireplot(position,width,outerdiameter,innerdiameter,handle);
            springpos=[wheelbase/2 springpositionyr springpositionz];
            airspringplot(obj,springpos,handle);
           airtankplot(obj,[airspringpositionyfront airspringpositionyright 0],handle);
            position=[-wheelbase/2 -(vehiclewidth-width)/2 (outerdiameter/2-groundclearance)];
            obj.tireplot(position,width,outerdiameter,innerdiameter,handle);
            springpos=[-wheelbase/2 springpositionyr springpositionz];
            airspringplot(obj,springpos,handle);
            airtankplot(obj,[airspringpositionyback airspringpositionyright 0],handle);
        end
        function tireplot(obj,position,width,outerdiameter,innerdiameter,handle)
            
            colr = [0.1 0.1 0.1];
            alph = 1;
            %% Calculate vertices at origin
            
            vo1=[0.5*outerdiameter*cos(linspace(0,pi,40));-0.5*width+0*(linspace(0,pi,40));0.5*outerdiameter*sin(linspace(0,pi,40))];
            vo2=[0.5*outerdiameter*cos(linspace(0,pi,40));0.5*width+0*(linspace(0,pi,40));0.5*outerdiameter*sin(linspace(0,pi,40))];
            vo3=[0.5*outerdiameter*cos(linspace(pi,2*pi,40));-0.5*width+0*(linspace(0,pi,40));0.5*outerdiameter*sin(linspace(pi,2*pi,40))];
            vo4=[0.5*outerdiameter*cos(linspace(pi,2*pi,40));0.5*width+0*(linspace(0,pi,40));0.5*outerdiameter*sin(linspace(pi,2*pi,40))];
            vi1=[0.5*innerdiameter*cos(linspace(0,pi,40));-0.5*width+0*(linspace(0,pi,40));0.5*innerdiameter*sin(linspace(0,pi,40))];
            vi2=[0.5*innerdiameter*cos(linspace(0,pi,40));0.5*width+0*(linspace(0,pi,40));0.5*innerdiameter*sin(linspace(0,pi,40))];
            vi3=[0.5*innerdiameter*cos(linspace(pi,2*pi,40));-0.5*width+0*(linspace(0,pi,40));0.5*innerdiameter*sin(linspace(pi,2*pi,40))];
            vi4=[0.5*innerdiameter*cos(linspace(pi,2*pi,40));0.5*width+0*(linspace(0,pi,40));0.5*innerdiameter*sin(linspace(pi,2*pi,40))];
            V1=obj.translate([vo1,flip(vo2,2)],position);
            V2=obj.translate([vo3,flip(vo4,2)],position);
            V3=obj.translate([vi1,flip(vi2,2)],position);
            V4=obj.translate([vi3,flip(vi4,2)],position);
            V5=obj.translate([vo1,flip(vi1,2)],position);
            V6=obj.translate([vo2,flip(vi2,2)],position);
            V7=obj.translate([vo3,flip(vi3,2)],position);
            V8=obj.translate([vo4,flip(vi4,2)],position);
            
            patch(handle,'Faces',[1:80],'Vertices',V1','FaceColor', colr,'FaceAlpha',alph)
            patch(handle,'Faces',[1:80],'Vertices',V2','FaceColor', colr,'FaceAlpha',alph)
            patch(handle,'Faces',[1:80],'Vertices',V3','FaceColor', colr,'FaceAlpha',alph)
            patch(handle,'Faces',[1:80],'Vertices',V4','FaceColor', colr,'FaceAlpha',alph)
            patch(handle,'Faces',[1:80],'Vertices',V5','FaceColor', colr,'FaceAlpha',alph)
            patch(handle,'Faces',[1:80],'Vertices',V6','FaceColor', colr,'FaceAlpha',alph)
            patch(handle,'Faces',[1:80],'Vertices',V7','FaceColor', colr,'FaceAlpha',alph)
            patch(handle,'Faces',[1:80],'Vertices',V8','FaceColor', colr,'FaceAlpha',alph)
        %    view(handle,[140 30])
            
            
        end
        function airspringplot(obj,position,handle)
            diameter=obj.airspringdiameter;
            length=obj.airspringheight;
            orient=[0 0 pi/2];
           % position(2)=position(2) + ((position(2)>0)*-(50+diameter))+ ((position(2)<0)*(50+diameter));
            vo1=[0.5*diameter*cos(linspace(0,pi,40));-0.5*length+0*(linspace(0,pi,40));0.5*diameter*sin(linspace(0,pi,40))];
            vo2=[0.5*diameter*cos(linspace(0,pi,40));0.5*length+0*(linspace(0,pi,40));0.5*diameter*sin(linspace(0,pi,40))];
            vo3=[0.5*diameter*cos(linspace(pi,2*pi,40));-0.5*length+0*(linspace(0,pi,40));0.5*diameter*sin(linspace(pi,2*pi,40))];
            vo4=[0.5*diameter*cos(linspace(pi,2*pi,40));0.5*length+0*(linspace(0,pi,40));0.5*diameter*sin(linspace(pi,2*pi,40))];
            V1=obj.rotate([vo1,flip(vo2,2)]',orient)';
            V2=obj.rotate([vo3,flip(vo4,2)]',orient)';
            V3=obj.rotate([vo1,vo3]',orient)';
            V4=obj.rotate([vo2,vo4]',orient)';
            
            V1=obj.translate(V1,position);
            V2=obj.translate(V2,position);
            V3=obj.translate(V3,position);
            V4=obj.translate(V4,position);
            
            colr=[0 1 1];
            patch(handle,'Faces',[1:80],'Vertices',V1','FaceColor', colr)
            patch(handle,'Faces',[1:80],'Vertices',V2','FaceColor', colr)
            patch(handle,'Faces',[1:80],'Vertices',V3','FaceColor', colr)
            patch(handle,'Faces',[1:80],'Vertices',V4','FaceColor', colr)
            
        end
        function airtankplot(obj,position,handle)
             diameter=obj.airtankdiameter;
             length=obj.airtanklength;
             position(1)=(position(1)>0 )* (position(1)-length/2)+(position(1)<0) * (position(1)+length/2);
             position(2)=(position(2)>0 )*(position(2)-diameter/2-10)+(position(2)<0) * (position(2)+diameter/2+10);
             position(3)=diameter/2+5;
             orient=[pi/2 0 0];
        vo1=[0.5*diameter*cos(linspace(0,pi,40));-0.5*length+0*(linspace(0,pi,40));0.5*diameter*sin(linspace(0,pi,40))];
        vo2=[0.5*diameter*cos(linspace(0,pi,40));0.5*length+0*(linspace(0,pi,40));0.5*diameter*sin(linspace(0,pi,40))];
        vo3=[0.5*diameter*cos(linspace(pi,2*pi,40));-0.5*length+0*(linspace(0,pi,40));0.5*diameter*sin(linspace(pi,2*pi,40))];
        vo4=[0.5*diameter*cos(linspace(pi,2*pi,40));0.5*length+0*(linspace(0,pi,40));0.5*diameter*sin(linspace(pi,2*pi,40))];
         V1=obj.rotate([vo1,flip(vo2,2)]',orient)';
            V2=obj.rotate([vo3,flip(vo4,2)]',orient)';
            V3=obj.rotate([vo1,vo3]',orient)';
            V4=obj.rotate([vo2,vo4]',orient)';
            
            V1=obj.translate(V1,position);
            V2=obj.translate(V2,position);
            V3=obj.translate(V3,position);
            V4=obj.translate(V4,position);
       
        colr=[0 0 0.4];
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
        
    end
end


