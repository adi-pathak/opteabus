classdef depot
    properties
        lines
        passengers
        numberofterminals
        terminals
        depotcoordinates=[103.966986,1.36885]; %-loyang ;%woodlands?[103.786499996390,1.43762790090815];
        timetable
        blocks
        energy
        plot=0;
        vehicle
        gCO2_km
        gCO2_passengerkm
        TCO
        TCO_passengerkm
        fleetsize
        vehicledemand
        occupancy
        occupancy_passengers
        dailyvehiclekm
        dailypassengerkm
        dialogbar
        arcs
        property
        heuristics=1;
    end
    
    methods
        function obj=depot(services,depot_parameters,vehicle_parameters,driving_cycle,dialogbar)
            %% Create a new vehicle concept using the vehicle parameters
            vehicle=vehicleconcept(vehicle_parameters);
            %   if vehicle.Battery.capacity*1000/vehicle.Grossvehiclemass~=17.135325
            if vehicle.constraints~=-1
                if ~isempty(dialogbar)
                    dialogbar=waitbar(0.1,dialogbar,'Dervived Vehicle Concept');
                end
                %% Evaluate the timetable and passenger occupancy in the vehicles
                if ~isempty(dialogbar)
                    dialogbar=waitbar(0.25,dialogbar,'Starting Route Simulation');
                end
                for i=1:size(services,1)
                    service{i,1}=route(services{i,1},[],vehicle); %route(services{i,1},services{i,2},vehicle);
                end
                obj.lines=[service{:}];
                service=obj.lines;
                
                %% find ther terminal ids for each line
                terminals=[];
                for i=1:size(service,2)
                    terminals=[terminals;(service(i).demand.direction1([1 end],[2 14 16 17]))];
                end
                [a,b]=unique(terminals(:,1));
                terminals=terminals(sortrows(b),:);
                obj.terminals=terminals;
                obj.numberofterminals=size(terminals,1);
                obj.timetable=merge(obj,service,terminals); % combine timetable of each service
                
                % mean passenger occupancy of all routes
                obj.occupancy=mean([service.meanoccupancy])/vehicle.Passengercapacity*100;
                obj.occupancy_passengers=mean([service.meanoccupancy]);
                obj.dailypassengerkm=sum([service.passengerkm]);
                passengers={service.totalpassengers}';
                obj.passengers=[cell2mat(passengers)];
                
                % simulate the energyconsumption and longitudinal performance
                % after obtaining average passenger occupancy
                if ~isempty(dialogbar)
                    dialogbar=waitbar(0.45,dialogbar,'Completed Longitudinal Vehicle Simulation');
                end
                if isempty(driving_cycle)  % if no driving cycle specified use default
                    driving_cycle= load(strcat(pwd,filesep(),'Inputs',filesep(),...
                        'DrivingCycles',filesep(),'drivingcycle_brt.mat'));
                    driving_cycle= cell2mat(struct2cell(driving_cycle));
                end
                [vehicle.Energyconsumption,vehicle.Properties.range, ...
                    vehicle.Properties.vmax,vehicle.Properties.amax,...
                    vehicle.Properties.gradeability,vehicle.Energythroughput]=EnergyConsumption(vehicle,driving_cycle,obj.occupancy);
                if isnan(vehicle.Energyconsumption)
                    return
                end
                %% calculate fleet size and find optimal vehicle duties/blocks
                if ~isempty(dialogbar)
                    dialogbar=waitbar(0.5,dialogbar,'Finding Optimal Vehicle Assignment');
                end
                if obj.heuristics~=1
                    %     obj.blocks=SDVSP(obj,obj.timetable,depot_parameters.DH.times/60,depot_parameters.DH.distances/1000);
                    [obj.blocks,gap]=SDVSP_depotcharging(obj,vehicle,obj.timetable,...
                        depot_parameters.DH.times/60,depot_parameters.DH.distances/1000,...
                        vehicle.Energyconsumption,depot_parameters.chargingpower);
                    %      [obj.arcs,obj.blocks]=SDVSP_charge(obj,vehicle,obj.timetable,...
                    %         depot_parameters.DH.times/60,depot_parameters.DH.distances/1000,vehicle.Energyconsumption,depot_parameters.chargingpower);
                    if gap<=10
                        obj.dailyvehiclekm=sum([obj.blocks{4,:}]);
                        dailykm=([obj.blocks{4,:}]);
                        obj.fleetsize=size(obj.blocks,2)-1;
                        meandailykm=obj.dailyvehiclekm/obj.fleetsize;
                        obj.energy=obj.dailyvehiclekm*vehicle.Energyconsumption;
                        
                        [vehicle.Battery.replacements,vehicle.TCO]=LCC(vehicle,dailykm,obj.blocks,obj);
                        obj.TCO=vehicle.TCO;
                        [vehicle.LCAemissions,vehicle.EOLemissions,vehicle.WTTemissions,vehicle.Distributionemissions...
                            vehicle.Productionemissions,obj.gCO2_km,obj.gCO2_passengerkm]=LCemissions(vehicle,(([obj.blocks{4,:}])),obj.dailypassengerkm,obj);
                        if ~isempty(dialogbar)
                            dialogbar=waitbar(0.75,dialogbar,'Calculating Life Cycle Emissions');
                            
                        end
                        vehicle=PropertyEvaluation(vehicle,obj);
                        obj.vehicle=vehicle;
                        %% get vehicle demand
                        x={obj.blocks{2,1:end-1}};
                        fun=@(x) x(1);
                        x_1=cellfun(fun,x);
                        fun=@(x) x(end);
                        x_2=cellfun(fun,x);
                        T=[x_1' ones(size(x_1,2),1);x_2' -1*ones(size(x_2,2),1)];
                        T=sortrows(T);
                        [vehicledemand(:,1),vehicledemand(:,2)]=stairs([T(1,1) ;T(:,1)], [0;cumsum(T(:,2))]);
                        vehicledemand(vehicledemand(:,1)>24,1)=vehicledemand(vehicledemand(:,1)>24,1)-24;
                        vehicledemand=sortrows(vehicledemand,1);
                        obj.vehicledemand=vehicledemand;
                        obj.dialogbar=dialogbar;
                        %% Property evaluation
                        obj.property=propertyevaluation(obj);
                        %  end
                    end
                else
                   [arcs,var,nodes,timetable]=MDVSP_TSN(obj,vehicle,obj.timetable,depot_parameters{1, 1},depot_parameters{1, 2},vehicle.Energyconsumption,150,1)
                    obj.blocks=dispatchsimfifo(obj,arcs,var,nodes,timetable,vehicle.Energyconsumption);
                     obj.dailyvehiclekm=sum([obj.blocks{5,:}]);
                        dailykm=([obj.blocks{5,:}]);
                        obj.fleetsize=size(obj.blocks,2);
                        meandailykm=obj.dailyvehiclekm/obj.fleetsize;
                        obj.energy=obj.dailyvehiclekm*vehicle.Energyconsumption;
                        
                        [vehicle.Battery.replacements,vehicle.TCO]=LCC(vehicle,dailykm,obj.blocks,obj);
                        obj.TCO=vehicle.TCO;
                        [vehicle.LCAemissions,vehicle.EOLemissions,vehicle.WTTemissions,vehicle.Distributionemissions...
                            vehicle.Productionemissions,obj.gCO2_km,obj.gCO2_passengerkm]=LCemissions(vehicle,dailykm,obj.dailypassengerkm,obj);
                        if ~isempty(dialogbar)
                            dialogbar=waitbar(0.75,dialogbar,'Calculating Life Cycle Emissions');
                            
                        end
                        vehicle=PropertyEvaluation(vehicle,obj);
                        obj.vehicle=vehicle;
                        obj.property=propertyevaluation(obj);
                       
                    
                end
                
            end
        end
        
        function [timetable]=merge(obj,services,terminals)
            timetable=[];
            b_1=string(table2array(terminals(:,1)));
            x=[1:size(terminals,1)]'; % index
            for i=1:size(services,2)
                %
                a_1=string(table2array((services(i).demand.direction1([1 end],2))));
                n=x(ismember(b_1,a_1));
                a=services(i).timetable(services(i).timetable(:,3)==1,:); % trip direction 1
                b=services(i).timetable(services(i).timetable(:,3)==-1,:); % trip direction 2
                a=[a(:,1:2) ones(size(a,1),1)*services(i).routelength a(:,3)]; % trip distance
                b=[b(:,1:2) ones(size(b,1),1)*services(i).routelength b(:,3)];
                a(:,5)=cumsum(a(:,4)); % trip number
                b(:,5)=cumsum(-1*b(:,4));
                %
                if ~isempty(a) & ~isempty(b)
                    timetable=[timetable;n(1)*ones(size(a,1),1) (n(2))*ones(size(a,1),1)   a i*ones(size(a,1),1)];
                    timetable=[timetable;(n(2))*ones(size(b,1),1) (n(1))*ones(size(b,1),1)   b i*ones(size(b,1),1)];
                elseif    ~isempty(a) & isempty(b)
                    timetable=[timetable;n(1)*ones(size(a,1),1) (n(2))*ones(size(a,1),1)   a i*ones(size(a,1),1)];
                    
                else
                    error
                end
                timetable=sortrows(timetable,3); % sort according to depart time
                
            end
            
        end
        function [blocks]=SDVSP(obj,timetable,DH,DH_dist)
            Econs=1; %kwh/km
            pulloutcost=150000;
            
            departl=[timetable(:,3) timetable(:,1)];
            arrivall=[timetable(:,4) timetable(:,2)];
            ui=timetable(:,end)*Econs; % Energy consumption of each trip
            
            
            N=1:size(timetable,1); % Set of nodes/trips sorted by travel/depart time
            % E=nchoosek(N,2); % Possible set of deadhead arcs between trips
            tic
            E=combinator(size(timetable,1),2,'c'); % faster than nchoosek
            toc
            tic
            for i=1:size(E,1)
                E(i,3)=departl(E(i,2))-(arrivall(E(i,1))+(((DH(arrivall(E(i,1),2),...
                    departl(E(i,2),2))/60)))); % check if trips in E can be joined- starting time -(arrival time +deadheading time)
                E(i,4)=DH_dist(timetable(E(i,1),2),timetable(E(i,2),1)); % deadheading distances between trips
            end
            toc
            DH_E=E(E(:,3)>=0  ,4);
            E=E(E(:,3)>=0 ,:);
            %             DH_E=DH_E(E(:,3)<2,:); % remove trip arcs if waiting time is over 6 hours
            %             E=E(E(:,3)<2,:); % remove trip arcs if waiting time is over 6 hours
            
            
            E=E(:,1:2);
            E1=E;
            DH1=DH_E;
            vij=Econs.*DH_E; %Energy consumed during deadheading
            E1cost=ui(E1(:,1))+vij;
            %%
            V=[N N(end)+1 N(end)+2]; % end depot start and end nodes
            E3=[N' V(end)*ones(size(N,2),1)]; % add arcs from trips to depot
            DH_depot_1=DH_dist(end,1:end-1); % from depot
            DH_depot_2=[DH_dist(1:end-1,end)]';
            rhoij=DH_depot_1(1,timetable(E3(:,1),1))*Econs';
            
            E2=[V(end-1)*ones(size(N,2),1) [N']]; % pull out arcs from depot
            DH2=DH_depot_2(1,timetable(E2(:,2),1));
            piij=DH2*Econs';
            
            E=[E; N' V(end)*ones(size(N,2),1)]; % add arcs to depot / pull in arcs
            E=[E; V(end-1)*ones(size(N,2),1) [N'] ]; % add arcs from depot to trips
            
            C=[vij;rhoij';piij']; % energy usage of each arc - does not include trip energy use
            % Ci=[vij;rhoij';batterycapacity-piij'];
            
            E_1=E1;
            
            CM=1:size(timetable,1);
            fixedcost=E1cost;
            
            
            variablecost=fixedcost;%variablecost(i,:);
            CN1=variablecost; % energy consumption for each arc in E2
            pull_outcosts=zeros(size(E2,1),1);
            DH2=DH_depot_1(1,timetable(E2(:,2)))';
            
            CN2=zeros(size(E2,1),1);
            CN2(1:end,:)=DH2*Econs; % deadheading energy consumption from depot to trips
            pull_outcosts(1:end,:)=pulloutcost+CN2; % cost of starting new vehicle + deadhead
            
            
            DH3=DH_depot_2(1,timetable(E3(:,1),2))';
            CN3=zeros(size(E3,1),1); % energy consumption when returning to depot
            CN3(1:end,:)=DH3.*Econs;  % energy consumption when returning to depot
            pull_incosts=zeros(size(E3,1),1);
            pull_incosts(1:end,:)=pull_incosts(1:end)+CN3; % cost of deadhead to depot
            
            CN=[CN1;CN2;CN3]; %energy consumed in each arc and trip
            DH_E=[DH1;DH2;DH3]; %deadheading
            
            E0=[E1 variablecost;E2 pull_outcosts;E3 pull_incosts]; %trip arcs with costs
            E0(:,4)=CN;
            E0(:,5)=DH_E;
            E0=sortrows(E0);
            Deadheading=E0(:,5);
            
            %% define optimization constraints
            
            %
            %             Aeqi=sparse(size(N,2),size(E,1));
            %             Aeqj=sparse(size(N,2),size(E,1));
            %             tic
            %             for i=1:size(N,2) % for each timetabled trip
            %                 whichIdxsi = (E0(:,1) == i); % find all dead head trips from timetable trip i
            %                 whichIdxsj = (E0(:,2) == i); % find all dead head trips to timetable trip i
            %                 Aeqi(i,whichIdxsi') = 1; % include in the constraint matrix
            %                 Aeqj(i,whichIdxsj') = 1; % include in the constraint matrix
            %             end
            %            toc
            %             Aeq=[Aeqi;Aeqj];
            Ce1=cell(size(N,2),2);
            Ce2=cell(size(N,2),2);
            tic
            for i=1:size(N,2) % for each timetabled trip
                colsi=find(E0(:,1) == i);
                colsj=find(E0(:,2) == i);
                Ce1{i,1}=colsi;
                Ce1{i,2}=ones(size(colsi,1),1)*i;
                Ce2{i,1}=colsj;
                Ce2{i,2}=ones(size(colsj,1),1)*i;
            end
            Ce1= cell2mat(Ce1);
            Ce1(:,3)=ones(size(Ce1,1),1);
            Ce2= cell2mat(Ce2);
            Ce2(:,3)=ones(size(Ce2,1),1);
            de1=find((E0(:,2) == size(V,2)));
            de1(:,2)=ones(size(de1,1),1)*((2*size(N,2))+1);
            de1(:,3)=ones(size(de1,1),1);
            de2=find((E0(:,1) == size(V,2)-1));
            de2(:,2)=1*ones(size(de2,1),1)*((2*size(N,2))+1);
            de2(:,3)=-1*ones(size(de2,1),1);
            de=[de1;de2];
            
            Aeq=sparse( [Ce1(:,2);size(N,2)+Ce2(:,2);de(:,2)], [Ce1(:,1);Ce2(:,1);de(:,1)], [Ce1(:,3);Ce2(:,3);de(:,3)],(2*size(N,2))+1,size(E,1));
            toc
            
            beq= ones( size(Aeq,1),1);
            beq(end)=0;
            
            %Aeq=[Aeq;[(E(:,2) == size(V,2))-(E(:,1) == size(V,2)-1)]']; % depot constraint
            %beq=[beq;0]; % depot constraint: no of vehicle pull out=vehicles pulling in
            intcon=(1:length(E0))'; % optimization variables are possible arcs
            lb=zeros(size(intcon,1),1);
            ub=ones(size(intcon,1),1);
            
            
            % define cost matrix
            cost=E0(:,3);
            
            %% optimize
            tic
            [x,fval]=intlinprog(cost,intcon,[],[],Aeq,beq,lb,ub);
            toc
            
            %% form chains/blocks
            z=E0(x.*E0(:,1)>0,:); % filter used arcs
            z_dh=Deadheading(x.*E0(:,1)>0,:); % filter used dead heading arcs
            blocks=cell(8,1);
            blocks{1}=0;
            n=1;
            j=1;
            while ~isempty(j)
                j=j(1);
                while j<=V(end)
                    if blocks{1,n}==0
                        blocks{1,n}= z(z(:,2)==j,1:2);
                        blocks{2,n}=[timetable(z(z(:,2)==j,2),3)-DH(end,timetable(z(z(:,2)==j,2),1))/60 timetable(z(z(:,2)==j,2),3:4)];
                        blocks{3,n}=z_dh(z(:,2)==j,:);
                        blocks{5,n}=[size(DH,1) timetable(z(z(:,2)==j,2),1:2)];
                        blocks{6,n}=[timetable(z(z(:,2)==j,2),3)-DH(end,timetable(z(z(:,2)==j,2),1))/60 timetable(z(z(:,2)==j,2),3)]; % deadheadtrips
                        blocks{7,n}=[timetable(z(z(:,2)==j,2),3:4)]; % timetabledtrips
                        z_dh(z(:,2)==j,:)=[];
                        z(z(:,2)==j,:)=[];
                        
                        j=z(z(:,1)==j,2);
                    elseif j==V(end)
                        blocks{3,n}=[blocks{3,n} z_dh(z(:,1)==blocks{1,n}(end),:)];
                        blocks{5,n}=[blocks{5,n} timetable(blocks{1,n}(end),2) size(DH,1) ];
                        blocks{6,n}=[blocks{6,n} blocks{2,n}(end) blocks{2,n}(end)+DH(timetable(blocks{1,n}(end),2),end)/60]; % deadheadtrips
                        z_dh(z(:,1)==blocks{1,n}(end),:)=[];
                        blocks{2,n}=[blocks{2,n}  timetable(blocks{1,n}(end),4)+DH(timetable(blocks{1,n}(end),2),end)/60];
                        z(z(:,1)==blocks{1,n}(end),:)=[];
                        blocks{1,n}=[blocks{1,n} V(end)];
                        j=j+1;
                    else
                        blocks{1,n}=[blocks{1,n} z(z(:,2)==j,2)];
                        if timetable(j,1)==blocks{5,n}(end) % if next trip from same stop
                            % add waiting time
                            blocks{8,n}=[blocks{8,n} blocks{2,n}(end) timetable(z(z(:,2)==j,2),3)]; % waitingtrips
                        else
                            %deadhead
                            blocks{6,n}=[blocks{6,n} blocks{2,n}(end) blocks{2,n}(end)+DH(timetable(blocks{1,n}(end-1),2),timetable(blocks{1,n}(end),1))/60];
                            % add waiting time
                            blocks{8,n}=[blocks{8,n} blocks{6,n}(end) timetable(z(z(:,2)==j,2),3)]; % waitingtrips
                            
                        end
                        blocks{2,n}=[blocks{2,n} timetable(z(z(:,2)==j,2),3:4)];
                        blocks{3,n}=[blocks{3,n} timetable(z(z(:,2)==j,2),5)+z_dh(z(:,2)==j,:)];
                        blocks{5,n}=[blocks{5,n} timetable(z(z(:,2)==j,2),1:2)];
                        blocks{7,n}=[blocks{7,n} timetable(z(z(:,2)==j,2),3:4)]; % timetabledtrips
                        z_dh(z(:,2)==j,:)=[];
                        z(z(:,2)==j,:)=[];
                        j=z(z(:,1)==j,2);
                    end
                end
                blocks{4,n}= sum(blocks{3,n});
                n=n+1;
                blocks{1,n}=0;
                j=sortrows(z(:,1));
                %  reshape(tb{7,1},2,6)'
            end
            
        end
        function [arcs,blocks]=SDVSP_charge(obj,vehicle,timetable,DH,DH_dist,Econs,chargingpower)
            tic
            uncapacitated=0;
            gurobi_solver=1;
            % Econs=0.59; %0.59; %kwh/km
            pulloutcost=150000;
            chargerpower=chargingpower; % 300 kW
            batterycapacity=vehicle.Battery.capacity  ; %120 kWh
            departl=[timetable(:,3) timetable(:,1)];
            arrivall=[timetable(:,4) timetable(:,2)];
            ui=timetable(:,5)*Econs; % Energy consumption of each trip
            Q=batterycapacity;
            
            N=1:size(timetable,1); % Set of nodes/trips sorted by travel/depart time
            % E=nchoosek(N,2); % Possible set of deadhead arcs between trips
            tic
            E=combinator(size(timetable,1),2,'c'); % faster than nchoosek
            toc
            tic
            for i=1:size(E,1)
                E(i,3)=departl(E(i,2))-(arrivall(E(i,1))+(((DH(arrivall(E(i,1),2),...
                    departl(E(i,2),2))/60)))); % check if trips in E can be joined- starting time -(arrival time +deadheading time)
                E(i,4)=DH_dist(timetable(E(i,1),2),timetable(E(i,2),1)); % deadheading distances between trips
                
            end
            toc
            DH_E=E(E(:,3)>=0  ,4);
            E=E(E(:,3)>=0 ,:);
            DH_E=DH_E(E(:,3)<2,:); % remove trip arcs if waiting time is over 6 hours
            E=E(E(:,3)<2,:); % remove trip arcs if waiting time is over 6 hours
            
            chargetime=E(:,3); %available time until next trip in hours
            E=E(:,1:2);
            E1=E;
            DH1=DH_E;
            vij=Econs.*DH_E; %Energy consumed during deadheading
            wij=chargerpower*chargetime;
            wij(wij>batterycapacity)=batterycapacity;
            E1cost=ui(E1(:,1))+vij;
            charge=ui(E1(:,1))+vij-wij;
            %%
            V=[N N(end)+1 N(end)+2]; % end depot start and end nodes
            E3=[N' V(end)*ones(size(N,2),1)]; % add arcs from trips to depot
            DH_depot_1=DH_dist(end,1:end-1); % from depot
            DH_depot_2=[DH_dist(1:end-1,end)]';
            rhoij=DH_depot_1(1,timetable(E3(:,1),1))*Econs';
            
            E2=[V(end-1)*ones(size(N,2),1) [N']]; % pull out arcs from depot
            DH2=DH_depot_2(1,timetable(E2(:,2),1));
            piij=DH2*Econs';
            
            E=[E; N' V(end)*ones(size(N,2),1)]; % add arcs to depot / pull in arcs
            E=[E; V(end-1)*ones(size(N,2),1) [N'] ]; % add arcs from depot to trips
            
            C=[vij;rhoij';piij']; % energy usage of each arc - does not include trip energy use
            % Ci=[vij;rhoij';batterycapacity-piij'];
            
            E_1=E1;
            
            CM=1:size(timetable,1);
            
            
            
            
            CN1=E1cost; % energy consumption for each arc in E2
            pull_outcosts=zeros(size(E2,1),1);
            DH2=DH_depot_1(1,timetable(E2(:,2)))';
            
            CN2=zeros(size(E2,1),1);
            CN2(1:end,:)=DH2*Econs; % deadheading energy consumption from depot to trips
            pull_outcosts(1:end,:)=pulloutcost+CN2; % cost of starting new vehicle + deadhead
            
            
            DH3=DH_depot_2(1,timetable(E3(:,1),2))';
            CN3=zeros(size(E3,1),1); % energy consumption when returning to depot
            CN3(1:end,:)=(DH3.*Econs)+ui;  % energy consumption when returning to depot
            pull_incosts=zeros(size(E3,1),1);
            pull_incosts(1:end,:)=pull_incosts(1:end)+CN3; % cost of deadhead to depot
            % CN3(1:end,:)=(DH3.*Econs)+ui-Q;  % energy consumption when returning to depot
            
            
            CN=[CN1;CN2;CN3]; %energy consumed in each arc and trip
            DH_E=[DH1;DH2;DH3]; %deadheading
            
            E0=[E1 E1cost;E2 pull_outcosts;E3 pull_incosts]; %trip arcs with costs
            E0(:,4)=[charge;CN2;CN3]; %CN;
            E0(:,5)=DH_E;
            E0(:,6)=CN; %energy discharged
            E0(1:size(wij,1),7)=wij;
            E0=sortrows(E0);
            
            Deadheading=E0(:,5);
            
            
            %%
            
            Ce1=cell(size(N,2),2);
            Ce2=cell(size(N,2),2);
            tic
            for i=1:size(N,2) % for each timetabled trip
                colsi=find(E0(:,1) == i);
                colsj=find(E0(:,2) == i);
                Ce1{i,1}=colsi;
                Ce1{i,2}=ones(size(colsi,1),1)*i;
                Ce2{i,1}=colsj;
                Ce2{i,2}=ones(size(colsj,1),1)*i;
            end
            Ce1= cell2mat(Ce1);
            Ce1(:,3)=ones(size(Ce1,1),1);
            Ce2= cell2mat(Ce2);
            Ce2(:,3)=ones(size(Ce2,1),1);
            de1=find((E0(:,2) == size(V,2)));
            de1(:,2)=ones(size(de1,1),1)*((2*size(N,2))+1);
            de1(:,3)=ones(size(de1,1),1);
            de2=find((E0(:,1) == size(V,2)-1));
            de2(:,2)=1*ones(size(de2,1),1)*((2*size(N,2))+1);
            de2(:,3)=-1*ones(size(de2,1),1);
            de=[de1;de2];
            
            Aeq=sparse([Ce1(:,2);size(N,2)+Ce2(:,2);de(:,2)], [Ce1(:,1);Ce2(:,1);de(:,1)], [Ce1(:,3);Ce2(:,3);de(:,3)],(2*size(N,2))+1,size(E,1));
            toc
            
            
            
            
            %%
            tic
            if uncapacitated~=1
                Ce0=cell(size(N,2),2);
                Ce01=cell(size(N,2),2);
                Ce1=cell(size(N,2),2);
                Ce2=cell(size(N,2),2);
                G=digraph(E0(:,1),E0(:,2));
                Q=batterycapacity;
                n=0;
                for i=1:size(N,2)
                    [idx,idy]=inedges(G,i);
                    [idx1,idy1]=outedges(G,i);
                    
                    
                    %               Ce1{i,1}=[idx; size(E,1)+i];
                    %               Ce1{i,2}=[ones(size(idx,1),1)*i;Ce1{i,2};i];
                    %               Ce1{i,3}=[E0(idx,4);-1];
                    %
                    Ce1{i,1}=[idx(end); size(E,1)+i];
                    Ce1{i,2}=[ones(size(idx(end),1),1)*i;i];
                    Ce1{i,3}=[E0(idx(end),4);1];
                    
                    
                    
                    
                    Ce2{i,1}=[idx1(end); size(E,1)+i];
                    Ce2{i,2}=[i;i];
                    Ce2{i,3}=[E0(idx1(end),4);-1];
                    
                    
                    for ii=1:size(idx1,1)
                        
                        
                        if idy1(ii)>N(end)
                            
                        else
                            n=n+1;
                            Ce0{n,1}=[idx1(ii);size(E,1)+idy1(ii);size(E,1)+i];
                            Ce0{n,2}=ones(size(Ce0{n,1},1),1)*n;
                            Ce0{n,3}=[E0(idx1(ii),4)+Q;1;-1];
                            
                            
                            Ce01{n,1}=[idx1(ii);size(E,1)+i];
                            Ce01{n,2}=ones(size(Ce01{n,1},1),1)*n;
                            Ce01{n,3}=[E0(idx1(ii),6);-1];
                        end
                        
                    end
                    
                end
                toc
                Ce1= cell2mat(Ce1);
                Ce2= cell2mat(Ce2);
                Ce0= cell2mat(Ce0);
                Ce01= cell2mat(Ce01);
                %  A=sparse([Ce0(:,2)],[Ce0(:,1)],[Ce0(:,3)]);
                %A=sparse([Ce0(:,2);n+Ce1(:,2)],[Ce0(:,1);Ce1(:,1)],[Ce0(:,3);Ce1(:,3)]);
                %  A=sparse([Ce0(:,2);n+Ce2(:,2)],[Ce0(:,1);Ce2(:,1)],[Ce0(:,3);Ce2(:,3)]);
                % A=sparse([Ce0(:,2);n+Ce2(:,2);n+size(N,2)+Ce1(:,2)],[Ce0(:,1);Ce2(:,1);Ce1(:,1)],[Ce0(:,3);Ce2(:,3);Ce1(:,3)]);
                A=sparse([Ce0(:,2);n+Ce2(:,2);n+size(N,2)+Ce1(:,2);n+size(N,2)+size(N,2)+Ce01(:,2)],[Ce0(:,1);Ce2(:,1);Ce1(:,1);Ce01(:,1)],[Ce0(:,3);Ce2(:,3);Ce1(:,3);Ce01(:,3)]);
                b=ones(n,1)*Q;
                b_1=zeros(size(N,2),1);
                % b=[b;b_1];
                b_2=ones(size(N,2),1)*Q;
                b_3=zeros(n,1);
                b=[b;b_1;b_2;b_3];
                Aeq(:,size(A,2))=0;
                intcon=(1:(length(E))'); % optimization variables are possible arcs
                
                intN=(1:length(N))';
                cost=E0(:,3);
                cost(size(A,2),:)=0;
                lb=zeros(size(intcon,2),1);
                ub=ones(size(intcon,2),1);
                lbu=zeros(size(intN,1),1);
                ubu=Q*ones(size(intN,1),1);
                % define cost matrix
                intcon=[intcon';intN];
                model.vtype=[repmat('B', size(E0,1), 1); repmat('C', size(N',1), 1)];
                lb=[lb;lbu];
                ub=[ub;ubu];
            else
                A=[];
                b=[];
                cost=E0(:,3);
                intcon=(1:(length(E))'); % optimization variables are possible arcs
                
                intN=(1:length(N))';
                lb=zeros(size(intcon,2),1);
                ub=ones(size(intcon,2),1);
                model.vtype=[repmat('B', size(E0,1), 1)];
            end
            beq= ones( size(Aeq,1)-1,1);
            beq=[beq;0];
            
            
            
            
            if gurobi_solver~=1
                
                tic
                [x,fval]=intlinprog(cost,intcon,A,b,Aeq,beq,lb,ub);
                toc
            else
                tic
                addpath('C:\gurobi902\win64\matlab')
                model.obj=cost;
                model.A=[A;Aeq]; % - not efficient/slow speed for large matrices
                model.sense=[repmat('<', size(A,1), 1);repmat('=', size(Aeq,1), 1)];
                model.rhs=[b;beq];
                model.modelsense='min';
                params.outputflag=1;
                params.TimeLimit = 350; % 300 seconds to find optimal
                model.lb    = lb;
                model.ub    = ub;
                
                r=gurobi(model,params);
                toc
                x=r.x;
            end
            
            arcs=E0;
            if uncapacitated~=1
                x1=x;
                x=x(1:end-(size(timetable,1)),:);
                c=x1((end-(size(timetable,1)-1):end));
            end
            %% form chains/blocks
            x=round(x);
            if isempty(x)
                disp('No solution found')
            else
                z=E0(x.*E0(:,1)>0,:); % filter used arcs
                %z_dh=Deadheading(x.*E(:,1)>0,:); % filter used dead heading arcs
                blocks=cell(8,1);
                blocks{1}=0;
                n=1;
                j=1;
                while ~isempty(j)
                    j=j(1);
                    while j<=V(end)
                        if blocks{1,n}==0
                            blocks{1,n}= z(z(:,2)==j,1:2); % trip joinings
                            blocks{2,n}=[timetable(z(z(:,2)==j,2),3)-DH(end,timetable(z(z(:,2)==j,2),1))/60 timetable(z(z(:,2)==j,2),3:4)]; % dispatch, departure and arrival times
                            blocks{3,n}= z(z(:,2)==j,5); % total distances (trip+deadhead)
                            blocks{5,n}=[size(DH,1) timetable(z(z(:,2)==j,2),1:2)]; %terminals
                            blocks{6,n}=[timetable(z(z(:,2)==j,2),3)-DH(end,timetable(z(z(:,2)==j,2),1))/60 timetable(z(z(:,2)==j,2),3)]; % deadheadtrips
                            blocks{7,n}=[timetable(z(z(:,2)==j,2),3:4)]; % timetabledtrips
                            blocks{10,n}=[ z(z(:,2)==j,6)]; % energy discharged
                            blocks{11,n}=[ z(z(:,2)==j,7)]; % energy charged
                            blocks{12,n}=[batterycapacity batterycapacity-blocks{10,n}(end)];
                            blocks{13,n}=[blocks{2,n}(1:end-1) ];
                            z(z(:,2)==j,:)=[];
                            j=z(z(:,1)==j,2);
                        elseif j==V(end)
                            blocks{3,n}=[blocks{3,n} z((z(:,1)==blocks{1,n}(end) & z(:,2)==j),5)];  % total distances (trip+deadhead)
                            blocks{10,n}=[blocks{10,n} E0(E0(:,1)==blocks{1,n}(end) & E0(:,2)==V(end),4)]; % energy discharged
                            z_end=E0(E0(:,1)==blocks{1,n}(end) & E0(:,2)==V(end),:);
                            blocks{12,n}=[blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)-z_end(:,6)];
                            blocks{13,n}=[blocks{13,n} timetable(z_end(:,1),3) timetable(blocks{1,n}(end),4)+DH(timetable(blocks{1,n}(end),2),end)/60];
                            % blocks{12,n}=[blocks{12,n} blocks{12,n}(end)-E0(E0(:,1)==blocks{1,n}(end) & E0(:,2)==V(end),4)];
                            blocks{5,n}=[blocks{5,n} timetable(blocks{1,n}(end),2) size(DH,1) ];%terminals
                            blocks{6,n}=[blocks{6,n} blocks{2,n}(end) blocks{2,n}(end)+DH(timetable(blocks{1,n}(end),2),end)/60]; % deadheadtrips
                            %  z_dh(z(:,1)==blocks{1,n}(end),:)=[];
                            blocks{2,n}=[blocks{2,n}  timetable(blocks{1,n}(end),4)+DH(timetable(blocks{1,n}(end),2),end)/60];  % dispatch, departure and arrival times
                            z(z(:,1)==blocks{1,n}(end),:)=[];
                            blocks{1,n}=[blocks{1,n} V(end)];
                            j=j+1;
                        else
                            blocks{1,n}=[blocks{1,n} z(z(:,2)==j,2)]; % trip joinings
                            blocks{10,n}=[blocks{10,n} z(z(:,2)==j,6)]; % energy discharged
                            blocks{3,n}=[blocks{3,n} timetable(z(z(:,2)==j,2),5)+z(z(:,2)==j,5)];  % total distances (trip+deadhead)
                            
                            
                            
                            if timetable(j,1)==blocks{5,n}(end) % if next trip from same stop
                                % add waiting time
                                blocks{8,n}=[blocks{8,n} blocks{2,n}(end) timetable(z(z(:,2)==j,2),3)]; % waitingtrips
                                
                                charge_block=blocks{12,n}(end)-z(z(:,2)==j,6)+z(z(:,2)==j,7);
                                if charge_block>0
                                    
                                    if charge_block>batterycapacity
                                        blocks{12,n}= [blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)-z(z(:,2)==j,6)]; % amount discharged
                                        e_charged=batterycapacity-(blocks{12,n}(end));
                                        blocks{12,n}= [blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)+e_charged];
                                        blocks{11,n}=[blocks{11,n} e_charged];
                                        charge_endtime=e_charged/chargerpower;
                                        if isnan(charge_endtime)
                                            charge_endtime=0;
                                        end
                                        blocks{13,n}= [blocks{13,n} timetable(z(z(:,2)==j,1),3:4)];
                                        blocks{13,n}= [blocks{13,n} blocks{13,n}(end) timetable(z(z(:,2)==j,1),4)+charge_endtime];
                                        
                                        
                                    else
                                        e_charged=z(z(:,2)==j,7);
                                        blocks{11,n}=[blocks{11,n} e_charged]; % energy charged
                                        %  blocks{12,n}=[blocks{12,n} blocks{12,n}(end)-z(z(:,2)==j,6) blocks{12,n}(end)-z(z(:,2)==j,6)+z(z(:,2)==j,7)];
                                        blocks{12,n}= [blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)-z(z(:,2)==j,6)]; % update discharge level
                                        blocks{12,n}= [blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)+e_charged]; %update charge level after trip i
                                        % calculate charge duration
                                        charge_endtime=e_charged/chargerpower;
                                        if isnan(charge_endtime)
                                            charge_endtime=0;
                                        end
                                        blocks{13,n}= [blocks{13,n} timetable(z(z(:,2)==j,1),3:4) timetable(z(z(:,2)==j,1),4) timetable(z(z(:,2)==j,1),4)+charge_endtime];
                                        % no waiting between charging and next trip blocks{13,n}= [blocks{13,n} blocks{13,n}(end) timetable(z(z(:,2)==j,2),3)];
                                    end
                                end
                                
                            else
                                %deadhead
                                blocks{6,n}=[blocks{6,n} blocks{2,n}(end) blocks{2,n}(end)+DH(timetable(blocks{1,n}(end-1),2),timetable(blocks{1,n}(end),1))/60];% deadheadtrips
                                % add waiting time
                                blocks{8,n}=[blocks{8,n} blocks{6,n}(end) timetable(z(z(:,2)==j,2),3)]; % waitingtrips
                                
                                charge_block=blocks{12,n}(end)-z(z(:,2)==j,6)+z(z(:,2)==j,7);
                                if charge_block>0
                                    
                                    if charge_block>batterycapacity
                                        blocks{12,n}= [blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)-z(z(:,2)==j,6)]; % amount discharged
                                        e_charged=batterycapacity-(blocks{12,n}(end));
                                        blocks{12,n}= [blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)+e_charged];
                                        blocks{11,n}=[blocks{11,n} e_charged];
                                        charge_endtime=e_charged/chargerpower;
                                        if isnan(charge_endtime)
                                            charge_endtime=0;
                                        end
                                        blocks{13,n}= [blocks{13,n} timetable(z(z(:,2)==j,1),3) blocks{6,n}(end)];
                                        blocks{13,n}= [blocks{13,n} blocks{13,n}(end) blocks{6,n}(end)+charge_endtime];
                                        
                                        
                                    else
                                        e_charged=z(z(:,2)==j,7);
                                        blocks{11,n}=[blocks{11,n} e_charged]; % energy charged
                                        %  blocks{12,n}=[blocks{12,n} blocks{12,n}(end)-z(z(:,2)==j,6) blocks{12,n}(end)-z(z(:,2)==j,6)+z(z(:,2)==j,7)];
                                        blocks{12,n}= [blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)-z(z(:,2)==j,6)]; % update discharge level
                                        blocks{12,n}= [blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)+e_charged]; %update charge level after trip i
                                        % calculate charge duration
                                        charge_endtime=e_charged/chargerpower;
                                        if isnan(charge_endtime)
                                            charge_endtime=0;
                                        end
                                        blocks{13,n}= [blocks{13,n} timetable(z(z(:,2)==j,1),3) blocks{6,n}(end) blocks{6,n}(end) blocks{6,n}(end)+charge_endtime];
                                        % no waiting between charging and next trip blocks{13,n}= [blocks{13,n} blocks{13,n}(end) timetable(z(z(:,2)==j,2),3)];
                                    end
                                end
                                
                                
                            end
                            blocks{2,n}=[blocks{2,n} timetable(z(z(:,2)==j,2),3:4)];
                            
                            blocks{5,n}=[blocks{5,n} timetable(z(z(:,2)==j,2),1:2)]; % deadheadtrips
                            blocks{7,n}=[blocks{7,n} timetable(z(z(:,2)==j,2),3:4)]; % timetabledtrips
                            % z_dh(z(:,2)==j,:)=[];
                            z(z(:,2)==j,:)=[];
                            j=z(z(:,1)==j,2);
                        end
                    end
                    blocks{4,n}= sum(blocks{3,n}); % total distance driven of vehicle
                    blocks{9,n}= (batterycapacity-cumsum(blocks{3,n}))/batterycapacity*100; %wrong add energy consumption
                    blocks{14,n}= (blocks{12,n})/batterycapacity; % total distance driven of vehicle
                    %    blocks{12,n}= batterycapacity-cumsum(blocks{10,n});%;
                    n=n+1;
                    blocks{1,n}=0;
                    j=sortrows(z(:,1));
                    %  reshape(tb{7,1},2,6)'
                end
            end
            toc
        end
        function [blocks,gap]=SDVSP_depotcharging(obj,vehicle,timetable,DH,DH_dist,Econs,chargingpower)
            gurobi_solver=1;
            pulloutcost=150000;
            chargerpower=chargingpower; % 300 kW
            chargercost=5000;
            batterycapacity=vehicle.Battery.capacity;
            departl=[timetable(:,3) timetable(:,1)];
            arrivall=[timetable(:,4) timetable(:,2)];
            ui=timetable(:,5)*Econs; % Energy consumption of each trip
            minimumchargetime=0.5;%.3*batterycapacity/chargerpower;%10/60; % minimum available time to charge 5 mins
            chargeraccesstime=5/60;
            chargerconnectiontime=3/60; % minimum time between charging 2 vehicles
            minimumlayovertime=5/60;
            N=1:size(timetable,1); % Set of nodes/trips sorted by travel/depart time
            %E=nchoosek(N,2); % Possible set of deadhead arcs between trips
            tic
            maxchargetime=batterycapacity/chargerpower;
            
            E=combinator(size(timetable,1),2,'c'); % faster than nchoosek
            toc
            tic
            for i=1:size(E,1)
                E(i,3)=departl(E(i,2))-(arrivall(E(i,1))+(((DH(arrivall(E(i,1),2),...
                    departl(E(i,2),2))/60)))); % check if trips in E can be joined- starting time -(arrival time +deadheading time)
                E(i,4)=DH_dist(timetable(E(i,1),2),timetable(E(i,2),1)); % deadheading distances between trips
                if E(i,3)>0 % if trips can be joined
                    timetodepot=DH(arrivall(E(i,1),2),end)/60; % from trip i to depot
                    disttodepot=DH_dist(arrivall(E(i,1),2),end); % from trip i to depot
                    timefromdepot=DH(end,departl(E(i,2),2))/60; % time from depot to trip j in E(i)
                    distfromdepot=DH_dist(end,departl(E(i,2),2)); % time from depot to trip j in E(i)
                    timeatdepot=departl(E(i,2),1)-(timetodepot+timefromdepot+arrivall(E(i,1),1));
                    E(i,5)=(timeatdepot)*(timeatdepot>0); % add the time at depot if greater than zero
                    E(i,6)=disttodepot+distfromdepot; % deadhead distances from trip i to depot + depot to connecting trip
                end
            end
            toc
            DH_E=E(E(:,3)>=0  ,4);
            E=E(E(:,3)>=0 ,:);
            CNE=E(E(:,5)>minimumchargetime+chargeraccesstime,:);
            DH_E=DH_E(E(:,3)<8,:); % remove trip arcs if waiting time is over 1 hours
            E=E(E(:,3)<8,:); % remove trip arcs if waiting time is over 1 hours
            
            %interlining
            %no off service trips
            %             DH_E([arrivall(E(:,1),2)-departl(E(:,2),2)]==0,:)=[];
            %             E([arrivall(E(:,1),2)-departl(E(:,2),2)]==0,:)=[]; %removes connecting trips that need a deadhead without depot
            %             CNE=E(1,:); % no going to depot
            %
            CNE(:,4)=CNE(:,6); %  update deadheading distance with a trip to from depot
            vijd=Econs.*CNE(:,6); %Energy consumed during deadheading
            wijd=chargerpower*CNE(:,5);
            wijd(wijd>batterycapacity)=batterycapacity;
            chargetime=E(:,3); %available time until next trip in hours
            E=E(:,1:2);
            E1=E;
            DH1=DH_E;
            vij=Econs.*DH_E; %Energy consumed during deadheading
            wij=0*chargerpower*chargetime; % no charging at terminals
            
            wij(wij>batterycapacity)=batterycapacity;
            E1cost=ui(E1(:,1))+vij;
            charge=ui(E1(:,1))+vij-wij;
            charge_d=ui(CNE(:,1))+vijd-wijd;
            cost_d=ui(CNE(:,1))+vijd;
            %%
            V=[N N(end)+1 N(end)+2]; % end depot start and end nodes
            E3=[N' V(end)*ones(size(N,2),1)]; % add arcs from trips to depot
            DH_depot_1=DH_dist(end,1:end-1); % from depot
            DH_depot_2=[DH_dist(1:end-1,end)]';
            rhoij=DH_depot_1(1,timetable(E3(:,1),1))*Econs';
            E2=[V(end-1)*ones(size(N,2),1) [N']]; % pull out arcs from depot
            DH2=DH_depot_2(1,timetable(E2(:,2),1));
            piij=DH2*Econs';
            E=[E; N' V(end)*ones(size(N,2),1)]; % add arcs to depot / pull in arcs
            E=[E; V(end-1)*ones(size(N,2),1) [N'] ]; % add arcs from depot to trips
            C=[vij;rhoij';piij']; % energy usage of each arc - does not include trip energy use
            % Ci=[vij;rhoij';batterycapacity-piij'];
            
            CN1=E1cost; % energy consumption for each arc in E2
            pull_outcosts=zeros(size(E2,1),1);
            DH2=DH_depot_1(1,timetable(E2(:,2)))';
            
            CN2=zeros(size(E2,1),1);
            CN2(1:end,:)=DH2*Econs; % deadheading energy consumption from depot to trips
            pull_outcosts(1:end,:)=pulloutcost+CN2; % cost of starting new vehicle + deadhead
            
            
            DH3=DH_depot_2(1,timetable(E3(:,1),2))';
            CN3=zeros(size(E3,1),1); % energy consumption when returning to depot
            CN3(1:end,:)=(DH3.*Econs)+ui;  % energy consumption when returning to depot
            pull_incosts=zeros(size(E3,1),1);
            pull_incosts(1:end,:)=pull_incosts(1:end)+CN3; % cost of deadhead to depot
            % CN3(1:end,:)=(DH3.*Econs)+ui-Q;  % energy consumption when returning to depot
            
            
            CN=[CN1;CN2;CN3;ui(CNE(:,1))+vijd]; %energy consumed in each arc and trip
            DH_E=[DH1;DH2;DH3;CNE(:,4)]; %deadheading distance
            
            E0=[E1 E1cost;E2 pull_outcosts;E3 pull_incosts;CNE(:,1:2) cost_d]; %trip arcs with costs
            E0(:,4)=[charge;CN2;CN3;charge_d]; %CN;
            E0(:,5)=DH_E;
            E0(:,6)=CN; %energy discharged
            E0(:,7)=[wij;pull_outcosts*0;pull_incosts*0;wijd];
            E0(:,8)=[wij*0;pull_outcosts*0;pull_incosts*0;wijd*0+1];
            E0=sortrows(E0);
            
            Deadheading=E0(:,5);
            
            
            %% Equality/Flow Constraints
            Ce1=cell(size(N,2),2);
            Ce2=cell(size(N,2),2);
            tic
            for i=1:size(N,2) % for each timetabled trip
                colsi=find(E0(:,1) == i);
                colsj=find(E0(:,2) == i);
                Ce1{i,1}=colsi;
                Ce1{i,2}=ones(size(colsi,1),1)*i;
                Ce2{i,1}=colsj;
                Ce2{i,2}=ones(size(colsj,1),1)*i;
            end
            Ce1= cell2mat(Ce1);
            Ce1(:,3)=ones(size(Ce1,1),1);
            Ce2= cell2mat(Ce2);
            Ce2(:,3)=ones(size(Ce2,1),1);
            de1=find((E0(:,2) == size(V,2)));
            de1(:,2)=ones(size(de1,1),1)*((2*size(N,2))+1);
            de1(:,3)=ones(size(de1,1),1);
            de2=find((E0(:,1) == size(V,2)-1));
            de2(:,2)=1*ones(size(de2,1),1)*((2*size(N,2))+1);
            de2(:,3)=-1*ones(size(de2,1),1);
            de=[de1;de2];
            
            Aeq=sparse([Ce1(:,2);size(N,2)+Ce2(:,2);de(:,2)], [Ce1(:,1);Ce2(:,1);de(:,1)], [Ce1(:,3);Ce2(:,3);de(:,3)],(2*size(N,2))+1,size(E0,1));
            toc
            
            %% Range/Charge Constraints
            
            Ce0=cell(size(N,2),2);
            CE0=cell(size(N,2),2);
            CE01=cell(size(N,2),2);
            Ce01=cell(size(N,2),2);
            Ce1=cell(size(N,2),2);
            Ce2=cell(size(N,2),2);
            G=digraph(E0(:,1),E0(:,2));
            Q=batterycapacity;
            n=0;
            j=1;
            for i=1:size(N,2)
                [idx,idy]=inedges(G,i); % incoming edges for trip i
                [idx1,idy1]=outedges(G,i); % outgoing edges from trip i
                Ce1{i,1}=[idx(end); size(E0,1)+i];  % arc from depot
                Ce1{i,2}=[ones(size(idx(end),1),1)*i;i]; %
                Ce1{i,3}=[E0(idx(end),4);1];
                Ce2{i,1}=[idx1(end); size(E0,1)+i];
                Ce2{i,2}=[i;i];
                Ce2{i,3}=[E0(idx1(end),4);-1];
                
                if size(idx1,1)~=1
                    i1=idx1(1:end-1);
                    i2=size(E0,1)+idy1(1:end-1);
                    i3=ones(size(i1,1),1)*size(E0,1)+i;
                    k=size(i1,1);
                    i21=ones(size(i1,1),1).*[j:j+k-1]';
                    i22=ones(size(i1,1),1).*[j:j+k-1]';
                    i23=ones(size(i1,1),1).*[j:j+k-1]';
                    i31=E0(idx1(1:end-1),4)+Q;
                    i32=ones(size(i1,1),1);
                    i33=-1*ones(size(i1,1),1);
                    i41=E0(idx1(1:end-1),6);
                    i42=-1*ones(size(i1,1),1);
                    CE0(j:j+k-1,1)= num2cell([i1 i2 i3]',1)';
                    CE0(j:j+k-1,2)= num2cell([i21 i22 i23]',1)';
                    CE0(j:j+k-1,3)= num2cell([i31 i32 i33]',1)';
                    CE01(j:j+k-1,1)= num2cell([i1  i3]',1)';
                    CE01(j:j+k-1,2)= num2cell([i21  i22]',1)';
                    CE01(j:j+k-1,3)= num2cell([i41  i42]',1)';
                    j=j+k;
                end
            end
            n=size(CE0,1);
            Ce0=CE0;
            Ce01=CE01;
            Ce1= cell2mat(Ce1);
            Ce2= cell2mat(Ce2);
            Ce0= cell2mat(Ce0);
            Ce01= cell2mat(Ce01);
            A=sparse([Ce0(:,2);n+Ce2(:,2);n+size(N,2)+Ce1(:,2);n+size(N,2)+size(N,2)+Ce01(:,2)],[Ce0(:,1);Ce2(:,1);Ce1(:,1);Ce01(:,1)],[Ce0(:,3);Ce2(:,3);Ce1(:,3);Ce01(:,3)]);
            b=ones(n,1)*Q;
            b_1=zeros(size(N,2),1);
            % b=[b;b_1];
            b_2=ones(size(N,2),1)*Q;
            b_3=zeros(n,1);
            b=[b;b_1;b_2;b_3];
            Aeq(:,size(A,2))=0;
            intcon=(1:(length(E0))'); % optimization variables are possible arcs
            intN=(1:length(N))';
            cost=E0(:,3);
            cost(size(A,2),:)=0;
            lb=zeros(size(intcon,2),1);
            ub=ones(size(intcon,2),1);
            lbu=zeros(size(intN,1),1);
            ubu=Q*ones(size(intN,1),1);
            % define cost matrix
            intcon=[intcon';intN];
            model.vtype=[repmat('B', size(E0,1), 1); repmat('C', size(N',1), 1)];
            lb=[lb;lbu];
            ub=[ub;ubu];
            
            
            beq= ones( size(Aeq,1)-1,1);
            beq=[beq;0];
            
            
            
            
            if gurobi_solver~=1
                
                tic
                [x,fval]=intlinprog(cost,intcon,A,b,Aeq,beq,lb,ub);
                toc
            else
                tic
                addpath('C:\gurobi902\win64\matlab')
                model.obj=cost;
                model.A=[A;Aeq]; % - not efficient/slow speed for large matrices
                model.sense=[repmat('<', size(A,1), 1);repmat('=', size(Aeq,1), 1)];
                model.rhs=[b;beq];
                model.modelsense='min';
                params.outputflag=1;
                params.TimeLimit = 250; % 300 seconds to find optimal
                model.lb    = lb;
                model.ub    = ub;
                
                r=gurobi(model,params);
                toc
                x=r.x;
            end
            gap=r.mipgap*100;
            arcs=E0;
            x1=x;
            x=x(1:end-(size(timetable,1)),:);
            c=x1((end-(size(timetable,1)-1):end));
            %% form chains/blocks
            x=round(x);
            if isempty(x)
                disp('No solution found')
            else
                z=E0(x.*E0(:,1)>0,:); % filter used arcs
                %z_dh=Deadheading(x.*E(:,1)>0,:); % filter used dead heading arcs
                blocks=cell(8,1);
                blocks{1}=0;
                n=1;
                j=1;
                while ~isempty(j)
                    j=j(1);
                    while j<=V(end)
                        if blocks{1,n}==0
                            blocks{1,n}= z(z(:,2)==j,1:2); % trip joinings
                            blocks{2,n}=[timetable(z(z(:,2)==j,2),3)-DH(end,timetable(z(z(:,2)==j,2),1))/60 timetable(z(z(:,2)==j,2),3:4)]; % dispatch, departure and arrival times
                            blocks{3,n}= z(z(:,2)==j,5); % total distances (trip+deadhead)
                            blocks{5,n}=[size(DH,1) timetable(z(z(:,2)==j,2),1:2)]; %terminals
                            blocks{6,n}=[timetable(z(z(:,2)==j,2),3)-DH(end,timetable(z(z(:,2)==j,2),1))/60 timetable(z(z(:,2)==j,2),3)]; % deadheadtrips
                            blocks{7,n}=[timetable(z(z(:,2)==j,2),3:4)]; % timetabledtrips
                            blocks{10,n}=[ z(z(:,2)==j,6)]; % energy discharged
                            blocks{11,n}=[ z(z(:,2)==j,7)]; % energy charged
                            blocks{12,n}=[batterycapacity batterycapacity-blocks{10,n}(end)];
                            blocks{13,n}=[blocks{2,n}(1:end-1) ];
                            z(z(:,2)==j,:)=[];
                            j=z(z(:,1)==j,2);
                        elseif j==V(end)
                            blocks{3,n}=[blocks{3,n} z((z(:,1)==blocks{1,n}(end) & z(:,2)==j),5)];  % total distances (trip+deadhead)
                            blocks{10,n}=[blocks{10,n} E0(E0(:,1)==blocks{1,n}(end) & E0(:,2)==V(end),4)]; % energy discharged
                            z_end=E0(E0(:,1)==blocks{1,n}(end) & E0(:,2)==V(end),:);
                            blocks{12,n}=[blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)-z_end(:,6)];
                            blocks{13,n}=[blocks{13,n} timetable(z_end(:,1),3) timetable(blocks{1,n}(end),4)+DH(timetable(blocks{1,n}(end),2),end)/60];
                            % blocks{12,n}=[blocks{12,n} blocks{12,n}(end)-E0(E0(:,1)==blocks{1,n}(end) & E0(:,2)==V(end),4)];
                            blocks{5,n}=[blocks{5,n} timetable(blocks{1,n}(end),2) size(DH,1) ];%terminals
                            blocks{6,n}=[blocks{6,n} blocks{2,n}(end) blocks{2,n}(end)+DH(timetable(blocks{1,n}(end),2),end)/60]; % deadheadtrips
                            %  z_dh(z(:,1)==blocks{1,n}(end),:)=[];
                            blocks{2,n}=[blocks{2,n}  timetable(blocks{1,n}(end),4)+DH(timetable(blocks{1,n}(end),2),end)/60];  % dispatch, departure and arrival times
                            z(z(:,1)==blocks{1,n}(end),:)=[];
                            blocks{1,n}=[blocks{1,n} V(end)];
                            j=j+1;
                        else
                            blocks{1,n}=[blocks{1,n} z(z(:,2)==j,2)]; % trip joinings
                            blocks{10,n}=[blocks{10,n} z(z(:,2)==j,6)]; % energy discharged
                            blocks{3,n}=[blocks{3,n} timetable(z(z(:,2)==j,2),5)+z(z(:,2)==j,5)];  % total distances (trip+deadhead)
                            
                            
                            
                            if timetable(j,1)==blocks{5,n}(end) % if next trip from same stop
                                % add waiting time
                                blocks{8,n}=[blocks{8,n} blocks{2,n}(end) timetable(z(z(:,2)==j,2),3)]; % waitingtrips
                                
                                charge_block=blocks{12,n}(end)-z(z(:,2)==j,6)+z(z(:,2)==j,7);
                                if charge_block>0
                                    
                                    if charge_block>batterycapacity
                                        blocks{12,n}= [blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)-z(z(:,2)==j,6)]; % amount discharged
                                        e_charged=batterycapacity-(blocks{12,n}(end));
                                        blocks{12,n}= [blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)+e_charged];
                                        blocks{11,n}=[blocks{11,n} e_charged];
                                        charge_endtime=e_charged/chargerpower;
                                        if isnan(charge_endtime)
                                            charge_endtime=0;
                                        end
                                        blocks{13,n}= [blocks{13,n} timetable(z(z(:,2)==j,1),3:4)];
                                        blocks{13,n}= [blocks{13,n} blocks{13,n}(end) timetable(z(z(:,2)==j,1),4)+charge_endtime];
                                        
                                        
                                    else
                                        e_charged=z(z(:,2)==j,7);
                                        blocks{11,n}=[blocks{11,n} e_charged]; % energy charged
                                        %  blocks{12,n}=[blocks{12,n} blocks{12,n}(end)-z(z(:,2)==j,6) blocks{12,n}(end)-z(z(:,2)==j,6)+z(z(:,2)==j,7)];
                                        blocks{12,n}= [blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)-z(z(:,2)==j,6)]; % update discharge level
                                        blocks{12,n}= [blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)+e_charged]; %update charge level after trip i
                                        % calculate charge duration
                                        charge_endtime=e_charged/chargerpower;
                                        if isnan(charge_endtime)
                                            charge_endtime=0;
                                        end
                                        blocks{13,n}= [blocks{13,n} timetable(z(z(:,2)==j,1),3:4) timetable(z(z(:,2)==j,1),4) timetable(z(z(:,2)==j,1),4)+charge_endtime];
                                        % no waiting between charging and next trip blocks{13,n}= [blocks{13,n} blocks{13,n}(end) timetable(z(z(:,2)==j,2),3)];
                                    end
                                end
                                
                            else
                                %deadhead
                                blocks{6,n}=[blocks{6,n} blocks{2,n}(end) blocks{2,n}(end)+DH(timetable(blocks{1,n}(end-1),2),timetable(blocks{1,n}(end),1))/60];% deadheadtrips
                                % add waiting time
                                blocks{8,n}=[blocks{8,n} blocks{6,n}(end) timetable(z(z(:,2)==j,2),3)]; % waitingtrips
                                
                                charge_block=blocks{12,n}(end)-z(z(:,2)==j,6)+z(z(:,2)==j,7);
                                if charge_block>0
                                    
                                    if charge_block>batterycapacity
                                        blocks{12,n}= [blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)-z(z(:,2)==j,6)]; % amount discharged
                                        e_charged=batterycapacity-(blocks{12,n}(end));
                                        blocks{12,n}= [blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)+e_charged];
                                        blocks{11,n}=[blocks{11,n} e_charged];
                                        charge_endtime=e_charged/chargerpower;
                                        if isnan(charge_endtime)
                                            charge_endtime=0;
                                        end
                                        blocks{13,n}= [blocks{13,n} timetable(z(z(:,2)==j,1),3) blocks{6,n}(end)];
                                        blocks{13,n}= [blocks{13,n} blocks{13,n}(end) blocks{6,n}(end)+charge_endtime];
                                        
                                        
                                    else
                                        e_charged=z(z(:,2)==j,7);
                                        blocks{11,n}=[blocks{11,n} e_charged]; % energy charged
                                        %  blocks{12,n}=[blocks{12,n} blocks{12,n}(end)-z(z(:,2)==j,6) blocks{12,n}(end)-z(z(:,2)==j,6)+z(z(:,2)==j,7)];
                                        blocks{12,n}= [blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)-z(z(:,2)==j,6)]; % update discharge level
                                        blocks{12,n}= [blocks{12,n} blocks{12,n}(end) blocks{12,n}(end)+e_charged]; %update charge level after trip i
                                        % calculate charge duration
                                        charge_endtime=e_charged/chargerpower;
                                        if isnan(charge_endtime)
                                            charge_endtime=0;
                                        end
                                        blocks{13,n}= [blocks{13,n} timetable(z(z(:,2)==j,1),3) blocks{6,n}(end) blocks{6,n}(end) blocks{6,n}(end)+charge_endtime];
                                        % no waiting between charging and next trip blocks{13,n}= [blocks{13,n} blocks{13,n}(end) timetable(z(z(:,2)==j,2),3)];
                                    end
                                end
                                
                                
                            end
                            blocks{2,n}=[blocks{2,n} timetable(z(z(:,2)==j,2),3:4)];
                            
                            blocks{5,n}=[blocks{5,n} timetable(z(z(:,2)==j,2),1:2)]; % deadheadtrips
                            blocks{7,n}=[blocks{7,n} timetable(z(z(:,2)==j,2),3:4)]; % timetabledtrips
                            % z_dh(z(:,2)==j,:)=[];
                            z(z(:,2)==j,:)=[];
                            j=z(z(:,1)==j,2);
                        end
                    end
                    blocks{4,n}= sum(blocks{3,n}); % total distance driven of vehicle
                    blocks{9,n}= (batterycapacity-cumsum(blocks{3,n}))/batterycapacity*100; %wrong add energy consumption
                    blocks{14,n}= (blocks{12,n})/batterycapacity; % total distance driven of vehicle
                    %    blocks{12,n}= batterycapacity-cumsum(blocks{10,n});%;
                    n=n+1;
                    blocks{1,n}=0;
                    j=sortrows(z(:,1));
                    %  reshape(tb{7,1},2,6)'
                end
            end
            toc
            
        end
         function [arcs,var,nodes,timetable]=MDVSP_TSN(obj,vehicle,timetable,DH,DH_dist,Econs,chargingpower,numdepots)
             
             % timetable [origin-stop destination-stop start-time end-time depot]
             D=numdepots; % list of depots
             % Create pull in/pull out arcs
             S=unique([timetable(:,1);timetable(:,2)]); %unique stops
             numterminals=size(S,1); %number of terminals
             timetable(:,6)=1:size(timetable,1);
             %Econs=1.2; % energy consumption
             %%
             % for each station i add pull-in arc from depot
             PI=[];%pull in arcs
             PO=[];%pull out arcs
             DH_arcs=[]; % dead head arcs
             DH_firstarcs=[]; %dead head arcs first match (TSN - Kliewer 2006)
             DH_firstarcsn=[]; %dead head arcs first match (TSN - Kliewer 2006)
             DH_latestfirstarcs=[]; %dead head arcs latest first match (TSN - Kliewer 2006)
             W_arcs=[]; %waiting arcs
             W_arcsn=[]; %waiting arcs
             CI=[]; %charging station in arcs
             CO=[]; %charging station pull out arc
             N=[]; % no of nodes
             N1=[1:size(timetable,1)]'; % nodes with outgoing trips
             N1(:,2:4)=[timetable(:,[1 3 6])];
             N2=[size(timetable,1)+1:size(timetable,1)*2]';
             N2(:,2:4)=[timetable(:,[2 4 6])];
             numdepots=length(D);
             DHarcs=[];
             W_arcn=[];
             dhc=1;
             triparcs=[];
             N=[N1;N2];
             D0_1=[];
             D0_2=[];
             D0_3=[];
             D0_4=[];
             D0_5=[];
             Dend=[];
             minbuses=[.412]; % minimum number of vehicles constrain
             %%
             
             
             for d=1:length(D)
                 trip_arcs=[N1(:,1) N2(:,1) timetable(:,3:6)];        % trip arcs
                 
                 for z=1:size(trip_arcs,1) % for each trip define pull in/out arcs
                     PO=[PO;D(d) trip_arcs(z,1) trip_arcs(z,3)-DH(numterminals+d,N1(trip_arcs(z,1),2))./60 trip_arcs(z,3)  ...
                         DH_dist(numterminals+d,N1(trip_arcs(z,1),2)) ones(size(trip_arcs(z,6),1),1).*d ];
                     PI=[PI;trip_arcs(z,2) D(d) trip_arcs(z,4) trip_arcs(z,4)+DH(N1(trip_arcs(z,1),2),numterminals+d)./60 ...
                         DH_dist(N1(trip_arcs(z,1),2),numterminals+d) ones(size(trip_arcs(z,6),1),1).*d ];
                 end
                 
                 for k=1:length(S)
                     
                     trips_k=timetable(timetable(:,2)==S(k),:); %arrival at stop k
                     for i=1:size((trips_k),1) %for each arriving trip at stop k
                         atimek=trips_k(i,4); % arrival time of trip i at k
                         %find the first connecting trip of trip i at all stops
                         for l=1:length(S) %find all the possible deadheading trips with stops l
                             if l~=k
                                 trips_l=timetable(timetable(:,1)==S(l),:); % find outgoing trips from stop l
                                 connectiontime=atimek+DH(k,l)./60; % connection time
                                 dh_distance=DH_dist(k,l);
                                 compatible_trips=trips_l(((trips_l(:,3)-connectiontime)>=0),:); % check if compatible in time
                                 if ~isempty(compatible_trips)
                                     %                        DH_firstarcsn=[DH_firstarcsn;N2(trips_k(i,6),1)  N1(compatible_trips(1,6),1) connectiontime compatible_trips(1,3) ...
                                     %                          dh_distance    ];
                                     
                                     if dh_distance<DH_dist(k,numterminals+d)+DH_dist(numterminals+d,l) % check if a return trip to the depot is 'cheaper'
                                         
                                         DH_firstarcsn=[DH_firstarcsn;N2(trips_k(i,6),1)  N1(compatible_trips(1,6),1) trips_k(i,4) connectiontime ...
                                             dh_distance    ];
                                     else
                                         disp('going to the depot yeaaahh')
                                     end
                                 end
                             end
                         end
                         
                     end
                     Wn=sortrows(N(N(:,2)==k,:),3); % find all nodes from terminal k (both departure and arrivals)
                     Wn=[Wn(1:end-1,1) Wn(2:end,1) Wn(1:end-1,3) Wn(2:end,3)];
                     W_arcsn=[W_arcsn;Wn];
                 end
                 
                 if isempty(Dend)
                     PO(:,1)=N2(end,1)+1:N2(end,1)+size(N1,1);
                     PI(:,2)=PO(end,1)+1:PO(end,1)+size(N1,1);
                 else
                     PO(:,1)=Dend(:,2)+1:Dend(:,2)+size(N1,1);
                     PI(:,2)=PO(end,1)+1:PO(end,1)+size(N1,1);
                 end
                 Dend=PI(end,:);
                 Dn=sortrows([PO;PI(:,[2,1,4,3,5,6])],3);
                 Dn=[Dn(1:end-1,1) Dn(2:end,1) Dn(1:end-1,3) Dn(2:end,3) ...
                     zeros(size([PO;PI],1)-1,1)  ones(size([PO;PI],1)-1,1).*d];
                 D1=[max(Dn(:,2))+1 Dn(1,1) 0 Dn(1,3) 0 d];
                 Dend=[Dn(end,2) max(Dn(:,2))+1 0 Dn(end,4) 0 d];
                 D0_1=[D0_1;PO ];
                 D0_2=[D0_2;PI ];
                 D0_3=[D0_3;Dn ];
                 D0_4=[D0_4;D1 ];
                 D0_5=[D0_5;Dend ];
                 
                 DHarcs=[DHarcs;DH_firstarcsn ones(size(DH_firstarcsn,1),1).*d];
                 DH_firstarcsn=[];
                 W_arcn=[W_arcn; W_arcsn zeros(size(W_arcsn,1),1) ones(size(W_arcsn,1),1).*d];
                 W_arcsn=[];
                 triparcs=[triparcs;trip_arcs(:,1:5) ones(size(trip_arcs(:,5),1),1).*d];
                 trip_arcs=[];
                 PO=[];
                 PI=[];
             end
             %% add costs
             PO=D0_1;
             PI=D0_2;
             Dn=D0_3;
             D1=D0_4;
             Dend=D0_5;
             PO(:,7)=(PO(:,5))*Econs;  % fixed cost of starting a new vehicle
             %PO(:,1)=N2(end,1)+PO(:,1);
             PO(:,8)=-1;
             %PI(:,2)=N2(end,1)+PI(:,2);
             PI(:,8)=-2;
             PI(:,7)=(PI(:,5))*Econs*0.85; % hourly cost of driving
             Dn(:,7)=(Dn(:,4)-Dn(:,3))*0;  % waiting cost at depot
             Dn(:,8)=0;
             D1(:,7)=50000;
             Dend(:,7)=50;
             D1(:,8)=-5;
             Dend(:,8)=5;
             % D1(:,5)=D1(:,3);
             % Dend(:,5)=Dend(:,3);
             triparcs(:,7)=(triparcs(:,5))*Econs; % hourly cost of driving
             triparcs(:,8)=1;
             DHarcs(:,7)=(DHarcs(:,5))*Econs*.7; % fixed cost of deadheading
             W_arcn(:,7)=(W_arcn(:,4)-W_arcn(:,3))*05 + .5; %no cost of waiting
             DHarcs(:,8)=2; % id of deadhead/wait/trip/pullout
             W_arcn(:,8)=0;
             
             E=[D1;Dn;Dend;PO;PI;triparcs;DHarcs;W_arcn]; % all arcs
             %E=sortrows(E,3);
             E(:,9)=1:size(E,1);
             %% Add constraints
             N=[N1;N2];
             Ce1=cell(size(N1,2),2);
             Ce2=cell(size(N1,2),2);
             Ce3=[];
             n=1;
             tic
             for i=1:size(N,1) % for each timetabled trip
                 for j=1:numdepots
                     colsi=E(E(:,1) == N(i,1) & E(:,6)==j,9); % find outgoing arcs from node i for depot j
                     colsj=E(E(:,2) == N(i,1) & E(:,6)==j ,9); % find incoming arcs to node i
                     Ce1{n,1}=colsi;
                     Ce1{n,2}=ones(size(colsi,1),1)*n;
                     Ce2{n,1}=colsj;
                     Ce2{n,2}=ones(size(colsj,1),1)*n;
                     n=n+1;
                 end
                 if i<=size(timetable,1)
                     id=find(E(:,1)==i & E(:,2)==N2(i,1));
                     Ce3=[Ce3;id ones(size(id,1),1)*i];
                 end
             end
             
             Ce1= cell2mat(Ce1);
             Ce1(:,3)=-ones(size(Ce1,1),1);
             Ce2= cell2mat(Ce2);
             Ce2(:,3)=ones(size(Ce2,1),1);
             Ce3(:,2)=Ce3(:,2)+n-1;
             Ce3(:,3)=ones(size(Ce3,1),1);
             
             de1=zeros(size(N1,2)*numdepots,3);
             de2=zeros(size(N2,2)*numdepots,3);
             n=1;
             Did1= unique(PO(:,1));
             Did2= unique(PI(:,2));
             %% depot constraints
             n1=1;
             x=1;
             
             % add minumum bus constraint
             for j=1:numdepots
                 for i=Dn(1,1):D1(end,1)
                     
                     colsi=E(E(:,1) == i & E(:,6)==j,9);
                     colsj=E(E(:,2) == i & E(:,6)==j,9);
                     de1(n:n+size(colsi,1)-1,:)=[colsi  ones(size(colsi,1),1).*x ones(size(colsi,1),1)];
                     de2(n1:n1+size(colsj,1)-1,:)=[colsj  ones(size(colsj,1),1).*x -ones(size(colsj,1),1)];
                     n=n+size(colsi,1);
                     n1=n1+size(colsj,1);
                     x=x+1;
                 end
                 colsi=E(E(:,8)==-5,9); % depot -vehicle pullout arcs
                 Aer(j,:)=colsi;
             end
             de=[de1;de2];
             de(:,2)=de(:,2)+Ce3(end,2);
             Aeq=sparse([Ce1(:,2);Ce2(:,2);Ce3(:,2);de(:,2)],... % rows [out going arcs;incoming arcs; depot arcs]
                 [Ce1(:,1);Ce2(:,1);Ce3(:,1);de(:,1)], ... % columns / variable ids
                 [Ce1(:,3);Ce2(:,3);Ce3(:,3);de(:,3)],... % value
                 de(end,2),size(E,1)); % matrix size rowsxcolumns
             
             A=sparse(numdepots,Aer,-1,1,size(E,1));
             b=[-minbuses];
             cost=E(:,7);
             intcon=(1:(length(E))'); % optimization variables are possible arcs
             
             intN=(1:length(N))';
             lb=zeros(size(intcon,2),1);
             ub=ones(size(intcon,2),1);
             
             
             beq= zeros( size(Aeq,1),1);
             beq(Ce3(:,2),1)=1;%size(timetable,1);
             model.vtype=[repmat('I', size(E,1), 1)];
             model.obj=cost;
             model.A=[A;Aeq]; % - not efficient/slow speed for large matrices
             model.sense=[repmat('<', size(A,1), 1);repmat('=', size(Aeq,1), 1)];
             model.rhs=[b;beq];
             model.modelsense='min';
             params.outputflag=1;
             model.lb    = lb;
             model.ub    = ub*Inf;
               addpath('C:\gurobi901\win64\matlab')
               
             r=gurobi(model,params);
             x=r.x;
             z=E(round(x)>0,:);
             toc
             arcs=E;var=x;
             z1=z;
             nodes=N;
         end
         function vehicles=dispatchsimfifo(obj,arcs,var,nodes,timetable,Econs)
             %% This function decomposes the optimal flow using a heuristic by assigning
             % trips to buses that have the highest state of charge
             %%
             tic
             chargerpower=150; %charger power in kW
             timestep = 1; % timestep is 1 min
             usedarcs=arcs(var>0,:);
             x=var(var>0);
             usedarcs(:,end+1)=x;
             departures=sortrows([usedarcs(:,1) usedarcs(:,3)  usedarcs(:,6)  usedarcs(:,8)  usedarcs(:,9:10)],2); % nodeid, departime, depotnum, status, arcno
             arrivals=sortrows([usedarcs(:,2) usedarcs(:,4) usedarcs(:,6)  usedarcs(:,8)  usedarcs(:,9:10)],2); % nodeid, departime, depotnum, status, arcno
             terminals=unique(nodes(:,2));
             depots=unique(usedarcs(:,6));
             events=unique(sortrows([departures(:,2);arrivals(:,2)]));
             events=(events(2:end,:));
             %n_tripsteps = round(t_trip/timestep);
             energyconsumed = timestep*4; %Energy consumed in the Schedule.timestep [kWh]
             energycharged = timestep*chargerpower; %Charged energy per Schedule.timestep [kWh]
             vehicledemand=zeros(size(usedarcs(usedarcs(:,8)==-1,:),1)*2,2);
             nd=1; % counter for depot departures/arrivals
             batterycapacity=350; % 270 kwh
             %t_load = cell2mat(Demand(1,1));
             %demand_load = cell2mat(Demand(1,2));
             %[t_load, index] = unique(t_load);
             %t = departures(3,2):timestep: arrivals(end,2); %timevector
             t=floor(departures(3,2)*60):timestep:round((24+departures(3,2))*60)-1; % until next day first trip
             % t=sort([t';events]);
             % t=unique(t);
             if depots==1
                 numVehicles=var(1);
             end
             Nchargers=400;
             Usedchargers=0;
             % Econs=1.2;
             Econsempty=0.9*Econs;
             Charging = zeros(numVehicles,4); % charger number, available, vehiclenum, energy charged,
             Availablechargers = Nchargers;
             vehicles=cell(6,numVehicles);
             %Parameter for finding continuum
             threshold = -100; %Abort threshold for average energy drained at the end of the day
             iterationcounter = 0;
             continuous = false;
             flag=0; % flag if solution is infeasible
            % while not(continuous)
                 %Preallocate matrices
                 Ebat = zeros(numVehicles,length(t)); %Amount of energy in battery [kWh]
                 Driving = zeros(numVehicles,length(t)); %Driving = 1, not-driving = 0
                 Charging = zeros(numVehicles,length(t)); %Charging = 1, not-charging = 0
                 vehiclestarted = zeros(numVehicles,length(t));
                 vehiclesdone = zeros(numVehicles,length(t));
                 chargestarted = zeros(numVehicles,length(t));
                 chargedone = zeros(numVehicles,length(t));
                 %modulenode = zeros(numVehicles,length(t));
                 locations=zeros(terminals(end)+depots,numVehicles); %,length(t)
                 queue=zeros(numVehicles,6,terminals(end)+depots);% arrival time, vehicle no, dist travelled, energy discharge, bat capacity, soc
                 chargingqueue=zeros(Nchargers,6); %vehicle,  charge start time, endtime, energy charged
                 % Initialize depots with all vehicles
                 for i=1:depots
                     locations(terminals(end)+i,:)=ones(numVehicles(1),1);
                     queue(:,:,terminals(end)+i)=[zeros(numVehicles,1) [1:numVehicles]' zeros(numVehicles,1) zeros(numVehicles,1)...
                         ones(numVehicles,1)*batterycapacity ones(numVehicles,1)];
                 end
                 
                 for i = 1:length(t)
                     % check if any departure/arrival/wait event occurs
                     currentevents=usedarcs(floor(usedarcs(:,3)*60)==t(i),:);
                     c01=[];
                     for j=1:size(currentevents,1) % some arcs are used more than once therefore duplicate events to assign correctly
                         % repeat only arcs that are not waiting at depot
                         if currentevents(j,1)<=nodes(end,1)
                             c01=[c01;repmat(currentevents(j,:),currentevents(j,10),1)];
                         elseif currentevents(j,8)==-1 & currentevents(j,1)>nodes(end,1)% duplicate pullout arcs if used more tahn once
                             c01=[c01;repmat(currentevents(j,:),currentevents(j,10),1)];
                         else
                             c01=currentevents(j,:);
                         end
                         if j==size(currentevents,1)
                             c01=sortrows(c01,8);
                         end
                     end
                     currentevents=c01;
                     if ~isempty(currentevents)
                         currentevents=sortrows(currentevents,3);
                     end
                     
                     %for each event
                     for j = 1:size(currentevents,1)
                         
                         if currentevents(j,8)==-1 % if the event is to pull out of depot
                             %available = find(Driving(:,i)==0); %look for available modules
                             available= queue(queue(:,2,terminals(end)+currentevents(j,6))~=0,:,...
                                 terminals(end)+currentevents(j,6));% find available vehicles at location/depot
                             if isempty(available) % No idling/waiting vehicles
                                 % Find vehicle with the highest SOC at the charging
                                 % station
                                 [~,Index]=max(chargingqueue(:,6)); % find available charger used least
                                 vehicle = chargingqueue(Index,:); % select vehicle
                                 vehicle(1)=t(i)/60; % return to depot after disconnection
                                 % add back to depot queue
                                 queue(vehicle(1,2),:,end)=vehicle; % remove from queue at current stop
                                 available=vehicle;
                                 chargingqueue(Index,:)=0;
                                 Availablechargers=Availablechargers+1;
                                 Usedchargers=Usedchargers-1;
                             end
                             [~,Index]=max(available(:,4)); % find vehicle with max charge levels
                             vehicle = available(Index,2);
                             vehicle= queue(vehicle,:,end);
                             queue(vehicle(2),:,terminals(end)+currentevents(j,6))=0;
                             vehicle(1)=currentevents(j,4);
                             vehicle(3)=currentevents(j,5)+vehicle(3); % distance travelled
                             vehicle(4)=vehicle(4)-currentevents(j,5)*Econsempty;
                             vehicle(6)=(vehicle(5)+vehicle(4))/vehicle(5); %SOC
                             if vehicle(6)<0
                                 flag=1;
                             end
                             % update queue at destination
                             queue(vehicle(2),:,nodes(currentevents(j,2),2))=vehicle;
                             % update vehicle
                             vehicles{1,vehicle(2)}=[vehicles{1,vehicle(2)} currentevents(j,1:2)]; % update used arcs
                             vehicles{2,vehicle(2)}=[vehicles{2,vehicle(2)} currentevents(j,3:4)]; % update time
                             vehicles{3,vehicle(2)}=[vehicles{3,vehicle(2)} -1]; % update status
                             vehicles{4,vehicle(2)}=[vehicles{4,vehicle(2)} currentevents(j,5)]; % update km
                             vehicles{5,vehicle(2)}=sum(vehicles{4,vehicle(2)}); % total km
                             vehicles{6,vehicle(2)}=[vehicles{6,vehicle(2)}-currentevents(j,5)*1.2]; % soc
                             vehicledemand(nd,:)=[t(i) 1];
                             nd=nd+1;
                         elseif currentevents(j,8)==-2 % if the event is a pull in
                             terminal= nodes(currentevents(j,1),2); % terminal of trip departure
                             % find vehicles available at terminal
                             available=queue(queue(:,2,terminal)~=0,:,terminal);
                             available=available(available(:,1)<=currentevents(j,3),:); % vehicles that have arrived before departure time of event j
                             if isempty(available)
                                 2
                             end
                             [~,Index]=min(available(:,1)); % xx find vehicle with min charge levels xx % send vehicle with the lowest SOC back to depot
                             vehicle = available(Index,:); % select vehicle
                             vehicle=queue(vehicle(1,2),:,terminal);
                             queue(vehicle(1,2),:,terminal)=0; % remove from queue at current stop
                             vehicle(1)=currentevents(j,4);
                             vehicle(3)=currentevents(j,5)+vehicle(3); % distance travelled
                             vehicle(4)=vehicle(4)-currentevents(j,5)*Econsempty+0;
                             vehicle(6)=(vehicle(5)+vehicle(4))/vehicle(5); %SOC
                             if vehicle(6)<0
                                 flag=1;
                             end
                             % update queue at destination
                             queue(vehicle(2),:,terminals(end)+currentevents(j,6))=vehicle;
                             %
                             distance=vehicle(1,3)+currentevents(j,5);  % total distance at next stop
                             charge=vehicle(1,4)-currentevents(j,5)*Econs+0; % charge level at next stop
                             
                             
                             vehicles{1,vehicle(1,2)}=[vehicles{1,vehicle(1,2)} currentevents(j,1:2)]; % update used arcs
                             vehicles{2,vehicle(1,2)}=[vehicles{2,vehicle(1,2)} currentevents(j,3:4)]; % update time
                             vehicles{3,vehicle(1,2)}=[vehicles{3,vehicle(1,2)} -2]; % update status
                             vehicles{4,vehicle(1,2)}=[vehicles{4,vehicle(1,2)} currentevents(j,5) ]; % update km
                             vehicles{5,vehicle(1,2)}=sum(vehicles{4,vehicle(1,2)}); % total km
                             vehicles{6,vehicle(1,2)}=[vehicles{6,vehicle(1,2)} charge]; % soc
                             
                             vehicledemand(nd,:)=[t(i) -1];
                             nd=nd+1;
                             
                         elseif currentevents(j,8)==1 % if the event is start a timetabled trip
                             terminal= nodes(currentevents(j,1),2); % terminal of trip departure
                             % find vehicles available at terminal
                             available=queue(queue(:,2,terminal)~=0,:,terminal);
                             available=available(available(:,1)<=currentevents(j,3),:); % vehicles that have arrived before departure time of event j
                             [~,Index]=min(available(:,1)); % find vehicle waiting the longest
                             vehicle = available(Index,:); % select vehicle
                             vehicle=queue(vehicle(1,2),:,terminal);
                             queue(vehicle(1,2),:,terminal)=0; % remove from queue at current stop
                             % update queue at destination
                             distance=vehicle(1,3)+currentevents(j,5);  % total distance at next stop
                             charge=vehicle(1,4)-currentevents(j,5)*Econs; % charge level at next stop
                             vehicle(1)=currentevents(j,4);
                             vehicle(3)=distance; % distance travelled
                             vehicle(4)=charge;
                             vehicle(6)=(vehicle(5)+vehicle(4))/vehicle(5); %SOC
                             if vehicle(6)<0
                                 flag=1;
                             end
                             
                             queue(vehicle(1,2),:,nodes(currentevents(j,2),2))=vehicle;
                             endind= find(t==currentevents(j,4)); % endtime of event
                             Driving(vehicle(1,2),i:endind) = 1; %Start this module on a trip
                             laststarted(vehicle(1,2)) = i; %Mark point for checking which modules are done driving
                             vehiclestarted(vehicle(1,2),i:endind) = 1; %Mark point for plotting markers
                             vehicles{1,vehicle(1,2)}=[vehicles{1,vehicle(1,2)} currentevents(j,1:2)]; % update used arcs
                             vehicles{2,vehicle(1,2)}=[vehicles{2,vehicle(1,2)} currentevents(j,3:4)]; % update time
                             vehicles{3,vehicle(1,2)}=[vehicles{3,vehicle(1,2)} 1]; % update status
                             vehicles{4,vehicle(1,2)}=[vehicles{4,vehicle(1,2)} currentevents(j,5)]; % update km
                             vehicles{5,vehicle(1,2)}=sum(vehicles{4,vehicle(1,2)}); % total km
                             vehicles{6,vehicle(1,2)}=[vehicles{6,vehicle(1,2)} charge]; % soc
                             
                         elseif currentevents(j,8)==2 % if the event is to deadhead
                             terminal= nodes(currentevents(j,1),2); % terminal of trip departure
                             % find vehicles available at terminal
                             available=queue(queue(:,2,terminal)~=0,:,terminal);
                             available=available(available(:,1)<=currentevents(j,3),:); % vehicles that have arrived before departure time of event j
                             [~,Index]=min(available(:,4)); % find vehicle that arrived first
                             vehicle = available(Index,:); % select vehicle
                             vehicle=queue(vehicle(1,2),:,terminal);
                             queue(vehicle(1,2),:,terminal)=0; % remove from queue at current stop
                             % update queue at destination
                             distance=vehicle(1,3)+currentevents(j,5);  % total distance at next stop
                             charge=vehicle(1,4)-currentevents(j,5)*Econsempty; % charge level at next stop
                             vehicle(1)=currentevents(j,4);
                             vehicle(3)=distance; % distance travelled
                             vehicle(4)=charge;
                             vehicle(6)=(vehicle(5)+vehicle(4))/vehicle(5); %SOC
                             if vehicle(6)<0
                                 flag=1;
                             end
                             
                             queue(vehicle(1,2),:,nodes(currentevents(j,2),2))=vehicle;
                             endind= find(t==currentevents(j,4)); % endtime of event
                             Driving(vehicle(1,2),i:endind) = 1; %Start this module on a trip
                             laststarted(vehicle(1,2)) = i; %Mark point for checking which modules are done driving
                             vehiclestarted(vehicle(1,2),i:endind) = 1; %Mark point for plotting markers
                             vehicles{1,vehicle(1,2)}=[vehicles{1,vehicle(1,2)} currentevents(j,1:2)]; % update used arcs
                             vehicles{2,vehicle(1,2)}=[vehicles{2,vehicle(1,2)} currentevents(j,3:4)]; % update time
                             vehicles{3,vehicle(1,2)}=[vehicles{3,vehicle(1,2)} 2]; % update status
                             vehicles{4,vehicle(1,2)}=[vehicles{4,vehicle(1,2)} currentevents(j,5)]; % update km
                             vehicles{5,vehicle(1,2)}=sum(vehicles{4,vehicle(1,2)}); % total km
                             vehicles{6,vehicle(1,2)}=[vehicles{6,vehicle(1,2)} charge]; % soc
                             
                         elseif currentevents(j,8)==0 % if the event is to keep waiting
                             
                             if currentevents(j,1)<=nodes(end,1) % if the waiting event is not at the depot
                                 
                                 terminal= nodes(currentevents(j,1),2); % terminal of trip departure
                                 % find vehicles available at terminal
                                 available=queue(queue(:,2,terminal)~=0,:,terminal);
                                 available=available(available(:,1)<=currentevents(j,3),:); % vehicles that have arrived before departure time of event j
                                 if size(available,1)>1
                                     2;
                                 end
                                 [~,Index]=max(available(:,4)); % find vehicle that has the least soc
                                 vehicle = available(Index,:); % select vehicle
                                 vehicle=queue(vehicle(1,2),:,terminal);
                                 queue(vehicle(1,2),:,terminal)=0; % remove from queue at current stop
                                 % update queue at destination
                                 distance=vehicle(1,3)+currentevents(j,5);  % total distance at next stop
                                 charge=vehicle(1,4)-currentevents(j,5)*Econs; % charge level at next stop
                                 % vehicle(1)=currentevents(j,4);
                                 vehicle(3)=distance; % distance travelled
                                 vehicle(4)=charge;
                                 vehicle(6)=(vehicle(5)+vehicle(4))/vehicle(5); %SOC
                                 queue(vehicle(1,2),:,nodes(currentevents(j,2),2))=vehicle;
                                 if vehicle(6)<0
                                     flag=1;
                                 end
                                 
                                 % if no charging:
                                 vehicles{1,vehicle(1,2)}=[vehicles{1,vehicle(1,2)} currentevents(j,1:2)]; % update used arcs
                                 vehicles{2,vehicle(1,2)}=[vehicles{2,vehicle(1,2)} currentevents(j,3:4)]; % update time
                                 vehicles{3,vehicle(1,2)}=[vehicles{3,vehicle(1,2)} 3]; % update status
                                 vehicles{4,vehicle(1,2)}=[vehicles{4,vehicle(1,2)} currentevents(j,5)]; % update km
                                 vehicles{5,vehicle(1,2)}=sum(vehicles{4,vehicle(1,2)}); % total km
                             end
                             
                             
                         elseif currentevents(j,8)==5 % if the event is to end service
                             []; % do nothing - make sure all vehicles return to depot
                         end
                         
                     end
                     
                     %% assign vehicles to available chargers
                     availablevehicles=queue(queue(:,2,end)~=0,:,end);
                     availablevehicles=availablevehicles(availablevehicles(:,6)...
                         <1,:);
                     if Availablechargers>0 & size(availablevehicles,1)>0 % if there is a charger available and vehicle needs charging
                         if size(availablevehicles,1)==1
                             
                             Availablechargers=Availablechargers-1; % connect vehicle to charger
                             Usedchargers=Nchargers-Availablechargers;
                             [~,Index]=min(chargingqueue(:,1)); % find available charger used least
                             % find vehicle with lowest SOC to connect to charger
                             vehicle=availablevehicles;
                             queue(vehicle(2),:,end)=0; %remove vehicle from queue at depot
                             chargingqueue(Index,1)=t(i)/60 + 0; % charging start time instantaneous connection
                             chargingqueue(Index,2)=vehicle(2); % vehicle num
                             chargingqueue(Index,3)=vehicle(3); % vehicle distance travelled
                             chargingqueue(Index,4)=vehicle(4); % current charge level
                             chargingqueue(Index,5)=vehicle(5); % battery capacity level
                             chargingqueue(Index,6)=vehicle(6); % soc level
                             
                             
                             
                         else % if more than one vehicle needs to be charged
                             
                             if size(availablevehicles,1)<=Availablechargers
                                 for ii=1:size(availablevehicles,1)
                                     Availablechargers=Availablechargers-1; % connect vehicle to charger
                                     Usedchargers=Nchargers-Availablechargers;
                                     [~,Index]=min(chargingqueue(:,1)); % find available charger
                                     % all vehicles in the queue can be connected to a charger
                                     vehicle=availablevehicles(ii,:);
                                     queue(vehicle(2),:,end)=0; %remove vehicle from queue at depot
                                     chargingqueue(Index,1)=t(i)/60+ 0; % charging start time instantaneous connection
                                     chargingqueue(Index,2)=vehicle(2); % vehicle num
                                     chargingqueue(Index,3)=vehicle(3); % vehicle distance travelled
                                     chargingqueue(Index,4)=vehicle(4); % current charge level
                                     chargingqueue(Index,5)=vehicle(5); % battery capacity level
                                     chargingqueue(Index,6)=vehicle(6); % soc level
                                     
                                 end
                                 
                             else % if more vehicles need to be charged than available chargers
                                 % find vehicles withh the lowest soc
                                 for ii=1:Availablechargers % for each available charger connect vehicle with the lowest SOC
                                     Availablechargers=Availablechargers-1; % connect vehicle to charger
                                     Usedchargers=Nchargers-Availablechargers;
                                     [~,Index]=min(chargingqueue(:,1)); % find available charger
                                     % find vehicle with lowest SOC to connect to charger
                                     [~,Index1]=min(availablevehicles(:,4)); % xx find vehicle with min charge levels xx % send vehicle with the lowest SOC back to depot
                                     vehicle = availablevehicles(Index1,:); % select vehicle
                                     queue(vehicle(2),:,end)=0; %remove vehicle from queue at depot
                                     chargingqueue(Index,1)=t(i)/60+ 0; % charging start time instantaneous connection
                                     chargingqueue(Index,2)=vehicle(2); % vehicle num
                                     chargingqueue(Index,3)=vehicle(3); % vehicle distance travelled
                                     chargingqueue(Index,4)=vehicle(4); % current charge level
                                     chargingqueue(Index,5)=vehicle(5); % battery capacity level
                                     chargingqueue(Index,6)=vehicle(6); % soc level
                                     
                                 end
                             end
                         end
                     end
                     % update amount charged to the vehicles until the next time step
                     if i<length(t-1) & Usedchargers>=0
                         %  find time difference until next time step
                         t_step=t(i+1)-t(i);
                         energycharged=t_step*chargerpower/60; % kwh in 1 min
                         chargedemand(i)=energycharged*Usedchargers;
                         chargingqueue(chargingqueue(:,2)>0,4)=...
                             chargingqueue(chargingqueue(:,2)>0,4)+energycharged;
                         chargingqueue(:,6)=(chargingqueue(:,5)+chargingqueue(:,4))./...
                             chargingqueue(:,5); % update SOC
                         % stop charging vehicles that have completed
                         donecharging=chargingqueue(chargingqueue(:,6)>=1);
                         
                         if size(donecharging,1)>0
                             for ii=1:size(donecharging)
                                 [~,Index]=max(chargingqueue(:,6)); % find available charger used least
                                 vehicle = chargingqueue(Index,:); % select vehicle
                                 vehicle(1)=t(i+1)+5/60; % return to depot after disconnection
                                 vehicle(4)=0;
                                 vehicle(6)=1;
                                 % add back to depot queue
                                 queue(vehicle(1,2),:,end)=vehicle; % remove from queue at current stop
                                 chargingqueue(Index,:)=0;
                                 Availablechargers=Availablechargers+1;
                                 Usedchargers=Usedchargers-1;
                             end
                         end
                         
                     end
                     currentevents=[];
                     if flag==1
                         disp('No feasible solution')
                         break
                     end
                 end
                 toc
                 %Connect free modules that are in need of charge to available
                 %charger(s)
                 [x,y]=stairs(vehicledemand(:,1),cumsum(vehicledemand(:,2)));
                 vehicledemand=[[x(1) 0];x y];
                 
             end
        
             
        
        function property=propertyevaluation(obj)
            
            property.serviceperformance.waitingtime=mean([obj.lines.meanwaitingtime]);
            
            property.serviceperformance.dwellingtime= mean([obj.lines.meandwellingtime]); % not including stops with zero dwelling time
            
            property.serviceperformance.traveltime=mean([obj.lines.meaninvehicletime])+mean([obj.lines.meanwaitingtime;]);
            property.serviceperformance.seatavailability=mean(cell2mat({obj.lines.seatavailabilityinmin}'));
            property.serviceperformance.meanoccupancy=obj.occupancy;
            passengers=cell2mat({(obj.lines.totalpassengers)}');
            property.serviceperformance.missedboardings=sum(passengers(:,6))/size(passengers,1)*100;
            property.accessibility.wheelchairzones=obj.vehicle.Interior.wheelchairzones;
            property.accessibility.entryheight=150;%obj.vehicle.Body.groundclearance;
            property.accessibility.doorratio=obj.vehicle.Passengercapacity/(obj.vehicle.Interior.numberofdoors*1.2);
            property.accessibility.saturationflow=20;
            property.accessibility.decks=1;
            property.accessibility.bicycles=0;
            property.comfort.seatsize=(obj.vehicle.Interior.seatwidth/1000)*obj.vehicle.Interior.seatpitch/1000; %-correct
            property.comfort.legroom=obj.vehicle.Interior.seatpitch;
            property.comfort.standingarea=1/obj.vehicle.Interior.standingspace;
            property.comfort.temperature=24;
            property.comfort.windowsize=1;
            property.comfort.headroom=250;
            property.comfort.nvh=1;
            property.comfort.ridecomfort=1;
            property.functionality.information=0;
            property.functionality.payment=0;
            property.functionality.storagespace=0;
            property.functionality.usablearea=0;
            property.functionality.poweroutlet=0;
            property.luxury.seats=0;
            property.luxury.privacyseating=0;
            property.luxury.privatestanding=0;
            property.safety.rollstability=0;
            property.safety.cameras=0;
            property.safety.seatbelts=0;
            property.safety.handles=0;
            property.safety.illumination=0;
            property.longitudinaldynamics.speed=obj.vehicle.Properties.TopSpeed.Value*3.6;
            property.longitudinaldynamics.acceleration=obj.vehicle.Properties.Acceleration.Value;
            property.longitudinaldynamics.grade=obj.vehicle.Properties.Gradeability.Value;
            property.costs.TCO=obj.TCO.passengerkm;
            property.costs.acquisition=obj.TCO.year0/obj.TCO.total;
            property.costs.operatio=1-obj.TCO.year0/obj.TCO.total;
            property.environment.emissions=obj.gCO2_passengerkm;
            property.environment.usephaseemissions   =   obj.vehicle.Productionemissions.CO2/sum(obj.vehicle.LCAemissions.CO2);
            property.environment.productionemissions=1-property.environment.usephaseemissions;
            
        end
        
    end
    
    
    
end
