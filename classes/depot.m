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
             [vehicle.Energyconsumption,vehicle.Properties.range, ...
                vehicle.Properties.vmax,vehicle.Properties.amax,...
                vehicle.Properties.gradeability,vehicle.Energythroughput]=EnergyConsumption(vehicle,driving_cycle);
            
            %% calculate fleet size and find optimal vehicle duties/blocks
            if ~isempty(dialogbar)
                dialogbar=waitbar(0.5,dialogbar,'Finding Optimal Vehicle Assignment');
            end
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
