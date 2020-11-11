classdef route
    properties
        id
        vehicleconcept
        routelength
        topology % single direction/bidirectional/loop
        busstops
        departures
        timetable
        meaninvehicletime
        interstoptime
        interstopspeed
        OD %origin-destination
        demand
        passengersonboard
        headwaypolicy
        desiredoccupancypolicy
        passengercapacity %persons/h/direction
        totalpassengers
        passengers
        averagevehiclespeed
        meanwaitingtime %in mins
        meandwellingtime %in secs
        missedboardings
        comfort
        accessibility
        stopcoordinates
        frequency
        headway
        passengerkm
        meanoccupancy %- passengers
        trajectories
        simulationdata
        volume
        capacity
        volumecapacityratio
        dwellingtime
        seatavailabilityinmin
    end
    
    methods
        function obj=route(service,serviceno,vehicle)
            obj.id=serviceno;
            obj.routelength=service.direction1.distance(end);
            obj.topology=1; % single direction/bidirectional/loop
            obj.busstops=size(service.stopcoords,1);
            obj.vehicleconcept=vehicle;
            obj.stopcoordinates=service.stopcoords;
      
            
            obj.demand=service;
   
            obj=generatetimetable(obj,vehicle.Passengercapacity);
            obj=passengerarrival(obj);
            obj=departure_simulation(obj);
            %obj=simulatefleet(obj);
            
        end
        function obj=passengerarrival(obj) %arrival of passengers at bus stop
            
            Tmax=60;
            if isfield(obj.demand,'direction_2')>0 % generate timetable for direction 2
                OD=obj.demand.direction_2.OD;
                boardingpassengers=cell(1,size(OD(:,:,:),1)); %create cell array boardingpassengers
                for i=5:size(OD,3) % for each hour from 5 am to 24 am/ last 
                    ODH=OD(:,:,i);
                    for j=1:size(ODH,1) %for each boarding stop
                        
                        lambda=sum(ODH(j,:)); %find average arrival rate and destination of passengers
                        if lambda~=0
                            passengerdestination=[ODH(j,:)' (1:size(ODH,2))']; % find destination stops and their frequency
                            passengerdestination=repelem(passengerdestination(:,2),passengerdestination(:,1)); % create matrix with frequency of stops
                            passengerarrival=poissonod(obj,lambda/60,lambda,Tmax,passengerdestination);
                            passengerarrival(:,1)=passengerarrival(:,1)+i;
                            boardingpassengers{j}=[boardingpassengers{j};passengerarrival];
                        end
                    end
                    
                end
                for i=1:size(ODH,1)
                    boardingpassengers{i}=sortrows(boardingpassengers{i});
                end
                obj.demand.direction_2.boardingpassengers=boardingpassengers;
                obj.demand.direction_2.inter_boardingpassengers=boardingpassengers;
            end
            
            if isfield(obj.demand,'direction1')>0 % generate timetable for direction 1
                boardingdemand=obj.demand.direction_1.OD;
                boardingpassengers=cell(1,size(boardingdemand,1));
                OD=boardingdemand;
                for i=5:size(OD,3)
                    ODH=OD(:,:,i);
                    for j=1:size(ODH,1)
                        
                        lambda=sum(ODH(j,:)); %find average arrival rate
                        if lambda~=0
                            passengerdestination=[ODH(j,:)' (1:size(ODH,2))']; % find destination stops and their frequency
                            passengerdestination=repelem(passengerdestination(:,2),passengerdestination(:,1)); % create matrix with frequency of stops
                            passengerarrival=poissonod(obj,lambda/60,lambda,Tmax,passengerdestination);
                            passengerarrival(:,1)=passengerarrival(:,1)+i;
                            x=boardingpassengers{j};
                            boardingpassengers{j}=[x;passengerarrival];
                        end
                    end
                end
                for i=1:size(ODH,1)
                    boardingpassengers{i}=sortrows(boardingpassengers{i});
                end
                obj.demand.direction_1.boardingpassengers=boardingpassengers;
                obj.demand.direction_1.inter_boardingpassengers=boardingpassengers;
            end
            
            
            
        end
        function obj=generatetimetable(obj,vehiclecapacity)
            
            if isfield(obj.demand,'direction_2')>0 % generate timetable for direction 2
                passengerload=obj.demand.direction_2.occupancy_passengers;
                stopdist=[obj.demand.direction_2.interstopdistance];
                try
                    frequency_policy=obj.demand.direction_2.policy*01;
                catch
                    frequency_policy=zeros(1,24);
                end
                
                firstdeparture=datevec(cell2mat(table2array(obj.demand.direction2(1,6))),'HHMM');
                firstdeparture=firstdeparture(4)+(firstdeparture(5)/60);
                lastdeparture=datevec(cell2mat(table2array(obj.demand.direction2(1,7))),'HHMM');
                lastdeparture=lastdeparture(4)+(lastdeparture(5)/60);
                [obj.departures(:,1),obj.headway(:,1)]=getdepartures(obj,passengerload,stopdist,frequency_policy,vehiclecapacity,firstdeparture,lastdeparture,0*size(obj.demand.direction_2.occupancy,1));
                obj.departures(:,2)=-1;
            end
            if isfield(obj.demand,'direction_1')>0 % generate timetable for direction 1
                passengerload=obj.demand.direction_1.occupancy_passengers;
                stopdist=[obj.demand.direction_1.interstopdistance];
                try
                    frequency_policy=obj.demand.direction_1.policy;
                catch
                    frequency_policy=zeros(1,24);
                end
                
                firstdeparture=datevec(cell2mat(table2array(obj.demand.direction1(1,6))),'HHMM');
                firstdeparture=firstdeparture(4)+(firstdeparture(5)/60);
                lastdeparture=datevec(cell2mat(table2array(obj.demand.direction1(1,7))),'HHMM');
                lastdeparture=lastdeparture(4)+(lastdeparture(5)/60);
                [departures(:,1),obj.headway(:,2)]=getdepartures(obj,passengerload,stopdist,frequency_policy,vehiclecapacity,firstdeparture,lastdeparture,0*size(obj.demand.direction_1.occupancy,1));
                departures(:,2)=1;
                obj.departures=[obj.departures; departures];
                obj.departures= sortrows(obj.departures);
                
            end
            
            
        end
        function [departures,headway]=getdepartures(obj,passengerload,stopdist,policy,vehiclecapacity,firstbus,lastbus,predefineddepartures)
            %% Rework of method 4 - ceder: point check vs ride check methods
            
            for i=1:size(passengerload,2) % find max load observed in each time period
                maxpass(i)=max(passengerload(:,i));
                % alighting.maxpass(i)=max(alighting.pass(:,i));
            end
            A=stopdist'*passengerload;
            L=obj.routelength; % route length
            Pmj=maxpass;
            % policy for hourly bus frequency
            Fmj=policy*1;
            %Fmj=0*ceil((1:24)/24);%ceil([0 0 0 0 5 60/9 60/9 60/9 60/9 60/10.5 ...
            %60/10.5 60/10.5 60/10.5 60/10.5 60/10.5 60/10.5 60/10.5 60/12.5 60/12.5 60/12.5 ...
            %60/13.5 60/13.5 0 0]);
            %ceil([0 0 0 0 0 0 60/5 60/5 60/5 60/5 ...
            %    60/5 60/5 60/5 60/5 60/5 60/5 60/5 60/5 60/5 60/5 ...
            %   60/14.5 60/14.5 60/14.5 60/14.5]);%%1*ceil((1:24)/24);%[1 2]; %[3 3 3 3 3];%
            c=vehiclecapacity;
            do=0.55*c*ceil((1:24)/24);%50*ceil((1:24)/24); %[0 0 0 0 0 60 70	70	60	50	50	50	50	50	60	70	70	60	60	50	50	50	50	50];
            %ceil(c*0.30)*ceil((1:24)/24);%[50 60]; %[ceil(c*0.55) 50 50 50 50]; desired occupancy
            b1=0.1; % allowable portion of the route length to exceed desired occupancy
            
            %% frequency determination method 4 - Ceder
            F_1=Pmj./do; % method 1
            F_2=max([Pmj./do;Fmj]); % method 2
            %F=max([((A./(do*L)).*((A./(do*L))>Fmj));(((A./(do*L))<Fmj).*(Pmj/c)) ;((((A/(do*L))<Fmj).*((Pmj/c)<Fmj)).*Fmj)]); % ceder pg 77
            F_3=max(([(A./(do*L)); ((Pmj/c)); Fmj])); %method 3
            F=F_2;
          
            %% timetable generation
            if predefineddepartures>0
                s=(predefineddepartures-1)/sum(F);
                F=F*s; %scale frequency based on predetermined departures
            end
            if lastbus>1 
            lasthour=floor(lastbus*(lastbus>1)+(lastbus+24)*(lastbus<1));
            else
                lasthour=24;
                lastbus=(lastbus*(lastbus>1)+(lastbus+24)*(lastbus<1));
            end
              
            n=1; % cumulative vehicle frequency
            x=0;
            c=1; 
            dt=0;
            if ~isempty(firstbus)
            ii= floor(firstbus);
            Depart(1)=firstbus;
            n=n+1;
            
             while x<1
                    
                    if F(ii)==0 break % if frequency of the ith hour is not zero , continue
                    end
                    m=(1+F(ii)-1)/(ii+1-firstbus); %change of slope due to non zero start
                    x=(n-c)/m; % x is change in time for next departure, n=Fx+c
                    x=firstbus-ii+x;
                    if x>1 break
                    end
                    dt(n)=x;
                    Depart(n)=ii+dt(n);
                    n=n+1;
             end
             c=c+F(ii);
             ii=ii+1;
             x=0;
              
             for i=ii:lasthour-1
                
                 while x<1
                    
                     if F(i)==0 break % if frequency of the ith hour is not zero , continue
                     end
                     x=(n-c)/F(i);
                     if x>1 break
                     end
                     dt(n)=x;
                     Depart(n)=(i)+dt(n);
                     n=n+1;
                 end
                 x=0;
                 c=c+F(i); % slope
             end
            ii=lasthour;
            while x<1
            m=(c+F(ii)-c)/(lastbus-ii);
            x=(n-c)/m;
            x= x*(x<lastbus-ii)+ (1-(x<=lastbus-ii));
            if x<1 
                dt(n)=x;
                Depart(n)=(ii)+dt(n);
                n=n+1;
                if x==0 break
                end
            end
            end     
             
            else
                ii=1;
                
            for i=ii:lasthour-1 %size(passengerload,2)
               
                while x<1
                    
                    if F(i)==0 break % if frequency of the ith hour is not zero , continue
                    end
                    x=(n-c)/F(i);
                    
                    if x>1 break
                    end
                    dt(n)=x;
                    Depart(n)=(i)+dt(n);
                    n=n+1;
                end
                x=0;
                c=c+F(i); % slope
                
            end
            
            ii=lasthour;
            while x<1
            m=(c+F(ii)-c)/(lastbus-ii);
            x=(n-c)/m;
            x= x*(x<lastbus-ii)+ (1-(x<=lastbus-ii));
            if x<1 
                dt(n)=x;
                Depart(n)=(ii)+dt(n);
                n=n+1;
                if x==0 break
                end
            end
            end   
            end
            Depart(end)= lastbus;
            departures=Depart;
            frequency=F;
            headway=round(60./F);
        end
        function passengerarrival = poissonod(obj,lambda,eventnum,Tmax,passengerdestination)
            % lambda=280;      % arrival rate in hr
            % Tmax=1;         % maximum time in hr
            
            w(1) = 0.0;
            w(2:eventnum+1) = - log ( rand ( eventnum, 1 ) ) / lambda;
            %
            %  The time til event I is the sum of the waiting times 0 through I.
            
            t(1:eventnum+1) = cumsum ( w(1:eventnum+1) );
            passengerarrival(:,1)=t(2:end)./60;
            passengerarrival(:,2)=datasample(passengerdestination,eventnum,'Replace',false);
            passengerarrival(:,3)=0; % for missed boarding count
            return
            
        end
        function [boardingpassengers,arrivaltime,occupancy,passengers,arrivaltime_stops,dwellt,boarding,alighting]=pax_exchange(obj,demand,departure)
            S=0; % S=0 if low floor, S=1 otherwise
            D=0; % D=0 if single deck D=1 otherwise
            vehiclecapacity=obj.vehicleconcept.Passengercapacity;
            seats=obj.vehicleconcept.Interior.numberseats;
            interstop_distance=[demand.interstopdistance*1000; 0];
            stopdistance=cumsum(interstop_distance);
            interstop_time=demand.interstoptime;
            meaninterstop_time=demand.meaninterstoptime;
            arrivaltime=departure; % time of first boarding/alighting
            boardingpassengers=demand.inter_boardingpassengers;
            %missedboardings=[];
            %passengercount=zeros(1,size(interstop_distance,1)); %passenger boarding count
            %pbdt=zeros(size(Depart,2),demand.stops(end)); %passenger boarding time
            occupancy=zeros(1,size(interstop_distance,1)); %occupancy of each vehicle
            dwellt=zeros(1,size(interstop_distance,1)); %dwelling time matrix
           % at=zeros(1,size(interstop_distance,1)); %travel time
          %  alightingpassengers=zeros(6,size(interstop_distance,1)); %passengers alighting at each stop
            boarding_passengers=cell(1,size(interstop_distance,1));
            alighting_passengers=cell(1,size(interstop_distance,1));
            alightingpassengers=[0 0 0 0 0 0 0];
            %bp=[];
            %ap=[];
            z=[];
            distancecovered=0;
            seatavailability=seats;
            passengers=[];
            if size(boardingpassengers,2)~=size(boarding_passengers,2)
             boardingpassengers{1,size(boarding_passengers,2)}=[];
            end
            for k1=1:size(boarding_passengers,2)-1 %at each stop
              
                % check if passenger is alighting & remove passengers from
                % vehicle..
                if ~isempty(alightingpassengers(alightingpassengers(:,3)==k1)) % find alighting passengers
                    occupancy(k1)=occupancy(k1)-size((alightingpassengers(alightingpassengers(:,3)==k1)),1); % update vehicle occupancy/remove alighting pax
                   % size((alightingpassengers(alightingpassengers(:,3)==k1)),1)
                    if seatavailability==0 % if no seats available
                        if size(alightingpassengers((alightingpassengers(:,3)==k1 & alightingpassengers(:,6)~=0)),1)>1
                            idx=find(alightingpassengers(:,6)==0); % find standing pax wanting to sit
                            if size(idx,1)~=1 & ~isempty(idx)
                                [~,b]=sortrows(alightingpassengers(idx,3)-alightingpassengers(idx,4));
                                idx=idx(b);
                                a= size((alightingpassengers(alightingpassengers(:,3)==k1)),1);
                                if a>size(idx,1)
                                    
                                    idx=idx(2:end);
                                    
                                else
                                    idx=idx(end-a+1:end);
                                end
                                if idx(1)==1
                                    idx=idx(2:end);
                                end
                                alightingpassengers(idx,6)=1-((distancecovered(k1)-alightingpassengers(idx,7))./stopdistance(alightingpassengers(idx,3)));
                                seatavailability=seatavailability+a-size(idx,1);
                            else
                                seatavailability=seatavailability+size((alightingpassengers(alightingpassengers(:,3)==k1)),1);
                            end
                        elseif size(alightingpassengers((alightingpassengers(:,3)==k1 & alightingpassengers(:,6)~=0)),1)==1
                            idx=find(alightingpassengers(:,6)==0);
                            if idx==1
                                seatavailability=seatavailability+ size((alightingpassengers(alightingpassengers(:,3)==k1)),1);
                            else
                            [~,b]=max(alightingpassengers(idx,3)-alightingpassengers(idx,4));
                            alightingpassengers(idx(b),6)=1-(distancecovered(k1)-alightingpassengers(idx(b),7))/sum([interstop_distance(1:alightingpassengers(idx(b),3))]);
                            seatavailability=0;
                            end
                        else % all pax have seats
                            
                            seatavailability=0;
                            
                        end
                        
                    else
                    seatavailability=seatavailability+size((alightingpassengers(alightingpassengers(:,3)==k1)),1);%size(alightingpassengers((alightingpassengers(:,3)==k1 & alightingpassengers(:,6)==1)),1);
                        
                    end
                    
                    z(2,k1)=  size((alightingpassengers(alightingpassengers(:,3)==k1)),1);
                    alighting_passengers{k1}=[(ones(size((alightingpassengers(alightingpassengers(:,3)==k1)),1),1)*arrivaltime(k1)) (alightingpassengers((alightingpassengers(:,3)==k1),:))  distancecovered(k1)*ones(size((alightingpassengers(alightingpassengers(:,3)==k1)),1),1) ]; %vehicle arrival time passenger arrival time stop ];
                    passengers=[passengers; alighting_passengers{k1}]; 
                   % ap=[ap; alightingpassengers((alightingpassengers(:,3)==k1),:)];
                    alightingpassengers((alightingpassengers(:,3)==k1),:)=[];
                  
                end
                
                %%
                if ~isempty(boardingpassengers{1, k1}) % if there are passengers boarding
                    % check if there are waiting passengers when vehicle arrives
                    passengercount=boardingpassengers{1, k1}((((boardingpassengers{1, k1}(:,1))<arrivaltime(k1))),:);      %number of people at each stop before departure
                    
                    if sum(passengercount(:,1))>0%if passengers waiting at stop start boarding
                        if (occupancy(k1)+size(passengercount,1))<=vehiclecapacity % if number of passenger less than vehicle capacity
                            boarding_passengers(k1)={passengercount};
                            %bp=[bp;passengercount];
                            %passengers=[passengers; [passengercount k1*ones(size(passengercount,1),1) arrivaltime(end)*ones(size(passengercount,1),1) arrivaltime(end)-passengercount(:,1)]];
                            z(1,k1)=size(passengercount,1);
                            occupancy(k1)=occupancy(k1)+size(passengercount,1);
                            % find the alighting stops of the boarding passengers
                            
                            if seatavailability-size(passengercount,1)>0
                               seatingpassengers=passengercount(:,1)*0+1;
                               seatavailability=seatavailability-size(passengercount,1);
                            else
                                seatingpassengers=passengercount(:,1)*0;
                                seatingpassengers(1:seatavailability,:)=1;
                                seatavailability=0;
                            end
                            alightingpassengers=[alightingpassengers;[arrivaltime(k1).*ones(size(passengercount,1),1) passengercount(:,1:2) k1*ones(size(passengercount,1),1) passengercount(:,3) seatingpassengers distancecovered(k1)*ones(size(passengercount,1),1)] ]; %vehicle arrival time passenger arrival time stop
                            % remove the boarded passengers from the waiting
                            % passengers
                            boardingpassengers{1, k1}((((boardingpassengers{1, k1}(:,1))<arrivaltime(k1) ) ),:)=[];
                            
                        else % if waiting passengers exceed capacity, find the number of pax that can be allowed to board
                            ab=vehiclecapacity-(occupancy(k1)); % max passengers allowed to board
                            ab=(ab*(ab<=size(passengercount,1))+(size(passengercount,1))*(ab>size(passengercount,1))); %ensures correct pax board the stop
                            occupancy(k1)=ab+occupancy(k1); % add pax to the vehicle
                            ff=find(passengercount(:,1)); % find index of all available boarding passengers
                            boardingpassengers{1, k1}(ab+1:size(passengercount,1),3)= cell2mat({boardingpassengers{1, k1}(ab+1:size(passengercount,1),3)})+1; % indicate if passenger missed boarding
                         %   missedboardings=[missedboardings;boardingpassengers{1, k1}(ab+1:size(passengercount,1),:)];
                            ff = ff(1:ab); % limit the number of pax to board
                            boarding_passengers(k1)={boardingpassengers{1, k1}(ff,:)};
                           % bp=[bp;boardingpassengers{1, k1}(ff,:)];
                            %passengers=[passengers; [passengercount k1*ones(size(passengercount,1),1) arrivaltime(end)-passengercount(:,1)]];
                            z(1,k1)=size(boardingpassengers{1, k1}(ff,:),1);
                            seatingpassengers=(boardingpassengers{1, k1}(ff,1))*0;
                            alightingpassengers=[alightingpassengers;[arrivaltime(k1).*ones(size(boardingpassengers{1, k1}(ff,:),1),1) boardingpassengers{1, k1}(ff,1:2)  k1*ones(size(boardingpassengers{1, k1}(ff,:),1),1) boardingpassengers{1, k1}(ff,3) seatingpassengers distancecovered(k1)*ones(size(boardingpassengers{1, k1}(ff,:),1),1)]];
                            boardingpassengers{1, k1}(ff,:)=[]; % remove the passengers boarded from bus stop
                            
                        end
                    end
                    
                    
                    
                end
                
                
                % driving time to next stop in min
                try
                    if ~isnan(interstop_time(k1,k1+1, floor(arrivaltime(k1)))) % time to next stop in min
                        time_interstop=(interstop_time(k1,k1+1, floor(arrivaltime(k1))))/60;
                        
                    elseif ~isnan(meaninterstop_time(k1,k1+1))
                        time_interstop=(meaninterstop_time(k1,k1+1))/60;
                    else
                        time_interstop=interstop_distance(k1)/(15*1000); %assuming 15kmph average speed
                    end
                catch
                    time_interstop=interstop_distance(k1)/(15*1000); %assuming 15kmph average speed
                    
                end
                
                % boarding time
                B=(size(boarding_passengers{k1},1)-1);
                boarding(k1)=size(boarding_passengers{k1},1);
                if B<0
                    B=0;
                end
                dwellt(k1)=B*(1.951+(-0.017*(B))+(0.047*D)+(0.156*S)+(0.340*(occupancy(k1)/vehiclecapacity)));
                
                % alighting time
                A=(size(alighting_passengers{k1},1)-1); % alighting number
                alighting(k1)=(size(alighting_passengers{k1},1));
                if A<0
                    A=0;
                end
                dwellt(k1)=dwellt(k1)+(A*(1.691+(-0.014*A)+(0.217*D)+(0.016*S)+(-0.082*(occupancy(k1)/vehiclecapacity))));
                arrivaltime(k1+1)=arrivaltime(k1)+(dwellt(k1)/3600)+(time_interstop); %arrival time to next stop
                distancecovered(k1+1)=distancecovered(k1)+interstop_distance(k1);
                seating(k1)=seatavailability;
                arrivaltime(k1+1)=arrivaltime(k1)+(dwellt(k1)/3600)+(time_interstop); %arrival time to next stop
                timeinterstop(k1)=time_interstop*60;
                if k1<size(occupancy,2)
                    occupancy(k1+1)=occupancy(k1);
                end
                
            end
            k1=k1+1;
            
            if ~isempty(alightingpassengers(alightingpassengers==k1))
                occupancy(k1)=occupancy(k1-1)-size((alightingpassengers(alightingpassengers==k1)),1);
                z(2,k1)=  size((alightingpassengers(alightingpassengers==k1)),1);
                alighting_passengers{k1}=[(ones(size((alightingpassengers(alightingpassengers==k1)),1),1)*arrivaltime(k1)) (alightingpassengers((alightingpassengers(:,3)==k1),:)) distancecovered(k1)*ones(size((alightingpassengers(alightingpassengers(:,3)==k1)),1),1)];
                %ap=[ap; alightingpassengers((alightingpassengers(:,3)==k1),:)];
                passengers=[passengers; alighting_passengers{k1}];
                alightingpassengers((alightingpassengers(:,3)==k1),:)=[];
                boarding(k1)=size(boarding_passengers{k1},1);
                alighting(k1)=(size(alighting_passengers{k1},1));
                %passengers=[passengers; alighting_passengers{k1}];
                %                           dwellt(k1)=dwellt(k1)+(size(alighting_passengers{k1},1))*(1.691+(-0.014*((size(alighting_passengers{k1},1))-1))+(0.217*0)+(-0.082*(occupancy(k1)/vehiclecapacity)));
            else
                occupancy(k1)=0; 
                boarding(k1)=size(boarding_passengers{k1},1);
                alighting(k1)=(size(alighting_passengers{k1},1));
            end
            z=cumsum(diff(z,1)*-1);
            arrivaltime_stops=arrivaltime;
            arrivaltime=arrivaltime(end);
        end
        function obj=departure_simulation(obj)
            occupancy_1=[];
            occupancy_2=[];
            occupancy__1=[];
            occupancy__2=[];
            passengers_2=[];
            passengers_1=[];
            timetable=[];
            stoparrivaltime_1=[];
            stoparrivaltime_2=[];
            dwellingtime_1=[];
            dwellingtime_2=[];
            volume_1=zeros(24,size(obj.demand.direction_1.distance,1));
            capacity_1=zeros(24,size(obj.demand.direction_1.distance,1));
            volume_2=zeros(24,size(obj.demand.direction_2.distance,1));
            capacity_2=zeros(24,size(obj.demand.direction_2.distance,1));
            passengercapacity=obj.vehicleconcept.Passengercapacity;
            stops_1=size(obj.demand.direction1.distance,1);
            stops_2=size(obj.demand.direction2.distance,1);
            stops_coords1=obj.stopcoordinates(1:stops_1,:);
            stops_coords2=obj.stopcoordinates(stops_1+1:end,:);
            num_departures_1=0;
            num_departures_2=0;
            numseats=obj.vehicleconcept.Interior.numberseats;
            
            for i=1:length(obj.departures)  % for each departure
                
                if obj.departures(i,2)==1; %
                    [obj.demand.direction_1.inter_boardingpassengers,arrivaltime,occupancy1,passengers1,arrivaltime_stops,dwellt,boarding1,alighting1]=pax_exchange(obj,obj.demand.direction_1,obj.departures(i)); %find arrival time to terminal & passenger occupancy for the trip
                    % arrivaltime(n)=sum(obj.interstoptime)/60+obj.departures(i); %calculate time of arrival at terminal
                    num_departures_1=num_departures_1+1;
                    arrivaltimet(i)=arrivaltime;
                    passengers_1=[passengers_1;passengers1];
                    occupancy_1=[occupancy_1;occupancy1];
                    % stoparrivaltime_1=[stoparrivaltime_1;arrivaltime_stops];
                    passengerload(occupancy1<=numseats)=1; % seat availability
                    passengerload(occupancy1>numseats & occupancy1<passengercapacity)=0.5; %standing availability
                    passengerload(occupancy1==passengercapacity)=0; % bus full
                    
                    
                    route_data_1=[arrivaltime_stops' [1:stops_1]' dwellt' boarding1' alighting1' occupancy1' passengerload' stops_coords1 num_departures_1.*ones(stops_1,1)]; %arrival time at stop, stop number, dwellingtime, boarding num, alighting num, occupancy, loading, coordinates
                    routedata_1{:,num_departures_1}=route_data_1;
                    trajectory_time_1(:,num_departures_1)=arrivaltime_stops';
                    trajectory_stops_1(:,num_departures_1)=[1:stops_1]';
                    dwellingtime_1=[dwellingtime_1;dwellt];
                    timetable=[timetable;[obj.departures(i,1) arrivaltime obj.departures(i,2)] ];
                    arrivaltime_stops(arrivaltime_stops>24)=arrivaltime_stops(arrivaltime_stops>24)-24;
                    arrivaltime_stops=floor(arrivaltime_stops);
                    arrivaltime_stops(arrivaltime_stops==0)=24;
                    volume_1(unique((arrivaltime_stops)),:)=volume_1(unique((arrivaltime_stops)),:)+((arrivaltime_stops)==[unique((arrivaltime_stops))]').*occupancy1;
                    capacity_1(unique((arrivaltime_stops)),:)=capacity_1(unique((arrivaltime_stops)),:)+((arrivaltime_stops)==[unique((arrivaltime_stops))]').*passengercapacity;
                    passengerload=[];
                else
                    [obj.demand.direction_2.inter_boardingpassengers,arrivaltime,occupancy2,passengers2,arrivaltime_stops,dwellt,boarding2,alighting2]=pax_exchange(obj,obj.demand.direction_2,obj.departures(i));
                    %arrivaltime(n)=sum(obj.interstoptime)/60+obj.departures(i); %calculate time of arrival
                    num_departures_2=num_departures_2+1;
                    arrivaltimed(i)=arrivaltime;
                    occupancy_2=[occupancy_2;occupancy2];
                    passengers_2=[passengers_2;passengers2];
                    %stoparrivaltime_2=[stoparrivaltime_2;arrivaltime_stops];
                    
                    passengerload(occupancy2<=numseats)=1; % seat availability
                    passengerload(occupancy2>numseats & occupancy2<passengercapacity)=0.5; %standing availability
                    passengerload(occupancy2==passengercapacity)=0; % bus full
                    
                    route_data_2=[arrivaltime_stops' [1:stops_2]' dwellt' boarding2' alighting2' occupancy2' passengerload' stops_coords2  num_departures_2.*ones(stops_2,1)];
                    routedata_2{:,num_departures_2}=route_data_2;
                    trajectory_time_2(:,num_departures_2)=arrivaltime_stops';
                    trajectory_stops_2(:,num_departures_2)=[1:stops_2]';
                    dwellingtime_2=[dwellingtime_2;dwellt];
                    timetable=[timetable;[obj.departures(i,1) arrivaltime obj.departures(i,2)] ];
                    arrivaltime_stops(arrivaltime_stops>24)=arrivaltime_stops(arrivaltime_stops>24)-24;
                    arrivaltime_stops=floor(arrivaltime_stops);
                    arrivaltime_stops(arrivaltime_stops==0)=24;
                    volume_2(unique((arrivaltime_stops)),:)=volume_2(unique((arrivaltime_stops)),:)+((arrivaltime_stops)==[unique((arrivaltime_stops))]').*occupancy2;
                    capacity_2(unique((arrivaltime_stops)),:)=capacity_2(unique((arrivaltime_stops)),:)+((arrivaltime_stops)==[unique((arrivaltime_stops))]').*passengercapacity;
                    passengerload=[];
                end
            end
            obj.timetable=timetable;
            obj.missedboardings= sum(passengers_1(:,6))+sum(passengers_2(:,6));
            obj.meaninvehicletime=mean([mean((passengers_1(:,1)-passengers_1(:,2))*60);mean((passengers_2(:,1)-passengers_2(:,2))*60)]);
            obj.meanwaitingtime=mean([mean((passengers_1(:,2)-passengers_1(:,3))*60);mean((passengers_2(:,2)-passengers_2(:,3))*60)]);
            obj.passengerkm=sum([(passengers_1(:,9)-passengers_1(:,8))/1000;(passengers_2(:,9)-passengers_2(:,8))/1000]);
            obj.passengersonboard.direction1=occupancy_1;
            obj.passengersonboard.direction2=occupancy_2;
            obj.meanoccupancy=mean([mean(mean(occupancy_2));mean(mean(occupancy_1))]);
            obj.passengers.d1=passengers_1;
            obj.passengers.d2=passengers_2;
            obj.totalpassengers=[passengers_1 ;passengers_2];
            obj.seatavailabilityinmin=((obj.totalpassengers(:,7).*(obj.totalpassengers(:,1)-obj.totalpassengers(:,2)))*60);
            obj.trajectories.d1.time=trajectory_time_1;
            obj.trajectories.d1.space=trajectory_stops_1;
            obj.trajectories.d2.time=trajectory_time_2;
            obj.trajectories.d2.space=trajectory_stops_2;
            obj.simulationdata.d1=routedata_1;
            obj.simulationdata.d2=routedata_2;
            obj.dwellingtime.d1=dwellingtime_1;
            obj.dwellingtime.d2=dwellingtime_2;
            if ~isempty(dwellingtime_2)
            obj.meandwellingtime=mean([mean(dwellingtime_1(dwellingtime_1~=0)) mean(dwellingtime_2(dwellingtime_2~=0))]);
            else
              obj.meandwellingtime=mean(dwellingtime_1(dwellingtime_1~=0));  
            end
            obj.volume.d1=volume_1;
            obj.volume.d2=volume_2;
            obj.capacity.d1=capacity_1;
            obj.capacity.d2=capacity_2;
            obj.volumecapacityratio.d1=round(volume_1./capacity_1*100);
            obj.volumecapacityratio.d1(isnan(obj.volumecapacityratio.d1))=0;
            obj.volumecapacityratio.d2=round(volume_2./capacity_2*100);
            obj.volumecapacityratio.d2(isnan(obj.volumecapacityratio.d2))=0;
        end
        
    end
    
end