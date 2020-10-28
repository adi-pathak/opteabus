            uncapacitated=0;
            gurobi_solver=0;
            Econs=0.59; %0.59; %kwh/km
            pulloutcost=150000;
            chargerpower=0; % 300 kW
            batterycapacity=140; %120 kWh
            departl=[timetable(:,3) timetable(:,1)];
            arrivall=[timetable(:,4) timetable(:,2)];
            ui=timetable(:,end)*Econs; % Energy consumption of each trip
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
          if uncapacitated~=1
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
              
              if size(idx1,1)~=1
              i1=idx1(1:end-1);
              i2=size(E,1)+idy1(1:end-1);
              i3=ones(size(i1,1),1)*size(E,1)+i;
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
              
%               tic
%               for ii=1:size(idx1,1)
%                 
%                   
%                   if idy1(ii)>N(end)
%                       
%                   else
%                       n=n+1;
%                       Ce0{n,1}=[idx1(ii);size(E,1)+idy1(ii);size(E,1)+i];
%                       Ce0{n,2}=ones(size(Ce0{n,1},1),1)*n;
%                       Ce0{n,3}=[E0(idx1(ii),4)+Q;1;-1];
%                      
%                       
%                     
%                       Ce01{n,1}=[idx1(ii);size(E,1)+i];
%                       Ce01{n,2}=ones(size(Ce01{n,1},1),1)*n;
%                       Ce01{n,3}=[E0(idx1(ii),6);-1];
%                   end
% 
%              end
%                 toc
            end
            n=size(CE0,1);
            Ce0=CE0;
            Ce01=CE01;
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
               addpath('C:\gurobi901\win64\matlab')
               model.obj=cost;
               model.A=[A;Aeq]; % - not efficient/slow speed for large matrices
               model.sense=[repmat('<', size(A,1), 1);repmat('=', size(Aeq,1), 1)];
               model.rhs=[b;beq];
               model.modelsense='min';
               params.outputflag=1;
               model.lb    = lb;
               model.ub    = ub;
              
               r=gurobi(model,params);
               toc
               x=r.x;
            end
            
            
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
                %    blocks{12,n}= batterycapacity-cumsum(blocks{10,n});%;
                    n=n+1;
                    blocks{1,n}=0;
                    j=sortrows(z(:,1));
                    %  reshape(tb{7,1},2,6)'
                end
            end
            