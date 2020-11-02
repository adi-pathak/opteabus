classdef AEVtoolkit
    %%
    % This class contains all the functions used to plot and evaluate
    % vehicle concepts in the GUI - AEV.mlapp
    %%
    properties
        GUI;
        slope=[-0.222 -0.0667 -0.1333 0.1333 -0.0889 -2 0 0; % Service Performance
            4.0000   -0.0444   -0.1333    1.5556   -4.0000    1.6444         0         0; % Accessibility
            0.5000    0.0444   -0.6222    0.4444    0.5000    0.5000    0.0889    0.5000; % Comfort
            0.5 0.5 0.5 0.5 0.5 0 0 0; % Functionality
            0.5 0.5 0.5 0 0 0 0 0; % Luxury
            0.5 4 4 0.5 0.5 0 0 0; % Safety
            0.1333    3.6444 0.5 0 0 0 0 0; % Lngitudinal Dynamics
            0.5 0.5 0.5 0 0 0 0 0; % Costs
            0.5 0.5 0.5 0 0 0 0 0; % Environment
            ];
        inflexion=[15 51.75 28.65 22.75 40 0.1 0 0; % Service Performance
            0.44 110.35 38.75 2 2.0 0.05 0 0; % Accessibility
            0.5000 930.0000    7.6782   19.4828    0.5000    0.5000  449.6552    0.5000; % Comfort
            0.5 0.5 0.5 0.5 0.5 0 0 0; % Functionality
            0.5 0.5 0.5 0 0 0 0 0; % Luxury
            0.5 0.5 0.3793 0.5 0.5 0 0 0; % Safety
            65.2874 1.0000 0.5000 0  0         0         0         0; % Lngitudinal Dynamics
            0.5 0.5 0.5 0 0 0 0 0; % Costs
            0.5 0.5 0.5 0 0 0 0 0; % Environment
            ];
        xl=[0 0 0 0 0 0 0 0; % Service Performance
            0 0 10 0.5 1 0 0 0; % Accessibility
            0 660 2 0 0 0 240 0; % Comfort
            0 0 0 0 0 0 0 0; % Functionality
            0 0 0 0 0 0 0 0; % Luxury
            0 0 0 0 0 0 0 0; % Safety
            40 0 0 0 0 0 0 0; % Lngitudinal Dynamics
            0 0 0 0 0 0 0 0; % Costs
            0 0 0 0 0 0 0 0; % Environment
            ];
        xu=[60 120 60 60 100 5 0 0; % Service Performance
            2 300 60 4 2 2 0 0; % Accessibility
            0.5 1200 12 25 0.5 0.5 1000 0.5; % Comfort
            0.5 0.5 0.5 0.5 0.5 0 0 0; % Functionality
            0.5 0.5 0.5 0 0 0 0 0; % Luxury
            0.5 1 1 0.5 0.5 0 0 0; % Safety
            120 2 0.5 0 0 0 0 0; % Lngitudinal Dynamics
            0.5 1 1 0 0 0 0 0; % Costs
            0.5 1 1 0 0 0 0 0; % Environment
            ];
        xname={'Waiting Time in min' 'Mean Dwelling Time in sec' ...
            'Travel Time in min' 'Mean Seat Availability in min' ...
            'Mean Occupancy in %' 'Missed Boardings' '' ''; % Service Performance
            'Wheelchair Zones' 'Entry Height in mm' 'Rated Capacity to Door Ratio'...
            'Saturation Flow in Pass/s' 'Number of Decks' 'Bicycle Zones' '' ''; % Accessibility
            'Seat Size in m2' 'Legroom in mm' 'Standing Area in P/m2' ...
            'Temperature in deg C' '' 'Window Area in m2' 'Headroom in mm' 'NVH'; % Comfort
            '' 0.5 0.5 0.5 0.5 0 0 0; % Functionality
            0.5 0.5 0.5 0 0 0 0 0; % Luxury
            0.5 0.5 0.5 0.5 0.5 0 0 0; % Safety
            'Top Speed in km/h' 'Acceleration' 0.5 0 0 0 0 0; % Longitudinal Dynamics
            'Total Cost of Ownership in $/passenger-km' ...
            0.5 0.5 0 0 0 0 0; % Costs
            'GHG Emissions in gCO2/passenger-km' 0.5 0.5 0 0 0 0 0; % Environment
            };
        catergorical_weights={[0.18204 0.17122 0.17122 0.17122 0.13308 0.17122];... % Service Performance - Waitingtime, dwelltime,traveltime,seatavail,crowdedness,missedboards
            [0.20031 0.15728 0.17293 0.17293 0.14397 0.15258 ]; ... % Wheelchairzones,lowentry,ndoors,doorwidth,numdecks,bicycles
              [0.111619/2 0.111619/2 0.109344/2 0.15548 0.29309 0.087995 0.109344/2 0.119474 ];... %seatcomf,legroom,standcomf,climcomf,ridecomf,visib,headroom,nvh
             [0.31267 0.35167 0.335655/3 0.335655/3 0.335655/3];... %information,payment,storagespace,usablearea,poweroutlets (luggage,stoller,usbpower)
             [0.3381 .3327 .32916];... %luxury,privacyseats,privatestandingspace
             [0.27716 0.26703 0.1915285/2 0.2771639 0.1915285/2];... %rollsafety,cameras,seatbelts.supporthandles,illumination
             [0.33 0.33 0.33];... %speed,acceleration,gradability
             [0.33 0.33 0.33];... 
             [0.33 0.33 0.33];}
         overallweights=[0.34662359;0.04511067;0.16792956;0.09271263;...
             0.05938627;0.11912941;0.16910785;0;0];
    end
    
    methods
        
        function plotProperties(obj,handle,numberofproperties,targets,LevelNum,PropertyNames)
            n=numberofproperties; % number of properties
            %%m=comparators; % number of simultaneous comparisons
            s=5;
            R=targets-s;
            R=[R R(:,1)];
            R(R<0)=0;
            % Plot property targets
            [Theta,M]=meshgrid(2*pi/n*[0:n]+pi/n,ones(1,size(R,1)));
            X=R.*sin(Theta);
            Y=R.*cos(Theta);
            cla(handle)
            
            
            fillPlot=1;
            if fillPlot==1
                color=lines;
                %color=color(randi([1 256],1),:);
                for i=1:size(X,1)
                    h1(i)= fill(handle,X(i,:)',Y(i,:)',color(i,:),'FaceAlpha',0.5);
                    hold(handle,'on')
                    h(i)=plot(handle,X(i,:)',Y(i,:)','LineWidth',2,'Color',color(i,:)); % plot lines
                    
                end
            else
                h=plot(handle,X',Y','LineWidth',2); % plot lines
            end
            
            % draw axis
            MAXAXIS=10-s;
            axis(handle,[-MAXAXIS MAXAXIS -MAXAXIS MAXAXIS]);
            axis(handle,'equal')
            AxisR=linspace(0,MAXAXIS,LevelNum);
            
            % Axis Labels
            text(handle,AxisR(1)*sin(pi/n)*1.15,AxisR(1)*cos(pi/n)*1.15,num2str(AxisR(1),2),'FontSize',12)
            %text(AxisR(end)*sin(pi/n)*1.0,AxisR(end)*cos(pi/n)*1.0,num2str(AxisR(end),2),'FontSize',12,'FontWeight','bold')
            
            
            
            [M,AxisR]=meshgrid(ones(1,n),AxisR);
            AxisR=[AxisR AxisR(:,1)];
            [AxisTheta,M]=meshgrid(2*pi/n*[0:n]+pi/n,ones(1,size(AxisR,1)));
            AxisX=AxisR.*sin(AxisTheta);
            AxisY=AxisR.*cos(AxisTheta);
            hold(handle,'on')
            plot(handle,AxisX,AxisY,'k')%,':k')
            plot(handle,AxisX',AxisY','k')%,':k')
            hold(handle,'on')
            
            
            % plot names of properties
            LableTheta=2*pi/n*[0:n-1]+pi/n;
            LableR=MAXAXIS;
            LableX=LableR.*sin(LableTheta);
            LableY=LableR.*cos(LableTheta);
            Lable=PropertyNames;
            
            for i=1:n
                if ~sum(strcmpi({'' },Lable(i)))
                    text(handle,LableX(i)*1.35, LableY(i)*1.35,sprintf(string(Lable(i))), 'FontSize',14,'HorizontalAlignment','center','Rotation',0)
                    %text(handle,LableX(i)*1.3, LableY(i)*1.3,sprintf(string(Lable(i))), 'FontSize',12,'FontWeight','bold','HorizontalAlignment','center','Rotation',0)
                    
                end
            end
            % Add axis Labels
            %             for i=1:size(AxisR,2)-1
            %                 text(handle,LableX(i)*1.05, LableY(i)*1.05,num2str(AxisR(end),2),'FontSize',12,'HorizontalAlignment','center','Rotation',0)
            %             end
            xlim(handle,[-MAXAXIS-3 MAXAXIS+3]) % change axis limits for better visibility
            ylim(handle,[-MAXAXIS-3 MAXAXIS+3])
            
        end
        function plotValueFunction(obj,handle,slope,inflexion,xl,xu,xlabelname)
            g0=4;
            gu=10;
            S=slope;
            W=inflexion;
            x=linspace(xl,xu,100);
            y=((g0-gu)./(1+(exp(S*(x-W)))))+gu;
            plot(handle,x,y)
            xlabel(handle,xlabelname)
            ylabel(handle,'Objective Value')
        end
        function plotTCO(ax,Results)
            hold(ax,'on')
            figure(1)
            
            barData=[];
            xTicker=[];
            if size(Results,1)==1
                barData=[Results.TCO.EOL/1000000 Results.TCO.year0/1000000 ...
                    Results.TCO.cleaning/1000000 Results.TCO.maintenance/1000000 ...
                    Results.TCO.energy/1000000 Results.TCO.personnel/1000000 ...
                    Results.TCO.batteryreplacement/1000000 Results.TCO.insurance/1000000 Results.TCO.roadtax/1000000; nan(1,9)];
                xTicker={'Concept 1'};
            else
                for i=1:size(Results,1)
                    barData=[barData;Results(i).TCO.EOL/1000000 Results(i).TCO.year0/1000000 ...
                        Results(i).TCO.cleaning/1000000 Results(i).TCO.maintenance/1000000 ...
                        Results(i).TCO.energy/1000000 Results(i).TCO.personnel/1000000 ...
                        Results(i).TCO.batteryreplacement/1000000 Results(i).TCO.insurance/1000000 Results(i).TCO.roadtax/1000000];
                    xTicker=[xTicker strcat('Concept ',string(i))];
                end
            end
            J=subplot(1,2,1)
            bar(barData,'stacked')
            legend('EOL','Acquisition','Cleaning','Maintenance','Energy','Personnel','Battery Replacement','Insurance','Road Tax')
            %xlabel(ax,'Vehicle')
            ylabel('Total Costs in Million SGD')
            colormap(parula(5))
            yt = get(J,'YTick');
            set(J, 'YTick', yt,'XTick',[1:i], 'XTickLabel', xTicker)
            textPosY=([barData(:,1)*0 cumsum(barData(:,2:end-1)')']+barData(:,2:end)./2);
            textPosY=[barData(:,1)/2 textPosY]';
            
            textPosX=(textPosY'.*0+[1:i]')';
            textData=barData';
            td=textData;
            j=(abs(barData')<0.01);
            textPosX(j)=[]
            textPosY(j)=[];
            textData(j)=[];
            textData=round(textData,2,'decimal');
            
            T=text(textPosX(:),textPosY(:),string(textData(:)),'HorizontalAlignment','center')
            for i=1:size(T,1)
                T(i).FontSize = 12;
            end
            
            
            %%
            
            barData=[];
            xTicker=[];
            
            if size(Results,1)==1
                n=Results.fleetsize;
                barData=[n*Results.TCO.EOL/1000000 n*Results.TCO.year0/1000000 ...
                    n*Results.TCO.cleaning/1000000 n*Results.TCO.maintenance/1000000 ...
                    n*Results.TCO.energy/1000000 n*Results.TCO.personnel/1000000 ...
                    n*Results.TCO.batteryreplacement/1000000 n*Results.TCO.insurance/1000000 n*Results.TCO.roadtax/1000000; nan(1,9)];
                xTicker={'Concept 1'};
            else
                for i=1:size(Results,1)
                    n=Results(i).fleetsize;
                    barData=[barData;n*Results(i).TCO.EOL/1000000 n*Results(i).TCO.year0/1000000 ...
                        n*Results(i).TCO.cleaning/1000000 n*Results(i).TCO.maintenance/1000000 ...
                        n*Results(i).TCO.energy/1000000 n*Results(i).TCO.personnel/1000000 ...
                        n*Results(i).TCO.batteryreplacement/1000000 n*Results(i).TCO.insurance/1000000 n*Results(i).TCO.roadtax/1000000];
                    xTicker=[xTicker strcat('Concept ',string(i))];
                end
            end
            J2=subplot(1,2,2)
            bar(barData,'stacked')
            legend(J2,'EOL','Acquisition','Cleaning','Maintenance','Energy','Personnel','Battery Replacement','Insurance','Road Tax')
            %xlabel(ax,'Vehicle')
            ylabel(J2,'Total Fleet Costs in Million SGD')
            %colormap(ax2,parula(5))
            yt = get(J2, 'YTick');
            set(J2, 'YTick', yt,'XTick',[1:i], 'XTickLabel', xTicker)
            textPosY=([barData(:,1)*0 cumsum(barData(:,2:end-1)')']+barData(:,2:end)./2);
            textPosY=[barData(:,1)/2 textPosY]';
            
            textPosX=(textPosY'.*0+[1:i]')';
            textData=barData';
            td=textData;
            j=(abs(barData')<0.01);
            textPosX(j)=[]
            textPosY(j)=[];
            textData(j)=[];
            textData=round(textData,2,'decimal');
            
            T=text(J2,textPosX(:),textPosY(:),string(textData(:)),'HorizontalAlignment','center')
            for i=1:size(T,1)
                T(i).FontSize = 12;
            end
            
        end
        function plotLCA(ax,Results)
            %% Vehicle level
            %% fleet level
            
            xTicker=[];
            barData=[];
            barDataY=[];
            barDataEOL=[];
            figure
            if size(Results,1)==1
                production=Results.vehicle.Productionemissions.CO2;
                distribution=Results.vehicle.Distributionemissions.CO2;
                EOL=Results.vehicle.EOLemissions.CO2;
                WTT=mean(Results.vehicle.WTTemissions);
                TTW=(Results.vehicle.TTWemissions);
                if isempty(TTW)
                    TTW=0;
                end
                n=1;
                barData=[[ production/1000 distribution/1000 WTT/1000 TTW/1000;0 0 0 0 ]];
                barDataY=[1 ;2];
                
                xTicker={'Concept 1'};
            else
                for i=1:size(Results,1)
                    n=1;%Results(i).fleetsize;
                    production=n*Results(i).vehicle.Productionemissions.CO2;
                    distribution=n*Results(i).vehicle.Distributionemissions.CO2;
                    EOL=n*Results(i).vehicle.EOLemissions.CO2;
                    WTT=n*mean(Results(i).vehicle.WTTemissions);
                    TTW=n*(Results(i).vehicle.TTWemissions);
                    if isempty(TTW)
                        TTW=0;
                    end
                    barData=[barData;[ production/1000 distribution/1000 WTT/1000 TTW/1000]];
                    barDataY=[barDataY;i];
                    barDataEOL=[barDataEOL;EOL/1000];
                    xTicker=[xTicker strcat('Concept ',string(i))];
                end
            end
            J1=subplot(1,2,1)
            barh(barDataY,barData,'stacked')
            hold on
            barh(barDataY,barDataEOL,'stacked')
            barData=[barDataEOL barData];
            yt = get(J1, 'YTick');
            set(J1, 'YTick', yt,'YTickLabel', xTicker)
            xlabel('Thousand kgCO2-eq')
            %legend('Production','Distribution','Well-to-Tank','Tank-to-Wheel','End-Of-Life')
            textPosX=([barData(:,1)*0 cumsum(barData(:,2:end-1)')']+barData(:,2:end)./2);
            textPosX=[barData(:,1)/2 textPosX];
            
            textPosY=(textPosX.*0+[1:i]');
            textData=barData;
            td=textData;
            j=(abs(barData)<0.01);
            textPosX(j)=[]
            textPosY(j)=[];
            textData(j)=[];
            textData=round(textData,2,'significant');
            
            T=text(textPosX(:),textPosY(:),string(textData(:)),'HorizontalAlignment','center')
            for i=1:size(T,1)
                T(i).FontSize = 12;
            end
            
            %% fleet level
            
            xTicker=[];
            barData=[];
            barDataY=[];
            barDataEOL=[];
            if size(Results,1)==1
                production=Results.vehicle.Productionemissions.CO2;
                distribution=Results.vehicle.Distributionemissions.CO2;
                EOL=Results.vehicle.EOLemissions.CO2;
                WTT=mean(Results.vehicle.WTTemissions);
                TTW=(Results.vehicle.TTWemissions);
                if isempty(TTW)
                    TTW=0;
                end
                n=1;
                barData=[[ production/1000 distribution/1000 WTT/1000 TTW/1000;0 0 0 0 ]];
                barDataY=[1 ;2];
                
                xTicker={'Concept 1'};
            else
                for i=1:size(Results,1)
                    n=Results(i).fleetsize;
                    production=n*Results(i).vehicle.Productionemissions.CO2;
                    distribution=n*Results(i).vehicle.Distributionemissions.CO2;
                    EOL=n*Results(i).vehicle.EOLemissions.CO2;
                    WTT=n*mean(Results(i).vehicle.WTTemissions);
                    TTW=n*(Results(i).vehicle.TTWemissions);
                    if isempty(TTW)
                        TTW=0;
                    end
                    barData=[barData;[ production/1000000 distribution/1000000 WTT/1000000 TTW/1000000]];
                    barDataY=[barDataY;i];
                    barDataEOL=[barDataEOL;EOL/1000000];
                    xTicker=[xTicker strcat('Concept ',string(i))];
                end
            end
            J2=subplot(1,2,2)
            barh(barDataY,barData,'stacked')
            hold on
            barh(barDataY,barDataEOL,'stacked')
            yt = get(J2, 'YTick');
            set(J2, 'YTick', yt,'YTickLabel', xTicker)
            barData=[barDataEOL barData];
            legend('Production','Distribution','Well-to-Tank','Tank-to-Wheel','End-Of-Life')
            xlabel('Million kgCO2-eq')
            textPosX=([barData(:,1)*0 cumsum(barData(:,2:end-1)')']+barData(:,2:end)./2);
            textPosX=[barData(:,1)/2 textPosX];
            
            textPosY=(textPosX.*0+[1:i]');
            textData=barData;
            td=textData;
            j=(abs(barData)<0.01);
            textPosX(j)=[]
            textPosY(j)=[];
            textData(j)=[];
            textData=round(textData,2,'significant');
            
            T=text(textPosX(:),textPosY(:),string(textData(:)),'HorizontalAlignment','center')
            for i=1:size(T,1)
                T(i).FontSize = 12;
            end
            
        end
        function [simdata,stopdata]=getsimdata(obj,results)
            
            numvehicles=results.fleetsize;
            vehicledata=cell(2,numvehicles);
            timetable=results.timetable;
            blocks=results.blocks;
            simdata=[];
            for i=1:numvehicles
                trips=timetable(blocks{1, i}(2:end-1)',:);
                vehicledata{1,i}=trips;
                SOC=reshape(blocks{14,i}',2,[])'; % find the pairs of charging/discharging
                SOC_dc=SOC(SOC(:,1)-SOC(:,2)>0,:);
                SOC_dc=SOC_dc(2:end,:); % ignore first and last trips to depot
                for j=1:size(trips,1)
                    
                    data=results.lines(1, trips(j,8)).simulationdata; %find route select data
                    dir=trips(j,6); % find direction 1/-1
                    tripno=trips(j,7);
                    if dir==1
                        SOC_trip=linspace(SOC_dc(j,1),SOC_dc(j,2),size(data.d1{tripno},1))';
                        vehicledata{2,i}=[vehicledata{2,i};[data.d1{tripno} ones(size(data.d1{tripno},1),1)*i SOC_trip]];
                    else
                        SOC_trip=linspace(SOC_dc(j,1),SOC_dc(j,2),size(data.d2{tripno},1))';
                        vehicledata{2,i}=[vehicledata{2,i};[data.d2{tripno} ones(size(data.d2{tripno},1),1)*i SOC_trip]];
                    end
                end
                simdata=[simdata;vehicledata{2,i}];
            end
            simdata=sortrows(simdata);
            numstops=sum([results.lines.busstops]); % no of total stops
            stopdata=cell(1,numstops); % initialise cell matrix - for boarding waiting alighting coord
            routestops=[results.lines.busstops]; %num stops per route
            numlines=size(results.lines,2);
            a=0;
            k=1;
            for i=1:numlines
                coords=results.lines(i).stopcoordinates;
                waiting=[ results.lines(1, i).totalpassengers(:,2)...
                    results.lines(1, i).totalpassengers(:,3) ...
                    results.lines(1, i).totalpassengers(:,5)+a ...
                    coords(results.lines(1, i).totalpassengers(:,5),:)...
                    ones(size(results.lines(1, i).totalpassengers(:,3),1),1)*i]; %arrival time stop coords
                for j=1:coords(:,1)
                    a1=waiting(waiting(:,3)== j,:);
                    if ~isempty(a1)
                        events=sortrows([a1(:,2) ones(size(a1,1),1);a1(:,1) -ones(size(a1,1),1);]);
                        events(:,3)=cumsum(events(:,2));
                        events(:,4)=a1(1,4);
                        events(:,5)=a1(1,5);
                        stopdata{1,a+j}=events;
                    end
                end
            end
        end
        function plotVehiclesfunction(obj,simdata,vehicles,frame,handle,nroutes,stopdata)
            numvehicles=size(vehicles,1);
            frame=ceil(frame);
            
            delete(handle.Children(1:end-1))%(1+nroutes)))
            for j=1:numvehicles
                h(j)=scatter(handle,[],[]); % vehicles
                k(j)=text(handle); % text
            end
            
            l=text(handle,104,1.2,'Time');
            %F(size(simdata,1)) = struct('cdata',[],'colormap',[]);
            i=frame;
            % find the last state of all vehicles
            data=ones(numvehicles,6);
            time=simdata(frame,1);
          %  plotstop(obj,stopdata,handle,time)
            j=simdata(frame,11); % vehicle id
            for k1=1:numvehicles
                if vehicles(k1)~=j
                    tempdata=simdata(simdata(:,1)<time & simdata(:,11)==vehicles(k1),:);
                    if ~isempty(tempdata)
                    data(k1,:)=[tempdata(end,11) tempdata(end,8) tempdata(end,9) tempdata(end,6) tempdata(end,7) ceil(tempdata(end,12)*100)];
                    end
                else
                    tempdata=simdata(frame,:);
                    if ~isempty(tempdata)
                    data(k1,:)=[tempdata(end,11) tempdata(end,8) tempdata(end,9) tempdata(end,6) tempdata(end,7) ceil(tempdata(end,12)*100)];
                    end
                end
            end
            for i=1:size(data,1)
                n=data(i,1); %numvehicle
                hold(handle,'on')
                %  hold on
                % h=routedata{2,n};
                h(i).XData=data(i,2);
                h(i).YData=data(i,3);
                % h.XData=routedata(i,8);
                %h.YData=routedata(i,9);
                if data(i,5)==1
                    %  if routedata(i,7)==1
                    color=	[0.466 .674 .188];
                elseif data(i,5)==0.5
                    %elseif routedata(i,7)==0.5
                    color=	[0.929 .694 .125];
                else
                    color='r';
                end
                %routedata{2,n}.MarkerFaceColor=color;
                %routedata{2,n}.MarkerEdgeColor=color;
                h(i).MarkerFaceColor=color;
                h(i).MarkerEdgeColor=color;
                %scatter(handle,routedata(i,6),routedata(i,7),30,'filled','b');
                %text(routedata(i,6), routedata(i,7),(cellstr(['  ' num2str(routedata(i,4)) ' Passengers'])))
                k(i).String=cellstr([' Vehicle' num2str(n)  newline  ' ' num2str(data(i,4)) ' Pax' ...
                    newline 'SOC_:' num2str(data(i,6)) '%']);%(cellstr(['  ' num2str(data(i,4)) ' Passengers']));
                k(i).Position=[data(i,2), data(i,3),0];
                % pause(0.2)
                [hr,min,sec]=hms(hours(time));
                if hr>=24
                    hr=0;
                end
                
                l.String=cellstr(['Time: ' datestr(duration([hr min sec]))]);
                %routedata{3,n}.String=(cellstr(['  ' num2str(data(i,6)) ' Passengers']));
                %routedata{3,n}.Position=[data(i,8), data(i,9),0];
                % pause(0.000000000000005)
                %drawnow
               
            end
            
        end
        function plotstop(obj,simdata,handle,time)
            
              for j=1:size(simdata,2)
               
                if ~isempty(simdata{j})
                id=find((simdata{j}(:,1)<time));
                if ~isempty(id)
                                  
                persons=(simdata{j}(id(end),3));
                pos=(simdata{j}(id(end),4:5));
                persons=(persons==0)*.01+persons*(persons~=0);
                h(j)=scatter(handle,pos(:,1),pos(:,2),persons*3,...
                'b'); % vehicles
               % k(j)=text(handle); % text
                end
                end
              end
              
              
              
        end
        function [value,property]=sumproperty(obj,property)
            % return the weightes sum of the properties
            property=struct2cell(property);
            g0=4;
            gu=10;
            S=obj.slope;
            W=obj.inflexion;
            for i=1:size(property,1)
                x=struct2array(property{i});
                sz=length(x);
                y{i}=((g0-gu)./(1+( (exp(S(i,1:sz).*(x-W(i,1:sz)))))))+gu;
                property{i,2}=y{i};
                 score(i)=sum(y{i}.*obj.catergorical_weights{i});
            end
            value=sum(cell2mat(y))/length(cell2mat(y)); %% add weights
            
        end
        function result=optimizeConcept(obj,app,OptVar)
            options=nsgaopt();
            options.Name='AEV Optimization'; % options name
            options.numObj=2;                % number of objectives
            options.popsize=500;             % population size
            options.maxGen=20;              % maximum generations
            options.numVar=size(OptVar,1);                % number of design variables
            options.numCons=1;               % number of constraints
            options.lb=cell2mat(OptVar(:,2))';              % lower bound
            options.ub=cell2mat(OptVar(:,3))';                % upper bound
            options.vartype=cell2mat(OptVar(:,4))';           % Variable type: 1float/2 integer
            options.nameObj={'TCO in SGD/passenger-km'...
                ,'Property Fulfillment'};
            options.nameVar=OptVar(:,1)';
            options.useParallel='yes';
            options.poolsize=18;              % Number of parallel workers
            options.objfun=@objectiveFunction;    % Objective function
            options.vectorized = 'yes';
            %  options.mutation = {'gaussian',0.5,0.2};
            options.plotInterval = 2;            % Plot interval
          %  options.initfun={@initpop,'populations.txt',1};
            options.plotInterval=1;
            result=nsga2(options,obj,app.Services,app.depotparameters,app.drivingcycle,app.Opt,app.OptimisationStateTable);
            
        end
        function [fitness,constraints]=objectiveFun(x,obj,app)
            y = [0,0];
            cons = [0,0];
            
            y(1) = x(1);
            y(2) = (1+x(2)) / x(1);
            
            % calculate the constraint violations
            c = x(2) + 9*x(1) - 6;
            if(c<0)
                cons(1) = abs(c);
            end
            
            c = -x(2) + 9*x(1) - 1;
            if(c<0)
                cons(2) = abs(c);
            end
            fitness=y;
            constraints=cons;
        end
        function [fitness,constraints]=objectiveFunction(vehicleparameters,obj,services,depotparameters,drivingcycle,plothandle,tablehandle)
            fitness=[0,0];
            try
                vehicleparameters(12)=0;
                VC=depot(services,depotparameters,vehicleparameters,drivingcycle,[]);
                [fitness(2),VC.property]=sumproperty(obj,VC.property);
                fitness(2)=-fitness(2);
                fitness(1)=VC.TCO.passengerkm;
                constraints=0;
                dlmwrite('res.txt',[vehicleparameters(14) VC.fleetsize ...
                   VC.gCO2_passengerkm VC.TCO.passengerkm VC.vehicle.Energyconsumption VC.dailyvehiclekm ],'-append');

            catch
                constraints=010;
                fitness(1)=4;
                fitness(2)=2;
            end
            %% Add property eval, constraints
            
        end
        function plotPopulation(obj,data, gen,handle,tabhandle)
            % Function: plotPopulation(handles, gen)
            % Description: Plot population with the first two objective values.
            % Parameters:
            %   gen : the generation ID
            %
            %         LSSSSWC, NWPU
            %    Revision: 1.3  Data: 2011-07-26
            %*************************************************************************
            
            cla(handle)
            
            %*************************************************************************
            % Initialize data
            %*************************************************************************
            pop     = data.pops(gen, :);
            obj     = vertcat(pop.obj);
            popH    = data.pops(1:gen, :);
            objH     = vertcat(popH.obj);
            numObj  = data.opt.numObj;
            strObj1 = 'objective 1';
            strObj2 = 'objective 2';
            strObj3 = 'objective 3';
            
            %  When the result is read from file, there is no 'opt' structure.
            if(isfield(data, 'opt'))
                opt = data.opt;
                maxGen  = opt.maxGen;
                if( ~isempty(opt.nameObj) )
                    strObj1 = ['obj 1 : ', opt.nameObj{1}];
                    strObj2 = ['obj 2 : ', opt.nameObj{2}];
                    if(numObj == 3)
                        strObj3 = ['obj 3 : ', opt.nameObj{3}];
                    end
                end
            else
                maxGen = size(data.pops, 1);
            end
            
            % Determin if reference points exist
            refPoints = [];
            refPlotStyle = {'kd', ...
                'LineWidth', 1,...
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', 'g',...
                'MarkerSize',6};
            if( isfield(data, 'opt') && ~isempty(data.opt.refPoints) )
                refPoints = data.opt.refPoints;
            end
            
            
            %*************************************************************************
            % Plot population with different methods for every "numObj" number
            %*************************************************************************
            if(numObj == 2)
                scatter(handle,objH(:,1),-objH(:,2),'filled', 'markerFaceAlpha',0.3,'MarkerFaceColor',[128 193 219]./255);
                hold(handle,'on')
                plot(handle,obj(:,1), -obj(:,2), 'ok');
                xlabel(handle,strObj1, 'interpreter', 'none');
                ylabel(handle,strObj2, 'interpreter', 'none');
                drawnow
                % plot reference points
                if(~isempty(refPoints))
                    hold on
                    plot(refPoints(:, 1), refPoints(:, 2), refPlotStyle{:});
                end
            elseif(numObj == 3)
                [az, el] = view;    % Save the last view angles
                plot3(obj(:,1), obj(:,2), obj(:,3), 'ob');
                view(az,el);        % Restore the last view angles
                xlabel(strObj1, 'interpreter', 'none');
                ylabel(strObj2, 'interpreter', 'none');
                zlabel(strObj3, 'interpreter', 'none');
                
                % plot reference points
                if(~isempty(refPoints))
                    hold on
                    plot3(refPoints(:, 1), refPoints(:, 2), refPoints(:, 3), refPlotStyle{:});
                end
            else
                plot(handle,obj', 'b-');
                
                % set(gca, 'XGrid', 'on');
                xlabel('Objective number');
                ylabel('Objective value');
                
                % plot reference points
                if(~isempty(refPoints))
                    hold on
                    refPlotStyle{1} = 'gd-';
                    plot(refPoints', refPlotStyle{:});
                end
            end
            strTitle = sprintf('Generation %d / %d', gen, maxGen);
            %             if(handles.bLoadFromFile == 1)
            %                 strTitle = sprintf('%s\nLoad from : %s', strTitle, handles.strPopFile);
            %             end
            title(handle,strTitle, 'interpreter', 'none');
            
            %% update table
            tabhandle.Data=table(fieldnames(data.states),struct2cell(data.states(gen,:)));
            
        end
    end
end