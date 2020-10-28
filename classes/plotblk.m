n=1;
            block_i=blocks(:,1);
            [a,b]=(sortrows(([(reshape(block_i{6},2,[])) (reshape(block_i{7},2,[])) (reshape(block_i{8},2,[]))])',1));
            state=[(reshape(block_i{6},2,[]))*0+1 (reshape(block_i{7},2,[]))*0+2 (reshape(block_i{8},2,[]))*0+3];%state
            state=[0 state(1,b)];
            %y=zeros(n+1,size(state,2));
            y=zeros(size(block_i,2),size(state,2));
            y_1=[block_i{6}(1) diff(sortrows(([(reshape(block_i{6},2,[])) (reshape(block_i{7},2,[])) (reshape(block_i{8},2,[]))])',1)')];
            %y_1= duration(hours(y_1),'format','hh:mm')
            y(n,:)=y_1;
            x=1:size(y,1);
            x=[x;nan];
            y=[y;nan(1,size(y,2))];
            
            H=barh(x,y,'stacked');
            set(H(1),'visible','off')
            set(H(find(state==2)),'facecolor','blue')
            set(H(find(state==3)),'facecolor','red')
            set(H(find(state==1)),'facecolor','green')
            %l=legend(H,'Deadheading','Service Trip','Waiting');
            %app.VehicleBlocksAxes.Color=[0.941 0.941 0.941];
            xlabel('Time')
            vn=1;
            ylabel(strcat('Vehicle ',vn))
            hold on
            yyaxis right
            plot(block_i{13, 1} ,block_i{12, 1}/batterycapacity,'LineWidth',2,'Color','k')