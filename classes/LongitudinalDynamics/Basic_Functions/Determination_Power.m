%


function [P_EL_Motor, DATA] = Determination_Power (V, Motor)

NumofLines=length(V(:,1));

DATA=zeros(2,length(V(1,:)));
P_EL_Motor=zeros(NumofLines-1,length(V(1,:)));

[xd,yd]=ndgrid(Motor.eff_n_axis,Motor.eff_T_axis);
% Motor.eff(1:160,:)=(1/0.9025);
% Motor.eff(161:end,:)=(0.9025);
GRIDINT=griddedInterpolant(xd,yd, Motor.eff','linear','none' );
% for i=1:length(V(1,:))
%     if V(1,i)>0
%         EFF(i)=0.9025;
%     else
%         EFF(i)=1/0.9025;
%     end
% end
EFF=GRIDINT(V(2,:).*60, V(1,:));
%EFF(isnan(EFF))=1;
if max(V(1,:))>max(Motor.eff_T_axis)
    disp(strcat('Max Motor Torque exceeded by ',string(max(V(1,:))-max(Motor.eff_T_axis))))
end
P_EL_Motor(1,:)=V(1,:).*V(2,:).*2.*pi./EFF;

if NumofLines>2
    P_EL_Motor(2:NumofLines-1,:)=V(3:NumofLines,:);
end

DATA(:,:)=V(1:2,:); % DATA-Matrix mit Moment und Drehzahl zu den einzelnen Leistungen zur Nachvollziehung

end