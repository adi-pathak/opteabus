classdef propertyevaluation
   
    
    
    
end

function PF = PropertiesEvaluation (Weighting,Basisfahrzeug)
% Vehicle property evaluation
% Daniel Telschow
% 8th March
% IMVS - TUMCreate Singapore
MinReq=Basisfahrzeug.UInput.Target;
MinReq=[5 5 5 10 5 5];
%% Everyday Practicality
% 1. Road damage
AL = [300 890 4450 6230 8000 8900 13340];
RD = [0 0.0003 0.118 0.399 1 1.4 7.9];
ESAL = interp1(AL, RD, Basisfahrzeug.chassis.AxleLoad);
Phy.rd = [9e9 10 7.5 0];
Not.rd = [0 0 5 10];
Prop(1) = interp1(Phy.rd, Not.rd, ESAL);
Target(1) = MinReq(1);

% 2. Turning radius
if (Basisfahrzeug.chassis.numaxles == 2) && (Basisfahrzeug.chassis.RearWheelSteering == 1) % 4-WheelSteering
    tr = (Basisfahrzeug.dimensions.wheelbase/2)^2 + ((Basisfahrzeug.dimensions.wheelbase)*...
        ((cotd(46)+cotd(35))/2))^2;
else
    tr = Basisfahrzeug.dimensions.wheelbase/sind(22);
end
Phy.tr = [9e9 18.99 19 13 0];
Not.tr =[0 0 5 10 10];
Prop(2) = interp1(Phy.tr, Not.tr, tr);
Target(2)=MinReq(1);

% 3. Road space usage
Phy.rsu = [9e9 6 5.2 2.6 0];
Not.rsu = [0 0 5 10 10];
rsu = (Basisfahrzeug.dimensions.width *...
    Basisfahrzeug.dimensions.length)/ (Basisfahrzeug.interior.passengers * 0.17);
Prop(3) = interp1(Phy.rsu, Not.rsu, rsu);
Target(3)=MinReq(1);

%% Design
% 4. Aspect ratio
Phy.ar = [0 0.57 0.75 9e9];
Not.ar = [0 5 10 10];
ar = Basisfahrzeug.dimensions.width/Basisfahrzeug.dimensions.height;
Prop(4) = interp1(Phy.ar, Not.ar, ar);
Target(4)=MinReq(2);

% 5. Wheelbase ratio
Phy.wr = [0 0.4 0.6 9e9];
Not.wr = [0 5 10 10];
wr = Basisfahrzeug.dimensions.wheelbase/Basisfahrzeug.dimensions.length;
Prop(5) = interp1(Phy.wr, Not.wr, wr);
Target(5)=MinReq(2);

%% Longitudinal Dynamics
% 6. Top speed
Phy.vmax = [0 59 60 90 9e9];
Not.vmax = [0 0 5 10 10];
Prop(6) = interp1(Phy.vmax, Not.vmax, Basisfahrzeug.dynamics.vmax);
Target(6)=MinReq(3);

% 7. Max. Acceleration
Phy.amax = [0 0.6 1.5 1.75 9e9];
Not.amax = [0 5 10 10 10];
Prop(7) = interp1(Phy.amax, Not.amax, Basisfahrzeug.dynamics.amax);
Target(7)=MinReq(3);

% 8. Gradeability
Phy.grade = [0 9.9 10 20 9e9];
Not.grade = [0 0 5 10 10];
Prop(8) = interp1(Phy.grade, Not.grade, Basisfahrzeug.dynamics.grade);
Target(8)=MinReq(3);

%% Comfort
% 9. Seating-Standing ratio
Prop(9) = Basisfahrzeug.interior.CrowdingComfort;
Target(9)=MinReq(4);

% 10. Headspace
Phy.hs = [0 0.09 0.1 0.3 9e9];
Not.hs = [0 0 5 10 10];
hs=Basisfahrzeug.dimensions.height-(Basisfahrzeug.numdecks*1.8);
Prop(10) = interp1(Phy.hs, Not.hs, hs);
Target(10)=MinReq(4);

% 11. Seat width
Phy.sw = [0 0.39 0.4 0.45 9e9];
Not.sw = [0 0 5 10 10];
sw = Basisfahrzeug.interior.seatwidth;
Prop(11) = interp1(Phy.sw, Not.sw, sw);
Target(11)=MinReq(4);

% 12. Seat pitch
Phy.sp = [0 0.64 0.65 0.71 0.76 0.81 0.91 9e9];
Not.sp = [0 0 5 5.8 7.2 8.8 10 10];
sp = Basisfahrzeug.interior.pitch;
Prop(12) = interp1(Phy.sp, Not.sp, sp);
Target(12)=MinReq(4);

%% Availability
% 13. Peak hour frequency
Phy.phf = [0 3 6 9e9];
Not.phf = [0 5 10 10];
phf = Basisfahrzeug.infrastructure.MaxFrequency;
Prop(13) = interp1(Phy.phf, Not.phf, phf);
Target(13)=MinReq(5);

% 14. Off-Peak hour frequency
Phy.ohf = [0 0.3 1 9e9]; 
Not.ohf = [0 5 10 10];
ohf = Basisfahrzeug.infrastructure.MinFrequency;
Prop(14) = interp1(Phy.ohf, Not.ohf, ohf);
Target(14)=MinReq(5);

%% Accessibility
% 15. Dwelling time
Phy.dwt = [9e9 61 60 20 0];
Not.dwt = [0 0 5 10 10];
Prop(15) = interp1(Phy.dwt, Not.dwt, Basisfahrzeug.interior.dwellingtime);
Target(15)=MinReq(6);

% 16. Number of doors
% Phy.nd = [0 1 4 9e9];
% Not.nd = [0 5 10 10];
% Prop(16) = interp1(Phy.nd, Not.nd, Basisfahrzeug.interior.numdoors);
% Target(16)=MinReq(6);

% 17. Door width
% Phy.dw = [0 0.64 0.65 1.2 9e9];
% Not.dw = [0 0 5 10 10];
% Prop(17) = interp1(Phy.nd, Not.nd, Basisfahrzeug.interior.doorwidth);
% Target(17)=MinReq(6);

% 18. Gangway width
% Phy.gan = [0 0.44 0.45 0.55 9e9];
% Not.gan = [0 0 5 10 10];
% Prop(18) = interp1(Phy.gan, Not.gan,Basisfahrzeug.interior.gangwidth);
% Target(18)=MinReq(6);

% 19. WCU zones
Phy.wcu = [0 1 2 9e9];
Not.wcu = [0 5 10 10];
wcu = Basisfahrzeug.interior.standing/5;
if wcu > 4
    wcu = 4;
end
if Basisfahrzeug.interior.layout==1
    wcu = 1;
end

Prop(16) = interp1(Phy.wcu, Not.wcu, wcu);
Target(16)=MinReq(6);

%% Comfort properties adjustment for No Seating
if Basisfahrzeug.interior.layout==4
    Prop(11)=5; % Seat width = min
    Prop(12)=5; % Seat pitch = min
%     Prop(18)=5; % Gangway width = min
end

%% Property Evaluation
n=length(Prop);
diff=zeros(n,1);
for j=1:n
    diff(j)=(Prop(j)-Target(j));
end

Ind_neg = find(diff<0);
diff(Ind_neg) = -1*((-1*diff(Ind_neg)).^Weighting.diffpotenz)/...
    (5^Weighting.diffpotenz);
Ind_pos = find(diff>=0);
diff(Ind_pos) = Weighting.verh_ue_u*(diff(Ind_pos).^Weighting.diffpotenz)/...
    (5^Weighting.diffpotenz);

%% Group evaluation
% Everyday Pracicality
ep = ((diff(1) * Weighting.rd^Weighting.gewpotenz) +...
    (diff(2) * Weighting.tr^Weighting.gewpotenz) +...
    (diff(3) * Weighting.rsu^Weighting.gewpotenz))/...
    (Weighting.rd^Weighting.gewpotenz +...
    Weighting.tr^Weighting.gewpotenz +...
    Weighting.rsu^Weighting.gewpotenz);

% Design
des = ((diff(4) * Weighting.ar^Weighting.gewpotenz) +...
    (diff(5) * Weighting.wr^Weighting.gewpotenz))/...
    (Weighting.ar^Weighting.gewpotenz + Weighting.wr^Weighting.gewpotenz);

% Longitudinal Dynamics
long = ((diff(6) * Weighting.vmax^Weighting.gewpotenz) +...
    (diff(7) * Weighting.amax^Weighting.gewpotenz) +...
    (diff(8) * Weighting.grade^Weighting.gewpotenz))/...
    (Weighting.vmax^Weighting.gewpotenz +...
    Weighting.amax^Weighting.gewpotenz +...
    Weighting.grade^Weighting.gewpotenz);

% Comfort
com = ((diff(9) * Weighting.ssr^Weighting.gewpotenz) + ...
    (diff(10) * Weighting.hs^Weighting.gewpotenz) + ...
    (diff(11) * Weighting.sw^Weighting.gewpotenz) + ...
    (diff(12) * Weighting.sp^Weighting.gewpotenz))/ ...
    (Weighting.ssr^Weighting.gewpotenz + ...
    Weighting.hs^Weighting.gewpotenz + ...
    Weighting.sw^Weighting.gewpotenz + ...
    Weighting.sp^Weighting.gewpotenz);

% Availability
av = ((diff(13) * Weighting.phf^Weighting.gewpotenz) +...
    (diff(14) * Weighting.ohf^Weighting.gewpotenz))/...
    (Weighting.phf^Weighting.gewpotenz + Weighting.ohf^Weighting.gewpotenz);

% Accessibility
ac = ((diff(15) * Weighting.dwt^Weighting.gewpotenz) + ...
    (diff(16) * Weighting.wcu^Weighting.gewpotenz))/ ...
    (Weighting.dwt^Weighting.gewpotenz + Weighting.wcu^Weighting.gewpotenz);

% Prop1 = (Prop(1) + Prop(2) + Prop(3))/3;
% Prop2 = (Prop(4) + Prop(5))/2;
% Prop3 = (Prop(6) + Prop(7) + Prop(8))/3;
% Prop4 = (Prop(9) + Prop(10) + Prop(11) + Prop(12))/4;
% Prop5 = (Prop(13) + Prop(14))/2;
% Prop6 = (Prop(15) + Prop(16) + Prop(17) + Prop(18) + Prop(19))/5;
% Basisfahrzeug.Properties=[Prop1 Prop2 Prop3 Prop4 Prop5 Prop6];

%% Property fulfillment
PF = 100 + ((ep * Weighting.ep^Weighting.gewpotenz +...
    des * Weighting.des^Weighting.gewpotenz +...
    long * Weighting.long^Weighting.gewpotenz +...
    com * Weighting.com^Weighting.gewpotenz +...
    av * Weighting.av^Weighting.gewpotenz +...
    ac * Weighting.ac^Weighting.gewpotenz) / ...
    (Weighting.ep^Weighting.gewpotenz +...
    Weighting.des^Weighting.gewpotenz +...
    Weighting.long^Weighting.gewpotenz +...
    Weighting.com^Weighting.gewpotenz +...
    Weighting.av^Weighting.gewpotenz +...
    Weighting.ac^Weighting.gewpotenz)) * 100;

if Basisfahrzeug.dynamics.amax==0
    PF=0;
elseif Basisfahrzeug.dynamics.grade==0
    PF=0;
end
end