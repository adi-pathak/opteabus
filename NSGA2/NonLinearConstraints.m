function [FunOut] = NonLinearConstraints (Basisfahrzeug)
Err=zeros(6,1);

    % Constraint no.2
if (Basisfahrzeug.rigid==1) && (Basisfahrzeug.numdecks == 2)
    Err(2)=1;
end    
    % Constraint no.3
if (Basisfahrzeug.dimensions.length > 12) && (Basisfahrzeug.chassis.numaxles==2)
    Err(3)=1;
end        
    % Constraint no.4
if (Basisfahrzeug.rigid == 1) && (Basisfahrzeug.chassis.numaxles == 2)
    Err(4)=1;
end    
    % Constraint no.5
if (Basisfahrzeug.dimensions.length-Basisfahrzeug.dimensions.wheelbase) < 0.6
    Err(5)=1;
end 
    % Constraint no.5
if (Basisfahrzeug.numdecks == 2) && (Basisfahrzeug.interior.layout==4) 
    Err(6)=1;
end
FunOut=sum(Err);