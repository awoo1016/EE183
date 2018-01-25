%1-DOF variables, each column in q corresponds to the one variable for each
%joint. The values are in degrees.
q = [0 0 0 0];
T05 = myFK(q);
%position vector of end effector for IK. The values are in cm.
p = [0; 38; 0];
%gets the approximate angles closest to the given position vector p
IK = myIK(p);

function FK = myFK(q)
%Transformation Matrices
T01 = [cosd(q(1)) -sind(q(1)) 0 0; sind(q(1)) cosd(q(1)) 0 0; 0 0 1 0; 0 0 0 1];

T12 = [cosd(q(2)) -sind(q(2)) 0 30; sind(q(2)) cosd(q(2)) 0 0; 0 0 1 0; 0 0 0 1];

T23 = [cosd(q(3)) -sind(q(3)) 0 0; sind(q(3)) cosd(q(3)) 0 30; 0 0 1 0; 0 0 0 1];

T34 = [cosd(q(4)) -sind(q(4)) 0 0; 0 0 -1 8; sind(q(4)) cosd(q(4)) 0 0; 0 0 0 1];

T45 = [0 -1 0 -4; 0 0 1 0; -1 0 0 0; 0 0 0 1];
%matrix multiply transformation matrices to get end effector position
FK= T01*T12*T23*T34*T45;
end

function IK = myIK(p)
%Check if specified position is in operation space
if (sqrt(p(1)^2 + p(2)^2 + p(3)^2)>68)||p(3)>4
    IK = Inf(4,4);
    return;
end
%starting point is at rest position
q = [0 0 0 0];
%compute FK
FK = myFK(q);
x = [FK(13); FK(14); FK(15)];
%Jacobian calculated by using 1 degree increments for q values
J = [-0.6672 -0.6626 -0.139 0.006; 0.4480 -0.0756 -0.071 0; 0 0 0 -0.0698];
%calculate psuedo inverse Jacobian
Jinv = (J.'*J)\J.';
%get approximate q values. Accuracy can be increased by calculating
%Jacobian with smaller difference angles.
while isequal(round(p-x),[0; 0; 0])~=1
    dx = (p-x)/100;
    dq = Jinv*dx;
    q = q + dq;
    FK = myFK(q);
    x = [FK(13); FK(14); FK(15)];
end
%return final transformation matrix
IK = FK;
end
