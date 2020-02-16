function [ V_f_new,IA,axEst ] = reflect_PointsOnAxis( F_f, V_f, V_t, V_p,V_new)
%reflect_PointsOnAxis reflect the points based on the axis of moment of
%inertial
%   Detailed explanation goes here
axEst = calculateExternalFrame(V_f, V_p, V_t) ;
[CM,IA,SV] = computeInertialAxes(F_f,V_f);
IA(:,3) = sign(dot(axEst.X,IA(:,3))) * IA(:,3)';
IA(:,1) = sign(dot(axEst.Y,IA(:,1))) * IA(:,1)';
IA(:,2) = sign(dot(axEst.Z,IA(:,2))) * IA(:,2)';
[a] = CM - (IA(:,1))' *100;
[b] = CM + (IA(:,1))' *400;
[c] = CM - (IA(:,2))' *100;
[d] = CM + (IA(:,2))' *400;
[e] = CM - (IA(:,3))' *100;
[f] = CM + (IA(:,3))' *400;
plotsurf(V_f,F_f);hold on; plot3([a(1);b(1)],[a(2);b(2)],[a(3);b(3)],'r');hold on;plot3([c(1);d(1)],[c(2);d(2)],[c(3);d(3)],'g');
plot3([e(1);f(1)],[e(2);f(2)],[e(3);f(3)],'b');
A = zeros(4);
A(1:3,1:3) = [IA(:,3)'; IA(:,1)';IA(:,2)'];
A(4,1:3) = CM;
A(4,4) = 1;
tformO2A = affine3d(A);
if(nargin>4)
    x = V_new(:,1);
    y = V_new(:,2);
    z = V_new(:,3);
    V_f_new = V_new;
    [u,v,w] = transformPointsInverse(tformO2A,x,y,z);

    B = [1 0 0 0;0 1 0 0; 0 0 -1 0;0 0 0 1];
    tformReflect = affine3d(B);
    [uReflect,vReflect,wReflect] = transformPointsForward(tformReflect,u,v,w);

    [uN,vN,wN] = transformPointsForward(tformO2A,uReflect,vReflect,wReflect);

    V_f_new(:,1) = uN;
    V_f_new(:,2) = vN;
    V_f_new(:,3) = wN;
else
    x = V_f(:,1);
    y = V_f(:,2);
    z = V_f(:,3);
    V_f_new = V_f;
    [u,v,w] = transformPointsInverse(tformO2A,x,y,z);

    B = [1 0 0 0;0 1 0 0; 0 0 -1 0;0 0 0 1];
    tformReflect = affine3d(B);
    [uReflect,vReflect,wReflect] = transformPointsForward(tformReflect,u,v,w);

    [uN,vN,wN] = transformPointsForward(tformO2A,uReflect,vReflect,wReflect);

    V_f_new(:,1) = uN;
    V_f_new(:,2) = vN;
    V_f_new(:,3) = wN;
end


end

