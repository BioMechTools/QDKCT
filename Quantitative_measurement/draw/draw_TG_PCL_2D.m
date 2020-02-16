function [] = draw_TG_PCL_2D(Landmark_TGPoint,Landmark_RefLine,Landmark_contour,manual_TGPoint)
%draw_TG_PCL_2D is aimed to compare the manual TG and automatic TG
%   Detailed explanation goes here
PxxInterp = Landmark_contour.PxxInterp;
PyyInterp = Landmark_contour.PyyInterp;
TGPoint = Landmark_TGPoint.TGPoint ;
RefLinePoint11 = Landmark_TGPoint.RefLinePoint11;
RefLinePoint22 = Landmark_TGPoint.RefLinePoint22;
testPoints1 = [RefLinePoint11(1) RefLinePoint22(1)];
testPoints2 = [RefLinePoint11(2) RefLinePoint22(2)];
LinePoint1 = Landmark_RefLine.LinePoint_Medial;
LinePoint2 = Landmark_RefLine.LinePoint_Lateral;

figure;plot(PxxInterp,PyyInterp,'b*',testPoints1,testPoints2,'r-');
hold on;
plot(TGPoint(1),TGPoint(2),'r*','markersize',10);
plot(manual_TGPoint(1:3,2),manual_TGPoint(1:3,3),'g*');
plot(manual_TGPoint(4,2),manual_TGPoint(4,3),'k*');
plot(manual_TGPoint(5,2),manual_TGPoint(5,3),'m*');
plot(LinePoint2(1),LinePoint2(2),'c*','markersize',10);
plot(LinePoint1(1),LinePoint1(2),'y*','markersize',10);
hold off;

end

