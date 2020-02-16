function [] = draw_RefLine_TG(Landmark_RefLine, Landmark_TGPoint,Vertices_surface,str_figure,n_subj,str_side)
%draw_mesh_coord Summary of this function goes here
%   Detailed explanation goes here

figure; plot(Vertices_surface.PxxInterp,Vertices_surface.PyyInterp,'r-',Vertices_surface.sorted_pointset(:,1),Vertices_surface.sorted_pointset(:,2),'b*')
axis equal;
testPoints1 = [Landmark_RefLine.LinePoint1(1) Landmark_RefLine.LinePoint2(1)];
testPoints2 = [Landmark_RefLine.LinePoint1(2) Landmark_RefLine.LinePoint2(2)];
figure; plot(Vertices_surface.PxxInterp,Vertices_surface.PyyInterp,'b*',testPoints1,testPoints2,'k-','LineWidth',2);
axis equal;
hold on;
% SimulatedLinePoints11 = [Landmark_TGPoint.TGPoint(1) Landmark_TGPoint.RefLinePoint11(1)];
% SimulatedLinePoints12 = [Landmark_TGPoint.TGPoint(2) Landmark_TGPoint.RefLinePoint11(2)];
% plot(SimulatedLinePoints11,SimulatedLinePoints12,'g-','LineWidth',2);
xmin = min(Vertices_surface.PxxInterp);
xmax = max(Vertices_surface.PxxInterp);
f1 = polyval(Landmark_TGPoint.p1,xmin:xmax);
f2 = polyval(Landmark_TGPoint.p2,xmin:xmax);
plot(xmin:xmax,f1,'g-','LineWidth',2);
plot(xmin:xmax,f2,'g-','LineWidth',2);
% SimulatedLinePoints21 = [Landmark_TGPoint.TGPoint(1) Landmark_TGPoint.RefLinePoint22(1)];
% SimulatedLinePoints22 = [Landmark_TGPoint.TGPoint(2) Landmark_TGPoint.RefLinePoint22(2)];
% plot(SimulatedLinePoints21,SimulatedLinePoints22,'g-','LineWidth',2);

plot(Landmark_TGPoint.TGPoint(1),Landmark_TGPoint.TGPoint(2),'ro','LineWidth',2,'MarkerFaceColor',[0.5,0.5,0.5]);
SimulatedLinePointsX = [Landmark_TGPoint.RefLinePoint11(1) Landmark_TGPoint.RefLinePoint22(1)];
SimulatedLinePointsY = [Landmark_TGPoint.RefLinePoint11(2) Landmark_TGPoint.RefLinePoint22(2)];
plot(SimulatedLinePointsX,SimulatedLinePointsY,'rs','LineWidth',2,'MarkerFaceColor',[0.5,0.5,0.5]);

plot(testPoints1,testPoints2,'ks','LineWidth',2,'MarkerFaceColor',[0.5,0.5,0.5]);

hold off;

str_file = [str_figure num2str(n_subj) 'TG_RL_' str_side '.png'];
saveas(gcf,str_file);
end

