function [Landmark_RefLine, Landmark_TGPoint,Landmark_contour,Vertices_surface] = get_TG_PCL(rotM,femCoords,V_csL,cellV,secI,str_Side,n_target_ind)
%get_Ref_TG is aimed to get the reference line from posterior condyle and
%the TG point
%   Detailed explanation goes here
%   The points on the reference line, point 2 is lateral and point 1 is
%   medial
% rotM = eye(3) * rotM' ;  
% [CA,sectionCenters,CMonCAaxis,V_csL,cellV,secI] = computeCrossSectionAreaAlongAxis_new(V,rotM,stepSize_perc) ;
% rotM =rotM;
% rotM = eye(3) * R_axes' ; 
%rotate vertices in the relevant frame
% manual_TGPoint = (rotM * manual_TGPoint')' ;

maxSliceVertices = cellV{n_target_ind};
[ sorted_pointset ] = sort_2Dpoints( maxSliceVertices );
pointsetTemp = zeros(size(sorted_pointset,1)+10,2);
pointsetTemp(1:5,:) = sorted_pointset(end-4:end,:);
pointsetTemp(end-4:end,:) = sorted_pointset(1:5,:);
pointsetTemp(6:end-5,:) = sorted_pointset;
Pxx = pointsetTemp(:,1);%smooth(pointsetTemp(:,1));
Pyy =pointsetTemp(:,2);%smooth(pointsetTemp(:,2));
numLength = length(sorted_pointset);
x = 1:1:numLength+10;
xx = 1:0.25:numLength+10;
PxxInterp = spline(x,Pxx,xx);
PyyInterp = spline(x,Pyy,xx);
PxxInterp = smooth(PxxInterp);
PxxInterp = PxxInterp(16:end-20)';
PyyInterp = smooth(PyyInterp);
PyyInterp = PyyInterp(16:end-20)';
maxSliceVertices = [PxxInterp' PyyInterp'];
[ sorted_pointset ] = sort_2Dpoints( maxSliceVertices );
% figure; plot(PxxInterp,PyyInterp,'r-',sorted_pointset(:,1),sorted_pointset(:,2),'b*')
Vertices_surface.PxxInterp = sorted_pointset(:,1)';
Vertices_surface.PyyInterp = sorted_pointset(:,2)';
Vertices_surface.sorted_pointset = sorted_pointset;

theseV = sorted_pointset;
[k] = convhull(theseV(:,1),theseV(:,2));
[k_list] = convhull_list(k,length(PxxInterp));
kDiff = cell2mat(cellfun(@length,k_list,'uni',false));%abs(k(2:end)-k(1:end-1));
% if(abs(k(2)-k(1)) > abs(k(end)-k(end-1)))
%     kDiff(1) = abs(k(2)-size(theseV,1));
% else
%     kDiff(end) = abs(k(end-1)-size(theseV,1));
% end

[maxDiff,maxDiffInd] = max(kDiff);
indlist = k_list{maxDiffInd};
LinePoint1 =  theseV(indlist(1),:);
LinePoint2 =  theseV(indlist(end),:);
%% 2D to 3D
V1 = V_csL(:,1) ;
targetV1 = V1(secI==n_target_ind,:);
PD_pos = mean(targetV1(:));
LinePoint1_3d = [PD_pos, LinePoint1(1), LinePoint1(2)];
LinePoint2_3d = [PD_pos, LinePoint2(1), LinePoint2(2)];
LinePoint1_3d = rotM * LinePoint1_3d';
LinePoint2_3d = rotM * LinePoint2_3d';
testPoints1 = [LinePoint1(1) LinePoint2(1)];
testPoints2 = [LinePoint1(2) LinePoint2(2)];
% figure;plot(PxxInterp,PyyInterp,'b*',testPoints1,testPoints2,'r-');
%figure;plot(PxxInterp,PyyInterp,'b*',theseV(k,1),theseV(k,2),'r-');
%%% determine which point is lateral
directionRL = LinePoint2_3d - LinePoint1_3d;
directionRL = directionRL/norm(directionRL);
signFlag =sign(dot(femCoords.ML(1:2),directionRL(1:2)));
directionRL2D = (LinePoint2-LinePoint1);
directionRL2D = directionRL2D/norm(directionRL2D);
if strcmp(str_Side,'R')
    if(signFlag<0)
        temp = LinePoint2_3d;
        LinePoint2_3d = LinePoint1_3d;
        LinePoint1_3d = temp;
        directionRL2D = -directionRL2D;
        temp = LinePoint2;
        LinePoint2 = LinePoint1;
        LinePoint1 = temp;
    end
else
    if(signFlag>0)
        temp = LinePoint2_3d;
        LinePoint2_3d = LinePoint1_3d;
        LinePoint1_3d = temp;
        directionRL2D = -directionRL2D;
        temp = LinePoint2;
        LinePoint2 = LinePoint1;
        LinePoint1 = temp;
    end
end
Landmark_RefLine.LinePoint_Lateral_3d = LinePoint2_3d;
Landmark_RefLine.LinePoint_Medial_3d = LinePoint1_3d;
Landmark_RefLine.LinePoint_Lateral = LinePoint2;
Landmark_RefLine.LinePoint_Medial = LinePoint1;

% [ TGPoint, RefLinePoint11, RefLinePoint22,TGPoint3D,p1,p2] = getTG_Femur( theseV, k,PD_pos,rotM );

[ TGPoint, RefLinePoint11, RefLinePoint22,TGPoint3D,p1,p2] = getTG(directionRL2D,theseV,k,k_list,rotM,PD_pos,str_Side);

% figure;plot(PxxInterp,PyyInterp,'b*',testPoints1,testPoints2,'r-');
% hold on;
% plot(TGPoint(1),TGPoint(2),'r*','markersize',10);
% plot(manual_TGPoint(1),manual_TGPoint(2),'g*');
% plot(LinePoint2(1),LinePoint2(2),'m*','markersize',10);
% plot(LinePoint1(1),LinePoint1(2),'k*','markersize',10);
% hold off;
Landmark_contour.PxxInterp = PxxInterp;
Landmark_contour.PyyInterp = PyyInterp;
Landmark_TGPoint.TGPoint = TGPoint;
Landmark_TGPoint.RefLinePoint11 = RefLinePoint11;
Landmark_TGPoint.RefLinePoint22 = RefLinePoint22;
Landmark_TGPoint.TGPoint3D = TGPoint3D;
Landmark_TGPoint.p1 = p1;
Landmark_TGPoint.p2 = p2;
end

function [TGPoint1,LinePoint1, LinePoint2,TGPoint3D,P1_new,P2_new] = getTG(directionRL,theseV,k,k_list,rotM,PD_pos,str_Side)

[LinePointList1,LinePointList2] = deal(zeros(length(k),2));

LinePointList1(1:end,:) =  theseV(k(1:end),:);
LinePointList2(1:end-1,:) =  theseV(k(2:end),:);
LinePointList2(end,:) =  theseV(k(1),:);
directionList = LinePointList2-LinePointList1;
%acos(dot(a,b)/(norm(a)*norm(a)))
projLenList = abs(dot(directionList,repmat(directionRL,length(directionList),1),2));
[len_proj,IndsLen] = sort(projLenList,'descend');
%%% angle
% directionListCan = directionList(IndsLen(2:4),:);
% u = directionListCan;
% v = repmat(directionRL,length(directionListCan),1);
% normU = sqrt(sum(directionListCan.*directionListCan,2));
% CosTheta = dot(u,v,2)./(normU*norm(directionRL));
% ThetaInDegrees = abs(acosd(CosTheta));
% [v,IndsDegree] = sort(ThetaInDegrees,'ascend');

if len_proj(2)-len_proj(3) < 4
    if str_Side == 'R'
        LinePoint2_1 = LinePointList1(IndsLen(2),:);
        LinePoint2_2 = LinePointList2(IndsLen(2)+1,:);
        center2 = LinePoint2_1+LinePoint2_2;
        LinePoint3_1 = LinePointList1(IndsLen(3),:);
        LinePoint3_2 = LinePointList2(IndsLen(3)+1,:);
        center3 = LinePoint3_1+LinePoint3_2;
        if center3(1)>center2(1)
            Inds = IndsLen(3);
        else
            Inds = IndsLen(2);
        end
    else
        LinePoint2_1 = LinePointList1(IndsLen(2),:);
        LinePoint2_2 = LinePointList2(IndsLen(2)+1,:);
        center2 = LinePoint2_1+LinePoint2_2;
        LinePoint3_1 = LinePointList1(IndsLen(3),:);
        LinePoint3_2 = LinePointList2(IndsLen(3)+1,:);
        center3 = LinePoint3_1+LinePoint3_2;
        if center3(1)<center2(1)
            Inds = IndsLen(3);
        else
            Inds = IndsLen(2);
        end
    end
else
    Inds = IndsLen(2);
end
LinePoint1 = LinePointList1(Inds,:);
LinePoint2 = LinePointList2(Inds+1,:);

% if k(Inds)>k(Inds+1)
% %     selectVertices = PointSet(convexInd(sortDiffInd):size(PointSet,1),:);
% %     selectVertices = [selectVertices PointSet(2:convexInd(sortDiffInd(1)))];
%     selectVertices = theseV(k(Inds+1):k(Inds),:);
% else
%     selectVertices = theseV(k(Inds):k(Inds+1),:);
% end
selectVertices = theseV(k_list{Inds},:);
errorSumList = [];
totalNum = size(selectVertices,1)-2;
P1List = zeros(totalNum,2);
P2List = zeros(totalNum,2);
for num = 2:size(selectVertices,1)-1
    line1Vertices = selectVertices(1:num,:);
    p1 = polyfit(line1Vertices(:,1),line1Vertices(:,2),1);
    P1List(num-1,:) = p1; 
    y1 = polyval(p1,line1Vertices(:,1));  
    line2Vertices = selectVertices(num+1:end,:);
    p2 = polyfit(line2Vertices(:,1),line2Vertices(:,2),1);
    P2List(num-1,:) = p2; 
    y2 = polyval(p2,line2Vertices(:,1));  
    diff1 = abs(y1-line1Vertices(:,2));
    diff2 = abs(y2-line2Vertices(:,2));
    errorSum = sum(diff1(:)) + sum(diff2(:));
    errorSumList = [errorSumList errorSum];
end
[minError,minErrorInd] = min(errorSumList(:));
P1_new = P1List(minErrorInd,:);
P2_new = P2List(minErrorInd,:);
TGPoint1 =  selectVertices(minErrorInd+1,:);
TGPoint3D = [PD_pos,TGPoint1(1),TGPoint1(2)];
TGPoint3D = rotM * TGPoint3D';
end

