function [ TTPoint3D ] = getTT(  R_axes,tibia_coordinate,node,face)%, manual_TT)
%getTT is aimed to get the TT point 
%   Detailed explanation goes here
AP = tibia_coordinate(:,3);
%%% cut the vertices
rotM = eye(3) * R_axes; 
nodeT = (rotM\(node'))';    % rotate the node to the provided coordinate
AP = rotM\AP;   % rotate AP to the provided coordinate
[faceT1, VcsLSelection] = cutMeshByPercent(face,nodeT,[60 100],1,1);    % remove the bottom of the shaft for the irregular shape of the tibia
% manual_TT = (rotM * TT')' ;
theseV = VcsLSelection(:,2:3);
[k] = convhull(theseV(:,1),theseV(:,2));
disttemp = (( theseV(k(2:end),1)-theseV(k(1:end-1),1)).^2 + (theseV(k(2:end),2)-theseV(k(1:end-1),2)).^2).^0.5;
[~,k_inds] = sort(disttemp,'descend');
cen_candidate = zeros(4,2);
cen_pos = zeros(4,1);
for num = 1:4
    cen_candidate(num,:) = (theseV(k(k_inds(num)),:)+theseV(k(k_inds(num)+1),:))/2;
    cen_pos(num) = AP(2:3,1)'*cen_candidate(num,:)';
end
[~,cen_inds] = sort(cen_pos,'descend');
%%% calculate the pixel more anterior
%%% point 1
if theseV(k(k_inds(cen_inds(1))),2)>theseV(k(k_inds(cen_inds(1))+1),2)
    k_point1 =k_inds(cen_inds(1));
else
    k_point1 = k_inds(cen_inds(1))+1;
end
%%% point 2
if theseV(k(k_inds(cen_inds(2))),2)>theseV(k(k_inds(cen_inds(2))+1),2)
    k_point2 = k_inds(cen_inds(2));
else
    k_point2 = k_inds(cen_inds(2))+1;
end
if k_point1>k_point2
    temp = k_point1;
    k_point1 = k_point2;
    k_point2 = temp;
end
n_band1 = length(k_point1:k_point2);
n_band2 = length([1:k_point1 k_point2:length(k)]);
if n_band1<n_band2
    target_range = k(k_point1:k_point2);
else
    target_range = [k(k_point2:length(k)-1)' k(1:k_point1)'];
end

x = theseV(target_range,1)';
if theseV(target_range(1),1) > theseV(target_range(end),1)
    xx = theseV(target_range(end),1):0.5:theseV(target_range(1),1);
else
    xx = theseV(target_range(1),1):0.5:theseV(target_range(end),1);
end
PxxInterp = spline(x,theseV(target_range,2)',xx);

n_targetTT = round(length(PxxInterp)/2);
% n_k_TT = PxxInterp(n_targetTT);
disttemp = (( theseV(target_range,1)-xx(n_targetTT)).^2 + (theseV(target_range,2)-PxxInterp(n_targetTT)).^2).^0.5;
[~, n_k_TT] = min(disttemp);
% figure;
% plot(VcsLSelection(:,2),VcsLSelection(:,3),'b.','markersize',3);
% hold on;
% plot(theseV(k,1),theseV(k,2),'k-');
% plot(manual_TT(2),manual_TT(3),'r*','markersize',3);
% plot(theseV(k(k_point1),1),theseV(k(k_point1),2),'g*');
% plot(theseV(k(k_point2),1),theseV(k(k_point2),2),'g*');
% plot(theseV(target_range(n_k_TT),1),theseV(target_range(n_k_TT),2),'k*');
% % plot(xx(n_targetTT),PxxInterp(n_targetTT),'k*');
% hold off;
% nodeT = (rotM*(node'))';
TTPoint3D = (rotM*(VcsLSelection(target_range(n_k_TT),:)'))';

end