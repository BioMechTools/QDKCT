function [ vn ] = vertexnormals( F,V )
%VERTEXNORMALS(F,V) Calculate vertex normals for a traingular mesh.
%   The vertex normal is defined as the mean normal of all adjacent faces
%--Input
%   F     -     Faces (M by 3 integer)
%   V     -     Vertices (N by 3)
%--Output
%   vn    -     Vertex normals. The ith value corresponds to the ith Vertex

a = V(F(:,1),:);
b = V(F(:,2),:);
c = V(F(:,3),:);
fn = cross((b-a),(c-a)); % un-normalized face normals

vn = zeros(size(V)); % preallocate for vertex normals
vn(F(:,1),:) = vn(F(:,1),:) + fn; % add the normals for the first point in each triangle
vn(F(:,2),:) = vn(F(:,2),:) + fn; % add the normals for the second point in each triangle
vn(F(:,3),:) = vn(F(:,3),:) + fn; % add the normals for the third point in each triangle

vn = vn./repmat(sqrt(sum(vn.^2,2)),[1,3]); % normalized vertex normals

end

