% function [body, edge, node] = tetrahedmesh(Nxyz)
C = imgC1{1}(193:320, 193:320, :) - 1050;
C = C(1:4:end, 1:4:end, 1:4:end);
[Ny,Nx,Nz] = size(C);

% [Ny,Nx,Nz] = tac([64 64 64]);

% I know the mesh is in size (Ny,Nx,Nz), NOT (Nx, Ny, Nz)!
[X, Y, Z] = meshgrid(1:Nx, 1:Ny, 1:Nz);
XYZ = [X(:) Y(:) Z(:)];

% C = sqrt((X-32).^2 + (Y-32).^2 + (Z-32).^2) - 10;

[bodynodeindex, nodeshift, edgeorder, p3] = meshconsts(Nx, Ny, Nz);

Nnode = Nx*Ny*Nz;
Nbody = (Nx-1)*(Ny-1)*(Nz-1)*6;
Nupnode = (Nx*2-1)*(Ny*2-1)*(Nz*2-1);
Nedge = Nupnode-Nnode;

nodeindex = reshape(1:Nnode, Ny, Nx, Nz);
upnodeindex = reshape(1:Nupnode, Ny*2-1, Nx*2-1, Nz*2-1);
upindex2edgenode = upnodeindex.*0 + 1;
upindex2edgenode(1:2:end,1:2:end,1:2:end) = 0;
upindexsedge = upindex2edgenode==1;
upindex2edgenode(upindexsedge) = 1:Nedge;
upindex2edgenode(1:2:end,1:2:end,1:2:end) = nodeindex;

nodeindex2body = nodeindex(1:end-1,1:end-1,1:end-1);
nodeindex2up = upnodeindex(1:2:end,1:2:end,1:2:end);

body = struct();
body.nodeindex = repelem(nodeindex2body(:), 6, 4) + repmat(nodeshift(bodynodeindex), Nbody/6, 1);
body.edgeindex = upindex2edgenode(sum(reshape(nodeindex2up(body.nodeindex(:, edgeorder(:))),Nbody, 6, 2), 3)./2);

edge.nodeindex = zeros(Nedge, 2);

for ii = 2:8
    edge_ii = upindex2edgenode(p3(ii,1)+1 : 2 : end, p3(ii,2)+1 : 2 : end, p3(ii,3)+1 : 2 : end);
    node0_ii = upindex2edgenode(1 : 2 : end-p3(ii,1)*2, 1 : 2 : end-p3(ii,2)*2, 1 : 2 : end-p3(ii,3)*2);
    node1_ii = upindex2edgenode(1+p3(ii,1)*2 : 2 : end, 1+p3(ii,2)*2 : 2 : end, 1+p3(ii,3)*2 : 2 : end);
    edge.nodeindex(edge_ii, :) = [node0_ii(:) node1_ii(:)];
end

% end

% % check nodes & edges
% for jj = 1:Nupnode
%     index_jj = upindex2edgenode(jj);
%     if ~upindexsedge(jj)
%         % node
%         XYZ_jj = XYZ(index_jj, :);
%     else
%         % edge
%         XYZ_jj = mean(XYZ(edge.nodeindex(index_jj, :)', :), 1);
%     end
%     XYZindex = XYZ_jj(:,2)*2-1 + (XYZ_jj(:,1)*2-2)*(Ny*2-1) + (XYZ_jj(:,3)*2-2)*(Nx*2-1)*(Ny*2-1);
%     if XYZindex~=jj
%         error('check node/edge %d!', jj);
%     end
% end
% % check body
% for jj = 1:Nbody
%     if any(reshape(edge.nodeindex(body.edgeindex(jj, :)', :), 1, []) ~= body.nodeindex(jj,edgeorder(:)))
%         error('check body %d!', jj);
%     end
% end



% C test
% C12 = C(edge.nodeindex);
% alpha = fillmC12(:,1)./(C12(:,1)-C12(:,2));
% s = alpha>=0 | alpha<=1;
s = C>=0;
s_body = ~all(s(body.nodeindex), 2) & any(s(body.nodeindex), 2);
s_edge = s(edge.nodeindex(:,1)) ~= s(edge.nodeindex(:,2));
Nsedge = sum(s_edge);
Nsbody = sum(s_body);
edge.sedgeindex = zeros(Nedge, 1);
edge.sedgeindex(s_edge) = 1:Nsedge;

alpha = C(edge.nodeindex(s_edge, 1))./(C(edge.nodeindex(s_edge, 1))-C(edge.nodeindex(s_edge, 2)));
R = XYZ(edge.nodeindex(s_edge, 1), :).*(1-alpha) + XYZ(edge.nodeindex(s_edge, 2), :).*alpha;

sbody_edge = body.edgeindex(s_body, :)';
sedge_sbody = s_edge(sbody_edge')';
s_t1 = sum(sedge_sbody, 1)==3;
s_t2 = sum(sedge_sbody, 1)==4;
Ntri = sum(s_t1) + sum(s_t2)*2;
T1 = reshape(sbody_edge(sedge_sbody & s_t1), 3, [])';
T2 = reshape(sbody_edge(sedge_sbody & s_t2), 4, [])';
T2 = [T2(:, 1:3); T2(:, 2:4)];
T = edge.sedgeindex([T1; T2]);
Nt = length(T);

figure;
h=trisurf(T, R(:,1),R(:,2),R(:,3), R(:,1).*0+1, 'FaceAlpha', 0.8); axis equal
shading interp
lightangle(-45,30)
% h.FaceLighting = 'gouraud';
h.FaceLighting = 'flat';
% '''

% hsc3 = scatter3(X(:)+rand(Nnode,1)-0.5, Y(:)+rand(Nnode,1)-0.5, Z(:)+rand(Nnode,1)-0.5, 5, repmat(1-C2,1,3)); axis equal;

% V = T
V = reshape(R(T(:), :), Nt, 3, 3);
V = V - mean(V, 2);
% V1 = R(T(:,2),:) - R(T(:,1),:);
% V2 = R(T(:,3),:) - R(T(:,2),:);
% V3 = R(T(:,1),:) - R(T(:,3),:);
R2 = R.*0;
W = zeros(Nsedge, 1);
a = 1.0;
for ii = 1:Nt
    R2(T(ii,:),:) = R2(T(ii,:),:) - squeeze(V(ii, :, :)).*a;
    W(T(ii,:)) = W(T(ii,:)) + 1;
end
R2 = R2./W + R;

figure;
h2=trisurf(T, R2(:,1),R2(:,2),R2(:,3), R2(:,1).*0+1, 'FaceAlpha', 0.8); axis equal
shading interp
lightangle(-45,30)
% h2.FaceLighting = 'gouraud';
h2.FaceLighting = 'flat';

function [bodynodeindex, nodeshift, edgeorder, p3] = meshconsts(Nx, Ny, Nz)
bodynodeindex = [1 2 4 8;
                 1 2 6 8;
                 1 5 6 8;
                 1 5 7 8;
                 1 3 7 8;
                 1 3 4 8];

nodeshift = [0  Ny  1  Ny+1  Ny*Nx  Ny*(Nx+1)  Ny*Nx+1  Ny*(Nx+1)+1];

edgeorder = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

s8 = (0:7)';
p3 = [mod(s8,2) floor(mod(s8, 4)./2)  floor(s8./4) ];

end