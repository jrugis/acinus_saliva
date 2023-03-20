
% Read in the mesh file from John Rugis, and modify it a bit to suit my own
% code. Change notation, and do all the precalculation of triangle areas,
% apical and basal triangles and nodes, etc.

clear all
close all
clc

load('mod_basal1data_smoothed_mesh.mat')
cell_no = 7;
p = p{cell_no};
tets = tets{cell_no};
tetvol = tets_volume{cell_no};
surftrilist = triangles{cell_no};
dist_to_apical = dist_ap_p{cell_no};        % the distance of NODES to the apical membrane
dist_to_basal_tets = dist_ba{cell_no};      % the distance of TETS to the BASAL membrane (not the basolateral)
np = size(p,1);
ntets = size(tets,1);


trisurf(surftrilist,p(:,1),p(:,2),p(:,3),'FaceColor','cyan','FaceAlpha',0.1)
hold on
duma = find(dist_to_apical(surftrilist(:,1))<0.5);
apicaltrilist = surftrilist(duma,:);
dumb = find(dist_to_apical(surftrilist(:,1))>0.5);
basaltrilist =  surftrilist(dumb,:);
trisurf(apicaltrilist,p(:,1),p(:,2),p(:,3),'FaceColor','red')
trisurf(basaltrilist,p(:,1),p(:,2),p(:,3),'FaceColor','blue')         
hold off

% Calculate the total apical and basal surface areas, which we use in the
% main program. Also precalculate the area of each surface triangle, and
% the total volume of the cell.

for i=1:size(apicaltrilist,1)
    pts = surftrilist(apicaltrilist(i),1:3);
    P1 = p(pts(1),:); P2 = p(pts(2),:); P3 = p(pts(3),:);  
    U1 = P2-P1;
    U2 = P3-P1;
    TT = [1 1 1; U1; U2];
    apicaltriarea(i) = 0.5*abs(det(TT));
end
apicalarea=sum(apicaltriarea);

for i=1:size(basaltrilist,1)
    pts = surftrilist(basaltrilist(i),1:3);
    P1 = p(pts(1),:); P2 = p(pts(2),:); P3 = p(pts(3),:);  
    U1 = P2-P1;
    U2 = P3-P1;
    TT = [1 1 1; U1; U2];
    basaltriarea(i) = 0.5*abs(det(TT));
end
basalarea=sum(basaltriarea);


% Finally, for each node precalculate the (approximate) shortest distance
% to the basal surface. Note this is NOT the distance to the basolateral
% surface, only to the basal membrane.

dist_to_basal = zeros(np,1);
for i=1:ntets
dist_to_basal(tets(i,:))  = dist_to_basal_tets(i);
end 

volume = sum(tetvol);
 
% To avoid confusion I don't just save everything. I only save the
% information related to the particular cell under consideration. This
% avoids the subsequent code having access to the variables for all the
% cells, which can create problems if you're not careful.

save('real_cell_7_mesh.mat','cell_no','p','tets','tetvol','surftrilist','apicaltrilist','basaltrilist','dist_to_apical','dist_to_basal',...
     'apicaltriarea','basaltriarea','np','ntets','apicalarea','basalarea','volume')


