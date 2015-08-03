% Preparation
clearvars; 
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));
addpath(path,genpath([pwd '../MeshLabfiles/fixed_for_crack']));

%  Set Parameters for later
options.FeatureType = 'ConfMax';
options.NumDensityPnts = 150;
options.AngleIncrement = 0.05;
options.NumFeatureMatch = 4;
options.GaussMinMatch = 'on';

% Prepare meshes
% Load tooth that has been fixed for crack
G = Mesh('off','../MeshLabfiles/fixed_for_crack/aj05_ffc.off');

% Convert to .mat files
[G.Aux.Area,G.Aux.Center] = G.Centralize('ScaleArea');
G.ComputeMidEdgeUniformization(options); 
G.Nf = G.ComputeFaceNormals;
G.Nv = G.F2V'*G.Nf';
G.Nv = G.Nv'*diag(1./sqrt(sum((G.Nv').^2,1)));
G.Aux.LB = G.ComputeCotanLaplacian; 

save('aj05_ffc.mat', 'G'); %these can be used forever

% delete isolated vertices
ConfMax = G.V(:,G.Aux.ConfMaxInds);
GaussMax = G.V(:,G.Aux.GaussMaxInds);
GaussMin = G.V(:,G.Aux.GaussMinInds);
ADMax = G.V(:,G.Aux.ADMaxInds);

dVInds = G.DeleteIsolatedVertex();
if ~isempty(dVInds) %% recompute uniformization
    [G.Aux.Area,G.Aux.Center] = G.Centralize('ScaleArea');
    [G.BV, G.BE] = G.FindBoundaries();
    ConfMax = (ConfMax-repmat(G.Aux.Center,1,size(ConfMax,2)))/sqrt(G.Aux.Area);
    GaussMax = (GaussMax-repmat(G.Aux.Center,1,size(GaussMax,2)))/sqrt(G.Aux.Area);
    GaussMin = (GaussMin-repmat(G.Aux.Center,1,size(GaussMin,2)))/sqrt(G.Aux.Area);
    ADMax = (ADMax-repmat(G.Aux.Center,1,size(ADMax,2)))/sqrt(G.Aux.Area);

    G.ComputeMidEdgeUniformization(options);
    G.Nf = G.ComputeFaceNormals;
    G.Nv = G.F2V'*G.Nf';
    G.Nv = G.Nv'*diag(1./sqrt(sum((G.Nv').^2,1)));
    G.Aux.LB = G.ComputeCotanLaplacian;

    TREE = kdtree_build(G.V');
    G.Aux.ConfMaxInds = kdtree_nearest_neighbor(TREE, ConfMax');
    G.Aux.GaussMaxInds = kdtree_nearest_neighbor(TREE, GaussMax');
    G.Aux.GaussMinInds = kdtree_nearest_neighbor(TREE, GaussMin');
    G.Aux.ADMaxInds = kdtree_nearest_neighbor(TREE, ADMax');
end
save('aj05_ffc.mat','G');
