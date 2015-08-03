%% Computing the cPdist matrix between teeth subsamples
% first generate the subsamples (in this case fix sample size)
%% Preparation
% clear vars;
% close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));
addpath(path,genpath([pwd '/PNAS/meshes']));

%%  Set Parameters for later
options.FeatureType = 'ConfMax';
options.NumDensityPnts = 100;
options.AngleIncrement = 0.05;
options.NumFeatureMatch = 4;
options.GaussMinMatch = 'on';

%% load mesh to be subsampled
load('Lemuridae_fixed_for_crack.mat');
% load('x23_fixed_for_crack.mat'); Lemuridae{1} = G;
% load('x08_fixed_for_crack.mat'); Lemuridae{6} = G;
% load('x09_fixed_for_crack.mat'); Lemuridae{9} = G;

NumPts = 1000;
ConfMaxInds = cell(1,10);
GaussMaxInds = cell(1,10);
GaussMinInds = cell(1,10);
ADMaxInds = cell(1,10);
startingpts = cell(1,10);

for j = 1:10
    ConfMaxInds{j} = Lemuridae{j}.Aux.ConfMaxInds;
    GaussMaxInds{j} = Lemuridae{j}.Aux.GaussMaxInds;
%     GaussMinInds{j} = Lemuridae{j}.Aux.GaussMinInds;
%     ADMaxInds{j} = Lemuridae{j}.Aux.ADMaxInds;
    startingpts{j} = [ConfMaxInds{j}' GaussMaxInds{j}']; % want to start fastmarching with major features of tooth
end

% preallocations for speed
Gs = cell(1,10);
D = cell(1,10);
subV = cell(1,10);
subF = cell(1,10);

for j = 1:10 
    subV{j} = Lemuridae{j}.GeodesicFarthestPointSampling(NumPts,startingpts{j});
    [subV{j},subF{j}] = Lemuridae{j}.PerformGeodesicDelauneyTriangulation(subV{j},[]); % Create new mesh for subsample
    Gs{j} = Mesh('VF',subV{j},subF{j});
    [Gs{j}.Aux.Area, Gs{j}.Aux.Center] = Gs{j}.Centralize('ScaleArea');
    [~,TriAreas] = Gs{j}.ComputeSurfaceArea;
    Gs{j}.Aux.VertArea = (TriAreas'*Gs{j}.F2V)/3;
    D{j} = Gs{j}.ComputeCPMS().ComputeUniformization(struct('method','LSCM','boundary_conditions','disc'));
    Gs{j}.Aux.UniformizationV = D{j}.V;
%%%     Estimate conformal factor using area distortion
    [~,AG] = Gs{j}.ComputeSurfaceArea;
    [~,DG] = D{j}.ComputeSurfaceArea;
    Gs{j}.Aux.Conf = (AG./DG)'*D{j}.F2V/3;
%%%     set the G.Aux.Conf to 0 on the boundary
    [Gs{j}.BV,Gs{j}.BE] = Gs{j}.FindBoundaries();
    Gs{j}.Aux.Conf(Gs{j}.BV) = 0;
%%%     keep features on the original fine mesh for the downsampled mesh
    TREE = kdtree_build(Gs{j}.V');
    Gs{j}.Aux.ConfMaxInds = kdtree_nearest_neighbor(TREE, Lemuridae{j}.V(:,Lemuridae{j}.Aux.ConfMaxInds)');
%%%     similarly, find nearest neighbors of
    %%%% G.Aux.ADMaxInds, G.Aux.GaussMaxInds, G.Aux.GaussMinInds
    Gs{j}.Aux.ADMaxInds = kdtree_nearest_neighbor(TREE, Lemuridae{j}.V(:,Lemuridae{j}.Aux.ADMaxInds)');
    Gs{j}.Aux.GaussMaxInds = kdtree_nearest_neighbor(TREE, Lemuridae{j}.V(:,Lemuridae{j}.Aux.GaussMaxInds)');
    Gs{j}.Aux.GaussMinInds = kdtree_nearest_neighbor(TREE, Lemuridae{j}.V(:,Lemuridae{j}.Aux.GaussMinInds)');
%%%     finding density points
    minds = [Gs{j}.Aux.GaussMaxInds;Gs{j}.Aux.GaussMinInds;Gs{j}.Aux.ConfMaxInds];
    minds = unique(minds);
    Gs{j}.Aux.DensityPnts = Gs{j}.GeodesicFarthestPointSampling(Gs{j}.nV,minds);
end
save('LemTeeth_ConfGauss_1000.mat','Gs');

%% Delete isolated vertex
for j=1:10
    ConfMax = Gs{j}.V(:,Gs{j}.Aux.ConfMaxInds);
    GaussMax = Gs{j}.V(:,Gs{j}.Aux.GaussMaxInds);
    GaussMin = Gs{j}.V(:,Gs{j}.Aux.GaussMinInds);
    ADMax = Gs{j}.V(:,Gs{j}.Aux.ADMaxInds);
    
    dVInds = Gs{j}.DeleteIsolatedVertex();
    if ~isempty(dVInds) %% recompute uniformization
        [Gs{j}.Aux.Area,Gs{j}.Aux.Center] = Gs{j}.Centralize('ScaleArea');
        [Gs{j}.BV, Gs{j}.BE] = FindBoundaries(Gs{j});
        ConfMax = (ConfMax-repmat(Gs{j}.Aux.Center,1,size(ConfMax,2)))/sqrt(Gs{j}.Aux.Area);
        GaussMax = (GaussMax-repmat(Gs{j}.Aux.Center,1,size(GaussMax,2)))/sqrt(Gs{j}.Aux.Area);
        GaussMin = (GaussMin-repmat(Gs{j}.Aux.Center,1,size(GaussMin,2)))/sqrt(Gs{j}.Aux.Area);
        ADMax = (ADMax-repmat(Gs{j}.Aux.Center,1,size(ADMax,2)))/sqrt(Gs{j}.Aux.Area);
        
        Gs{j}.ComputeMidEdgeUniformization(options);
        Gs{j}.Nf = Gs{j}.ComputeFaceNormals;
        Gs{j}.Nv = Gs{j}.F2V'*Gs{j}.Nf';
        Gs{j}.Nv = Gs{j}.Nv'*diag(1./sqrt(sum((Gs{j}.Nv').^2,1)));
        Gs{j}.Aux.LB = Gs{j}.ComputeCotanLaplacian;
        
        TREE = kdtree_build(Gs{j}.V');
        Gs{j}.Aux.ConfMaxInds = kdtree_nearest_neighbor(TREE, ConfMax');
        Gs{j}.Aux.GaussMaxInds = kdtree_nearest_neighbor(TREE, GaussMax');
        Gs{j}.Aux.GaussMinInds = kdtree_nearest_neighbor(TREE, GaussMin');
        Gs{j}.Aux.ADMaxInds = kdtree_nearest_neighbor(TREE, ADMax');
    end
end
save('LemTeeth_ConfGauss_1000.mat','Gs');

%% Create cP distance matrix
cPDist = zeros(10);
rslt = cell(10);
for i = 1:10
    for j = 1:10
        try
            rslt{i,j} = ComputeContinuousProcrustes4(Gs{i},Gs{j},options); 
            cPDist(i,j) = rslt{i,j}.cPdist;
        catch
            rslt{i,j} = ComputeContinuousProcrustes4(Gs{j},Gs{i},options); 
            cPDist(i,j) = rslt{i,j}.cPdist;
        end
    end
end

save('LemTeeth_ConfGauss_1000_Rslt.mat','cPDist');

%% Compute distance from each subsampled mesh to each full mesh
CPrslt = cell(10);
CPDself = zeros(10); 
for i = 1:10
    for j = 1:10
        CPrslt{i,j} = ComputeContinuousProcrustes4(Gs{i},Lemuridae{j},options); %each row represents a single subsample, each column a single full mesh
        CPDself(i,j) = CPrslt{i,j}.cPdist;
    end
end
save('LemTeeth_ConfGauss_1000_CPDself.mat','CPDself');