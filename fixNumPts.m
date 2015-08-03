%% purpose of this script is to fix submeshes that have fewer than the number of pts they should have w/o resampling all 10
%% Preparation
% clear vars;
% close all;
%% Preparation
% clear vars;
% close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));
% addpath(path,genpath([pwd '/PNAS/meshes']));

%%  Set Parameters for later
options.FeatureType = 'ConfMax';
options.NumDensityPnts = 100;
options.AnGleIncrement = 0.05;
options.NumFeatureMatch = 4;
options.GaussMinMatch = 'on';

%% load in mesh to check for insufficient number of vertices
NumPts = 600;
load('Saimiri.mat'); G = Saimiri{2};
load('aj01_600.mat'); %Sai

for j = 1:15
    if Sai{j}.nV < 600
        i = randi(G.nV); 
        % Get subsampled mesh
        subV = G.GeodesicFarthestPointSampling(NumPts,i);
        [subV,subF] = G.PerformGeodesicDelauneyTriangulation(subV,[]); % Create new mesh for subsample
        Sai{j} = Mesh('VF',subV,subF);
        % compute uniformization for subsampled mesh
        [Sai{j}.Aux.Area, Sai{j}.Aux.Center] = Sai{j}.Centralize('ScaleArea');
        [~,TriAreas] = Sai{j}.ComputeSurfaceArea;
        Sai{j}.Aux.VertArea = (TriAreas'*Sai{j}.F2V)/3;
        D = Sai{j}.ComputeCPMS().ComputeUniformization(struct('method','LSCM','boundary_conditions','disc'));
        Sai{j}.Aux.UniformizationV = D.V; 
    %%%     Estimate conformal factor using area distortion
        [~,AG] = Sai{j}.ComputeSurfaceArea;
        [~,DG] = D.ComputeSurfaceArea;
        Sai{j}.Aux.Conf = (AG./DG)'*D.F2V/3;
    %%%     set the G.Aux.Conf to 0 on the boundary
        [Sai{j}.BV,Sai{j}.BE] = FindBoundaries(Sai{j});
        Sai{j}.Aux.Conf(Sai{j}.BV) = 0;
    %%%     keep features on the original fine mesh for the downsampled mesh
        TREE = kdtree_build(Sai{j}.V');
        Sai{j}.Aux.ConfMaxInds = kdtree_nearest_neighbor(TREE, G.V(:,G.Aux.ConfMaxInds)');

        %% similarly, find nearest neighbors of
        %%%% G.Aux.ADMaxInds, G.Aux.GaussMaxInds, G.Aux.GaussMinInds
        Sai{j}.Aux.ADMaxInds = kdtree_nearest_neighbor(TREE, G.V(:,G.Aux.ADMaxInds)');
        Sai{j}.Aux.GaussMaxInds = kdtree_nearest_neighbor(TREE, G.V(:,G.Aux.GaussMaxInds)');
        Sai{j}.Aux.GaussMinInds = kdtree_nearest_neighbor(TREE, G.V(:,G.Aux.GaussMinInds)');

        %%% finding density points
        minds = [Sai{j}.Aux.GaussMaxInds;Sai{j}.Aux.GaussMinInds;Sai{j}.Aux.ConfMaxInds];
        minds = unique(minds);
        Sai{j}.Aux.DensityPnts = Sai{j}.GeodesicFarthestPointSampling(Sai{j}.nV,minds);
    end
end
save('aj01_600.mat','Sai');

%% delete isolated vertices
for j=1:10
    ConfMax = Sai{j}.V(:,Sai{j}.Aux.ConfMaxInds);
    GaussMax = Sai{j}.V(:,Sai{j}.Aux.GaussMaxInds);
    GaussMin = Sai{j}.V(:,Sai{j}.Aux.GaussMinInds);
    ADMax = Sai{j}.V(:,Sai{j}.Aux.ADMaxInds);
    
    dVInds = Sai{j}.DeleteIsolatedVertex();
    if ~isempty(dVInds) %% recompute uniformization
        [Sai{j}.Aux.Area,Sai{j}.Aux.Center] = Sai{j}.Centralize('ScaleArea');
        [Sai{j}.BV, Sai{j}.BE] = Sai{j}.FindBoundaries();
        ConfMax = (ConfMax-repmat(Sai{j}.Aux.Center,1,size(ConfMax,2)))/sqrt(Sai{j}.Aux.Area);
        GaussMax = (GaussMax-repmat(Sai{j}.Aux.Center,1,size(GaussMax,2)))/sqrt(Sai{j}.Aux.Area);
        GaussMin = (GaussMin-repmat(Sai{j}.Aux.Center,1,size(GaussMin,2)))/sqrt(Sai{j}.Aux.Area);
        ADMax = (ADMax-repmat(Sai{j}.Aux.Center,1,size(ADMax,2)))/sqrt(Sai{j}.Aux.Area);
        
        Sai{j}.ComputeMidEdgeUniformization(options);
        Sai{j}.Nf = Sai{j}.ComputeFaceNormals;
        Sai{j}.Nv = Sai{j}.F2V'*Sai{j}.Nf';
        Sai{j}.Nv = Sai{j}.Nv'*diag(1./sqrt(sum((Sai{j}.Nv').^2,1)));
        Sai{j}.Aux.LB = Sai{j}.ComputeCotanLaplacian;
        
        TREE = kdtree_build(Sai{j}.V');
        Sai{j}.Aux.ConfMaxInds = kdtree_nearest_neighbor(TREE, ConfMax');
        Sai{j}.Aux.GaussMaxInds = kdtree_nearest_neighbor(TREE, GaussMax');
        Sai{j}.Aux.GaussMinInds = kdtree_nearest_neighbor(TREE, GaussMin');
        Sai{j}.Aux.ADMaxInds = kdtree_nearest_neighbor(TREE, ADMax');
    end
end

save('aj01_600.mat','Sai');

%% cp distance
cPDist = zeros(10);
rslt = cell(10);
for i = 1:10
    for j = 1:10
        try
            rslt{i,j} = ComputeContinuousProcrustes4(Sai{i},Sai{j},options); 
            cPDist(i,j) = rslt{i,j}.cPdist;
        catch
            rslt{i,j} = ComputeContinuousProcrustes4(Sai{j},Sai{i},options); 
            cPDist(i,j) = rslt{i,j}.cPdist;
        end
    end
end

save('aj01_600Rslt.mat','cPDist');

%% Compare the subsamples to the original mesh
CPrslt = cell(1,11);
CPDself = zeros(1,11); 
CPrslt{1,1} = ComputeContinuousProcrustes4(G,G,options);
CPDself(1,1) = CPrslt{1,1}.cPdist;
for i = 1:10
    CPrslt{1,i+1} = ComputeContinuousProcrustes4(G,Sai{i},options);
    CPDself(1,i+1) = CPrslt{1,i+1}.cPdist;
end
save('aj01_600_CPDself.mat','CPDself');