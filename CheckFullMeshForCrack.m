%% purpose of this script is to identify teeth that are likely cracked. Problems with cracked teeth arise when
% ConfMax points are placed on both sides of a crack, even though their is
% truly only one feature. This script will delete one of the two.

%% preparation
path(pathdef);
addpath(path,genpath([pwd '/utils/']));
addpath(path,genpath([pwd '/PNAS/meshes']));

options.FeatureType = 'ConfMax';
options.NumDensityPnts = 100;
options.AngleIncrement = 0.05;
options.NumFeatureMatch = 4;
options.GaussMinMatch = 'on';

%% load mesh and prepare variables for use
% load('Lemuridae.mat'); G = Lemuridae{1}; %now have variable G
load('Lemuridae.mat'); G = Lemuridae{1};
NumConfMax = length(G.Aux.ConfMaxInds);
f = 0;

%% find distances between ConfMaxInds; check for distance less than 0.1
for j = 1:NumConfMax
    for k = 1:NumConfMax
        ConfMaxDist = sqrt(sum((G.V(:,G.Aux.ConfMaxInds(k)) - G.V(:,G.Aux.ConfMaxInds(j))).^2));
        if ConfMaxDist == 0 
            continue
        elseif ConfMaxDist < 0.1
            disp(ConfMaxDist);
            f = j;
            break
        else
            continue
        end
    end
end

%% delete 1 confmax point if Euclidean distance between any two confmax points is less than 0.1
if f ~= 0 
    disp('Tooth is likely cracked!');
    G.Aux.ConfMaxInds(f) = [];
    save('x23_fixed_for_crack_deliv.mat','G');

    NumPts = 400;
    % preallocations for speed
    i = zeros(1,10);
    Gs = cell(1,10);
    D = cell(1,10);
    subV = cell(1,10);
    subF = cell(1,10);
    rng('shuffle'); % want a different random starting point for every iteration
    for j = 1:10 
        i(j) = randi(G.nV);
        subV{j} = G.GeodesicFarthestPointSampling(NumPts,i(j));
        [subV{j},subF{j}] = G.PerformGeodesicDelauneyTriangulation(subV{j},[]); % Create new mesh for subsample
        Gs{j} = Mesh('VF',subV{j},subF{j});
        [Gs{j}.Aux.Area, Gs{j}.Aux.Center] = Gs{j}.Centralize('ScaleArea');
        [~,TriAreas] = Gs{j}.ComputeSurfaceArea;
        Gs{j}.Aux.VertArea = (TriAreas'*Gs{j}.F2V)/3;
        D{j} = Gs{j}.ComputeCPMS().ComputeUniformization(struct('method','LSCM','boundary_conditions','disc'));
        Gs{j}.Aux.UniformizationV = D{j}.V;
        [~,AG] = Gs{j}.ComputeSurfaceArea;
        [~,DG] = D{j}.ComputeSurfaceArea;
        Gs{j}.Aux.Conf = (AG./DG)'*D{j}.F2V/3;
    %%%     set the G.Aux.Conf to 0 on the boundary
        [Gs{j}.BV,Gs{j}.BE] = Gs{j}.FindBoundaries();
        Gs{j}.Aux.Conf(Gs{j}.BV) = 0;
    %%%     keep features on the original fine mesh for the downsampled mesh
        TREE = kdtree_build(Gs{j}.V');
        Gs{j}.Aux.ConfMaxInds = kdtree_nearest_neighbor(TREE, G.V(:,G.Aux.ConfMaxInds)');
        %% similarly, find nearest neighbors of
        %%%% G.Aux.ADMaxInds, G.Aux.GaussMaxInds, G.Aux.GaussMinInds
        Gs{j}.Aux.ADMaxInds = kdtree_nearest_neighbor(TREE, G.V(:,G.Aux.ADMaxInds)');
        Gs{j}.Aux.GaussMaxInds = kdtree_nearest_neighbor(TREE, G.V(:,G.Aux.GaussMaxInds)');
        Gs{j}.Aux.GaussMinInds = kdtree_nearest_neighbor(TREE, G.V(:,G.Aux.GaussMinInds)');
        %%% finding density points
        minds = [Gs{j}.Aux.GaussMaxInds;Gs{j}.Aux.GaussMinInds;Gs{j}.Aux.ConfMaxInds];
        minds = unique(minds);
        Gs{j}.Aux.DensityPnts = Gs{j}.GeodesicFarthestPointSampling(Gs{j}.nV,minds);
    end    
    save('x23_400_from_fixed_full_deliv.mat','Gs');
    %% compute cPdist matrix
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
    save('x23_400_from_fixed_full_deliv_Rslt.mat','cPDist');
    %% Compare the subsamples to the original mesh
    CPrslt = cell(1,11);
    CPDself = zeros(1,11); 
    CPrslt{1,1} = ComputeContinuousProcrustes4(G,G,options);
    CPDself(1,1) = CPrslt{1,1}.cPdist;
    for i = 1:10
        CPrslt{1,i+1} = ComputeContinuousProcrustes4(G,Gs{i},options);
        CPDself(1,i+1) = CPrslt{1,i+1}.cPdist;
    end
    save('x23_400_from_fixed_full_deliv_CPDself.mat','CPDself');
end

