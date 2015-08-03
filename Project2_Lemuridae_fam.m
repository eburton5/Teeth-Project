%% Goal of this script: calculate the cP distance between submeshes (size 1000 pts) from family Lemuridae
% Hapalemur genus: x23, v13, v09 (replaced with v12 as of 6/25/15)
% Eulemur genus: t12, t14, x08
% Lemur genus: x02, x07, x09, x17
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

%% Prepare meshes
% Lemuridae = cell(1,10);
% % Load Hapalemur genus
% Lemuridae{1} = load('./samples/PNAS/x23.mat'); Lemuridae{1} = Lemuridae{1}.G;
% Lemuridae{2} = load('./samples/PNAS/V13.mat'); Lemuridae{2} = Lemuridae{2}.G;
% Lemuridae{3} = load('v12_full.mat'); Lemuridae{3} = Lemuridae{3}.G;
% % Lemuridae{3} = load('./samples/PNAS/v09.mat'); Lemuridae{3} = Lemuridae{3}.G;
% % Load Eulemer genus
% Lemuridae{4} = load('./samples/PNAS/T12.mat'); Lemuridae{4} = Lemuridae{4}.G;
% Lemuridae{5} = load('./samples/PNAS/T14.mat'); Lemuridae{5} = Lemuridae{5}.G;
% Lemuridae{6} = load('./samples/PNAS/x08.mat'); Lemuridae{6} = Lemuridae{6}.G;
% %Load Lemur genus
% Lemuridae{7} = load('./samples/PNAS/x02.mat'); Lemuridae{7} = Lemuridae{7}.G;
% Lemuridae{8} = load('./samples/PNAS/x07.mat'); Lemuridae{8} = Lemuridae{8}.G;
% Lemuridae{9} = load('./samples/PNAS/x09.mat'); Lemuridae{9} = Lemuridae{9}.G;
% Lemuridae{10} = load('./samples/PNAS/x17.mat'); Lemuridae{10} = Lemuridae{10}.G;
% 
% % delete isolated vertices
% for j=1:10
%     ConfMax = Lemuridae{j}.V(:,Lemuridae{j}.Aux.ConfMaxInds);
%     GaussMax = Lemuridae{j}.V(:,Lemuridae{j}.Aux.GaussMaxInds);
%     GaussMin = Lemuridae{j}.V(:,Lemuridae{j}.Aux.GaussMinInds);
%     ADMax = Lemuridae{j}.V(:,Lemuridae{j}.Aux.ADMaxInds);
%     
%     dVInds = Lemuridae{j}.DeleteIsolatedVertex();
%     if ~isempty(dVInds) %% recompute uniformization
%         [Lemuridae{j}.Aux.Area,Lemuridae{j}.Aux.Center] = Lemuridae{j}.Centralize('ScaleArea');
%         [Lemuridae{j}.BV, Lemuridae{j}.BE] = Lemuridae{j}.FindBoundaries();
%         ConfMax = (ConfMax-repmat(Lemuridae{j}.Aux.Center,1,size(ConfMax,2)))/sqrt(Lemuridae{j}.Aux.Area);
%         GaussMax = (GaussMax-repmat(Lemuridae{j}.Aux.Center,1,size(GaussMax,2)))/sqrt(Lemuridae{j}.Aux.Area);
%         GaussMin = (GaussMin-repmat(Lemuridae{j}.Aux.Center,1,size(GaussMin,2)))/sqrt(Lemuridae{j}.Aux.Area);
%         ADMax = (ADMax-repmat(Lemuridae{j}.Aux.Center,1,size(ADMax,2)))/sqrt(Lemuridae{j}.Aux.Area);
%         
%         Lemuridae{j}.ComputeMidEdgeUniformization(options);
%         Lemuridae{j}.Nf = Lemuridae{j}.ComputeFaceNormals;
%         Lemuridae{j}.Nv = Lemuridae{j}.F2V'*Lemuridae{j}.Nf';
%         Lemuridae{j}.Nv = Lemuridae{j}.Nv'*diag(1./sqrt(sum((Lemuridae{j}.Nv').^2,1)));
%         Lemuridae{j}.Aux.LB = Lemuridae{j}.ComputeCotanLaplacian;
%         
%         TREE = kdtree_build(Lemuridae{j}.V');
%         Lemuridae{j}.Aux.ConfMaxInds = kdtree_nearest_neighbor(TREE, ConfMax');
%         Lemuridae{j}.Aux.GaussMaxInds = kdtree_nearest_neighbor(TREE, GaussMax');
%         Lemuridae{j}.Aux.GaussMinInds = kdtree_nearest_neighbor(TREE, GaussMin');
%         Lemuridae{j}.Aux.ADMaxInds = kdtree_nearest_neighbor(TREE, ADMax');
%     end
% end
% 
% save('Lemuridae_deliv.mat','Lemuridae');
load('Lemuridae_fixed_for_crack.mat');
NumPts = 1000;
rng('shuffle');

% preallocations for speed
i = zeros(1,10);
DLem = cell(1,10);
subV = cell(1,10);
subF = cell(1,10);
Lem = cell(1,10);

for j = 1:10 
    i(j) = randi(Lemuridae{j}.nV); % generate random starting point, will be different every time
    % get subsampled mesh
    subV{j} = Lemuridae{j}.GeodesicFarthestPointSampling(NumPts,i(j));
    [subV{j},subF{j}] = PerformGeodesicDelauneyTriangulation(Lemuridae{j},subV{j},[]); % Create new mesh for subsample
    Lem{j} = Mesh('VF',subV{j},subF{j});
    % compute uniformization for subsampled mesh
    [Lem{j}.Aux.Area, Lem{j}.Aux.Center] = Lem{j}.Centralize('ScaleArea');
    [~,TriAreas] = Lem{j}.ComputeSurfaceArea;
    Lem{j}.Aux.VertArea = (TriAreas'*Lem{j}.F2V)/3;
    DLem{j} = Lem{j}.ComputeCPMS().ComputeUniformization(struct('method','LSCM','boundary_conditions','disc'));
%     Lem{j}.ComputeMidEdgeUniformization(options); % beginning of midedge uniformization addition
%     Lem{j}.Nf = Lem{j}.ComputeFaceNormals;
%     Lem{j}.Nv = Lem{j}.F2V'*Lem{j}.Nf';
%     Lem{j}.Nv = Lem{j}.Nv'*diag(1./sqrt(sum((Lem{j}.Nv').^2,1)));
%     Lem{j}.Aux.LB = Lem{j}.ComputeCotanLaplacian; % end of midedge uniformization addition
    Lem{j}.Aux.UniformizationV = DLem{j}.V;
%%%     Estimate conformal factor using area distortion
    [~,AG] = Lem{j}.ComputeSurfaceArea;
    [~,DG] = DLem{j}.ComputeSurfaceArea;
    Lem{j}.Aux.Conf = (AG./DG)'*DLem{j}.F2V/3;
%%%     set the G.Aux.Conf to 0 on the boundary
    [Lem{j}.BV,Lem{j}.BE] = Lem{j}.FindBoundaries();
    Lem{j}.Aux.Conf(Lem{j}.BV) = 0;

%%%     keep features on the original fine mesh for the downsampled mesh
    TREE = kdtree_build(Lem{j}.V');
    Lem{j}.Aux.ConfMaxInds = kdtree_nearest_neighbor(TREE, Lemuridae{j}.V(:,Lemuridae{j}.Aux.ConfMaxInds)');
    %% similarly, find nearest neighbors of
    %%%% G.Aux.ADMaxInds, G.Aux.GaussMaxInds, G.Aux.GaussMinInds
    Lem{j}.Aux.ADMaxInds = kdtree_nearest_neighbor(TREE, Lemuridae{j}.V(:,Lemuridae{j}.Aux.ADMaxInds)');
    Lem{j}.Aux.GaussMaxInds = kdtree_nearest_neighbor(TREE, Lemuridae{j}.V(:,Lemuridae{j}.Aux.GaussMaxInds)');
    Lem{j}.Aux.GaussMinInds = kdtree_nearest_neighbor(TREE, Lemuridae{j}.V(:,Lemuridae{j}.Aux.GaussMinInds)');

    %%% finding density points
    minds = [Lem{j}.Aux.GaussMaxInds;Lem{j}.Aux.GaussMinInds;Lem{j}.Aux.ConfMaxInds];
    minds = unique(minds);
    Lem{j}.Aux.DensityPnts = Lem{j}.GeodesicFarthestPointSampling(Lem{j}.nV,minds);
end
save('Lem_fixed_for_crack_1000.mat','Lem');

%% delete isolated vertices
for j=1:10
    ConfMax = Lem{j}.V(:,Lem{j}.Aux.ConfMaxInds);
    GaussMax = Lem{j}.V(:,Lem{j}.Aux.GaussMaxInds);
    GaussMin = Lem{j}.V(:,Lem{j}.Aux.GaussMinInds);
    ADMax = Lem{j}.V(:,Lem{j}.Aux.ADMaxInds);
    
    dVInds = Lem{j}.DeleteIsolatedVertex();
    if ~isempty(dVInds) %% recompute uniformization
        [Lem{j}.Aux.Area,Lem{j}.Aux.Center] = Lem{j}.Centralize('ScaleArea');
        [Lem{j}.BV, Lem{j}.BE] = Lem{j}.FindBoundaries();
        ConfMax = (ConfMax-repmat(Lem{j}.Aux.Center,1,size(ConfMax,2)))/sqrt(Lem{j}.Aux.Area);
        GaussMax = (GaussMax-repmat(Lem{j}.Aux.Center,1,size(GaussMax,2)))/sqrt(Lem{j}.Aux.Area);
        GaussMin = (GaussMin-repmat(Lem{j}.Aux.Center,1,size(GaussMin,2)))/sqrt(Lem{j}.Aux.Area);
        ADMax = (ADMax-repmat(Lem{j}.Aux.Center,1,size(ADMax,2)))/sqrt(Lem{j}.Aux.Area);
        
        Lem{j}.ComputeMidEdgeUniformization(options);
        Lem{j}.Nf = Lem{j}.ComputeFaceNormals;
        Lem{j}.Nv = Lem{j}.F2V'*Lem{j}.Nf';
        Lem{j}.Nv = Lem{j}.Nv'*diag(1./sqrt(sum((Lem{j}.Nv').^2,1)));
        Lem{j}.Aux.LB = Lem{j}.ComputeCotanLaplacian;
        
        TREE = kdtree_build(Lem{j}.V');
        Lem{j}.Aux.ConfMaxInds = kdtree_nearest_neighbor(TREE, ConfMax');
        Lem{j}.Aux.GaussMaxInds = kdtree_nearest_neighbor(TREE, GaussMax');
        Lem{j}.Aux.GaussMinInds = kdtree_nearest_neighbor(TREE, GaussMin');
        Lem{j}.Aux.ADMaxInds = kdtree_nearest_neighbor(TREE, ADMax');
    end
end

save('Lem_fixed_for_crack_1000.mat','Lem');

%% load all 100 subsamples from all 10 teeth
% G = cell(1,100);
% cPDist_Lem_all = zeros(100);
% load('Lem_all_100_1000.mat'); load('Lem_all_100_1000_Rslt.mat'); % generates variables G and cPDist_Lem_all
% 
% cPDist_Lem_all(1:10,:) = zeros(10,100); cPDist_Lem_all(:,1:10) = zeros(100,10);
% cPDist_Lem_all(51:60,:) = zeros(10,100); cPDist_Lem_all(:,51:60) = zeros(100,10);
% cPDist_Lem_all(81:90,:) = zeros(10,100); cPDist_Lem_all(:,81:90) = zeros(100,10);
% 
% load('x23_1000_from_fixed_full.mat'); G(1,1:10) = Gs;
% load('x23_1000_from_fixed_full_Rslt.mat'); cPDist_Lem_all(1:10,1:10) = cPDist;
% load('x23_1000.mat'); G(1,1:10) = Gs;
% load('x23_1000Rslt.mat'); cPDist_Lem_all(1:10,1:10) = cPDist;
% load('v13_1000.mat'); G(1,11:20) = Gs;
% load('v13_1000Rslt.mat'); cPDist_Lem_all(11:20,11:20) = cPDist;
% load('v12_1000.mat'); G(1,21:30) = Gs;
% load('v12_1000Rslt.mat'); cPDist_Lem_all(21:30,21:30) = cPDist;
% load('t12_1000.mat'); G(1,31:40) = Gs;
% load('t12_1000Rslt.mat'); cPDist_Lem_all(31:40,31:40) = cPDist;
% load('t14_1000.mat'); G(1,41:50) = Gs;
% load('t14_1000Rslt.mat'); cPDist_Lem_all(41:50,41:50) = cPDist;
% load('x08_1000_from_fixed_full.mat'); G(1,51:60) = Gs;
% load('x08_1000_from_fixed_full_Rslt.mat'); cPDist_Lem_all(51:60,51:60) = cPDist;
% load('x08_1000.mat'); G(1,51:60) = Gs;
% load('x08_1000Rslt.mat'); cPDist_Lem_all(51:60,51:60) = cPDist;
% load('x02_1000.mat'); G(1,61:70) = Gs;
% load('x02_1000Rslt.mat'); cPDist_Lem_all(61:70,61:70) = cPDist;
% load('x07_1000.mat'); G(1,71:80) = Gs;
% load('x09_1000_from_fixed_full.mat'); G(1,81:90) = Gs;
% load('x09_1000_from_fixed_full_Rslt.mat'); cPDist_Lem_all(81:90,81:90) = cPDist;
% load('x09_1000.mat'); G(1,81:90) = Gs;
% load('x09_1000Rslt.mat'); cPDist_Lem_all(81:90,81:90) = cPDist;
% load('x17_1000.mat'); G(1,91:100) = Gs;
% load('x17_1000Rslt.mat'); cPDist_Lem_all(91:100,91:100) = cPDist;
% save('Lem_all_100_1000_fixed_for_crack.mat','G');

%% Create cP distance matrix
cPDist = zeros(10);
rslt = cell(10);
for i = 1:10
    for j = 1:10
%         if cPDist_Lem_all(i,j) == 0
            try
                rslt{i,j} = ComputeContinuousProcrustes4(Lem{i},Lem{j},options); 
                cPDist(i,j) = rslt{i,j}.cPdist;
            catch
                rslt{i,j} = ComputeContinuousProcrustes4(Lem{j},Lem{i},options); 
                cPDist(i,j) = rslt{i,j}.cPdist;
            end
%         end
    end
end

%% save data
save('Lem_fixed_for_crack_1000_Rslt.mat','cPDist');

%% Compute distance from each subsampled mesh to each full mesh
% load('Lem_all_100_1000_fixed_for_crack.mat');
% load('Lemuridae_fixed_for_crack.mat');
% CPrslt = cell(1,10);
% CPDself = zeros(10); 
% for i = 1:10
%     for j = 1:10
%         CPrslt{i,j} = ComputeContinuousProcrustes4(G{i},Lemuridae{j},options); %each row represents a single subsample, each column a single full mesh
%         CPDself(i,j) = CPrslt{i,j}.cPdist;
%     end
% end
% save('Lem_all_100_1000_fixed_for_crack_CPDself.mat','CPDself');