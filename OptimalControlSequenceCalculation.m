

% whether or not to show output
verbose = false;
verboseLevel = 2;

% Define model ------------------------------------------------------------
nx = 25; ny = 20; nz = 5;
G = cartGrid([nx, ny, nz]);
G = computeGeometry(G);
%rock.perm  = repmat(300*milli*darcy, [G.cells.num, 1]);
% Set up permeability based on K-indices
%rock.poro  = repmat(0.3, [G.cells.num, 1]);
[I, J, K] = gridLogicalIndices(G);

px       = 200*milli*darcy*ones(G.cells.num,1);
px(K==2) = 400*milli*darcy;
px(K==3) = 600*milli*darcy;
px(K==4) = 800*milli*darcy;
px(K==5) = 300*milli*darcy;

% Introduce anisotropy by setting K_x = 5*K_z.
perm = [px, px, 0.1*px];
rock = makeRock(G, perm, 0.3);

fluid = initSimpleFluid('mu' , [   1,  10] .* centi*poise     , ...
                        'rho', [1000, 700] .* kilogram/meter^3, ...
                        'n'  , [   2,   2]);
                    
fluid  = adjointFluidFields(fluid);                    

% Wells and initial rates -------------------------------------------------

totVol = sum(poreVolume(G, rock));
totTime = 720*day;

% one injection wells
c1 = 4*(nx*ny) + (2:15)'; c2= 1002:nx:1275; c3 = 3*(nx*ny) + (2:15)';
%nInj = 1;
%cellsWell1 =  1 : nx*ny : nx*ny*nz;
%radius     = .1;
W = struct([]);
W = addWell(W, G, rock, (500:nx*ny:2500),'Name', 'I', 'radius', 5*inch,...
'Type', 'rate', 'Val', 1.0/day(), 'comp_i', [1, 0], 'Sign', 1);

%cellsWell2 =  nx : nx*ny : nx*ny*nz;
%W = addWell(W, G, rock, cellsWell2, 'Type', 'rate', ...
 %           'Val', 1.0/day(), 'Radius', radius, 'comp_i', [1,0], 'sign', 1, 'name', 'I2');
disp('Well #1: '); display(W(1));

% Two production wells
%nProd = 2;
W = addWell(W, G, rock, c1, 'Name', 'P', 'radius', 5*inch, ...
'Type', 'bhp', 'Val', 1000*psia, 'comp_i', [0, 1], 'Sign', -1);
disp('Well #2: '); display(W(1));

W = addWell(W, G, rock, c2, 'Name', 'P', 'radius', 5*inch, ...
'Type', 'bhp', 'Val', 1000*psia, 'comp_i', [0, 1],'Sign', -1);
disp('Well #2: '); display(W(1));

W = addWell(W, G, rock, c3, 'Name', 'P', 'radius', 5*inch, ...
'Type', 'bhp', 'Val', 1000*psia, 'comp_i', [0, 1], 'Sign', -1);
disp('Well #3: '); display(W(1));

 %       disp('Well #4: '); display(W(4));
plotGrid(G, 'FaceColor', 'none'); view(3);
[ht, htxt, hs] = plotWell(G, W, 'radius', 0.5, 'height', 2, 'Color', 'k');
set(htxt, 'FontSize', 16);

figure,            
plotCellData(G, log(rock.perm(:,1)));
plotWell(G, W), view([1 1 1])
%% System components -------------------------------------------------------
S = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, ...
                     'InnerProduct', 'ip_tpf');
W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

% Initialize --------------------------------------------------------------
state = initResSol(G, 0.0);
state.wellSol = initWellSol(W, 0);

% Objective function ------------------------------------------------------
objectiveFunction = str2func('simpleNPV');

% Initialize schedule and controls ----------------------------------------
numSteps = 100;
schedule = initSchedule(W, 'NumSteps', numSteps, 'TotalTime', ...
                        totTime, 'Verbose', verbose);

% box constraints for each well [min rate, max rate]       
nInj = 1;
nProd = 3;
box = [repmat([0 inf], nInj, 1); repmat([-inf 0], nProd, 1)];
controls = initControls(schedule, 'ControllableWells', (1:numel(W)), ...
                                  'MinMax', box, ...
                                  'Verbose', verbose, ...
                                  'NumControlSteps', numSteps, ...
                                  'LinEqConst', {ones(1, numel(W)), 0}');                               

% Run optimization --------------------------------------------------------                              
[simRes, schedule, controls, out] = optimizeObjective(G, S, W, rock, ...
                                        fluid, state, schedule, ...
                                        controls, objectiveFunction, ...
                                        'gradTol',       1e-3, ...
                                        'objChangeTol',  5e-4, ...
                                        'VerboseLevel', verboseLevel);
                                    
