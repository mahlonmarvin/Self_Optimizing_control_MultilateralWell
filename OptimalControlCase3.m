

% whether or not to show output
verbose = false;
verboseLevel = 2;

% Define model ------------------------------------------------------------
nx = 20; ny = 20; nz = 5;
G = cartGrid([nx, ny, nz]);
G = computeGeometry(G);
rock.perm  = repmat(100*milli*darcy, [G.cells.num, 1]);
rock.poro  = repmat(0.3, [G.cells.num, 1]);

fluid = initSimpleFluid('mu' , [   1,  10] .* centi*poise     , ...
                        'rho', [1000, 700] .* kilogram/meter^3, ...
                        'n'  , [   1.5,   1.5]);
                    
fluid  = adjointFluidFields(fluid);                    

% Wells and initial rates -------------------------------------------------

totVol = sum(poreVolume(G, rock));
totTime = 730*day;

% Two injection wells
nInj = 2;
cellsWell1 =  1 : nx*ny : nx*ny*nz;
radius     = .1;
W = addWell([], G, rock, cellsWell1, 'Type', 'rate', ...
            'Val', 1.0/day(), 'Radius', radius, 'comp_i', [1,0], 'sign', 1, 'name', 'I');
disp('Well #1: '); display(W(1));

cellsWell2 =  nx : nx*ny : nx*ny*nz;
W = addWell(W, G, rock, cellsWell2, 'Type', 'rate', ...
            'Val', 1.0/day(), 'Radius', radius, 'comp_i', [1,0], 'sign', 1, 'name', 'I2');
disp('Well #2: '); display(W(2));

% Two production wells
nProd = 2;
cellsWell3 =  nx*ny : nx*ny : nx*ny*nz;
W = addWell(W, G, rock, cellsWell3, 'Type', 'rate', ...
            'Val', -1.0/day(), 'Radius', radius, 'comp_i', [1,1], 'sign', -1,'name', 'P');
disp('Well #3: '); display(W(3));

cellsWell4 =  nx*(ny-1)+1 : nx*ny : nx*ny*nz;
W = addWell(W, G, rock, cellsWell4, 'Type', 'rate', ...
            'Val', -1.0/day(), 'Radius', radius, 'comp_i', [1,1], 'sign', -1, 'name', 'P2');
disp('Well #4: '); display(W(4));

% System components -------------------------------------------------------
S = computeMimeticIP(G, rock, 'Type', 'comp_hybrid', 'Verbose', verbose, ...
                     'InnerProduct', 'ip_tpf');
W = assembleWellSystem(G, W, 'Type', 'comp_hybrid');

% Initialize --------------------------------------------------------------
state = initResSol(G, 0.0);
state.wellSol = initWellSol(W, 0);

% Objective function ------------------------------------------------------
objectiveFunction = str2func('simpleNPV');

% Initialize schedule and controls ----------------------------------------
numSteps = 730;
schedule = initSchedule(W, 'NumSteps', numSteps, 'TotalTime', ...
                        totTime, 'Verbose', verbose);

% box constraints for each well [min rate, max rate]       
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
                                    
