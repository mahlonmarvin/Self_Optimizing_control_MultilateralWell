
%% Define geometry and rock parameters
% Construct a Cartesian grid of size 20-by-20-by-5 cells, where each cell
% has dimension 1-by-1-by-1. Set the permeability $K$ to be homogeneous,
% isotropic and equal 100 mD and the porosity to be equal to 0.3.
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

figure,            
plotCellData(G, log(rock.perm(:,1))); view(3);

%% Define the two-phase fluid model
% The <matlab:help('initSimpleFluid') two-phase fluid model> has default values:
%
% * densities: [rho_w, rho_o] = [1000 700] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 10] cP.
fluid = initSimpleFluid('mu' , [   1,  10] .* centi*poise     , ...
                        'rho', [1000, 700] .* kilogram/meter^3, ...
                        'n'  , [   1.5,   1.5]);

%%
% The fluid model represented by the <matlab:help('fluid') fluid structure>
% is the two-phase incompressible counterpart to the fluid model of the
% Black Oil <matlab:help('pvt') 'pvt'> function.
%
s=linspace(0,1,20)'; kr=fluid.relperm(s);

%% Initialize and construct linear system
% Initialize solution structure with reservoir pressure equal 0 and initial
% water saturation equal 0.0 (reservoir is filled with oil). Compute the
% mimetic inner product from input grid and rock properties.
S  = computeMimeticIP(G, rock, 'Verbose', true);

%% Introduce wells
% <html>
% We will include four wells, all rate-controlled vertical wells. Wells are described
% using a Peacemann model, giving an extra set of equations that need to be
% assembled, see <a href="../../1ph/html/simpleWellExample.html#3"> "Using
% Peacemann well models"</a> for more details.
% </html>

u1=controls.well(1).values(1);
c1 = 4*(nx*ny) + (2:15)'; c2= 1002:nx:1275; c3 = 3*(nx*ny) + (2:15)';
%nInj = 1;
%cellsWell1 =  1 : nx*ny : nx*ny*nz;
%radius     = .1;
W = struct([]);
W = addWell(W, G, rock, (500:nx*ny:2500),'Name', 'I1', 'radius', 5*inch,...
'Type', 'rate', 'Val', 1.0/day(), 'Val', u1, 'comp_i', [1, 0], 'Sign', 1);
disp('Well #1: '); display(W(1));
%cellsWell1 =  1 : nx*ny : nx*ny*nz;
%radius     = .1;
%W = addWell([], G, rock, cellsWell1, 'Type', 'rate', ...
 %           'Val', u1, 'Radius', radius, 'comp_i', [1,0], 'sign', 1, 'name', 'I');
%disp('Well #1: '); display(W(1));

W = addWell(W, G, rock, c1, 'Name', 'P1', 'radius', 5*inch, ...
'Type', 'bhp', 'Val',  1000*psia, 'Val', -u1, 'comp_i', [0, 1], 'Sign', -1);
%disp('Well #2: '); display(W(2));

W = addWell(W, G, rock, c2, 'Name', 'P2', 'radius', 5*inch, ...
'Type', 'bhp', 'Val', 1000*psia, 'Val', -u1, 'comp_i', [0, 1],'Sign', -1);
%disp('Well #3: '); display(W(3));

W = addWell(W, G, rock, c3, 'Name', 'P3', 'radius', 5*inch, ...
'Type', 'bhp', 'Val', 1000*psia, 'Val', -u1, 'comp_i', [0, 1], 'Sign', -1);
%disp('Well #4: '); display(W(4));

% To check if the wells are placed as we wanted them, we plot them
plotGrid(G, 'FaceColor', 'none'); view(3);
[ht, htxt, hs] = plotWell(G, W, 'radius', 0.1, 'height', 2);
set(htxt, 'FontSize', 16);



%%
% Once the wells are added, we can generate the components of the linear
% system corresponding to the two wells and initialize the solution
% structure (with correct bhp)
%
rSol = initState(G, W, 0, [0, 1]);

%% Solve initial pressure in reservoir
% Solve linear system construced from S and W to obtain solution for flow
% and pressure in the reservoir and the wells.
gravity off
rSol = solveIncompFlow(rSol, G, S, fluid, 'wells', W);

%%
% Report initial state of reservoir
figure,
   plotCellData(G, convertTo(rSol.pressure(1:G.cells.num), barsa));
    view(3)

figure,
   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
   plotCellData(G, accumarray(cellNo, ...
      abs(convertTo(faceFlux2cellFlux(G, rSol.flux), meter^3/day))));
   title('Initial flux intensity'), view(3)

%% Transport loop
% We solve the two-phase system using a sequential splitting in which the
% pressure and fluxes are computed by solving the flow equation and then
% held fixed as the saturation is advanced according to the transport
% equation. This procedure is repeated for a given number of time steps
% (here we use 15 equally spaced time steps). The error introduced by this
% splitting of flow and transport can be reduced by iterating each time
% step until e.g., the residual is below a certain user-prescribed
% threshold (this is not done herein).
T      = 720*day();
dT     = T/100;
pv     = poreVolume(G,rock);

%%
% The transport equation will be solved by the single-point upstream method
% with either explicit or implicit time discretizations. Both schemes may
% use internal time steps to obtain a stable discretization. To represent
% the two solutions, we create new solution objects to be used by the
% solver with implicit transport step.
rISol = rSol;


mu=fluid.properties();


opr = @(m,q)sum(m(:,2).*q./sum(m,2));
wpr = @(m,q)sum(m(:,1).*q./sum(m,2));
wir = @(m,q)sum(m(:,1).*q./sum(m,2));

%Define prices in $/m^3 and discount factor
op=100/0.1589873; %100$/bbl
wp=10/0.1589873;
wi=10/0.1589873;
b=0;

%% Start the main loop
t  = 0; plotNo = 1; hi = 'Implicit: '; he = 'Explicit: ';
po2ocC1 = []; po3ocC1 = []; po4ocC1 = []; pw2ocC1 = []; Yw1 = []; pw3ocC1 = [];pw4ocC1 = [];
uw1 = []; NPV1oc=0; NPVocC1=[]; l=1;
uoc1C1=[];

while t < T,
   l
   u1=controls.well(1).values(l);
   W(1).val = u1;
   W(2).val = -u1; W(3).val = -u1; W(4).val = -u1;  
   
   rISol = implicitTransport(rISol, G, dT, rock, fluid, 'wells', W);


   % Update solution of pressure equation.
   rISol = solveIncompFlow(rISol, G, S, fluid, 'wells', W);
   
   mSol1=bsxfun(@rdivide, fluid.relperm(rISol.s(W(1).cells,:)), mu);
   mSol2=bsxfun(@rdivide, fluid.relperm(rISol.s(W(2).cells,:)), mu);
   mSol3=bsxfun(@rdivide, fluid.relperm(rISol.s(W(3).cells,:)), mu);
   mSol4=bsxfun(@rdivide, fluid.relperm(rISol.s(W(4).cells,:)), mu);

   % Measure water saturation in production cells in saturation

   OilPr2 = opr(mSol2, -rISol.wellSol(2).flux);            %#ok
   OilPr3 = opr(mSol3, -rISol.wellSol(3).flux);
   OilPr4 = opr(mSol4, -rISol.wellSol(4).flux);
   
   po2ocC1 = [po2ocC1; OilPr2];  
   po3ocC1 = [po3ocC1; OilPr3];  %#ok
   po4ocC1 = [po4ocC1; OilPr4];
   
  
   % Measure flux in production cells at every step (implicit solution)

   WatPr2 = wpr(mSol2, -rISol.wellSol(2).flux);
   WatPr3 = wpr(mSol3, -rISol.wellSol(3).flux); %#ok
   WatPr4 = wpr(mSol4, -rISol.wellSol(4).flux);
   
   pw2ocC1 = [pw2ocC1; WatPr2]; 
   pw3ocC1 = [pw3ocC1; WatPr3]; %#ok
   pw4ocC1 = [pw4ocC1; WatPr4];
   
   
   
   %Measure flux in inyection cell at every step 
   
   WatIn1 = wir(mSol1, rISol.wellSol(1).flux);            %#ok
   
   
   uw1 = [uw1; WatIn1];                 %#ok
   
   
   uoc1C1 = [uoc1C1; W(1).val];                 %#ok
   
   
   %Calculate the NVP value for each step
   
   NPV1oc=NPV1oc+dT*(-(WatIn1)*wi-(WatPr2+WatPr3+WatPr4)*wp+(OilPr2+OilPr3+OilPr4)*op)/(1+b)^t;
   NPVocC1=[NPVocC1; NPV1oc];

   % Increase time and continue 
   t = t + dT;
   l=l+1;

end


%%
% Plot the water and oil rates
n = size(po2ocC1,1);

figure
subplot(2,3,1),
   plot(1:n,3600*24*po2ocC1(:,1))
   title('Oil production in well P1');
subplot(2,3,2),
   plot(1:n,3600*24*po3ocC1(:,1))
   title('Oil production in well P2');
subplot(2,3,3),
   plot(1:n,3600*24*po4ocC1(:,1))
   title('Oil production in well P3');
subplot(2,3,4),
   plot(1:n,3600*24*pw2ocC1(:,1))
   title('Water production in well P1');
subplot(2,3,5),
   plot(1:n,3600*24*pw3ocC1(:,1))
   title('Water production in well P2');
subplot(2,3,6),
   plot(1:n,3600*24*pw4ocC1(:,1))
   title('Water production in well P3');

  
   figure,
   plot(1:n,3600*24*uoc1C1(:,1))
   title('Water injection in well I');
   
%Plot the NPV

figure
   plot(1:n,NPVocC1(:,1))
   title('Net Present Value');
   
   % total production
Tvec = (dT : dT : T)';
oilT1oc = cumtrapz ([0; Tvec], [zeros([1, 1]); po2ocC1]);
oilT2oc = cumtrapz ([0; Tvec], [zeros([1, 1]); po3ocC1]);
oilT3oc = cumtrapz ([0; Tvec], [zeros([1, 1]); po4ocC1]);
wat1oc = cumtrapz ([0; Tvec], [zeros([1, 1]); pw2ocC1]);
wat2oc = cumtrapz ([0; Tvec], [zeros([1, 1]); pw3ocC1]);
wat3oc = cumtrapz ([0; Tvec], [zeros([1, 1]); pw4ocC1]);
%%
uoc1C1 = 3600*24*uoc1C1;
po2ocC1 = 3600*24*po2ocC1;
po3ocC1 = 3600*24*po3ocC1;
po4ocC1 = 3600*24*po4ocC1;

pw2ocC1 = 3600*24*pw2ocC1;
pw3ocC1 = 3600*24*pw3ocC1;
pw4ocC1 = 3600*24*pw4ocC1;
displayEndOfDemoMessage(mfilename)

