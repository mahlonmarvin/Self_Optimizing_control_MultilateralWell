


%% Define geometry and rock parameters
% Construct a Cartesian grid of size 20-by-20-by-5 cells, where each cell
% has dimension 1-by-1-by-1. Set the permeability $K$ to be homogeneous,
% isotropic and equal 100 mD and the porosity to be equal to 0.3.

nx = 25; ny = 20; nz = 5;
G = processGRDECL(simpleGrdecl  ([nx,ny,nz], 0.3, 'physDims',[225, 22.5,1]));
G = computeGeometry(G);


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
%
fluid = initSimpleFluid('mu' , [   1,  10] .* centi*poise     , ...
                        'rho', [1000, 700] .* kilogram/meter^3, ...
                        'n'  , [   2,   2]);

%%
%
s=linspace(0,1,20)'; kr=fluid.relperm(s);

%% Initialize and construct linear system
S  = computeMimeticIP(G, rock, 'Verbose', true);

%% Introduce wells
% <html>
% We will include four wells, all rate-controlled vertical wells. Wells are described
% using a Peacemann model, giving an extra set of equations that need to be
% assembled, see <a href="../../1ph/html/simpleWellExample.html#3"> "Using
% Peacemann well models"</a> for more details.
% </html>

%u1=controls.well(1).values(1); u2=controls.well(2).values(1);

%u1 = 0.3848/day;           % initial points using BM/OC
% initial points using material balance
totTime = 720;
nInj = 1;
totVol = sum(poreVolume(G, rock));
u1 = (1/nInj)*totVol/totTime;
%u2 = u1;

%u1=controls.well(1).values(1);
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
%
rSol = initState(G, W, 0, [0, 1]);

%% Solve initial pressure in reservoir
% Solve linear system construced from S and W to obtain solution for flow
% and pressure in the reservoir and the wells.
gravity off
rSol = solveIncompFlow(rSol, G, S, fluid, 'wells', W);


%% Transport loop
% 
T      = 720;
dT     = T/100;
pv     = poreVolume(G,rock);

%%
%
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

%W(1).val=controls.well(1).values(1); W(2).val=controls.well(2).values(1); 
%W(3).val=-W(2).val; W(4).val=-W(1).val;     
rISol = initState(G, W, 0, [0, 1]);   
gravity off;
rISol = solveIncompFlow(rISol, G, S, fluid, 'wells', W);


%% Start the main loop
t  = 0; plotNo = 1; hi = 'Implicit: '; he = 'Explicit: ';
po2socC1 = []; po3socC1 = []; po4socC1 = []; pw2socC1 = []; Yw1 = []; pw3socC1 = []; pw4socC1 = [];
usoc1C1 = []; NPV1=0; NPVsocC1=[]; l = 0;

%u1 = controls.well(1).values(1); u2 = controls.well(2).values(1);
%u1 = 0.3848/day;
%u2 = 0.3848/day;

while t < T,
   
   l = l + 1; 
   
   
   rISol = implicitTransport(rISol, G, dT, rock, fluid, 'wells', W);
   W(1).val = u1; W(2).val = -u1;
   W(3).val = -u1; W(4).val = -u1;
   
   mSol1=bsxfun(@rdivide, fluid.relperm(rISol.s(W(1).cells,:)), mu);
   mSol2=bsxfun(@rdivide, fluid.relperm(rISol.s(W(2).cells,:)), mu);
   mSol3=bsxfun(@rdivide, fluid.relperm(rISol.s(W(3).cells,:)), mu);
   mSol4=bsxfun(@rdivide, fluid.relperm(rISol.s(W(4).cells,:)), mu);
   
   OilPr2 = opr(mSol2, -rISol.wellSol(2).flux);           
   OilPr3 = opr(mSol3, -rISol.wellSol(3).flux);
   OilPr4 = opr(mSol4, -rISol.wellSol(4).flux);
   WatPr2 = wpr(mSol2, -rISol.wellSol(2).flux);            
   WatPr3 = wpr(mSol3, -rISol.wellSol(3).flux);
   WatPr4 = wpr(mSol4, -rISol.wellSol(4).flux);
   WatIn1 = wir(mSol1, rISol.wellSol(1).flux);           
   
   
   po2socC1 = [po2socC1; OilPr2];                 %#ok
   po3socC1 = [po3socC1; OilPr3];
   po4socC1 = [po4socC1; OilPr4];
   
   pw2socC1 = [pw2socC1; WatPr2];                 %#ok
   pw3socC1 = [pw3socC1; WatPr3];
   pw4socC1 = [pw4socC1; WatPr4];
   
   usoc1C1 = [usoc1C1; W(1).val];                 %#ok
   
   
   NPV1=NPV1+dT*(-(WatIn1)*wi-(WatPr2+WatPr3+WatPr4)*wp+(OilPr2+OilPr3+OilPr4)*op)/(1+b)^t;
   NPVsocC1=[NPVsocC1; NPV1];
   
   % Update solution of pressure equation.
   rISol = solveIncompFlow(rISol, G, S, fluid, 'wells', W);
   
   %Define next value of the CVs, using the feedback control law if there
   %are all ready two available past histories
   %
   if l > 2 % number of history check
   
   ucv1 = -(Theta(7))^-1*(Theta(1)*po4socC1(l)+Theta(2)*pw4socC1(l)+Theta(3)*po4socC1(l-1)+Theta(4)*pw4socC1(l-1)+Theta(5)*po4socC1(l-2)+Theta(6)*pw4socC1(l-2)); 
   ucv2 = -(Theta(14))^-1*(Theta(8)*po3socC1(l)+Theta(9)*pw3socC1(l)+Theta(10)*po3socC1(l-1)+Theta(11)*pw3socC1(l-1)+Theta(12)*po3socC1(l-2)+Theta(13)*pw3socC1(l-2)); 
   ucv3 = -(Theta(21))^-1*(Theta(15)*po2socC1(l)+Theta(16)*pw2socC1(l)+Theta(17)*po2socC1(l-1)+Theta(18)*pw2socC1(l-1)+Theta(19)*po2socC1(l-2)+Theta(20)*pw2socC1(l-2));   
   
   ucv1 = u1;
   ucv2 = u1;
   ucv3 = u1;
   
   W(1).val = u1; W(2).val = -u1;
   W(3).val = -u1; W(4).val = -u1;
   
   gradient1 = Theta(1)*po4socC1(l)+Theta(2)*pw4socC1(l)+Theta(3)*po4socC1(l-1)+Theta(4)*pw4socC1(l-1)+Theta(5)*po4socC1(l-2)+Theta(6)*pw4socC1(l-2)+ Theta(7)*ucv1
   gradient2 = Theta(8)*po3socC1(l)+Theta(9)*pw3socC1(l)+Theta(10)*po3socC1(l-1)+Theta(11)*pw3socC1(l-1)+Theta(12)*po3socC1(l-2)+Theta(13)*pw3socC1(l-2)+ Theta(14)*ucv2 
   gradient3 = Theta(15)*po2socC1(l)+Theta(16)*pw2socC1(l)+Theta(17)*po2socC1(l-1)+Theta(18)*pw2socC1(l-1)+Theta(19)*po2socC1(l-2)+Theta(20)*pw2socC1(l-2)+ Theta(21)*ucv3 
   
   Un=[ucv1 ucv2 ucv3]

   end
   %}

   t = t + dT;
end


%%
% Plot the water and oil rates
n = size(po2socC1,1);

figure
subplot(2,3,1),
   plot(1:n,3600*24*po2socC1(:,1))
   title('Oil production in well P');
subplot(2,3,2),
   plot(1:n,3600*24*po3socC1(:,1))
   title('Oil production in well P');
subplot(2,3,3),
   plot(1:n,3600*24*po4socC1(:,1))
   title('Oil production in well P1');
subplot(2,3,4),
   plot(1:n,3600*24*pw2socC1(:,1))
   title('Water production in well P');
subplot(2,3,5),
   plot(1:n,3600*24*pw3socC1(:,1))
   title('Water production in well P');
subplot(2,3,6),
   plot(1:n,3600*24*pw4socC1(:,1))
   title('Water production in well P1');

   figure,
   plot(1:n,3600*24*usoc1C1(:,1))
   title('Water injection in well I');
   
%Plot the NPV

figure
   plot(1:n,NPVsocC1(:,1))
   title('Net Present Value');
   
% total production
Tvec = (dT : dT : T)';
oilT1 = cumtrapz ([0; Tvec], [zeros([1, 1]); po2socC1]);
oilT2 = cumtrapz ([0; Tvec], [zeros([1, 1]); po3socC1]);
oilT3 = cumtrapz ([0; Tvec], [zeros([1, 1]); po4socC1]);
wat1 = cumtrapz ([0; Tvec], [zeros([1, 1]); pw2socC1]);
wat2 = cumtrapz ([0; Tvec], [zeros([1, 1]); pw3socC1]);
wat3 = cumtrapz ([0; Tvec], [zeros([1, 1]); pw4socC1]);
%%
displayEndOfDemoMessage(mfilename)
