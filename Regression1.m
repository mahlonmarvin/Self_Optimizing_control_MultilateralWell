ln=2; nT=50;
dJ=[];
dU1=[];
TotNPV=NPVTraj(100, :);


while ln<=nT,
    
    dJ=[dJ, TotNPV(ln)-TotNPV(ln-1)];
    dU1=[dU1, U1Traj(:,ln)-U1Traj(:,ln-1)];

    
    ln=ln+1;
end
    dJ=dJ';
    dU1=dU1';
   
    
