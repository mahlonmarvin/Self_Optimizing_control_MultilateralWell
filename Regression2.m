ln=2; nT=50;
X=zeros(49,7); X2=zeros(49,7); X3=zeros(49,7); 
Y=dJ;
nS=100;
Yo1=OilTraj2'; Yo2=OilTraj3'; Yo3=OilTraj4';
Yw1=WProTraj2'; Yw2=WProTraj3'; Yw3=WProTraj4';
U1=U1Traj';
Theta=[]; Theta2=[];


while ln<=nT,
    sn=3; 
        
    while sn<=nS,
        
        X(ln-1,1)= X(ln-1,1)+ Yo1(ln-1, sn)*dU1(ln-1,sn);
        X(ln-1,2)= X(ln-1,2)+ Yw1(ln-1, sn)*dU1(ln-1,sn);
        X(ln-1,3)= X(ln-1,3)+ Yo1(ln-1, sn-1)*dU1(ln-1,sn);
        X(ln-1,4)= X(ln-1,4)+ Yw1(ln-1, sn-1)*dU1(ln-1,sn);
        X(ln-1,5)= X(ln-1,5)+ Yo1(ln-1, sn-2)*dU1(ln-1,sn);
        X(ln-1,6)= X(ln-1,6)+ Yw1(ln-1, sn-2)*dU1(ln-1,sn);
        X(ln-1,7)= X(ln-1,7)+ U1(ln-1, sn)*dU1(ln-1,sn);

       
        X2(ln-1,1)= X2(ln-1,1)+ Yo2(ln-1, sn)*dU1(ln-1,sn);
        X2(ln-1,2)= X2(ln-1,2)+ Yw2(ln-1, sn)*dU1(ln-1,sn);
        X2(ln-1,3)= X2(ln-1,3)+ Yo2(ln-1, sn-1)*dU1(ln-1,sn);
        X2(ln-1,4)= X2(ln-1,4)+ Yw2(ln-1, sn-1)*dU1(ln-1,sn);
        X2(ln-1,5)= X2(ln-1,5)+ Yo2(ln-1, sn-2)*dU1(ln-1,sn);
        X2(ln-1,6)= X2(ln-1,6)+ Yw2(ln-1, sn-2)*dU1(ln-1,sn);
        X2(ln-1,7)= X2(ln-1,7)+ U1(ln-1, sn)*dU1(ln-1,sn);
        
        
        X3(ln-1,1)= X3(ln-1,1)+ Yo3(ln-1, sn)*dU1(ln-1,sn);
        X3(ln-1,2)= X3(ln-1,2)+ Yw3(ln-1, sn)*dU1(ln-1,sn);
        X3(ln-1,3)= X3(ln-1,3)+ Yo3(ln-1, sn-1)*dU1(ln-1,sn);
        X3(ln-1,4)= X3(ln-1,4)+ Yw3(ln-1, sn-1)*dU1(ln-1,sn);
        X3(ln-1,5)= X3(ln-1,5)+ Yo3(ln-1, sn-2)*dU1(ln-1,sn);
        X3(ln-1,6)= X3(ln-1,6)+ Yw3(ln-1, sn-2)*dU1(ln-1,sn);
        X3(ln-1,7)= X3(ln-1,7)+ U1(ln-1, sn)*dU1(ln-1,sn);

        
        sn=sn+1;
        
    end
    
    ln=ln+1;
end

X = [X X2 X3];
[Theta,BINT,R,RINT,STATS]= regress(Y,X)

