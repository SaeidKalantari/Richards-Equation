%MATLAB code for 1-D Richard equation using h-based formulation, fully
%implicit, tangent approximation.
%Writen by Saeed Kalantari and Danial Amini Baghbadorani
%For class of "Numerical Methods for Water Engineering", instructor:
%Professor Ataie-Ashtiani, Spring 2018-2019
%Reference paper:
% % Shahraiyni, H. T., & Ataie-Ashtiani, B. (2011). Mathematical forms and
% % numerical schemes for the solution of unsaturated flow equations. Journal
% % of Irrigation and Drainage Engineering, 138(1), 63-72.


clc
clear

L = 40;
% dz = 0.4;
% dt = 0.5;
dz = 1;
dt = 10;
N = L/dz +1;
T = 360;

hh = zeros(N,3);
kk = zeros(N,3);
CC = zeros(N,3);
QQ = zeros(N,3);

k =@(h) (0.00944*1.175e6)/(1.175e6+abs(h)^4.74);

C =@(h) (-0.1e1) * 0.135246672e1 * 0.1000000e7 * (abs(h) ^ 0.296e1) * sign(h)...
    / (0.1000000e7 * 0.1611e1 + (abs(h) ^ 0.396e1)) ^ 2;

theta =@(h) (1.611e6*(0.287-0.075))/(1.611e6+abs(h)^(3.96))+0.075;


for i=1:N
   
    z = (i-1)*dz;
    hh(i,1) = -61.5;
    QQ(i,1) = theta(hh(i,1));
    kk(i,1) = k(hh(i,1));
    CC(i,1) = C(hh(i,1));
    
end

A = zeros(N,N);
B = zeros(N,1);
iterate = zeros(T/dt,1);
for tt=1:T/dt
    hh(:,2) = hh(:,1);
    
    for ii=1:N
        z = (ii-1)*dz;
        QQ(ii,1) = theta(hh(ii,1));
        kk(ii,1) = k(hh(ii,1));
        CC(ii,1) = C(hh(ii,1));  
        
        QQ(ii,2) = theta(hh(ii,2));
        kk(ii,2) = k(hh(ii,2));
        CC(ii,2) = C(hh(ii,2)); 
    end
    for jj=1:200% m=jj

        for ii=1:N
            z = (ii-1)*dz;
            QQ(ii,2) = theta(hh(ii,2));
            kk(ii,2) = k(hh(ii,2));
            CC(ii,2) = C(hh(ii,2));                
        end
        for i=2:N-1
            
            if i>2
                A(i,i-1) = -1/(2*dz^2)*(kk(i,2)+kk(i-1,2));
            end
            if i<N-1
                A(i,i+1) = -1/(2*dz^2)*(kk(i,2)+kk(i+1,2));
            end
            A(i,i) = CC(i,2)/dt + 1/(2*dz^2)*(kk(i+1,2) + 2*kk(i,2) + kk(i-1,2));
            B(i,1) = 1/(2*dz)*(-kk(i+1,2)+kk(i-1,2))+CC(i,2)*hh(i,1)/dt;
            if i == 2
                B(i,1) = B(i,1)-(-1/(2*dz^2)*(kk(i,2)+kk(i-1,2)))*(-20.7);
            end
            if i == N-1
                B(i,1) = B(i,1)-(-1/(2*dz^2)*(kk(i,2)+kk(i+1,2)))*(-61.5);
            end
            
        end
        
        A=A(2:N-1,2:N-1);
        B=B(2:N-1,1);
        
        hX = A\B;

        for i=2:N-1
        
           hh(i,3) = hX(i-1,1);
        end
        hh(1,3)=-20.7;
        hh(N,3)=-61.5;

        temp = 0;
        for i = 1:N

           temp = max(temp,abs(hh(i,2)-hh(i,3)));

        end
        hh(:,2) = hh(:,3);
        

        if temp < 0.01
            iterate(tt) = jj; 
           break; 

        end
    end
    hh(:,1)=hh(:,2);
end
sum(iterate)
z= 0:dz:40;
plot(z,hh(:,1))



