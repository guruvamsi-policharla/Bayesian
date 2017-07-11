function [Cpt, XIpt, E]=bayesRosslerLorenz(Cpr,XIpr,h,max_loops,eps,x,y,z,xx,yy,zz)
%infers the parameters and the noise within one block of data.


%---inputs---
%Cpr - prior vector of parameters
%XIpr - prior concentration matrix (i.e. inv(covariance) matrix)
%h - sampling rate (e.g. h=0.01)
%max_loops - number of loops to stop for the error for convergence
%eps - error/tolerance for convergence
%x,y,z,xx,yy,zz - input block of data

%---outputs---
%Cpt - posterior vector of parameters
%XIpt - poaterior concentration matrix (i.e. inv(covariance) matrix)
%E - noise matrix

%% ------------

%xS & yS - the midpoint time series
xS=(x(2:end)+x(1:end-1))/2; yS=(y(2:end)+y(1:end-1))/2;
zS=(z(2:end)+z(1:end-1))/2;
xxS=(xx(2:end)+xx(1:end-1))/2; yyS=(yy(2:end)+yy(1:end-1))/2;
zzS=(zz(2:end)+zz(1:end-1))/2;

%xT & yT - the derivative time series
xT=(x(2:end)-x(1:end-1))/h; yT=(y(2:end)-y(1:end-1))/h;
zT=(z(2:end)-z(1:end-1))/h;
xxT=(xx(2:end)-xx(1:end-1))/h; yyT=(yy(2:end)-yy(1:end-1))/h;
zzT=(zz(2:end)-zz(1:end-1))/h;
xyzT=[xT;yT;zT;xxT;yyT;zzT];  

L=6;   %number of dimensions 
M=60;  %number of base functions
K=M/L; %number of dimensions for the C vector
lastwarn(''); %reset the last warning message


%the following evaluations fill temporary variables to be used bellow in the main calculations 
p=calculateP(xS,yS,zS,xxS,yyS,zzS,K);
v1=calculateV(xS,yS,zS,xxS,yyS,zzS,K,1);
v2=calculateV(xS,yS,zS,xxS,yyS,zzS,K,2);
v3=calculateV(xS,yS,zS,xxS,yyS,zzS,K,3);
v4=calculateV(xS,yS,zS,xxS,yyS,zzS,K,4);
v5=calculateV(xS,yS,zS,xxS,yyS,zzS,K,5);
v6=calculateV(xS,yS,zS,xxS,yyS,zzS,K,6);


%loop to converge iteratively the result until desired precision is reached
C_old=Cpr; 
Cpt=Cpr; %initialize Cpt
for loop=1:max_loops
   
    [E]=calculateE(Cpt',xyzT,L,h,p);    
    [Cpt,XIpt]=calculateC(E,p,v1,v2,v3,v4,v5,v6,Cpr,XIpr,M,L,xyzT,h);
                   
	if(sum((C_old-Cpt).*(C_old-Cpt)./(Cpt.*Cpt))<eps)      
		return; 
    end
     C_old=Cpt;
end
display 'Warning: Desired precision not reached after max_loops.'  
%%

%----------additional functions: -------

%evaluates the base functions 
function [p]=calculateP(x,y,z,xx,yy,zz,K)   

    p(K,length(x))=0;
    
    p(1,:)=1;    
    p(2,:)=x;
    p(3,:)=y;
    p(4,:)=z;
    p(5,:)=x.*z;
    
    p(6,:)=xx;    
    p(7,:)=yy;
    p(8,:)=zz;
    p(9,:)=xx.*yy;
    p(10,:)=xx.*zz;      
    
    
%%

   
%evaluates the derivatives of base functions
function [v]=calculateV(x,y,z,xx,yy,zz,K,mr)
   v(K,length(x))=0; 
   
   if mr==1
        v(1,:)=0;    
        v(2,:)=1;
        v(3,:)=0;
        v(4,:)=0;
        v(5,:)=z;

        v(6,:)=0;    
        v(7,:)=0;
        v(8,:)=0;
        v(9,:)=0;
        v(10,:)=0; 
   else 
       if mr==2
        v(1,:)=0;    
        v(2,:)=0;
        v(3,:)=1;
        v(4,:)=0;
        v(5,:)=0;

        v(6,:)=0;    
        v(7,:)=0;
        v(8,:)=0;
        v(9,:)=0;
        v(10,:)=0; 
       else
           if mr==3
            v(1,:)=0;    
            v(2,:)=0;
            v(3,:)=0;
            v(4,:)=1;
            v(5,:)=x;

            v(6,:)=0;    
            v(7,:)=0;
            v(8,:)=0;
            v(9,:)=0;
            v(10,:)=0; 
           else
               if mr==4
                    v(1,:)=0;    
                    v(2,:)=0;
                    v(3,:)=0;
                    v(4,:)=0;
                    v(5,:)=0;

                    v(6,:)=1;    
                    v(7,:)=0;
                    v(8,:)=0;
                    v(9,:)=yy;
                    v(10,:)=zz;
               else
                   if mr==5
                        v(1,:)=0;    
                        v(2,:)=0;
                        v(3,:)=0;
                        v(4,:)=0;
                        v(5,:)=0;

                        v(6,:)=0;    
                        v(7,:)=1;
                        v(8,:)=0;
                        v(9,:)=xx;
                        v(10,:)=0; 
                   else
                       if mr==6
                            v(1,:)=0;    
                            v(2,:)=0;
                            v(3,:)=0;
                            v(4,:)=0;
                            v(5,:)=0;

                            v(6,:)=0;    
                            v(7,:)=0;
                            v(8,:)=1;
                            v(9,:)=0;
                            v(10,:)=xx; 
                       end
                   end
               end
               
           end
           
       end
   end
      
   

%%              

%evaluates the noise matrix
function [E]=calculateE(c,xyzT,L,h,p)
    E=zeros(L,L);
    
    E=E+(xyzT-c*p)*(xyzT-c*p)';  
    E=(h/length(xyzT(1,:)))*E;    
%%   
    

%final result - evaluates the parameters
function [Cpt,XIpt]=calculateC(E,p,v1,v2,v3,v4,v5,v6,Cpr,XIpr,M,L,xyzT,h)
    K=M/L;
    invr=inv(E);
        
    wr=lastwarn;
    if (~isempty(wr))&&((strcmp(wr(1:18),'Matrix is singular'))||(strcmp(wr(1:27),'Matrix is close to singular')))
        display('Singular matrix can lead to wrong and imprecise results. Please check the input parameters, signals or base functions. '); 
        error('Singular matrix.');
    end
    
    
    %evaluates the concentration matrix ------
    XIpt=zeros(M,M);
        
    XIpt(1:K,1:K)            =XIpr(1:K,1:K)+h*invr(1,1)*p*(p');
    XIpt(1:K,K+1:2*K)        =XIpr(1:K,K+1:2*K)+h*invr(1,2)*p*(p');
    XIpt(1:K,2*K+1:3*K)      =XIpr(1:K,2*K+1:3*K)+h*invr(1,3)*p*(p');
    XIpt(1:K,3*K+1:4*K)      =XIpr(1:K,3*K+1:4*K)+h*invr(1,4)*p*(p');
    XIpt(1:K,4*K+1:5*K)      =XIpr(1:K,4*K+1:5*K)+h*invr(1,5)*p*(p');
    XIpt(1:K,5*K+1:6*K)      =XIpr(1:K,5*K+1:6*K)+h*invr(1,6)*p*(p');
        
    XIpt(K+1:2*K,1:K)        =XIpr(K+1:2*K,1:K)+h*invr(2,1)*p*(p');
    XIpt(K+1:2*K,K+1:2*K)    =XIpr(K+1:2*K,K+1:2*K)+h*invr(2,2)*p*(p');
    XIpt(K+1:2*K,2*K+1:3*K)  =XIpr(K+1:2*K,2*K+1:3*K)+h*invr(2,3)*p*(p');
    XIpt(K+1:2*K,3*K+1:4*K)  =XIpr(K+1:2*K,3*K+1:4*K)+h*invr(2,4)*p*(p');
    XIpt(K+1:2*K,4*K+1:5*K)  =XIpr(K+1:2*K,4*K+1:5*K)+h*invr(2,5)*p*(p');
    XIpt(K+1:2*K,5*K+1:6*K)  =XIpr(K+1:2*K,5*K+1:6*K)+h*invr(2,6)*p*(p');
    
    XIpt(2*K+1:3*K,1:K)      =XIpr(2*K+1:3*K,1:K)+h*invr(3,1)*p*(p');
    XIpt(2*K+1:3*K,K+1:2*K)  =XIpr(2*K+1:3*K,K+1:2*K)+h*invr(3,2)*p*(p');
    XIpt(2*K+1:3*K,2*K+1:3*K)=XIpr(2*K+1:3*K,2*K+1:3*K)+h*invr(3,3)*p*(p');
    XIpt(2*K+1:3*K,3*K+1:4*K)=XIpr(2*K+1:3*K,3*K+1:4*K)+h*invr(3,4)*p*(p');
    XIpt(2*K+1:3*K,4*K+1:5*K)=XIpr(2*K+1:3*K,4*K+1:5*K)+h*invr(3,5)*p*(p');
    XIpt(2*K+1:3*K,5*K+1:6*K)=XIpr(2*K+1:3*K,5*K+1:6*K)+h*invr(3,6)*p*(p');
    
    XIpt(3*K+1:4*K,1:K)      =XIpr(3*K+1:4*K,1:K)+h*invr(4,1)*p*(p');
    XIpt(3*K+1:4*K,K+1:2*K)  =XIpr(3*K+1:4*K,K+1:2*K)+h*invr(4,2)*p*(p');
    XIpt(3*K+1:4*K,2*K+1:3*K)=XIpr(3*K+1:4*K,2*K+1:3*K)+h*invr(4,3)*p*(p');
    XIpt(3*K+1:4*K,3*K+1:4*K)=XIpr(3*K+1:4*K,3*K+1:4*K)+h*invr(4,4)*p*(p');
    XIpt(3*K+1:4*K,4*K+1:5*K)=XIpr(3*K+1:4*K,4*K+1:5*K)+h*invr(4,5)*p*(p');
    XIpt(3*K+1:4*K,5*K+1:6*K)=XIpr(3*K+1:4*K,5*K+1:6*K)+h*invr(4,6)*p*(p');
    
    XIpt(4*K+1:5*K,1:K)      =XIpr(4*K+1:5*K,1:K)+h*invr(5,1)*p*(p');
    XIpt(4*K+1:5*K,K+1:2*K)  =XIpr(4*K+1:5*K,K+1:2*K)+h*invr(5,2)*p*(p');
    XIpt(4*K+1:5*K,2*K+1:3*K)=XIpr(4*K+1:5*K,2*K+1:3*K)+h*invr(5,3)*p*(p');
    XIpt(4*K+1:5*K,3*K+1:4*K)=XIpr(4*K+1:5*K,3*K+1:4*K)+h*invr(5,4)*p*(p');
    XIpt(4*K+1:5*K,4*K+1:5*K)=XIpr(4*K+1:5*K,4*K+1:5*K)+h*invr(5,5)*p*(p');
    XIpt(4*K+1:5*K,5*K+1:6*K)=XIpr(4*K+1:5*K,5*K+1:6*K)+h*invr(5,6)*p*(p');
    
    XIpt(5*K+1:6*K,1:K)      =XIpr(5*K+1:6*K,1:K)+h*invr(6,1)*p*(p');
    XIpt(5*K+1:6*K,K+1:2*K)  =XIpr(5*K+1:6*K,K+1:2*K)+h*invr(6,2)*p*(p');
    XIpt(5*K+1:6*K,2*K+1:3*K)=XIpr(5*K+1:6*K,2*K+1:3*K)+h*invr(6,3)*p*(p');
    XIpt(5*K+1:6*K,3*K+1:4*K)=XIpr(5*K+1:6*K,3*K+1:4*K)+h*invr(6,4)*p*(p');
    XIpt(5*K+1:6*K,4*K+1:5*K)=XIpr(5*K+1:6*K,4*K+1:5*K)+h*invr(6,5)*p*(p');
    XIpt(5*K+1:6*K,5*K+1:6*K)=XIpr(5*K+1:6*K,5*K+1:6*K)+h*invr(6,6)*p*(p');
    
    
    %evaluates temp r ------
    r=zeros(K,L);
    ED=(E\xyzT); %same as: ED=(inv(E)*phiT); (but faster)

    r(:,1)=XIpr(1:K,1:K)*Cpr(:,1)+XIpr(1:K,K+1:2*K)*Cpr(:,2)+XIpr(1:K,2*K+1:3*K)*Cpr(:,3)+XIpr(1:K,3*K+1:4*K)*Cpr(:,4)+XIpr(1:K,4*K+1:5*K)*Cpr(:,5)+XIpr(1:K,5*K+1:6*K)*Cpr(:,6)    +h*(p*(ED(1,:)')-(1/2)*sum(v1,2));
    r(:,2)=XIpr(K+1:2*K,1:K)*Cpr(:,1)+XIpr(K+1:2*K,K+1:2*K)*Cpr(:,2)+XIpr(K+1:2*K,2*K+1:3*K)*Cpr(:,3)+XIpr(K+1:2*K,3*K+1:4*K)*Cpr(:,4)+XIpr(K+1:2*K,4*K+1:5*K)*Cpr(:,5)+XIpr(K+1:2*K,5*K+1:6*K)*Cpr(:,6)   +h*(p*(ED(2,:)')-(1/2)*sum(v2,2));
    r(:,3)=XIpr(2*K+1:3*K,1:K)*Cpr(:,1)+XIpr(2*K+1:3*K,K+1:2*K)*Cpr(:,2)+XIpr(2*K+1:3*K,2*K+1:3*K)*Cpr(:,3)+XIpr(2*K+1:3*K,3*K+1:4*K)*Cpr(:,4)+XIpr(2*K+1:3*K,4*K+1:5*K)*Cpr(:,5)+XIpr(2*K+1:3*K,5*K+1:6*K)*Cpr(:,6)   +h*(p*(ED(3,:)')-(1/2)*sum(v3,2));
    r(:,4)=XIpr(3*K+1:4*K,1:K)*Cpr(:,1)+XIpr(3*K+1:4*K,K+1:2*K)*Cpr(:,2)+XIpr(3*K+1:4*K,2*K+1:3*K)*Cpr(:,3)+XIpr(3*K+1:4*K,3*K+1:4*K)*Cpr(:,4)+XIpr(3*K+1:4*K,4*K+1:5*K)*Cpr(:,5)+XIpr(3*K+1:4*K,5*K+1:6*K)*Cpr(:,6)   +h*(p*(ED(4,:)')-(1/2)*sum(v4,2));
    r(:,5)=XIpr(4*K+1:5*K,1:K)*Cpr(:,1)+XIpr(4*K+1:5*K,K+1:2*K)*Cpr(:,2)+XIpr(4*K+1:5*K,2*K+1:3*K)*Cpr(:,3)+XIpr(4*K+1:5*K,3*K+1:4*K)*Cpr(:,4)+XIpr(4*K+1:5*K,4*K+1:5*K)*Cpr(:,5)+XIpr(4*K+1:5*K,5*K+1:6*K)*Cpr(:,6)   +h*(p*(ED(5,:)')-(1/2)*sum(v5,2));
    r(:,6)=XIpr(5*K+1:6*K,1:K)*Cpr(:,1)+XIpr(5*K+1:6*K,K+1:2*K)*Cpr(:,2)+XIpr(5*K+1:6*K,2*K+1:3*K)*Cpr(:,3)+XIpr(5*K+1:6*K,3*K+1:4*K)*Cpr(:,4)+XIpr(5*K+1:6*K,4*K+1:5*K)*Cpr(:,5)+XIpr(5*K+1:6*K,5*K+1:6*K)*Cpr(:,6)   +h*(p*(ED(6,:)')-(1/2)*sum(v6,2));
    


%final evaluation of parameters 
 C=(XIpt\vertcat(r(:,1),r(:,2),r(:,3),r(:,4),r(:,5),r(:,6)))'; 
 Cpt(:,1)=C(1:K); Cpt(:,2)=C(K+1:2*K); Cpt(:,3)=C(2*K+1:3*K); Cpt(:,4)=C(3*K+1:4*K); Cpt(:,5)=C(4*K+1:5*K); Cpt(:,6)=C(5*K+1:6*K);
%%


