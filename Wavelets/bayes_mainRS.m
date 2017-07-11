
function [tm,cc,e]=bayes_mainRS(x1,x2,x3,y1,y2,y3,win,h,ovr,pr,s)
%the main function for the inference, propagation and other evaluations

%---inputs---
%ph1,ph2 - phase time-series vectors
%win - window in seconds
%h - sampling step e.g. h=0.01
%ovr - overlaping of windows; ovr=1 is no overlap; ovr=0.75 will overlap
%      the last 1/4 of each window with the next window
%pr - propagation constant
%s - print progress status if s=1 
%bn - order of Fourier base function

%---outputs---
%tm - time vector for plotting 
%cc - inferred mean parameters
%e  - inferred noise

%example for default input parameters and call of 
%the function >> [tm,cc,e]=bayes_main(ph1,ph2,40,0.01,1,0.2,1,2);
%%

cc=[]; 
win=win/h;
w=ovr*win;
pw=win*h*pr;

M=60;
L=6;

Cpr=zeros(M/L,L);XIpr=zeros(M);


%set the right dimensions for the vectors
[m,n]=size(x1);
if m<n
    x1=x1';
    x2=x2';
    x3=x3';
    y1=y1';
    y2=y2';
    y3=y3';
end


%% do the main calculations for each window
for i=0:floor((length(x1)-win)/w)
    
    xx1=x1(i*w+1:i*w+win); xx2=x2(i*w+1:i*w+win); xx3=x3(i*w+1:i*w+win);
    yy1=y1(i*w+1:i*w+win); yy2=y2(i*w+1:i*w+win); yy3=y3(i*w+1:i*w+win);
    
    %-----bayesian inference for one window------
    [Cpt,XIpt,E]=bayesRosslerLorenz(Cpr,XIpr,h,200,0.0001,xx1',xx2',xx3',yy1',yy2',yy3');
   
   
     %the propagation for the next window
     [XIpr,Cpr] = Propagation_function_XIpt(Cpt,XIpt,pw); 
    
  
 
 e(i+1,:,:)=E; 
 cc(i+1,:)=Cpt(:);
   
 %display progress 
 if s 
   display(['processed so far: t= ' num2str((i+1)*w*h) 's /' num2str(length(x1)*h) 's ;']);
   %Cpt
   %E
 end
 
end

%time vector for plotting
tm = (win/2:w:length(x1)-win/2)*h;
%%

 

function [ XIpr,Cpr ] = Propagation_function_XIpt(Cpt,XIpt,p)
% Propagation function with covariance
% find the new prior for the next block

% The average is not supposed to change
Cpr=Cpt;   


% Prepare the diffusion matrix
% Set a value for a particular parameter
Inv_Diffusion=zeros(length(XIpt));
invXIpt=inv(XIpt);
for i=1:length(Cpt(:)) 
    Inv_Diffusion(i,i) = p*p* invXIpt(i,i);
end

% The gaussian of the posterior is convoluted with another 
% gaussian which express the diffusion of the parameter.
XIpr=inv(( invXIpt + Inv_Diffusion ));
%%

 function [ XIpr,Cpr ] = Propagation_function_Cpt(Cpt,XIpt,p)
% Propagation function with parameters
% find the new prior for the next block

% The average is not supposed to change
Cpr=Cpt;   


% Prepare the diffusion matrix
% Set a value for a particular parameter
Inv_Diffusion=zeros(length(XIpt));
invXIpt=inv(XIpt);
for i=1:length(Cpt(:)) 
    Inv_Diffusion(i,i) = p*p* (Cpt(i));
end

% The gaussian of the posterior is convoluted with another 
% gaussian which express the diffusion of the parameter.
XIpr=inv(( invXIpt + Inv_Diffusion ));
%%