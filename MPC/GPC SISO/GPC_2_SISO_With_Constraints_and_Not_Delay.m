%% CREATED BY: 
% Cesar Hernandez-Hernandez (PhD)
% e-mail: ceelch@gmail.com
%%
% THIS PROGRAM SHOWS AN EXAMPLE OF GENERALIZED PREDICTIVE CONTROL (GPC)
% WITH CONSTRAINTS AND WITHOUT DELAY

% TO SOLVE THIS PROBLEM WITH QUADPROG SHOULD HAVE THE FOLLOWING STRUCTURE 

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        min 0.5 u'Hu+b'u
%         u 
%
%         s.t Ru < = c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example:

% Continuous time transfer function

% num=[0.4 11.6]
% den=[1 2.231]

% sysc=tf(num,den)

%          0.4 s + 11.6
% G(s) =   -------------
%           s + 2.231
%

% Discrete time transfer function using 'c2d' and time
% of sampling of T = 0.1 seconds:

% sysd=c2d(sysc,0.1)

%          0.4 z + 0.6
% G(z) =   -------------
%             z - 0.8

%%
clc;
clear all;
close all; 
%%
disp('Example of an unconstrained GPC algorithm')
%% Continuous model
num=[0.4 11.6];
den=[1 2.231];

sysc=tf(num,den)
%% Discrete Model
T=0.1; % Sampling Time
sysd=c2d(sysc,T)
%%
B = [sysd.Numerator{1}(1) sysd.Numerator{1}(2)]; % Numerator
A = [sysd.Denominator{1}(1) sysd.Denominator{1}(2)];  % Denominator

% B = [0.4 0.6]; % Numerator
% A = [1 -0.8];  % Denominator
%%
na=length(A);
nb=length(B);

disp(sprintf('\n'));
disp('GPC model of the process:')
B
A

% THE CONTROL HORIZON MAY BE SMALLER THAN THE MAXIMUM PREDICTION HORIZON

N1=1; % Minimum prediction horizon
N2=3; % Maximum prediction horizon
Nu=3;  % Control horizon
lambda=0.3; % Control effort weight factor

disp(sprintf('\n'));
disp('Control Parameters:')
disp('[N1 ,N2 ,Nu]');
[N1 N2 Nu]
lambda

%% CONSTRAINTS IN THE CONTROL SIGNAL AND THE OUTPUT
Umin=0;Umax=1;
Ymin=0.4; Ymax=1.1;
%% The polynomials E and F are calculated
disp(sprintf('\n'));
disp('The polynomials E and F are calculated')

Ap=conv(A,[1 -1]); % A'(z)=A(z)*(1-z^{-1})

Dividendo=[1 zeros(1,length(Ap)-1)];
for j=1:N2
    [Eaux,Faux]=deconv(Dividendo,Ap);
    F(j,:)=Faux(2:end);
    Dividendo=[F(j,:) 0];
    E(j:N2,j)=ones(N2-j+1,1)*Eaux;
end    
      
F=F(1:end,:); % The matrices are taken from j=d+1.

E
F
%% The matrix G=Ej(z)*B(z) and G'(z) are obtained: 
for j=1:N2
    Gaux(j,:)=conv(E(j,:),B);
end

disp(sprintf('\n'));
disp('The matrix G=Ej(z)*B(z) and Gp(z) are obtained:')

Gaux

disp(sprintf('\n'));
disp('Matrix G is obtained from polynomials Gj')
disp('by removing the last element of each polynomial')

G=zeros(size(Gaux,1));
for i=1:size(Gaux,1)
    k=1;
    for j=i:-1:1
        G(i,k)=Gaux(i,j);
        k=k+1;
    end
end

G=G(1:N2,1:N2); % The effect of the delay is eliminated

G

disp(sprintf('\n'));
disp('The matrix Gp is obtained using the last elements of each polynomial Gj')
disp('forming a column vector:')

% The matrix Gp is obtained:
Gp=zeros(size(G,1),1);
for i=1:size(G,1)
    for j=1:1
        Gp(i,j)=Gaux(i,size(Gaux,2)-size(Gaux,1)+i+j-1);
    end        
end

Gp
%% The matrix H is obtained:
disp(sprintf('\n'));
disp('The matrix H is obtained:')
H = 2*(G'*G+lambda*eye(size(G'*G)))
%% THE MATRIX T IS CALCULATED
T=eye(Nu);

for i=1:Nu
    for j=1:Nu
        if (i>j)
            T(i,j)=1;
        end
    end
end

T
%% THE MATRIX R OF THE CONSTRAINTS IS CALCULATED
%R=[eye(Nu); -eye(Nu); T; -T; G; -G];
R=[T; -T; G; -G];
R
%% VECTOR 1
UNO=ones(N2,1);
UNO
%% THE EXAMPLE IS SIMULATED
tfinal=50; % SIMULATION TIME

% VARIABLE INICIALIZATION  
du=zeros(tfinal,1);
u=zeros(tfinal,1);
y=zeros(tfinal,1);

ref=ones(tfinal,1); % REFERENCE

refYmin= Ymin*ones(tfinal,1);  % TO PLOT THE CONSTRAINTS
refYmax= Ymax*ones(tfinal,1);  % 

inicio=Nu;
for k=inicio:tfinal-(N2-N1)    
    
    w=ref(k:k+N2-N1); % FUTURE REFERENCE. 
    
    % THIS PART TAKES DELTA u(t-1). ALWAYS GETTING A SCALAR
    dupasadas=[du(k-1:-1:k-1)]; 
    
    % THIS PART TAKES THE LAST 2 "y". ALWAYS GETTING A MATRIZ 2X1            
    ypasadas=[y(k:-1:k-na+1)]; 
    
    f = Gp*dupasadas+F*ypasadas; % FREE RESPONSE
    
    %du(k)=K1*(w-f);
            
    %%% GETTING THE VECTOR b TO BE USED IN THE QP PROBLEM
    b=2*G'*(f-w);
    
    %
    %ELEMENTS OF THE RIGHT OF INEQUALITY. 
    %
    
    %CONSTRAINTS IN THE CONTROL SIGNAL 
    
    T1_aux=UNO*Umax-UNO*u(k-1);
    T1_aux_2=-UNO*Umin+UNO*u(k-1);
    
    % CONSTRAINS IN THE OUTPUT
    G1_aux=UNO*Ymax-f;
    G1_aux_2=-UNO*Ymin+f;
    
    % MATRIX C
    
    c=[T1_aux; T1_aux_2; G1_aux; G1_aux_2];
    
    % THE QP IS CALCULATED    
    du_restricciones=quadprog(H,b,R,c,[],[],[],[],[],optimset('Display','off','LargeScale','off'));
    du(k)=du_restricciones(1,1); % JUST THE FIRST CONTROL SIGNAL IS CONSIDERED
    u(k)=du(k)+u(k-1); % CURRENT CONTROL SIGNAL
      
    for i=2:na
        y(k+1)=y(k+1)-A(1,i)*y(k-i+2);
    end
    
    for j=1:nb
        y(k+1)=y(k+1)+B(1,j)*u(k-(j-1));
    end     
end
%% Graphic Results
figure(1)
t=0:T:(tfinal-inicio-(N2-N1))*T;
subplot(2,1,1)
%stairs(t',ref(inicio:tfinal-(N2-N1)),'r:');
plot(t',ref(inicio:tfinal-(N2-N1)),'k:','LineWidth',2);
hold on;
plot(t',refYmin(inicio:tfinal-(N2-N1)),'g:','LineWidth',2);
hold on;
plot(t',refYmax(inicio:tfinal-(N2-N1)),'g:','LineWidth',2);
hold on;
plot(t',y(inicio:tfinal-(N2-N1)),'r','LineWidth',2);
title('System response')
xlabel('time')
ylabel('Output')
subplot(2,1,2)
plot(t',u(inicio:tfinal-(N2-N1)),'LineWidth',2);
hold on;
xlabel('time')
ylabel('Control signal')
