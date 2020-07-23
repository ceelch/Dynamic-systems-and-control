%% CREATED BY: 
% Cesar Hernandez-Hernandez (PhD)
% e-mail: ceelch@gmail.com
%%
% THIS SCRIPT SHOW AN EXAMPLE OF MUTIVARIBALE GPC WITHOUT CONSTRAINTS AND WITHOUT DELAYS 

% THIS EXAMPLE CAN BE FOUND IN PAG. 144. "MODEL PREDICTIVE CONTROL" BY CAMACHO AND BORDONS

% Example:

% Continuous time transfer function

% Transfer Matrix
 
% |Y1(s)|   [ 1/(1+0.7s)  5/(1+0.3s) ][U1(s)]
% |Y2(s)| = [ 1/(1+0.5s)  2/(1+0.4s) ][U2(s)]

% Discretizing the model with a sampling time of 0.03 minutes

% |y1(t)|   [ 0.0420z^-1/(1-0.9580z^-1)  0.4758z^-1/(1-0.9048z^-1) ][u1(t)]  
% |y2(t)| = [ 0.0582z^-1/(1-0.9418z^-1)  0.1445z^-1/(1-0.9277z^-1) ][u2(t)] 

% |y1(t)|   [ 0.0420/(1-0.9580z^-1)  0.4758/(1-0.9048z^-1) ][u1(t-1)]  
% |y2(t)| = [ 0.0582/(1-0.9418z^-1)  0.1445/(1-0.9277z^-1) ][u2(t-1)] 

% And of the structure A(z^-1)y(t)= B(z^-1)u(t-1)

% [(1-0.9580z^-1)(1-0.9048z^-1)                 0            ][y1(t)]   [ 0.0420(1-0.9048z^-1)     0.4758(1-0.9580z^-1)][u1(t-1)]  
% [          0                   (1-0.9418z^-1)(1-0.9277z^-1)][y2(t)] = [ 0.0582(1-0.9277z^-1)     0.1445(1-0.9418z^-1)][u2(t-1)] 

% this is

% [1-1.8628z^-1+0.8668z^-2                 0            ][y1(t)]   [ 0.0420-0.03637z^-1     0.4758-0.4558z^-1][u1(t-1)]  
% [          0                   1-1.8695z^-1+0.8737z^-2][y2(t)] = [ 0.0582-0.05545z^-1     0.1445-0.1361z^-1][u2(t-1)] 
%%
clc
clear all;
close all;
%%
format long
%%
disp('Ejemplo para un algoritmo GPC multivariable sin restricciones y sin retardo')

nOutputs=2;
mInputs=2;

T=0.03; % Sampling time

% THE CONTROL HORIZON MAY BE SMALLER THAN THE MAXIMUM PREDICTION HORIZON

N1=1; % Minimum prediction horizon
N2=3; % Maximum prediction horizon
Nu=2;  % Control horizon
lambda=0.05; % Control effort weight factor

q=0.05; % Coefficient that will multiply the control weight matrix Q
r=1; %Coefficient that will multiply the control weight matrix R 

disp(sprintf('\n'));
disp('Control Parameters:')
disp('[N1 ,N2 ,Nu]');
[N1 N2 Nu]
%lambda;

A=cell(mInputs,mInputs);
B=cell(mInputs,mInputs);

A{1,1}=[1 -1.8628 0.8667984];
A{1,2}=[0 0 0];
A{2,1}=[0 0 0];
A{2,2}=[1 -1.8695 0.873707786];

B{1,1}=[0.04020 -0.03637296];
B{1,2}=[0.4758 -0.4558164];
B{2,1}=[0.0582 -0.055451214];
B{2,2}=[0.1445 -0.1360901];
%% Tha matrix A tilde is calculated
At=cell(nOutputs,mInputs);

for i=1:nOutputs
    for j=1:mInputs
        At{i,j}=conv(A{i,j},[1,-1]);
    end
end

E1aux=cell(nOutputs,mInputs);

for i=1:nOutputs
    for j=1:mInputs
        if i==j
            E1aux{i,j}=[1 zeros(1,length(At{1,1})-1)];
        else
            E1aux{i,j}=[zeros(1,length(At{1,1}))];
        end
    end
end
    
F1aux=cell(nOutputs,mInputs);

for i=1:nOutputs
    for j=1:mInputs
        F1aux{i,j}=E1aux{i,j}-At{i,j};
    end
end
%% The Matrix E, F and Rj are calculated
    
E=cell(N2,1);
F=cell(N2,1);
Rj=cell(N2,1);

% E1 is calculated

for k=1:1
    
    for i=1:nOutputs
        for j=1:mInputs
            if i==j
                E{k,1}{i,j}=[1];
            else
                E{k}{i,j}=0;
            end
        end
    end
    
    % F1 is calculated
    
    for i=1:nOutputs
        for j=1:mInputs
            F{k,1}{i,j}=F1aux{i,j}(2:length(F1aux{i,i}));
        end
    end
    
    % R1 is calculated
    
    for i=1:nOutputs
        for j=1:mInputs
            Rj{k,1}{i,j}=F{k}{i,j}(k);
        end
    end    
end
  
for i=2:N2
        
    % E2 is calculated to E(N2)
    
    for j=1:nOutputs
        for k=1:mInputs            
            if j==k
            E{i,1}{j,k}=[E{i-1,1}{j,k} Rj{i-1,1}{j,k}(1)];
            else
              E{i,1}{j,k}=[zeros(1,length(E{i,1}{1,1}))];
            end
        end
    end
    
    % F2 is calculated to F(N2)
    
    for j=1:nOutputs
        for k=1:mInputs            
            if j==k               
               for m=2:length(At{j,j})
                   if m<length(At{j,j})
                    F{i,1}{j,k}(m-1)=F{i-1,1}{j,k}(m)-Rj{i-1,1}{j,k}*At{j,k}(m);
                   else
                       if m==length(At{j,j})
                         F{i,1}{j,k}(m-1)=-Rj{i-1,1}{j,k}*At{j,k}(m);
                       end
                   end
               end
               
               % This part fills the cells with zero value polynomials
               
            else
                F{i,1}{j,k}=[zeros(1,length(F{i,1}{1,1}))];
               
            end
        end
    end
    
    % Rj2 is calculated to Rj(N2)
        
    for j=1:nOutputs
        for k=1:mInputs
            Rj{i,1}{j,k}=F{i,1}{j,k}(1);
        end
    end    
end
      
% matrix F is re-ordered

matrixFaux=cell(N2,1);

% matrizFaux{1,1}{1,1}(1,1)=F{1,1}{1,1}(1);
% matrizFaux{1,1}{1,1}(1,2)=F{1,1}{1,2}(1);
for i=1:N2
    for j=1:N2
        for k=1:Nu
            for m=1:Nu
                matrixFaux{i,1}{1,j}(k,m)=F{i,1}{k,m}(j);   
            end
        end
    end
end
%%

% We calculate Ej(z^-1)B(z^-1) that we will call Gaux, to extract the matrix Gj and G' from there, 
% we will form the matrix G from the matrix Gj

 Gaux=cell(N2,1);
 
%  Gaux{1,1}{1,1}=conv(E{1,1}{1,1},B{1,1});
%  Gaux{1,1}{1,2}=conv(E{1,1}{1,1},B{1,2});
for i=1:N2
    for j=1:mInputs
        for k=1:mInputs
            Gaux{i,1}{j,k}=conv(E{i,1}{j,j},B{j,k});
        end
    end
end
 
 % The Gj matrices are obteined
 
 Gj=cell(N2,1);
 
%  Gj{1,1}{1,1}=[Gaux{1,1}{1,1}(1:1)];
%  Gj{1,1}{1,2}=[Gaux{1,1}{1,2}(1:1)];
for i=1:N2    
    for j=1:mInputs
        for k=1:mInputs
            Gj{i,1}{j,k}=[Gaux{i,1}{j,k}(1:i)];
        end
    end    
end
 
 % The G' matrix is calculated
 
 Gprima=cell(N2,1);
 
%  Gp{1,1}{1,1}=Gaux{1,1}{1,1}(2);
%  Gp{1,1}{1,2}=Gaux{1,1}{1,2}(2);
for i=1:N2
    for j=1:mInputs
        for k=1:mInputs
            Gprima{i,1}{j,k}=Gaux{i,1}{j,k}(i+1);
        end
    end
end
 %% The G matrix is calculated
 
 G=cell(N2,Nu);
 
%  G{1,1}{1,1}(1)=Gj{1,1}{1,1}(1);
%  G{1,1}{1,1}(2)=Gj{1,1}{1,2}(1);

% In this loop column 1 of G is calculated, that is, from G{1,1} to G{N2,1}

for i=1:N2    
    for j=1:nOutputs
        for k=1:mInputs            
            G{i,1}{j,1}(k)=Gj{i,1}{j,k}(i);
        end
    end
end

% We create the top diagonal of zeros
 
Gsup=cell(nOutputs,1);
  
for i=1:nOutputs
    Gsup{i,1}=[zeros(1:mInputs)];
end         
 
% In this loop the matrix G (cell type) is filled 
 
for i=1:N2    
    for j=1:Nu        
        if i==j
            G{i,j}=G{1,1};
        else
            if (i>j&&~(i==N2&&j==1))
                G{i,j}=G{2,1};
            else
                if (i==N2&&j==1)
                    G{i,j}=G{3,1};
                else
                    if i<j
                        G{i,j}=Gsup;
                    end
                end
            end
        end
    end    
end 
 %% THE WEIGHT Q AND R ARE CREATED
 
 Q=cell(Nu,Nu);
 
for i=1:Nu
    for j=1:Nu
        if i==j
            Q{i,j}=q*eye(nOutputs);
        else
            Q{i,j}=zeros(nOutputs,mInputs);
        end
    end
end

R=cell(N2,N2);

for i=1:N2
    for j=1:N2
        if i==j
            R{i,j}=r*eye(nOutputs);
        else
            R{i,j}=zeros(nOutputs,mInputs);
        end
    end
end 
 %% WE CONVERT TO MATRICES THE CELL TYPE MATRICES: F, Gprima, G, Q y R
 
 %We convert matrix F whose name was matrixFaux cell type to matrix and call it "matrixF"
  
 for i=1:N2
     for j=1:1
         matrixF{i,j}= cell2mat(matrixFaux{i,j});
     end
 end
 
matrixF=cell2mat(matrixF)
 
% We convert the matrix Gprima cell type to matrix and call it "matrixGp" 
for i=1:N2
    for j=1:1
        matrixGp{i,j}= cell2mat(Gprima{i,j});
    end
end
 
matrixGp=cell2mat(matrixGp);

%We convert matrix G cell type to matrix and call it "matrixG"
 
for i=1:N2
     for j=1:Nu
         matrixG{i,j}= cell2mat(G{i,j});
     end
end
 
matrixG=cell2mat(matrixG);

%We convert the matrix Q cell type to matrix and call it "matrixQ" 
matrixQ=cell2mat(Q);

%We convert the matrix R cell type to matrix and call it "matrixR"
matrixR=cell2mat(R);
%% FINALLY THE CONTROLLER K MATRIX IS OBTAINED, FROM WHICH WE WILL ONLY USE THE FIRST m ROWS
K=inv((matrixG'*matrixR*matrixG+matrixQ))*matrixG'*matrixR;
K1=K(1:mInputs,:);
%% WE GET THE REFERENCE MATRIX
tfinal=100; 

ref=ones(tfinal,1); % Reference

w=cell(1,mInputs);

for i=1:mInputs
    w{1,i}=ref;
end

% we transform the matrix cell reference to matrix

matrixW=cell2mat(w);
  
%W(:,1)=0.5*matrixW(:,1);
W(:,2)=0.3*matrixW(:,2);
  
for i=1:tfinal
    if i<=(tfinal/2)
        W(i,1)=0.5*matrixW(i,1);
    else
        W(i,1)=0.4*matrixW(i,1);
    end
end
%% WE INITIALIZE THE VARIABLES u, du, y
  
u_aux=zeros(tfinal,1); 
du_aux=zeros(tfinal,1);
y_aux=zeros(tfinal,1);

%u=cell(1,mEntradas);

for i=1:mInputs
    u{1,i}=u_aux;
    du{1,i}=du_aux;
    y{1,i}=y_aux;
end

% we transform the control matrices u, du and y from cell to matrices
 
U=cell2mat(u);
dU=cell2mat(du);
Y=cell2mat(y);

start=Nu+1;

for k=start:tfinal-(N2-N1)
    
    %  dupasadas=[dU(inicio-1:-1:inicio-1); dU(inicio-1:-1:inicio-1)];
    %
    %  ypasadas=[Y(inicio,1); Y(inicio,2); Y(inicio-1,1); Y(inicio-1,2); Y(inicio-2,1); Y(inicio-2,2)]
    %
    %  refW=[W(inicio,1); W(inicio,2); W(inicio+1,1); W(inicio+1,2);W(inicio+N2-N1,1); W(inicio+N2-N1,2)]; % REFERENCIA FUTURA.
    %
    %  f=matrizGp*dupasadas+matrizF*ypasadas
    %
    %  dU(inicio,1)=K1(1,:)*(refW-f);
    %  dU(inicio,2)=K1(2,:)*(refW-f);
    %
    %  U(inicio,1)=dU(inicio,1)+U(inicio-1,1);
    %  U(inicio,2)=dU(inicio,2)+U(inicio-1,2);
    %
    %  %Y(inicio+1,1)=Y(inicio+1,1)-(A{1,1}(2)*Y(inicio))
    %
    %  for i=2:length(A{1,1})
    %      for j=1:2
    %         Y(inicio+1,1)=Y(inicio+1,1)-(A{1,1}(i)*Y(inicio-i+2,j));
    %      end
    %  end
    %
    %  for i=1:length(B{1,1})
    %        Y(inicio+1)=Y(inicio+1)+(B{1,1}(i)*U(inicio-i+1,1));
    %  end
    %
    %  for i=1:length(B{1,2})
    %        Y(inicio+1)=Y(inicio+1)+(B{1,2}(i)*U(inicio-i+1,2));
    %  end
    
    dupasadas=[dU(k-1,1); dU(k-1,2)];
    ypasadas=[Y(k,1); Y(k,2); Y(k-1,1); Y(k-1,2); Y(k-2,1); Y(k-2,2)];
    refW=[W(k,1); W(k,2); W(k+1,1); W(k+1,2);W(k+N2-N1,1); W(k+N2-N1,2)]; % FUTURE REFERENCE.
    
    f=matrixGp*dupasadas+matrixF*ypasadas;
    
    dU(k,1)=K1(1,:)*(refW-f);
    dU(k,2)=K1(2,:)*(refW-f);
    
    U(k,1)=dU(k,1)+U(k-1,1);
    U(k,2)=dU(k,2)+U(k-1,2);
    
    % to calculate Y1 future
    
    for i=2:length(A{1,1})
        Y(k+1,1)=Y(k+1,1)-(A{1,1}(i)*Y(k-i+2,1));
    end
    
    for i=1:length(B{1,1})
        Y(k+1,1)=Y(k+1,1)+(B{1,1}(i)*U(k-i+1,1));
    end
    
    for i=1:length(B{1,2})
        Y(k+1,1)=Y(k+1,1)+(B{1,2}(i)*U(k-i+1,2));
    end
    
    % to calculate Y2 future
    
    for i=2:length(A{2,2})
        Y(k+1,2)=Y(k+1,2)-(A{1,1}(i)*Y(k-i+2,2));
    end
    
    for i=1:length(B{1,1})
        Y(k+1,2)=Y(k+1,2)+(B{2,1}(i)*U(k-i+1,1));
    end
    
    for i=1:length(B{1,2})
        Y(k+1,2)=Y(k+1,2)+(B{2,2}(i)*U(k-i+1,2));
    end
end
%% Graphic Results

figure(1)
%t=0:T:(tfinal-inicio-(N2-N1))*T;
t=0:(tfinal-start-(N2-N1));
subplot(2,1,1)
%stairs(t',ref(inicio:tfinal-(N2-N1)),'r:');
plot(t',W(start:tfinal-(N2-N1),1),'r:','LineWidth',2);
hold on;
plot(t',W(start:tfinal-(N2-N1),2),'g:','LineWidth',2);  
hold on;
plot(t',Y(start:tfinal-(N2-N1),1),'b','LineWidth',2);
hold on;
plot(t',Y(start:tfinal-(N2-N1),2),'k','LineWidth',2);
legend('ref1','ref2','Y1','Y2')
title('System response: Y1, Y2')
xlabel('time')
ylabel('Outputs')

subplot(2,1,2)
plot(t',U(start:tfinal-(N2-N1),1),'b','LineWidth',2);
hold on;
plot(t',U(start:tfinal-(N2-N1),2),'k','LineWidth',2);
legend('U1','U2')

title('U1, U2')
xlabel('time')
ylabel('Control Signals') 
 