%% Program: Árvore
%
%   @DESCRIÇÃO Atividade 4: Definir os modos de vibrações de um pórtico
%
%   @AUTOR Renan Miranda Portela @
%
%% Limpar memória e fechar janelas
clear 
close all
clc

%% Dados de entrada

n = 4; %quantidade de níveis
r = 0.7; %fator de escala
phi = [pi/9,-pi/9]; %ângulo de abertura
xb = [0,0]; %coordenadas de base
yb = [0,6.9];
E=210e9;               %módulo de young
rho=805e-3;                %densidade


h=figure;
 axes('Parent',h,'Color',[1,1,1]);

angulo=length(phi);
if length(r)==1
    rM=ones(1,angulo)*r;
elseif length(r)==angulo
    rM=r;
end

a=sqrt((xb(1)-xb(2))^2+(yb(1)-yb(2))^2);   % scale factor
b=1;
c=sqrt((1-xb(2))^2+(yb(2))^2);
alpha=acos((a^2+b^2-c^2)/(2*a*b));
k=(yb(2)-yb(1))/(xb(2)-xb(1));                       % define the tangent
d=(yb(1)*xb(2)-xb(1)*yb(2))/(xb(2)-xb(1));           % define the free member
 if yb(2)>=yb(1)
       ksi=alpha;
   else
       ksi=-alpha; 
 end
 
 %    auxiliary vectors for theta and rD (see below their definition)
    psi=zeros(n,angulo^n);
    coeficiente=ones(n,angulo^n);
  %   psi and ralt have Kantor array`s structur !!!
  
    for i=1:1:n
        z=1;
       for j=1:1:angulo^(i-1)
           for k=1:1:angulo
               for m=1:1:angulo^(n-i)
                   psi(i,z)=phi(k);               % define psi on the base of phi
                   coeficiente(i,z)=rM(k);               % coeficiente de redução do comprimento  
                   z=z+1;
               end
           end
       end
   end
  
   theta=zeros(n,angulo^n);  % vetor dos angulos entre cada galho do nível fractal em relação ao eixo das abcissas
   L=ones(n,angulo^n);      % comprimento do galho
   Diametro=ones(n,angulo^n);      % diametro do galho
   
   
   for i=1:1:angulo^n
       for j=1:1:n
           for k=1:1:j
           theta(j,i)=theta(j,i)+psi(k,i);        % angulo de abertura
           L(j,i)=L(j,i)*coeficiente(k,i);             % comprimento do galho
           Diametro(j,i)=Diametro(j,i)*coeficiente(k,i);             % comprimento do galho
            end
       end
   end
   
   theta=theta+ksi;
    
         
%% Matriz de coordenadas

   coord=ones(n+1,angulo^n,2);
% initial coordinates   
   coord(1,:,1)=ones(1,angulo^n)*xb(2);   % x-coordinate
   coord(1,:,2)=ones(1,angulo^n)*yb(2);   % y-coordinate
   for j=1:1:angulo^n
       for i=1:1:n
           
           % define following coordinates
           coord(i+1,j,1)=coord(i,j,1)+a*L(i,j)*cos(theta(i,j));
           coord(i+1,j,2)=coord(i,j,2)+a*L(i,j)*sin(theta(i,j));
       end
   end
     
 % By visualisation the method of inverse trace is used    
    tau=1;
      for i=1:angulo:angulo^n
          z=1;
           for k=1:1:angulo-1
                for j=z:1:n-1
                 line([coord(j,i,1),coord(j+1,i,1)],[coord(j,i,2),coord(j+1,i,2)],'Color',[0.8,0.5,0.5],'LineWidth',1.5);
                end
              z=z+1;
          end       
     end
 
 % -------------------------- end branches ----------------------------    
      for i=1:1:angulo^n
        line([coord(n,i,1),coord(n+1,i,1)],[coord(n,i,2),coord(n+1,i,2)],'Color',[0.3,0.75,0.4],'LineWidth',1);
      end   
 % -------------- Trunk -----------------         
    line([xb(1),xb(2)],[yb(1),yb(2)],'Color',[0,0,0],'LineWidth',2);
    
    coordx = coord(:,:,1)';
    coordy = coord(:,:,2)';
    
    A = unique([coordx(:,1),coordy(:,1)],'rows');
    B = unique([coordx(:,2),coordy(:,2)],'rows');
    C = unique([coordx(:,3),coordy(:,3)],'rows');
    D = unique([coordx(:,4),coordy(:,4)],'rows');
    coord = [A;B;C;D];
    
    if n >= 5
        E = unique([coordx(:,5),coordy(:,5)],'rows');
        coord = [A;B;C;D;E];
    end
    
    coord=[0,0;coord];
    
    for i = 1 : size(coord,1)
    coord(i,3) = i;
    end
    
    figure(2)
    scatter(coord(:,1),coord(:,2))
        
    %% matriz de incidencia
    
    d0 = 18e-2;
    d1 = d0*r;
    d2 = d1*r;
    d3 = d2*r;
    
    inci = [1 d0 coord(1,3) coord(2,3);
                 2 d1 coord(2,3) coord(3,3); %x = 2.41 y = 11
                 3 d1 coord(2,3) coord(4,3); %x = -2.41 y = 11
                 4 d2 coord(4,3) coord(8,3); %x = 5.34 y = 12
                 5 d2 coord(4,3) coord(7,3); %x = 2.41 y = 14
                 6 d2 coord(3,3) coord(6,3); %x = -2.41 y = 14
                 7 d2 coord(3,3) coord(5,3); %x = -5.34 y = 12
                 8 d3 coord(8,3) coord(16,3); %x = 7 y = 12
                 9 d3 coord(8,3) coord(15,3); %x = 6 y = 14
                 10 d3 coord(7,3) coord(14,3); %x = 3 y = 16
                 11 d3 coord(7,3) coord(13,3); %x = 1 y = 16
                 12 d3 coord(6,3) coord(12,3); %x = -1 y = 16
                 13 d3 coord(6,3) coord(11,3); %x = -3 y = 16
                 14 d3 coord(5,3) coord(10,3);
                 15 d3 coord(5,3) coord(9,3)];
             
     nnos=size(coord,1);                               % conta o número de nós
     nel=size(inci,1);                            % conta o número de elementos
     ngdl=3;                                        % número de graus lib p nó

% ------ condições de contorno ---------------------------------------------%
%
%   bc=[nó | grau de liberdade | valor]
%
%   Grau de liberdade 1 --> x
%   Grau de liberdade 2 --> y
%   Grau de liberdade 3 --> z
%   Grau de liberdade 4 --> ox
%   Grau de liberdade 5 --> oy
%   Grau de liberdade 6 --> oz


bc= [ 1 1 0;
      1 2 0;
      1 3 0];
  
% ------ carregamentos externos ------------------------------------------%
%
%   F=[nó | grau | valor]
%
%   Grau de liberdade 1 --> Fx
%   Grau de liberdade 2 --> Fy
%   Grau de liberdade 3 --> Fz
%   Grau de liberdade 4 --> Mx
%   Grau de liberdade 5 --> My
%   Grau de liberdade 6 --> Mz

    Load = [ 8 2 -1]   ; 
    
% ------ matriz identificação ------------------------------------------%

id=ones(ngdl,nnos);
nbc = size(bc,1); %numero de condições de contorno

for i=1:nbc
    id(bc(i,2),bc(i,1))=bc(i,3);
end

neq = 0; %número de equações

for i= 1:nnos
    for j = 1:ngdl
        if id(j,i)== 1
            neq = neq +1;
            id(j,i) = neq;
        end
    end
end

kg = zeros(neq,neq); % pré locação da matriz de rigidez global
mg = zeros(neq,neq); % pré locação da matriz de massa global
F = zeros(neq,1); %pré locação da matriz coluna de forças 

nloads = size(Load,1); %número de carregamentos

for i=1:nloads
    F(id(Load(i,2),Load(i,1)),1) = Load(i,3); 
end

for i = 1:nel
    no1 = inci(i,3); %nó 1 do elemento
    no2 = inci(i,4); %nó 2 do elemento
    
    x1 = coord(no1,1); %coordenada x do ponto 1
    x2 = coord(no2,1); %coordenada x do ponto 2
    
    y1 = coord(no1,2); %coordenada y do ponto 1
    y2 = coord(no2,2); %coordenada y do ponto 2
    
    L = sqrt((x2-x1)^2+(y2-y1)^2); %comprimento da barra
    
    a=pi*inci(i,2)^2/4;                         %área da barra
    I=pi*inci(i,2)^4/64;                        %Momento de Inércia 
    
    %% ângulo da barra em relação ao eixo x
    
    c=(x2-x1)/L;                                %cosseno
    s=(y2-y1)/L;                                %seno
    
    
    T=[    c     s    0   0   0   0;     %matriz T
          -s     c     0   0   0   0;
           0     0     1   0   0   0;
           0     0     0   c   s  0;
           0     0     0   -s   c   0;
           0     0     0   0   0   1];
       
   k=a/L; cc=I/L^3;
   ke=E*[k      0       0       -k  0           0;
         0      cc*12    cc*6*L   0   -cc*12       cc*6*L;
         0      cc*6*L   cc*4*L^2 0   -cc*6*L     cc*2*L^2;
         -k     0       0       k   0           0;
         0      cc*-12   -cc*6*L  0   12*cc        -6*L*cc;
         0      cc*6*L   cc*2*L^2 0   cc*-6*L      cc*4*L^2];  
     
     ke=T'*ke*T;
     
     ke = round(ke,4);
     
     kk=1/420;
    me=rho*a*L*[1/3    0        0       1/6     0           0;
                0      kk*156    kk*22*L  0       kk*54        -kk*13*L;
                0      kk*22*L   kk*4*L^2 0       kk*13*L      -kk*L^2;
                1/6    0        0       1/3     0           0;
                0      kk*54     kk*13*L  0       kk*156       kk*-22*L ;
                0      -kk*13*L  -kk*L^2  0       kk*-22*L     kk*4*L^2];
    me=T'*me*T;
    
    me = round(me,4);
    
    loc = [id(1,no1) id(2,no1) id(3,no1) id(1,no2) id(2,no2) id(3,no2)]; 
    
    for j = 1:6
        if loc(j) ~= 0
            for k = 1:6
                if loc(k) ~=0
                    kg(loc(j),loc(k))=kg(loc(j),loc(k))+ke(j,k);
                    mg(loc(j),loc(k))=mg(loc(j),loc(k))+me(j,k);
                end
            end
        end
    end
end    

u = kg\F; 

[theta,D]=eig(kg,mg); %theta = modo de vibração

% for i = 1:size(theta,1)
%     theta(:,i) = theta(:,i)/max(abs(theta(:,i)));    
% end

for i = 1:size(theta,1)
    theta(:,i) = theta(:,i)/max(abs(theta(:,i)));    
end

lambda = diag(D); %autovalores [K] - lambda*[M] = 0 
omega = sqrt(lambda); %frequência de vibração

d = zeros(size(coord,1),3);
 
for i=1:nnos
    if id(1,i)==0
        d(i,1)=0;
    else
        d(i,1) = u(id(1,i));
    end
    if id(2,i)==0
        d(i,2)=0;
    else
        d(i,2) = u(id(2,i));
    end
    if id(3,i)==0
        d(i,3)=0;
    else
        d(i,3) = u(id(3,i));
    end
end

coord2 = coord(:,1:2)+10^4*d(:,1:2);

xlabel('X','FontSize',10, 'Color',[.8 .8 .8]); 
ylabel('Y','FontSize',10, 'Color',[.8 .8 .8]);
xmin = min(coord(:,1));
xmax = max(coord(:,1));
ymin = min(coord(:,2));
ymax = max(coord(:,2));
Lx = xmax - xmin;
Ly = ymax - ymin;
xmin = xmin - 0.8*Lx;
xmax = xmax + 0.8*Lx;
ymin = ymin - 0.8*Ly;
ymax = ymax + 0.8*Ly;

figure('Position',get(0,'ScreenSize'))              
hold on

 for i=1:nel
    plot(coord2(inci(i,3:4),1),coord2(inci(i,3:4),2),'*-','Color','blue');
    plot(coord(inci(i,3:4),1),coord(inci(i,3:4),2),'-o','Color','black');
 end
 
u=zeros((size(theta,1)+3)/ngdl,size(theta,2));
v=zeros((size(theta,1)+3)/ngdl,size(theta,2));
k=1; m=1;i=1;

for z=1:size(theta,2)                    
    for m=1:size(id,2)           
        if id(1,m)==0
            u(k,z)=0;
            k=k+1;
        elseif id(1,m)~=0
            u(k,z)=theta(i,z);
            i=i+3;
            k=k+1;
        else
        end
    end
    k=1;i=1;
end

n=1;k=1;i=2;
for z=1:size(theta,2)
    for m=1:size(id,2)
        if id(2,m)==0
            v(n,z)=0;
            n=n+1;
        elseif id(2,m)~=0
            v(n,z)=theta(i,z);
            n=n+1;
            i=i+3;
        else
        end
    end
    
    n=1;
    k=1;
    i=2;
end

modosu=zeros(16,45);
modosv=zeros(16,45);
for i=1:16
    for j=1:45
        modosu(i,j)=coord(i,1)+u(i,j);
        modosv(i,j)=coord(i,2)+v(i,j);
    end
end

modos = 5;

for j=1:modos
    figure('Position',get(0,'ScreenSize'))      %ativar em caso de plots somente
    title('Modo de vibração')

    for i=1:nel
        
        hold on
        grid on
        grid minor
        axis equal
        plot(modosu(inci(i,3:4),j),modosv(inci(i,3:4),j),'*-','Color','blue');
        plot(coord(inci(i,3:4),1),coord(inci(i,3:4),2),'-o','Color','black');      

    end   
      Leg=cell(1,1);  
      Leg{1}=strcat('Modo  ', num2str(j));
      legend(Leg)

 end