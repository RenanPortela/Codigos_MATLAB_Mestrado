%% Program: Treliça
%
%   @DESCRIÇÃO Atividade 2: Definir os deslocamentos nodais de uma barra 
%   com seção constante, tensões internas e reações@
%
%   @AUTOR Renan Miranda Portela @
%
%% Limpar memória e fechar janelas
clear 
close all
clc

%% ------ Matriz de Coordenadas do sistema ---------------------------------%
%   coord=[nº | X | Y] => Matriz de coordenadas dos nós


coord = [1 0 0;
        2 4 0;
        3 8 0;
        4 12 0;
        5 4 3;
        6 8 3];
    
%-------Matriz de incidência ---------------------------------------------%
% inci=[nº | tab.mat | tab.geo | nó1 | nó2 ] => Matriz de
% incidência do elemento

inci = [1 1 1 1 2;
        2 1 1 2 3;
        3 1 1 3 4;
        4 1 1 1 5;
        5 1 1 2 5;
        6 1 1 3 5;
        7 1 1 3 6;
        8 1 1 4 6;
        9 1 1 5 6];
    
% ------ Tabela de Materiais ---------------------------------------------%
%
%  Tmat[ Material 1 | Material 2 | Material 3 ...]
%
%  coluna 1:Mod de elasticidade
%  coluna 2:Coeficiente de Poisson
  
tabmat=[210E9 0.33];

% ------ Tabela de Geometria ---------------------------------------------%
%
%  Tgeo[ Geometria 1 | Geometria 2 | Geometria 3 ...]
%
%  coluna 1: Área da seção transversal /m²

tabgeo = 0.0025 ;

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


bc=[1 1 0;
    1 2 0;
    4 2 0];

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

    Load=[3 2 -1200;
          6 1 400];
      
%%   tamanho das matrizes

nnos = size(coord,1); %numero de nós
nel = size(inci,1); %número de elementos
nF = size(Load,1); %numero de forças externas
neq = 0; %número de equações
ngdl = 2; %número de graus de liberdade

%% Cálculo das reações

M = 0;

    for i = 1:nF
        if Load(i,2) == 1
          M = M - Load(i,3)*coord(Load(i,1),3);
        elseif Load(i,2) == 2
          M = M + Load(i,3)*coord(Load(i,1),2);  
        end
    end

    R4 = -M/coord(4,2); %reação no apoio
    R1x = -R4; %reação no engaste no eixo x
    R1y = 0; %reação no engaste no eixo y

%reação no engaste
    for i = 1:nF
        if Load(i,2) == 2
            R1x = R1x - Load(i,3);
        else
            R1y = R1y - Load(i,3);
        end
    end

fprintf(' Reações da treliça \n')
fprintf('  R4      R1x     R1y \n')
fprintf(' %d N   %d N   %d N\n',R4,R1x,R1y)
fprintf('\n') %pular linha de comando

id=ones(ngdl,nnos);
nbc = size(bc,1); %numero de condições de contorno

for i=1:nbc
    id(bc(i,2),bc(i,1))=bc(i,3);
end

for i= 1:nnos
    for j = 1:ngdl
        if id(j,i)== 1
            neq = neq +1;
            id(j,i) = neq;
        end
    end
end

kg = zeros(neq,neq); % pré locação da matriz de rigidez global

for i = 1:nel
    no1 = inci(i,4); %nó 1 do elemento
    no2 = inci(i,5); %nó 2 do elemento
    
    x1 = coord(no1,2); %coordenada x do ponto 1
    x2 = coord(no2,2); %coordenada x do ponto 2
    
    y1 = coord(no1,3); %coordenada y do ponto 1
    y2 = coord(no2,3); %coordenada y do ponto 2
    
    L = sqrt((x2-x1)^2+(y2-y1)^2); %comprimento da barra
    
    mat = inci(i,2); %material do elemento
    geo = inci(i,3); %geometria do elemento
    E = tabmat(1,mat); %módulo de elasticidade do material
    A = tabgeo(1,mat); %área do elemento
    
    %% ângulo da barra em relação ao eixo x
    
    if (x2-x1)==0 
        if y2>y1
            beta = 2*atan(1);
        else
            beta = -2*atan(1);
        end
    else
        beta = atan((y2-y1)/(x2-x1));
    end
    
    %% matriz de rigidez do elemento
    
    c = cos(beta);
    s = sin(beta);
    ke = E*A/L*[c*c c*s -c*c -c*s;
           c*s s*s -c*s -s*s;
           -c*c -c*s c*c c*s;
           -c*s -s*s c*s s*s];
        
    loc=[id(1,no1) id(2,no1) id(1,no2) id(2,no2)]; %matriz de localização da matriz de rigidez local na global
    
    %% matriz de rigidez global
    
    ngdl = 4;
    
    for j = 1:ngdl
        if loc(j) ~= 0
            for k = 1:ngdl
                if loc(k) ~=0
                    kg(loc(j),loc(k))=kg(loc(j),loc(k))+ke(j,k); %matriz de rigidez global
                end
            end
        end
    end    
      
end

%% matriz coluna Força

F = zeros(neq,1); %pré locação da matriz coluna de forças 
nloads = size(Load,1); %número de carregamentos
        for i=1:nloads
            F(id(Load(i,2),Load(i,1)),1) = Load(i,3); 
        end

u = kg\F; 

[ d ] = matriz_deslocamento( id,nnos,u ); %matriz deslocamento 

fprintf(' Deslocamento em cada nó \n')
disp(d)
fprintf('\n') %pular linha de comando

stress = zeros(nel,1);

for i = 1:nel
    
    no1 = inci(i,4); %nó 1 do elemento
    no2 = inci(i,5); %nó 2 do elemento
    
    x1 = coord(no1,2); %coordenada x do ponto 1
    x2 = coord(no2,2); %coordenada x do ponto 2
    
    y1 = coord(no1,3); %coordenada y do ponto 1
    y2 = coord(no2,3); %coordenada y do ponto 2
    
    L = sqrt((x2-x1)^2+(y2-y1)^2); %comprimento da barra
    
    mat = inci(i,2); %material do elemento
    geo = inci(i,3); %geometria do elemento
    E = tabmat(1,mat); %módulo de elasticidade do material
    A = tabgeo(1,mat); %área do elemento
    
    %% ângulo da barra em relação ao eixo x
    
    if (x2-x1)==0 
        if y2>y1
            beta = 2*atan(1);
        else
            beta = -2*atan(1);
        end
    else
        beta = atan((y2-y1)/(x2-x1));
    end
    
    %% matriz de rigidez do elemento
    
    c = cos(beta);
    s = sin(beta);
    ke = E*A/L*[c*c c*s -c*c -c*s;
           c*s s*s -c*s -s*s;
           -c*c -c*s c*c c*s;
           -c*s -s*s c*s s*s];

     
    %% deslocamento da matriz local
    eldisp = zeros(4,1);

    nd1 = inci(i,4);
    nd2 = inci(i,5);

    eldisp(1) = d(nd1,1);
    eldisp(2) = d(nd1,2);
    eldisp(3) = d(nd2,1);
    eldisp(4) = d(nd2,2);
    
    elforce = ke * eldisp; %matriz força do elemento
     
    stress(i) = sqrt(elforce(1)^2+elforce(2)^2)/A;
    
    if ((x2-x1)*elforce(3))<0
        stress(i)=-stress(i);
    end
end 

%% calculo da tensão em cada barra da treliça

fprintf(' Tensão em cada barra \n')
disp(stress)


% xmin = min(coord(:,2));
% xmax = max(coord(:,2));
% ymin = min(coord(:,3));
% ymax = max(coord(:,3));
% Lx = xmax - xmin;
% Ly = ymax - ymin;
% xmin = xmin - 0.1*Lx;
% xmax = xmax + 0.1*Lx;
% ymin = ymin - 0.1*Ly;
% ymax = ymax + 0.1*Ly;



figure(1)
hold on
axis equal
plot(coord(:,2),coord(:,3),'b.')

L = 0;
for i=1:nel
    xp = [coord(inci(i,4),2) coord(inci(i,5),2)];
    yp = [coord(inci(i,4),3) coord(inci(i,5),3)];
    L1 = sqrt((xp(2)-xp(1))^2+(yp(2)-yp(1))^2);
    if L1>L
       L = L1; 
    end
    plot(xp,yp,'r')
    text((xp(2)+xp(1))/2,(yp(2)+yp(1))/2,num2str(inci(i,1)))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Final truss plot %%%%%%%%%%%


for i=1:nel
    xp = [coord(inci(i,4),2)+1e4*d(inci(i,4),1) coord(inci(i,5),2)+1e4*d(inci(i,5),1)];
    yp = [coord(inci(i,4),3)+1e4*d(inci(i,4),2) coord(inci(i,5),3)+1e4*d(inci(i,5),2)];
    plot(xp,yp,'--k')
end

title('Treliça')
ylabel('Deslocamento vertical')
xlabel('Deslocamento horizontal')


%% plotando as forcas

scale=0.0025;
for i=1:length(Load(:,1))
    if Load(i,2)==1
        quiver(coord(Load(i,1),2),coord(Load(i,1),3),Load(i,3)*scale,0,'LineWidth',2,'color','k');
        text(coord(Load(i,1),2)+Load(i,3)*scale/2,coord(Load(i,1),3)+ 0.3,num2str(Load(i,3)));
    else
        quiver(coord(Load(i,1),2),coord(Load(i,1),3),0,Load(i,3)*scale,'LineWidth',2,'color','k');
         text(coord(Load(i,1),2)+0.3,coord(Load(i,1),3)+Load(i,3)*scale/2,num2str(Load(i,3)));
    end
end 

%% plotando as BCs

for i=1:length(bc(:,1))
    if bc(i,2)==1
        quiver(coord(bc(i,1),2)-0.2,coord(bc(i,1),3),0,0,'Marker','>','MarkerSize',11,'MaxHeadSize',3,'color','k');
    else
        quiver(coord(bc(i,1),2),coord(bc(i,1),3)-0.2,0,0,'Marker','^','MarkerSize',11,'MaxHeadSize',3,'color','k');
    end
end 