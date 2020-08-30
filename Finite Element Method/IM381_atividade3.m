%% Program: Torção
%
%   @DESCRIÇÃO Atividade 3: Definir o modo de vibração de uma barra
%   livre-livre e engastada
%
%   @AUTOR Renan Miranda Portela @
%
%% Limpar memória e fechar janelas
clear all
close all
clc

%% Matriz de Coordenadas do sistema 
%   coord=[nº | X] => Matriz de coordenadas dos nós

node = 49; %número incial de nós
L = 1; %comprimento inicial

% while conv == 0

coord = zeros(node-1,2);

for i = 1:node
coord(i,1) = i; %número do nó
coord(i,2) = L/(node-1)*(i-1); %coordenada x
coord(i,3) = 0; %coordenada y
end


%-------Matriz de incidência ---------------------------------------------%
% inci=[nº | tab.mat | tab.geo | nó1 | nó2 ] => Matriz de
% incidência do elemento

inci = zeros(node-1,5); 

for i = 1:node-1
inci(i,1) = i; %número do elemento
inci(i,2) = 1; %material
inci(i,3) = 1; %geometria
inci(i,4) = i; %nó 1
inci(i,5) = i +1; %nó 2
end

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
%  coluna 1: Área m²
%  coluna 2: Densidade kg/m³

tabgeo = [0.1 7860];

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

% bc=[];
bc=[1 1 0];

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

    Load=[];

    %------------------Fim dos dados------------------------------------------%
%   tamanho dos vetores

nbc = size(bc,1); %numero de condições de contorno
nF = size(Load,1); %numero de forças externas
neq = 0; %número de equações
ngdl = 1; %número de graus de liberdade

%%  Montagem das matrizes globais de massa e rigidez 

%cria matriz id, travando inicialmente todos os nós
%em seguida destrava os graus de liberdade necessários

id=ones(1,node);

for i=1:nbc
    id(bc(i,2),bc(i,2))=bc(i,3);
end

for i= 1:node
    for j = 1:ngdl
        if id(j,i)== 1
            neq = neq +1;
            id(j,i) = neq;
        end
    end
end

kg = zeros(neq,neq); % pré locação da matriz de rigidez global
mg = zeros(neq,neq); % pré locação da matriz de massa global

for i = 1:node-1 %montagem das matrizes globais de massa e rigidez
    no1 = inci(i,4);
    no2 = inci(i,5);
    x1 = coord(no1,2);
    x2 = coord(no2,2);
    l = x2 - x1;
    mat = inci(i,2);
    geo = inci(i,3);
    E = tabmat(1,mat);
    a = tabgeo(1);
    rho = tabgeo(2);
    ke = 1/l*[1 -1; -1 1]; %matriz de rigidez do elemento
    me = 1*l/6*[2 1;1 2]; %matriz de rigidez do elemento
    loc = [id(1,no1),id(1,no2)]; 
    for j = 1:2
        if loc(j) ~= 0;
            for k = 1:2
                if loc(k) ~=0
                    kg(loc(j),loc(k))=kg(loc(j),loc(k))+ke(j,k);
                    mg(loc(j),loc(k))=mg(loc(j),loc(k))+me(j,k);
                end
            end
        end
    end
end

[theta,D]=eig(kg,mg); %theta = modo de vibração

for i = 1:size(theta,1)
    theta(:,i) = theta(:,i)/max(abs(theta(:,i)));    
end

lambda = diag(D); %autovalores [K] - lambda*[M] = 0 

poisson = tabmat(2);

G = E/(2*(1+poisson));
 
b = 1;%sqrt(G/rho);

omega = sqrt(lambda)*b; %frequência de vibração

coord_2 = coord; %coord_2 é a coordenada dos nós após a vibração

y_1 = zeros(node,1);

if isempty(bc)
    
    for i = 1 : node-1 %primeiro modo de vibração
        coord_2(i+1,3) = coord_2(i+1,3) + theta(i,2);
        y_1(i+1) = sin((pi*coord_2(i+1,2))/2);
    end
    
    figure(1)
    plot(coord_2(:,2),coord_2(:,3),'*-')
    
    title('Modos de vibração')
    xlabel('Comprimento da barra sob torção (m)') % x-axis label
    ylabel('Deslocamento vertical do nó') % x-axis label

    hold on
    
    coord_2 = coord; %coord_2 é a coordenada dos nós após a vibração

    if node >= 3

        y_2 = zeros(node,1);

        for i = 1 : node-1 %segundo modo de vibração
            coord_2(i+1,3) = coord_2(i+1,3) - theta(i,3);
            y_2(i+1) = sin((3*pi*coord_2(i+1,2))/2);
        end

        erro_2 = norm(abs(coord_2(:,3))-abs(y_2));

        plot(coord_2(:,2),coord_2(:,3),'k')

        coord_2 = coord; %coord_2 é a coordenada dos nós após a vibração
    end

    if node >= 4

        y_3 = zeros(node,1);

        for i = 1 : node-1 %terceiro modo de vibração
            coord_2(i+1,3) = coord_2(i+1,3) + theta(i,4);
            y_3(i+1) = sin((5*pi*coord_2(i+1,2))/2);
        end

        erro_3 = norm(abs(coord_2(:,3))-abs(y_3));

        plot(coord_2(:,2),coord_2(:,3),'r--')

        legend('Primeiro','Segundo','Terceiro')

        coord_2 = coord; %coord_2 é a coordenada dos nós após a vibração

    end
    
else

    for i = 1 : node-1 %primeiro modo de vibração
        coord_2(i+1,3) = coord_2(i+1,3) - theta(i,1);
        y_1(i+1) = sin((pi*coord_2(i+1,2))/2);
    end

    erro = norm(abs(coord_2(:,3))-abs(y_1));

    figure(1)
    plot(coord_2(:,2),coord_2(:,3),'*-')
    % hold on
    % plot(coord_2(:,2),y,'r*')

    title('Modos de vibração')
    xlabel('Comprimento da barra sob torção (m)') % x-axis label
    ylabel('Deslocamento vertical do nó') % x-axis label

    hold on


    coord_2 = coord; %coord_2 é a coordenada dos nós após a vibração

    if node >= 3

        y_2 = zeros(node,1);

        for i = 1 : node-1 %segundo modo de vibração
            coord_2(i+1,3) = coord_2(i+1,3) + theta(i,2);
            y_2(i+1) = sin((3*pi*coord_2(i+1,2))/2);
        end

        erro_2 = norm(abs(coord_2(:,3))-abs(y_2));

        plot(coord_2(:,2),coord_2(:,3),'k')

        coord_2 = coord; %coord_2 é a coordenada dos nós após a vibração
    end

    if node >= 4

        y_3 = zeros(node,1);

        for i = 1 : node-1 %terceiro modo de vibração
            coord_2(i+1,3) = coord_2(i+1,3) - theta(i,3);
            y_3(i+1) = sin((5*pi*coord_2(i+1,2))/2);
        end

        erro_3 = norm(abs(coord_2(:,3))-abs(y_3));

        plot(coord_2(:,2),coord_2(:,3),'r--')

        legend('Primeiro','Segundo','Terceiro')

        coord_2 = coord; %coord_2 é a coordenada dos nós após a vibração

    end
end

fprintf('\n\n******* Autovalores *******\n')
fprintf('       %f\n',lambda)
fprintf('\n\n******* Frequências de vibração *******\n')
fprintf('       %f\n',omega)
fprintf('\n\n******* 1º modo de vibração *******\n')
fprintf('       %f\n',theta(:,1))