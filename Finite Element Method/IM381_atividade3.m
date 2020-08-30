%% Program: Tor��o
%
%   @DESCRI��O Atividade 3: Definir o modo de vibra��o de uma barra
%   livre-livre e engastada
%
%   @AUTOR Renan Miranda Portela @
%
%% Limpar mem�ria e fechar janelas
clear all
close all
clc

%% Matriz de Coordenadas do sistema 
%   coord=[n� | X] => Matriz de coordenadas dos n�s

node = 49; %n�mero incial de n�s
L = 1; %comprimento inicial

% while conv == 0

coord = zeros(node-1,2);

for i = 1:node
coord(i,1) = i; %n�mero do n�
coord(i,2) = L/(node-1)*(i-1); %coordenada x
coord(i,3) = 0; %coordenada y
end


%-------Matriz de incid�ncia ---------------------------------------------%
% inci=[n� | tab.mat | tab.geo | n�1 | n�2 ] => Matriz de
% incid�ncia do elemento

inci = zeros(node-1,5); 

for i = 1:node-1
inci(i,1) = i; %n�mero do elemento
inci(i,2) = 1; %material
inci(i,3) = 1; %geometria
inci(i,4) = i; %n� 1
inci(i,5) = i +1; %n� 2
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
%  coluna 1: �rea m�
%  coluna 2: Densidade kg/m�

tabgeo = [0.1 7860];

% ------ condi��es de contorno ---------------------------------------------%
%
%   bc=[n� | grau de liberdade | valor]
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
%   F=[n� | grau | valor]
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

nbc = size(bc,1); %numero de condi��es de contorno
nF = size(Load,1); %numero de for�as externas
neq = 0; %n�mero de equa��es
ngdl = 1; %n�mero de graus de liberdade

%%  Montagem das matrizes globais de massa e rigidez 

%cria matriz id, travando inicialmente todos os n�s
%em seguida destrava os graus de liberdade necess�rios

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

kg = zeros(neq,neq); % pr� loca��o da matriz de rigidez global
mg = zeros(neq,neq); % pr� loca��o da matriz de massa global

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

[theta,D]=eig(kg,mg); %theta = modo de vibra��o

for i = 1:size(theta,1)
    theta(:,i) = theta(:,i)/max(abs(theta(:,i)));    
end

lambda = diag(D); %autovalores [K] - lambda*[M] = 0 

poisson = tabmat(2);

G = E/(2*(1+poisson));
 
b = 1;%sqrt(G/rho);

omega = sqrt(lambda)*b; %frequ�ncia de vibra��o

coord_2 = coord; %coord_2 � a coordenada dos n�s ap�s a vibra��o

y_1 = zeros(node,1);

if isempty(bc)
    
    for i = 1 : node-1 %primeiro modo de vibra��o
        coord_2(i+1,3) = coord_2(i+1,3) + theta(i,2);
        y_1(i+1) = sin((pi*coord_2(i+1,2))/2);
    end
    
    figure(1)
    plot(coord_2(:,2),coord_2(:,3),'*-')
    
    title('Modos de vibra��o')
    xlabel('Comprimento da barra sob tor��o (m)') % x-axis label
    ylabel('Deslocamento vertical do n�') % x-axis label

    hold on
    
    coord_2 = coord; %coord_2 � a coordenada dos n�s ap�s a vibra��o

    if node >= 3

        y_2 = zeros(node,1);

        for i = 1 : node-1 %segundo modo de vibra��o
            coord_2(i+1,3) = coord_2(i+1,3) - theta(i,3);
            y_2(i+1) = sin((3*pi*coord_2(i+1,2))/2);
        end

        erro_2 = norm(abs(coord_2(:,3))-abs(y_2));

        plot(coord_2(:,2),coord_2(:,3),'k')

        coord_2 = coord; %coord_2 � a coordenada dos n�s ap�s a vibra��o
    end

    if node >= 4

        y_3 = zeros(node,1);

        for i = 1 : node-1 %terceiro modo de vibra��o
            coord_2(i+1,3) = coord_2(i+1,3) + theta(i,4);
            y_3(i+1) = sin((5*pi*coord_2(i+1,2))/2);
        end

        erro_3 = norm(abs(coord_2(:,3))-abs(y_3));

        plot(coord_2(:,2),coord_2(:,3),'r--')

        legend('Primeiro','Segundo','Terceiro')

        coord_2 = coord; %coord_2 � a coordenada dos n�s ap�s a vibra��o

    end
    
else

    for i = 1 : node-1 %primeiro modo de vibra��o
        coord_2(i+1,3) = coord_2(i+1,3) - theta(i,1);
        y_1(i+1) = sin((pi*coord_2(i+1,2))/2);
    end

    erro = norm(abs(coord_2(:,3))-abs(y_1));

    figure(1)
    plot(coord_2(:,2),coord_2(:,3),'*-')
    % hold on
    % plot(coord_2(:,2),y,'r*')

    title('Modos de vibra��o')
    xlabel('Comprimento da barra sob tor��o (m)') % x-axis label
    ylabel('Deslocamento vertical do n�') % x-axis label

    hold on


    coord_2 = coord; %coord_2 � a coordenada dos n�s ap�s a vibra��o

    if node >= 3

        y_2 = zeros(node,1);

        for i = 1 : node-1 %segundo modo de vibra��o
            coord_2(i+1,3) = coord_2(i+1,3) + theta(i,2);
            y_2(i+1) = sin((3*pi*coord_2(i+1,2))/2);
        end

        erro_2 = norm(abs(coord_2(:,3))-abs(y_2));

        plot(coord_2(:,2),coord_2(:,3),'k')

        coord_2 = coord; %coord_2 � a coordenada dos n�s ap�s a vibra��o
    end

    if node >= 4

        y_3 = zeros(node,1);

        for i = 1 : node-1 %terceiro modo de vibra��o
            coord_2(i+1,3) = coord_2(i+1,3) - theta(i,3);
            y_3(i+1) = sin((5*pi*coord_2(i+1,2))/2);
        end

        erro_3 = norm(abs(coord_2(:,3))-abs(y_3));

        plot(coord_2(:,2),coord_2(:,3),'r--')

        legend('Primeiro','Segundo','Terceiro')

        coord_2 = coord; %coord_2 � a coordenada dos n�s ap�s a vibra��o

    end
end

fprintf('\n\n******* Autovalores *******\n')
fprintf('       %f\n',lambda)
fprintf('\n\n******* Frequ�ncias de vibra��o *******\n')
fprintf('       %f\n',omega)
fprintf('\n\n******* 1� modo de vibra��o *******\n')
fprintf('       %f\n',theta(:,1))