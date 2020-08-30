%% Program: Barra de seção variável
%
%   @DESCRIÇÃO Atividade 1: Definir o deslocamento horizontal de uma barra com seção variável@
%
%   @AUTOR Renan Miranda Portela @
%
%% Limpar memória e fechar janelas
clear all
close all
clc

%% Valores de entrada
tic
nel = 100; %numero de elementos
L = 1; %m

% ------ Matriz de Coordenadas do sistema ---------------------------------%
%   coord=[nº | X | Y | Z] => Matriz de coordenadas dos nós


coord = zeros(nel,2);

for i = 1:nel+1
coord(i,1) = i; %número do nó
coord(i,2) = (i-1)*L/nel; %coordenada x
end

%-------Matriz de incidência ---------------------------------------------%
% inci=[nº | tab.mat | tab.geo | nó1 | nó2 ] => Matriz de
% incidência do elemento

inci = zeros(nel,5); 

for i = 1:nel
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
%  linha 1: Momento de Inércia

tabgeo = 1 ;

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
    1 2 0];

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

    Load=[nel+1 1 1000];
    
%------------------Fim dos dados------------------------------------------%
%   tamanho dos vetores

nnos = size(coord,1); %numero de nós
nbc = size(bc,1); %numero de condições de contorno
nF = size(Load,1); %numero de forças externas
neq = 0; %número de equações
ngdl = 1; %número de graus de liberdade


%%  Calculando os deslocamentos 

%cria matriz id, travando inicialmente todos os nós
%em seguida destrava os graus de liberdade necessários

id=ones(1,nnos);

for i=1:nbc
    id(bc(i,2),bc(i,2))=bc(i,3);
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
l = L/nel; %tamanho do elemento

for i = 1:nel
    
    no1 = inci(i,4);
    no2 = inci(i,5);
    x1 = coord(no1,2);
    x2 = coord(no2,2);
    l = x2 - x1;
    mat = inci(i,2);
    geo = inci(i,3);
    E = tabmat(1,mat);
    a = area(x1,L);
    ke = E*a/l*[1 -1; -1 1]; %matriz de rigidez do elemento
    %disp(ke)
    loc = [id(1,no1),id(1,no2)]; 
    for j = 1:2
        if loc(j) ~= 0;
            for k = 1:2
                if loc(k) ~=0
                    kg(loc(j),loc(k))=kg(loc(j),loc(k))+ke(j,k);
                end
            end
        end
    end
end
%disp(kg)
%%       Força

F = zeros(neq,1); %pré locação da matriz coluna de forças 
nloads = size(Load,1); %número de carregamentos
for i=1:nloads
    F(id(Load(i,2),Load(i,1)),1) = Load(i,3); 
end

%% Deslocamento

u = kg\F;

%% Pós processamento
figure(2)
xmin = min(coord(:,2));
xmax = max(coord(:,2));
Lx = xmax - xmin;
xmin = xmin - 0.1*Lx;
xmax = xmax + 0.1*Lx;

axis([xmin xmax -1 1]);

desloc = zeros(nnos,1);

for i=1:nnos
    if id(1,i)==0
        desloc(i,1)=0;
    else
        desloc(i,1) = u(id(1,i));
    end
end

m = zeros(nnos,1);

%axis('off')
for i=1:nnos
    m(i) = desloc(i,1);
end

coord1 = zeros(nel,2);

for i = 1:nel+1
coord1(i,1) = i; %número do nó
coord1(i,2) = ((i-1)*L/nel)+m(i)*1e6; %coordenada x
end

y=zeros(nel+1,1);

scatter(coord(:,2),y)
hold on
scatter(coord1(:,2),y,'*')

%% Tensão em cada elemento de barra

sigma = zeros(nel,1);
deformation_matrix = zeros(nel,1);

for i = 1:nel
    node1 = inci(i,4);
    node2 = inci(i,5);
    l = L / nel ; % comprimento do elemento
    E = tabmat(1,1); %Módulo de Young/Elasticidade
    deformation_matrix(i,1) = (desloc(node2)-desloc(node1))/l; %Deformação de cada elemento de barra
    sigma(i,1) = E*(desloc(node2)-desloc(node1))/l; %Tensão em cada elemento de barra
end

%% Solução analítica

F = 1000; %força em N
L = 1; %comprimento em m
E = 210E9; %módulo de Young em N/m² ou Pa
A = 0.1; %área final da barra em m²
% u_ana deslocamento do nó em m

x = linspace(0,1,11);
u_ana = F*L/(2*E*A)*log(1-2*x/(3*L)); %solução analítica do deslocamento 

figure(3)
plot(x,abs(u_ana))

hold on
plot(coord(:,2),desloc)
xlabel('x  (m)')
ylabel('\epsilon  (m)')
%% Cálculo do erro no deslocamento

u_ana = zeros(nel+1,1);

for i = 1:nel+1
x = coord(i,2); %coordenada x
u_ana(i,1) = -F*L/(2*E*A)*log(1-2*x/(3*L));%valor exato do deslocamento
end

erro = norm(u_ana-desloc);
fprintf('\n\n***** Erro do cálculo do deslocamento *****\n')
fprintf('       Quantidade de elementos         Valor\n')
fprintf('                  %d                %d\n',round(nel),erro)
%% Cálculo das tensões em cada barra
figure(4)

sigma = zeros(nel,1);

for i = 1:nel
    l = coord(i+1,2)-coord(i,2);
    sigma(i,1) = E*(desloc(i+1,1)-desloc(i,1))/l;
end

coord1 = zeros(nel,1);

for i = 1:nel
    coord1(i,1) = coord(i+1,2);
end

plot(coord1(:,1),sigma,'r')
xlabel('x (m)')
ylabel('\sigma  (Pa)')

hold on
sigma_ana = zeros(nel,1);

x = linspace(0,1,nel);
for i = 1 : nel 
    l = coord(i+1,2)-coord(i,2);   
    sigma_ana(i,1) = E*(u_ana(i+1,1)-u_ana(i,1))/l; %solução analítica da tensão 
end
plot(x,sigma_ana,'--b')
%% Cálculo do erro na tensão axial

erro_sigma = norm(sigma_ana-sigma);
fprintf('\n\n***** Erro do cálculo da tensão *****\n')
fprintf('       Quantidade de elementos         Valor\n')
fprintf('                  %d                %d\n',round(nel),erro_sigma)
%% Imprimindo deslocamento do último nó

fprintf('\n\n***** Deslocamento do último nó *****\n')
fprintf('       Nó         Valor\n')
fprintf('       %d       %d\n',round(nel+1),u(nel,1))

toc