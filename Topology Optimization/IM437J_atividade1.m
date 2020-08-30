% UNIVERSIDADE ESTADUAL DE CAMPINAS
% FACULDADE DE ENGENHARIA MECANICA
% METODOS DE OTIMIZACAO TOPOLOGICA EVOLUCIONARIA - IM437 J
% ATIVIDADE 1
%
% DIPL. -ENG RENAN MIRANDA PORTELA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% cleaning

clc; clear; close all

%% coordinate matrix

nx = 33; % nodes number in x
ny = 17; % nodes number in y
Lx = 16; % length in x
Ly = 8; % length in y
ex = nx - 1; % element number in x
ey = ny - 1; % element number in y
a = Lx/ex; % element length
b = Ly/ey; % element hight

coordx = linspace(0,Lx,nx); % coordinate in x
coordy = linspace(0,Ly,ny); % coordinate in y
[X,Y] = meshgrid(coordx,coordy);
X = X';
Y = Y';
num = 1:nx*ny; % number of nodes
coord = [num(:),X(:),Y(:)]; % coordinate matrix

%% incidence matrix

inci = zeros(ex*ey,6); % incidence matrix pre-location
A = coord(:,1);
A = reshape(A,[nx,ny]);
k = 1;
for j = 1:ny-1
    for i = 1:nx-1
        if k > ex*ey/2
            mat = 1;
        else
            mat = 1;
        end
        inci(k,:)=[mat, A(i,j),A(i+1,j),A(i+1,j+1),A(i,j+1),k];
        k = k + 1;
    end
end

%% plot

figure()
patch('Faces',inci(:,2:5),'Vertices',coord(:,2:3),'FaceColor','blue')
axis equal
chr = int2str(coord(:,1));
text(coord(:,2),coord(:,3),chr)

%% boundary conditions matrix

bc(1,:) = [1,2,0];
bc(2,:) = [1,1,0];
bc(3,:) = [nx,2,0];
bc(4,:) = [nx,1,0];

%% load matrix

if mod(ex,2) == 0
    k = (nx+1)/2;
else
    k = nx/2;
end

load = [k 2 -10000];

figure()
scatter(coord(:,2),coord(:,3),'x')
axis equal
hold on
scatter(coord(k,2),coord(k,3),'o','red')
scatter(coord(bc(:,1),2),coord(bc(:,1),3),200,'^','red')
hold off

%% material matrix
            %E          nu
material = [210e9 7860 0.3;  %steel
            70e9 2700 0.27]; %aluminium
        
%% solver
        
nel = size(inci,1);    % element number
nnodes = size(coord,1);% nodes number

alldof = 1:nnodes*2; % degrees of freedom

kg = zeros(2*nnodes); % global stiffness matrix pre-location

F = zeros(2*nnodes,1); %load matrix pre-location 
nload = size(load,1); %load quantity

for i = 1:nload
    if load(i,2) == 1
        F(2*load(i,1)-1) = load(i,3);        
    else
        F(2*load(i,1)) = load(i,3);
    end    
end

for i = 1:nel
    no1 = inci(i,2); % first node element 
    no2 = inci(i,3); % second node element
    no3 = inci(i,4); % third node element
    no4 = inci(i,5); % fourth node element
    
    E = material(inci(i,1),1); % young's module
    
    nu = material(inci(i,1),3); % poisson module
    
    k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
   -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];

    ke = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                      k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                      k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                      k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                      k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                      k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                      k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                      k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
    
    loc = [no1*2-1 no1*2 no2*2-1 no2*2 no3*2-1 no3*2 no4*2-1 no4*2]; % localization vector
    
    kg(loc,loc) = kg(loc,loc) + ke; % global stiffness matrix assemble
end

freedof = alldof;

for k = 1 : size(bc,1)
    freedof(2*bc(k,1)-(2-bc(k,2))) = 0;
end

kg_aux = kg(logical(freedof),logical(freedof)); % column & rows elimination
F = F(logical(freedof),1); % column & rows elimination

u = kg_aux\F;

freedof_aux = freedof;
k = 1;

for i = 1:2*size(coord,1)
    if freedof(1,i) == 0
        freedof_aux(1,i) = 0;
    else
        freedof_aux(1,i) = k;
        k = k + 1;
    end
end

desloc = zeros(2*nnodes,1);

for i=1:2*nnodes
    if freedof_aux(1,i)==0
        desloc(i,1)=0;
    else
        desloc(i,1) = u(freedof_aux(1,i));
    end
end

u_aux = desloc;
desloc = reshape(desloc,2,[]);
desloc = desloc';

X = coord(:,2) + 1e7*desloc(:,1);
Y = coord(:,3) + 1e7*desloc(:,2);

coord_aux = [X(:),Y(:)];

%% stress

stress = zeros(nel,1);

Bt=[-0.5/a,0,0.5/a,0,0.5/a,0,-0.5/a,0;
    0,-0.5/b,0,-0.5/b,0,0.5/b,0,0.5/b;
    -0.5/b,-0.5/a,-0.5/b,0.5/a,0.5/b,0.5/a,0.5/b,-0.5/a];

for i = 1:nel
    no1 = inci(i,2); % first node element 
    no2 = inci(i,3); % second node element
    no3 = inci(i,4); % third node element
    no4 = inci(i,5); % fourth node element
    
    nu = material(inci(i,1),3); % poisson module
    D = E/(1-nu^2)*[1 nu 0;
                   nu 1 0;
                   0 0 (1-nu)/2];
    
    loc = [no1*2-1 no1*2 no2*2-1 no2*2 no3*2-1 no3*2 no4*2-1 no4*2]; % localization vector
    epi = u_aux(loc,1);
    
    T_el=D*Bt*epi; % stress within the element
    
    stress(i,1)=sqrt((T_el(1,1)^2)+(T_el(2,1)^2)-(T_el(1,1)*T_el(2,1))+3*(T_el(3,1)^2));
end

%% stress plot 

xi = zeros(nel,4);
yi = zeros(nel,4);

for i = 1:nel
    xi(i,:) = coord(inci(i,2:5),2);
    yi(i,:) = coord(inci(i,2:5),3);
end

figure()
patch(xi',yi',[stress';stress';stress';stress']);
colorbar
axis equal

%% plot

y_mean = zeros(nel,1);
x_mean = zeros(nel,1);

for i = 1 : nel
   y_mean(i,:) = mean(coord_aux(inci(i,2:5),2));
   x_mean(i,:) = mean(coord_aux(inci(i,2:5),1));
end

figure()
patch('Faces',inci(:,2:5),'Vertices',coord_aux(:,1:2),'FaceColor','none')
axis equal