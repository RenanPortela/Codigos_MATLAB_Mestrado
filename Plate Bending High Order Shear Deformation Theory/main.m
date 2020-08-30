%% CLEAN MEMORY
clear all; colordef white; clc; close all;

%% INPUT

load_type = 'sinusoidal'; % 'sinusoidal' // 'uniform'
typeBC = 'ssss';          % 'ssss' // 'cccc'
thickness = 0.01; 
material = 'isotropic';   % 'isotropic' // 'orthotropic' // 'laminate_3_layers' // 'laminate_4_layers'
nel = 16;

%% MESH GENERATION

mesh.nex = nel; 
mesh.ney = nel;
mesh.Lx = 1;   % length in X 
mesh.Ly = 1;   % length in Y
mesh.type = 'Q4';

[ mesh ] = create_mesh( mesh );

%% MATERIAL

switch material
    case 'isotropic'
        E = 10920;
        poisson = 0.3;
        
        
        C = [E/(1-poisson^2), E*poisson/(1-poisson^2), 0, 0, 0;...
             E*poisson/(1-poisson^2), E/(1-poisson^2), 0, 0, 0;...
             0, 0, E/2/(1+poisson),   0,   0;
             0, 0,   0, E/2/(1+poisson),   0;
             0, 0,   0,   0, E/2/(1+poisson)]; 
         
         geo = [1  thickness/2  -thickness/2  0];
         
    case 'orthotropic'
        E1 = 25e3; E2 = 1e3; 
        poisson12 = 0.25; poisson21 = poisson12*E2/E1;
        G12 = E2/2; G13 = E2/2; G23 = E2/5;
        
        C = [E1/(1-poisson12*poisson21), E1*poisson21/(1-poisson12*poisson21), 0, 0, 0;...
             E1*poisson21/(1-poisson12*poisson21), E2/(1-poisson12*poisson21), 0, 0, 0;...
             0, 0, G12,   0,   0;
             0, 0,   0, G13,   0;
             0, 0,   0,   0, G23]; 
         
         geo = [1  thickness/2  -thickness/2  0];
         
    case 'laminate_3_layers'
        E1 = 25e3; E2 = 1e3; 
        poisson12 = 0.25; poisson21 = poisson12*E2/E1;
        G12 = E2/2; G13 = E2/2; G23 = E2/5;
        
        C = [E1/(1-poisson12*poisson21), E1*poisson21/(1-poisson12*poisson21), 0, 0, 0;...
             E1*poisson21/(1-poisson12*poisson21), E2/(1-poisson12*poisson21), 0, 0, 0;...
             0, 0, G12,   0,   0;
             0, 0,   0, G13,   0;
             0, 0,   0,   0, G23]; 
         
        geo = [1  thickness/2  thickness/6  0;
               2  thickness/6 -thickness/6  pi/2;
               3 -thickness/6 -thickness/2  0];

    case 'laminate_4_layers'     
        E1 = 25e3; E2 = 1e3; 
        poisson12 = 0.25; poisson21 = poisson12*E2/E1;
        G12 = E2/2; G13 = E2/2; G23 = E2/5;
        
        C = [E1/(1-poisson12*poisson21), E1*poisson21/(1-poisson12*poisson21), 0, 0, 0;...
             E1*poisson21/(1-poisson12*poisson21), E2/(1-poisson12*poisson21), 0, 0, 0;...
             0, 0, G12,   0,   0;
             0, 0,   0, G13,   0;
             0, 0,   0,   0, G23];
         
         geo = [1  thickness/2  thickness/4  0;
                2  thickness/4            0  pi/2;
                3            0 -thickness/4  pi/2;
                4 -thickness/4 -thickness/2  0];
        
end

[KG] = laminate_stiffness_matrix_HSDT(C, geo, mesh);

%% LOAD 

% LOAD AMPLITUDE 
Load = 1e-3;

switch load_type
    case 'uniform'
        [FG] = load_vector_HSDT(mesh, Load);
    case 'sinusoidal'
        [FG] = load_sin_vector_HSDT(mesh, Load);
end

%% BOUNDARY CONDITIONS

dof.all = 1 : 7*mesh.nnodes;

[ fixedNode ] = boundary_condition_HSDT( mesh, typeBC );

[dof] = free_DOF_HSDT(dof, fixedNode, mesh);

KG_aux = KG(logical(dof.free), logical(dof.free));
FG_aux = FG(logical(dof.free), 1);


%% DISPLACEMENT

U = zeros(size(KG,1),1);
U(logical(dof.free),1) = sparse(KG_aux)\sparse(FG_aux);

figure('Name','Transversal displacement','NumberTitle','off');
plot3(mesh.coord(:,2),mesh.coord(:,3),U(2*mesh.nnodes+1:3*mesh.nnodes),'.') 
grid on; grid minor

xx = mesh.coord(:,2); 
node = find(xx == 0.5); 
xx = linspace(0,1,size(node,1));
displacement = U(2*mesh.nnodes+1:3*mesh.nnodes);
figure('Name','Profile x = 0','NumberTitle','off');
plot(xx, displacement(node))
grid on; grid minor;

format long 
fprintf('displacement_max = %4.6e \n', max(U(2*mesh.nnodes+1:3*mesh.nnodes)))

%% STRESS

Q  = C(1:3, 1:3); 
if strcmp(material,'laminate_3_layers') || strcmp(material,'laminate_4_layers')
    s = sin(geo(2,4));
    c = cos(geo(2,4));
    
    Tr = [c^2, s^2,    -2*s*c, 0, 0;
         s^2, c^2,     2*s*c, 0, 0;
         s*c,-s*c, c^2 - s^2, 0, 0;
           0,   0,         0, c,-s;
           0,   0,         0, s, c];
       
    Q_aux = Tr'*C*Tr;
    Qs = Q_aux(4:5, 4:5);
else
    Qs = C(4:5, 4:5);
end

[stress] = stress_HSDT(Q, Qs, mesh, thickness, U);
fprintf('Sxx = %4.4f\n',stress(1));
fprintf('Syy = %4.4f\n',stress(2));
fprintf('Sxy = %4.4f\n',stress(3));
fprintf('Sxz = %4.4f\n',stress(4));

e_n = FG'*U;
fprintf('||e|| = %4.5e\n',e_n)