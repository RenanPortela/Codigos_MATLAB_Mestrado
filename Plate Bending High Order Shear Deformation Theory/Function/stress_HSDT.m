function [stress] = stress_HSDT(Q, Qs, mesh, thickness, U)

c1 = - 4/3/thickness^2;
c2 = 3*c1;

% IN PLANE STRESSES
i = (mesh.nex^2)/2 - mesh.nex/2;
line = mesh.inci(i,:);
posxy = mesh.coord(line, 2:3);
loc = [line line + mesh.nnodes line + 2*mesh.nnodes line + 3*mesh.nnodes line + 4*mesh.nnodes line + 5*mesh.nnodes line + 6*mesh.nnodes];
ndof = length(line);

csi = 1; eta = 1;
switch mesh.type
    case 'Q4'
        [ ~, naturalDerivatives ] = quad4( csi, eta );
    case 'Q9'
        [ ~, naturalDerivatives ] = quad9( csi, eta );
end
% JACOBIAN
J = naturalDerivatives*posxy;

% JACOBIAN'S INVERSE
iJ = inv(J);

B_x = iJ(1,1)*naturalDerivatives(1,:) + iJ(1,2)*naturalDerivatives(2,:);
B_y = iJ(2,1)*naturalDerivatives(1,:) + iJ(2,2)*naturalDerivatives(2,:);

% [B] MATRIX 
B_0 = zeros(3,7*ndof);
B_0(1, 1 : ndof) = B_x;
B_0(2, ndof+1:2*ndof) = B_y;
B_0(3, 1 : ndof) = B_y;
B_0(3, ndof + 1 : 2*ndof) = B_x;

B_1 = zeros(3,7*ndof);
B_1(1, 4*ndof+1:5*ndof) = B_x;
B_1(2, 3*ndof+1:4*ndof) =-B_y;
B_1(3, 3*ndof+1:4*ndof) =-B_x;
B_1(3, 4*ndof+1:5*ndof) = B_y;

B_2 = zeros(3,7*ndof);
B_2(1, 4*ndof+1:5*ndof) = B_x;
B_2(1, 5*ndof+1:6*ndof) = B_x;
B_2(2, 3*ndof+1:4*ndof) =-B_y;
B_2(2, 6*ndof+1:7*ndof) = B_y;
B_2(3, 3*ndof+1:4*ndof) =-B_x;
B_2(3, 4*ndof+1:5*ndof) = B_y;
B_2(3, 5*ndof+1:6*ndof) = B_y;
B_2(3, 6*ndof+1:7*ndof) = B_x;

z = thickness/2;

d = U(loc,1);

epsilon = (B_0 + z.*B_1 + c1.*z.^3.*B_2)*d;
sigma = Q*epsilon;
stress(1) = sigma(1);
stress(2) = sigma(2);

i = 1;
line = mesh.inci(i,:);
posxy = mesh.coord(line, 2:3);
loc = [line line + mesh.nnodes line + 2*mesh.nnodes line + 3*mesh.nnodes line + 4*mesh.nnodes line + 5*mesh.nnodes line + 6*mesh.nnodes];
ndof = length(line);

csi = -1; eta = -1;
switch mesh.type
    case 'Q4'
        [ ~, naturalDerivatives ] = quad4( csi, eta );
    case 'Q9'
        [ ~, naturalDerivatives ] = quad9( csi, eta );
end
% JACOBIAN
J = naturalDerivatives*posxy;

% JACOBIAN'S INVERSE
iJ = inv(J);

B_x = iJ(1,1)*naturalDerivatives(1,:) + iJ(1,2)*naturalDerivatives(2,:);
B_y = iJ(2,1)*naturalDerivatives(1,:) + iJ(2,2)*naturalDerivatives(2,:);

% [B] MATRIX 
B_0 = zeros(3,7*ndof);
B_0(1, 1 : ndof) = B_x;
B_0(2, ndof+1:2*ndof) = B_y;
B_0(3, 1 : ndof) = B_y;
B_0(3, ndof + 1 : 2*ndof) = B_x;

B_1 = zeros(3,7*ndof);
B_1(1, 4*ndof+1:5*ndof) = B_x;
B_1(2, 3*ndof+1:4*ndof) =-B_y;
B_1(3, 3*ndof+1:4*ndof) =-B_x;
B_1(3, 4*ndof+1:5*ndof) = B_y;

B_2 = zeros(3,7*ndof);
B_2(1, 4*ndof+1:5*ndof) = B_x;
B_2(1, 5*ndof+1:6*ndof) = B_x;
B_2(2, 3*ndof+1:4*ndof) =-B_y;
B_2(2, 6*ndof+1:7*ndof) = B_y;
B_2(3, 3*ndof+1:4*ndof) =-B_x;
B_2(3, 4*ndof+1:5*ndof) = B_y;
B_2(3, 5*ndof+1:6*ndof) = B_y;
B_2(3, 6*ndof+1:7*ndof) = B_x;

z = thickness/2;

d = U(loc,1);

epsilon = (B_0 + z.*B_1 + c1.*z.^3.*B_2)*d;
sigma = Q*epsilon;
stress(3) = sigma(3);

% OUT PLANE STRESSES
i = (mesh.nex/2 - 1)*mesh.nex + 1;
line = mesh.inci(i,:);
posxy = mesh.coord(line, 2:3);
loc = [line line + mesh.nnodes line + 2*mesh.nnodes line + 3*mesh.nnodes line + 4*mesh.nnodes line + 5*mesh.nnodes line + 6*mesh.nnodes];
ndof = length(line);

csi = -1; eta = 1;
switch mesh.type
    case 'Q4'
        [ ShapeFunction, naturalDerivatives ] = quad4( csi, eta );
    case 'Q9'
        [ ShapeFunction, naturalDerivatives ] = quad9( csi, eta );
end
% JACOBIAN
J = naturalDerivatives*posxy;

% JACOBIAN'S INVERSE
iJ = inv(J);

B_x = iJ(1,1)*naturalDerivatives(1,:) + iJ(1,2)*naturalDerivatives(2,:);
B_y = iJ(2,1)*naturalDerivatives(1,:) + iJ(2,2)*naturalDerivatives(2,:);

B_3 = zeros(2,7*ndof);
B_3(1, 2*ndof+1:3*ndof) = B_x;
B_3(1, 4*ndof+1:5*ndof) = ShapeFunction;
B_3(2, 2*ndof+1:3*ndof) = B_y;
B_3(2, 3*ndof+1:4*ndof) =-ShapeFunction;

B_4 = zeros(2,7*ndof);
B_4(1, 4*ndof+1:5*ndof) = ShapeFunction;
B_4(1, 5*ndof+1:6*ndof) = ShapeFunction;
B_4(2, 3*ndof+1:4*ndof) =-ShapeFunction;
B_4(2, 6*ndof+1:7*ndof) = ShapeFunction;

z = 0;

d = U(loc,1);

epsilon = (B_3 + c2*z^2*B_4)*d;
sigma = Qs*epsilon;

stress(4) = sigma(1);
stress(5) = sigma(2);

end