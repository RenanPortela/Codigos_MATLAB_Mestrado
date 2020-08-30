function [KG] = laminate_stiffness_matrix_HSDT(C, geo, mesh)

A = zeros(3); B = zeros(3); D = zeros(3); E = zeros(3); F = zeros(3); 
H = zeros(3); G = zeros(2); S = zeros(2); T = zeros(2);

for i = 1 : size(geo,1)
    top = geo(i,2);
    bottom = geo(i,3);
    
    s = sin(geo(i,4));
    c = cos(geo(i,4));
    
    Tr = [c^2, s^2,    -2*s*c, 0, 0;
         s^2, c^2,     2*s*c, 0, 0;
         s*c,-s*c, c^2 - s^2, 0, 0;
           0,   0,         0, c,-s;
           0,   0,         0, s, c];
       
    Q_aux = Tr'*C*Tr;
    
    Q  = Q_aux(1:3, 1:3);
    
    Qs = Q_aux(4:5, 4:5);
    
    A = A + Q*(top-bottom);
    B = B + 1/2*Q*(top^2 - bottom^2);
    D = D + 1/3*Q*(top^3 - bottom^3);
    E = E + 1/4*Q*(top^4 - bottom^4);
    F = F + 1/5*Q*(top^5 - bottom^5);
    H = H + 1/7*Q*(top^7 - bottom^7);
    
    G = G + Qs*(top-bottom);
    S = S + 1/3*Qs*(top^3 - bottom^3);
    T = T + 1/5*Qs*(top^5 - bottom^5);
end

% GLOBAL STIFFNESS MATRIX PRE-ALOCATION
KG = zeros(7*mesh.nnodes);

thickness = geo(1,2) - geo(end,3);

c1 = - 4/3/thickness^2;
c2 = 3*c1;

% BENDING PART
switch mesh.type
    case 'Q4'
        [gaussWeights,gaussLocations] = gaussQuadrature('complete');
    case 'Q9'
        [gaussWeights,gaussLocations] = gaussQuadrature_Q9('complete');
end

% BENDING STIFFNESS MATRIX ASSEMBLING
for i = 1 : mesh.nel
    line = mesh.inci(i,:);
    posxy = mesh.coord(line, 2:3);
    loc = [line line + mesh.nnodes line + 2*mesh.nnodes line + 3*mesh.nnodes line + 4*mesh.nnodes line + 5*mesh.nnodes line + 6*mesh.nnodes];
    ndof = length(line);
    
    % INTEGRATION
    for j = 1 : size(gaussWeights,1)
        GaussPoint = gaussLocations(j,:);                                                     
        csi = GaussPoint(1);
        eta = GaussPoint(2);
        
        % SHAPE FUNCTIONS
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
        
        KG(loc, loc) = KG(loc, loc) + B_0'*A*B_0*gaussWeights(j)*det(J) + B_1'*D*B_1*gaussWeights(j)*det(J) + c1*B_2'*F*B_1*gaussWeights(j)*det(J) + c1*B_1'*F*B_2*gaussWeights(j)*det(J) + c1^2*B_2'*H*B_2*gaussWeights(j)*det(J);
    end
end

% SHEAR PART
switch mesh.type
    case 'Q4'
        [gaussWeights,gaussLocations] = gaussQuadrature('reduced');
    case 'Q9'
        [gaussWeights,gaussLocations] = gaussQuadrature_Q9('complete');
end

% SHEAR STIFFNESS MATRIX ASSEMBLING
for i = 1 : mesh.nel
    line = mesh.inci(i,:);
    posxy = mesh.coord(line, 2:3);
    loc = [line line + mesh.nnodes line + 2*mesh.nnodes line + 3*mesh.nnodes line + 4*mesh.nnodes line + 5*mesh.nnodes line + 6*mesh.nnodes];
    ndof = length(line);
    
    % INTEGRATION
    for j = 1 : size(gaussWeights,1)
        GaussPoint = gaussLocations(j,:);                                                     
        csi = GaussPoint(1);
        eta = GaussPoint(2);
        
        % SHAPE FUNCTIONS
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
        
        KG(loc, loc) = KG(loc, loc) + B_3'*G*B_3*gaussWeights(j)*det(J) + c2*B_4'*S*B_3*gaussWeights(j)*det(J) + c2*B_3'*S*B_4*gaussWeights(j)*det(J) + c2^2*B_4'*T*B_4*gaussWeights(j)*det(J);
    end
end

end