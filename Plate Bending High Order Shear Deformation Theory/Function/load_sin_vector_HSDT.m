function [FG] = load_sin_vector_HSDT(mesh, Load)

a = mesh.Lx;
b = mesh.Ly;

% LOAD VECTOR PRE-ALOCATION
FG = zeros(7*mesh.nnodes,1);

% GAUSSIAN QUADRATURE
switch mesh.type
    case 'Q4'
        [gaussWeights,gaussLocations] = gaussQuadrature('reduced');
    case 'Q9'
        [gaussWeights,gaussLocations] = gaussQuadrature_Q9('complete');
end

for i = 1 : mesh.nel
    % NODES IN EACH ELEMENT
    line = mesh.inci(i,:);
    posxy = mesh.coord(line, 2:3);
    coordx = posxy(:,1); coordy = posxy(:,2);
    
    for j = 1 : size(gaussWeights,1)
        GaussPoint = gaussLocations(j,:);                                                     
        csi = GaussPoint(1);
        eta = GaussPoint(2);       
        
        % SHAPE FUNCTIONS
        switch mesh.type
            case 'Q4'
                [ shapeFunction, naturalDerivatives ] = quad4( csi, eta );
            case 'Q9'
                [ shapeFunction, naturalDerivatives ] = quad9( csi, eta );
        end        
        
       % JACOBIAN
        J = naturalDerivatives*posxy; 
        
       % LOAD VECTOR
       FG(line + 2*mesh.nnodes) = FG(line + 2*mesh.nnodes) + shapeFunction'*Load*det(J)*gaussWeights(j).*sin(coordx*pi/a).*sin(coordy*pi/b);
    end
end

end