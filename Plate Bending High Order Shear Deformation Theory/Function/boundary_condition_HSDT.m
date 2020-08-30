function [ fixedNode ] = boundary_condition_HSDT( mesh, typeBC )
%PURPOSE: This function takes the nodes from the mesh, and the type of
%         boundary condition. It returns the fixed nodes in different
%         direction.
%VARIABLES:
%         fixedNode - fixed node structure
%         mesh      - mesh structure 
%         typeBC    - type of boundary conditions
%% -----------------------------------------------------------------------

% coordinate in x direction 
xx = mesh.coord(:,2);

% coordinate in y direction 
yy = mesh.coord(:,3);
  
fixedNode.BC = typeBC;

switch typeBC
    case 'cccc' % clamped laminate
        fixedNode.W =find(xx==max(mesh.coord(:,2))|...   % fixed node for a displacement in z direction                
                          xx==min(mesh.coord(:,2))|...
                          yy==min(mesh.coord(:,3))|...
                          yy==max(mesh.coord(:,3)));
       fixedNode.U     = fixedNode.W;  % fixed node for a displacement in x direction
       fixedNode.V     = fixedNode.W;% fixed node for a displacement in y direction
       fixedNode.phi_1 = fixedNode.W;  % fixed node for a displacement in y direction
       fixedNode.phi_2 = fixedNode.W;  % fixed node for a displacement in y direction                                  
       fixedNode.dWdx  = fixedNode.W;  % fixed node for a displacement in y direction                   
       fixedNode.dWdy  = fixedNode.W;  % fixed node for a displacement in y direction
                                  
        
    case 'fccc' % one free and three clamped
        fixedNode.W =find(xx==max(mesh.coord(:,2))|...   % fixed node for a displacement in z direction                
                          yy==min(mesh.coord(:,3))|...
                          yy==max(mesh.coord(:,3)));
        fixedNode.U = fixedNode.W;                       % fixed node for a displacement in x direction
        fixedNode.V = fixedNode.W;                       % fixed node for a displacement in y direction
    
    case 'ssss' % simply supported edges
        fixedNode.W =find(xx==max(mesh.coord(:,2))|...   % fixed node for a displacement in z direction                
                          xx==min(mesh.coord(:,2))|...
                          yy==min(mesh.coord(:,3))|...
                          yy==max(mesh.coord(:,3)));
        fixedNode.U = find(yy==max(mesh.coord(:,2))|...  % fixed node for a displacement in x direction
                           yy==min(mesh.coord(:,2)));
        fixedNode.V = find(xx==max(mesh.coord(:,2))|...  % fixed node for a displacement in y direction
                           xx==min(mesh.coord(:,2)));       
       fixedNode.phi_1 = find(xx==max(mesh.coord(:,2))|...  % fixed node for a displacement in y direction
                           xx==min(mesh.coord(:,2)));                   
       fixedNode.phi_2 = find(yy==max(mesh.coord(:,2))|...  % fixed node for a displacement in y direction
                           yy==min(mesh.coord(:,2)));                                    
       fixedNode.dWdx = find(yy==max(mesh.coord(:,2))|...  % fixed node for a displacement in y direction
                           yy==min(mesh.coord(:,2)));                   
       fixedNode.dWdy = find(xx==max(mesh.coord(:,2))|...  % fixed node for a displacement in y direction
                           xx==min(mesh.coord(:,2)));                               
    case 'fsss'
        fixedNode.W =find(xx==max(mesh.coord(:,2))|...   % fixed node for a displacement in z direction                
                          yy==min(mesh.coord(:,3))|...
                          yy==max(mesh.coord(:,3)));
        fixedNode.U = find(xx==max(mesh.coord(:,2)));    % fixed node for a displacement in x direction
        fixedNode.V = find(yy==max(mesh.coord(:,2))|...  % fixed node for a displacement in y direction
                           yy==min(mesh.coord(:,2)));                          
    case 'fsfs'
        fixedNode.W =find(yy==min(mesh.coord(:,3))|...   % fixed node for a displacement in z direction                
                          yy==max(mesh.coord(:,3)));
        fixedNode.U = [];                                % fixed node for a displacement in x direction
        fixedNode.V = find(yy==max(mesh.coord(:,2))|...  % fixed node for a displacement in y direction
                           yy==min(mesh.coord(:,2)));  
    case 'fcfc'
        fixedNode.W =find(yy==min(mesh.coord(:,3))|...   % fixed node for a displacement in z direction                
                          yy==max(mesh.coord(:,3)));
        fixedNode.U = fixedNode.W;                       % fixed node for a displacement in x direction                
        fixedNode.V = fixedNode.W;                       % fixed node for a displacement in y direction                  
end


end