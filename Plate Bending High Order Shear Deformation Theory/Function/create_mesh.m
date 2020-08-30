function [ mesh ] = create_mesh( mesh )
%PURPOSE: This function takes the number of elements, the length and the
%         type of element and returns the incidence matrix
%VARIABLES:
%         nex - number of elements in x direction
%         ney - number of elements in y direction
%         Lx - length in x direction
%         Ly - length in y direction
%         type - type of element ('Q4' or 'Q9')
%% -----------------------------------------------------------------------

switch mesh.type
    
    case 'Q4'
        mesh.nnx = mesh.nex + 1;  % number of nodes in x
        mesh.nny = mesh.ney + 1;  % number of nodes in y
        % coordinate matrix
        [ mesh ] = coordinate_matrix( mesh );
        % incidence matrix pre-location
        mesh.inci = zeros(mesh.nex*mesh.ney, 4); 
        % nodes reference pre-location
        node_aux = mesh.coord(:,1);    
        % reshaping matrix
        node_aux = reshape(node_aux,[mesh.nnx, mesh.nny]); 
        count = 1; % counter
        for j = 1 : mesh.ney
            for i = 1: mesh.nex
                % element-number first-node second-node third-node fourth-node
                mesh.inci(count,:) = [node_aux(i,j), node_aux(i+1,j)...
                                      node_aux(i+1,j+1), node_aux(i,j+1)];
                count = count + 1;
            end
        end
        
        mesh.nel = size(mesh.inci,1);               % number of elements
        mesh.nnodes = size(mesh.coord,1);           % number of nodes
        mesh.node_element = 4;                      % number of nodes by element
        
        dx = (mesh.coord(2,2) - mesh.coord(1,2))/10;
        dy = (mesh.coord(2,3) - mesh.coord(1,3))/10;
        
        figure('Name','Laminate','NumberTitle','off');
%         legend('Original layer')
        hold on
         for i = 1:size(mesh.inci,1)
             x = [mesh.coord(mesh.inci(i,1),2),mesh.coord(mesh.inci(i,2),2),mesh.coord(mesh.inci(i,3),2), mesh.coord(mesh.inci(i,4),2)];
             y = [mesh.coord(mesh.inci(i,1),3),mesh.coord(mesh.inci(i,2),3),mesh.coord(mesh.inci(i,3),3), mesh.coord(mesh.inci(i,4),3)];
             if mod(i,3) == 1
                patch(x,y,[1 1 1])
             elseif mod(i,3) == 2
                 patch(x,y,[0 1 1])
             else
                 patch(x,y,[1 1 0])
             end
         end
        scatter(mesh.coord(:,2), mesh.coord(:,3),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
        b = num2str(mesh.coord(:,1)); c = cellstr(b);
        text(mesh.coord(:,2,1) + dx, mesh.coord(:,3,1) + dy, c);
        grid on 
        grid minor
        view(3)
        
    case 'Q9'
        
        mesh.nnx = 2*mesh.nex + 1;  % number of nodes in x
        mesh.nny = 2*mesh.ney + 1;  % number of nodes in y 
        % coordinate matrix
        [ mesh ] = coordinate_matrix( mesh );
        % incidence matrix pre-location
        mesh.inci = zeros(mesh.nex*mesh.ney, 9); 
        % steps to write the incidence matrix
        step.x = 2; 
        step.y = 2 * mesh.nnx;
        step.pattern = [1 3 2*mesh.nnx+3 2*mesh.nnx+1 2 mesh.nnx+3 2*mesh.nnx+2 mesh.nnx+1 mesh.nnx+2];
        step.line = zeros(1,9);
        count = 1; % counter
        
        for i = 1 : mesh.ney
            for j = 1 : mesh.nex
                mesh.inci(count,:) = step.pattern + step.line;
                step.line = step.line + step.x;
                count = count + 1;
            end
            step.line = i*step.y;
        end
        
        mesh.nel = size(mesh.inci,1);               % number of elements
        mesh.nnodes = size(mesh.coord,1);           % number of nodes
        mesh.node_element = 9;                      % number of nodes by element
        
        figure('Name','Laminate','NumberTitle','off');
       
         for i = 1:size(mesh.inci,1)
             x = [mesh.coord(mesh.inci(i,1),2),mesh.coord(mesh.inci(i,2),2),mesh.coord(mesh.inci(i,3),2), mesh.coord(mesh.inci(i,4),2)];
             y = [mesh.coord(mesh.inci(i,1),3),mesh.coord(mesh.inci(i,2),3),mesh.coord(mesh.inci(i,3),3), mesh.coord(mesh.inci(i,4),3)];
             if mod(i,3) == 1
                patch(x,y,[1 1 1])
             elseif mod(i,3) == 2
                 patch(x,y,[0 1 1])
             else
                 patch(x,y,[1 1 0])
             end
         end
         hold on
        scatter(mesh.coord(:,2), mesh.coord(:,3), 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
        dx = (mesh.coord(2,2) - mesh.coord(1,2))/10;
        dy = (mesh.coord(2,3) - mesh.coord(1,3))/10;
        b = num2str(mesh.coord(:,1)); c = cellstr(b);
        text(mesh.coord(:,2,1) + dx, mesh.coord(:,3,1) + dy, c);
        view(3)         
        
end

end