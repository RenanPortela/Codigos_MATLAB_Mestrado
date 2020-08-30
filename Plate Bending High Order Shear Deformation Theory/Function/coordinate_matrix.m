function [ mesh ] = coordinate_matrix( mesh )
%PURPOSE: This function takes the number of nodes in x and y direction and
%         returns the coordinate matrix
%VARIABLES:
%         nnx - number of nodes in x direction
%         nny - number of nodes in y direction
%         Lx - length in x direction
%         Ly - length in y direction
%% -----------------------------------------------------------------------

if ( nargin > 1 )
    disp('Too many parameters specified for this function')
end

coord_x = linspace(0, mesh.Lx, mesh.nnx); % coordinate of the nodes in x direction
coord_y = linspace(0, mesh.Ly, mesh.nny); % coordinate of the nodes in y direction
[x,y] = meshgrid(coord_x, coord_y); x = x'; y = y'; num = 1:[mesh.nnx]*[mesh.nny];        % number of nodes
mesh.coord = [num(:), x(:), y(:)];    % coordinate matrix

end