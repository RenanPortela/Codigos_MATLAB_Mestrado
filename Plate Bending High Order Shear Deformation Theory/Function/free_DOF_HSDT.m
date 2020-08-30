function [dof] = free_DOF_HSDT(dof, fixedNode, mesh)
%PURPOSE: This function takes all degrees of freedom vector, the
%         displacement approximation vector, the fixed nodes structure,
%         geometric structure, and the mesh structure. It returns the free
%         degrees of freedom, that is the degrees of freedom reference of
%         the nodes which are not fixed.
%VARIABLES:
%         alldof    - all degrees of freedom
%         freedof   - free degrees of freedom
%         fixedNode - fixed node structure
%         mesh      - mesh structure 
%% -----------------------------------------------------------------------

% free degrees of freedom vector
dof.free = dof.all;

for i = 1 : size(fixedNode.U,1)
    node = fixedNode.U(i);
    dof.free(node) = 0;
end

for i = 1 : size(fixedNode.V,1)
    node = fixedNode.V(i) + mesh.nnodes;
    dof.free(node) = 0;
end

for i = 1 : size(fixedNode.W,1)
    node = fixedNode.W(i) + mesh.nnodes*2;
    dof.free(node) = 0;
end

for i = 1 : size(fixedNode.phi_1,1)
    node = fixedNode.phi_1(i) + mesh.nnodes*3;
    dof.free(node) = 0;
end

for i = 1 : size(fixedNode.phi_2,1)
    node = fixedNode.phi_2(i) + mesh.nnodes*4;
    dof.free(node) = 0;
end

for i = 1 : size(fixedNode.dWdx,1)
    node = fixedNode.dWdx(i) + mesh.nnodes*5;
    dof.free(node) = 0;
end

for i = 1 : size(fixedNode.dWdy,1)
    node = fixedNode.dWdy(i) + mesh.nnodes*6;
    dof.free(node) = 0;
end

end