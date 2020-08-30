function [ N, B ] = quad9( xi, eta )
%PURPOSE: This function takes the points of integration and returns the
%         values of the shape function and its derivative at a given point
%VARIABLES:
%         csi - point of integration 
%         eta - point of integration 
%         N - shape function
%         B - shape function's derivative 
%% -----------------------------------------------------------------------

N = zeros(1,9); % shape function pre-allocating 
B = zeros(2,9); % shape function's derivative pre-allocating 

%% SHAPE FUCNTION
N(1) = (xi*eta*(xi-1)*(eta-1))/4;
N(2) = (xi*eta*(xi+1)*(eta-1))/4;
N(3) = (xi*eta*(xi+1)*(eta+1))/4;
N(4) = (xi*eta*(xi-1)*(eta+1))/4;
N(5) = (-2*eta*(xi+1)*(xi-1)*(eta-1))/4;
N(6) = (-2*xi*(xi+1)*(eta+1)*(eta-1))/4;
N(7) = (-2*eta*(xi+1)*(xi-1)*(eta+1))/4;
N(8) = (-2*xi*(xi-1)*(eta+1)*(eta-1))/4;
N(9) = (4*(xi+1)*(xi-1)*(eta+1)*(eta-1))/4;

%% SHAPE FUNCTION'S DERIVATIVE

% DERIVATIVE IN CSI
B(1,1) = (eta*(2*xi-1)*(eta-1))/4;
B(1,2) = eta*(2*xi+1)*(eta-1)/4;
B(1,3) = eta*(2*xi+1)*(eta+1)/4;
B(1,4) = eta*(2*xi-1)*(eta+1)/4;
B(1,5) = -4*xi*eta*(eta-1)/4;
B(1,6) = -2*(2*xi+1)*(eta+1)*(eta-1)/4;
B(1,7) = -4*xi*eta*(eta+1)/4;
B(1,8) = -2*(2*xi-1)*(eta+1)*(eta-1)/4;
B(1,9) = 8*xi*(eta^2-1)/4;

% DERIVATIVE IN ETA
B(2,1) = xi*(xi-1)*(2*eta-1)/4;
B(2,2) = xi*(xi+1)*(2*eta-1)/4;
B(2,3) = xi*(xi+1)*(2*eta+1)/4;
B(2,4) = xi*(xi-1)*(2*eta+1)/4;
B(2,5) = -2*(xi+1)*(xi-1)*(2*eta-1)/4;
B(2,6) = -4*xi*eta*(xi+1)/4;
B(2,7) = -2*(xi+1)*(xi-1)*(2*eta+1)/4;
B(2,8) = -4*xi*eta*(xi-1)/4;
B(2,9) = 8*eta*(xi^2-1)/4;

end