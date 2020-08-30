function [ N, B ] = quad4( csi, eta )
%PURPOSE: This function takes the points of integration and returns the
%         values of the shape function and its derivative at a given point
%VARIABLES:
%         csi - point of integration (csi)
%         eta - point of integration (eta)
%         N - shape function
%         B - shape function's derivative 
%% -----------------------------------------------------------------------

N = zeros(1,4); % shape function pre-allocating 
B = zeros(2,4); % shape function's derivative pre-allocating 

%% SHAPE FUNCTION

N(1) = (1-csi)*(1-eta)/4;
N(2) = (1+csi)*(1-eta)/4;
N(3) = (1+csi)*(1+eta)/4;
N(4) = (1-csi)*(1+eta)/4;

%% SHAPE FUNCTION'S DERIVATIVE

% DERIVATIVE IN CSI
B(1,1) = (eta-1)/4;
B(1,2) = (1-eta)/4;
B(1,3) = (eta+1)/4;
B(1,4) =-(1+eta)/4;

% DERIVATIVE IN ETA
B(2,1) = (csi-1)/4;
B(2,2) =-(1+csi)/4;
B(2,3) = (1+csi)/4;
B(2,4) = (1-csi)/4;

end