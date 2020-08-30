function [w_col,locations_cols]=gaussQuadrature_Q9(option)
% Gauss quadrature for Q4 elements
% option 'complete' (2x2)
% option 'reduced'  (1x1)
% locations: Gauss point locations
% weights: Gauss point weights

switch option
    case 'complete'
        locations_cols=[ -sqrt(3/5) -sqrt(3/5)
                          sqrt(3/5) -sqrt(3/5)
                          sqrt(3/5)  sqrt(3/5)
                         -sqrt(3/5)  sqrt(3/5)
                             0      -sqrt(3/5)
                          sqrt(3/5)     0
                             0       sqrt(3/5)
                         -sqrt(3/5)     0
                             0          0     ];
                    
            w_col=[0.555555555555556*0.555555555555556; 0.555555555555556*0.888888888888889; 0.555555555555556*0.555555555555556; 0.555555555555556*0.888888888888889; 0.888888888888889*0.888888888888889; 0.555555555555556*0.555555555555556; 0.555555555555556*0.555555555555556; 0.555555555555556*0.888888888888889; 0.555555555555556*0.555555555555556 ];
    case 'reduced'
        locations_cols=[-0.577350269189626 -0.577350269189626
                        0.577350269189626 -0.577350269189626
                        0.577350269189626  0.577350269189626
                        -0.577350269189626  0.577350269189626];
        w_col=[ 1 1 1 1]; 
end