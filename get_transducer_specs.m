function [roc,diameters,frequency] = get_transducer_specs(transducer)
%GET_TRANSDUCER_SPECS Gets transducer radius of curvature and diameters
%   
% Usage:
%   [roc,diameters,frequency] = get_transducer_specs(transducer)
%   [roc,diameters,frequency] = get_transducer_specs('CTX500')
%   
% Inputs:
%   transducer: Options are 'CTX500', 'CTX560' or 'CTX250'
% 
% Author: Siti N. Yaakub 
arguments
    transducer  char
end

if strcmp(transducer, 'CTX500')
    roc = 63.2e-3;	% bowl radius of curvature [m]
    % aperture diameters of the elements given an inner, outer pairs [m]
    diameters = [0 1.28; 1.3 1.802; 1.82 2.19; 2.208 2.52].' * 0.0254;
    frequency = 500;
elseif strcmp(transducer, 'CTX560')
    roc = 63.2e-3;	% bowl radius of curvature [m]
    % aperture diameters of the elements given an inner, outer pairs [m]
    diameters = [0 1.28; 1.3 1.802; 1.82 2.19; 2.208 2.52].' * 0.0254;
    frequency = 560;
elseif strcmp(transducer, 'CTX250')
    roc = 63.2e-3;	% bowl radius of curvature [m]
    % aperture diameters of the elements given an inner, outer pairs [m]
    diameters = [0 1.274; 1.296 1.794; 1.816 2.186; 2.206 2.52].' * 0.0254;
    frequency = 250;
else
    error("Unknown transducer type! Options are 'CTX500' or 'CTX250'.")
end

end

