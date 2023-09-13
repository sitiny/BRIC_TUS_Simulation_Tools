function [roc,diameters,frequency] = get_transducer_specs(transducer)
%GET_TRANSDUCER_SPECS Gets transducer radius of curvature and diameters
%
% Usage:
%   [roc,diameters,frequency] = get_transducer_specs(transducer)
%   [roc,diameters,frequency] = get_transducer_specs('CTX500')
%
% Inputs:
%   transducer: Options are 'CTX500' or 'CTX250'
%
% Outputs:
%   roc:        Transducer radius of curvature [m]
%   diameters:  Aperture diameters of the transducer elements given as inner,
%               outer pairs [m]
%   frequency:  Transducer central frequency [Hz]
%
% Author: Siti N. Yaakub, University of Plymouth, 7 Sep 2022
arguments
    transducer  char
end

switch transducer
    case 'CTX500'
        roc = 63.2e-3;	% bowl radius of curvature [m]
        % aperture diameters of the elements given an inner, outer pairs [m]
        diameters = [0 1.28; 1.3 1.802; 1.82 2.19; 2.208 2.52].' * 0.0254;
        frequency = 500e3;
    case 'CTX250'
        roc = 63.2e-3;	% bowl radius of curvature [m]
        % aperture diameters of the elements given an inner, outer pairs [m]
        diameters = [0 1.274; 1.296 1.794; 1.816 2.186; 2.206 2.52].' * 0.0254;
        frequency = 250e3;
    otherwise
        error("Unknown transducer type! Options are 'CTX500' or 'CTX250'.")
end

end
