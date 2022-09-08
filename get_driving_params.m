function [pressure,phase] = get_driving_params(focus_depth,transducer,isppa)
%GET_DRIVING_PARAMS Gets transducer pressure and relative phase of each
%element for the Sonic Concepts CTX-series transducers.
%   Currently, this function only works for the CTX-500 transducer for
%   free-field acoustic simulation at 20 W/cm^2. You can edit the function
%   to use your own driving parameters .mat file.
%
% Usage:
%   [pressure,phase] = get_driving_params(focus_depth,transducer,isppa)
%   [pressure,phase] = get_driving_params(60,'CTX500',20)
%
% Inputs:
%   focus_depth:	Distance from transducer face to intended focus in mm, 
%                   rounded to the nearest integer.
%   transducer:     Options are 'CTX500' or 'CTX560'.
%   isppa:          Isppa in free-field in W/cm^2.
%
% Outputs:
%   pressure:   Source pressure in Pa.
%   phase:      4-element array of phases of each transducer element in 
%               degrees for the given focal depth.
% 
% Dependencies:
%   driving_params_Isppa20.mat
% 
% Author: Siti N. Yaakub, University of Plymouth, 7 Sep 2022

if strcmp(transducer,'CTX500')
    if isppa==20
        load('driving_params_Isppa20.mat', 'driving_params');
        idx = find(driving_params.dist == focus_depth);
        pressure = driving_params.amp{idx}; % source pressure [Pa]
        phase = driving_params.phase{idx,:}; % phase [rad]
    else
        warning(['No saved values for this Isppa, please provide your ' ...
            'own as optional Name-Value input arguments to the main ' ...
            'script if you intend to run the acoustic simulation.']);
    end
else
    warning(['No saved values for this transducer, please provide your ' ...
        'own as optional Name-Value input arguments to the main script ' ...
        'if you intend to run the acoustic simulation.']);
end

% pressure=91846;
% phase=[0	304.2	248.3	192.4];
end

