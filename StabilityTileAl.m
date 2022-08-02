function [nt,nx,iT] = StabilityTileAl(tmax,xmax,tol,Al_thickness_ratio, graphData, tileData)

%% StabilityTileAl function analysis
% This function gives appropriate nx and nt values to meet a tolerance set
% by the user

% Inputs
% tmax = the total time that the model is copmputed cover
% xmax = the maximum thickness of the tile
% tol = the model tolerance to a degree temperature
% The ratio of aluminium to the total thickness of the tile
% graph Data = the temperature and time data used in the model
% tileData = the properties of the tile needed to model the tile

% Outputs
% nt = an nt that is appropriate to the tolerance set by the user
% nx = an nx that is appropriate to the tolerance set by the user
% iT = the position of the transitional point of the aluminium and silica
% fibre


% Nx chosen for all modelling in the model
nx = 101;
iT = nx - round(nx * Al_thickness_ratio);
% Constant nx values used for stability analysis

constant_nx = 41;
iT_constant = constant_nx - round(constant_nx * Al_thickness_ratio);
% Initialise nt stability analysis
counter = 101;
terminate = 0;
stabilityRange = 50; % range of nt results to check for stability
i = 0;

while terminate ~= 1
    i = i + 1;
    % finding the inside temperature for a given nt
    [~, ~, u,~] = shuttleNeumann2mat(tmax, counter,iT_constant, xmax, constant_nx, 'crank nicolson',graphData, tileData);
    insideTemp(i) = u(end,end);
    
    % checking if the range analysed is stable
    if rem(i, stabilityRange) == 0
        newRange = insideTemp(i-stabilityRange+1:i);
        var = insideTemp(i) - insideTemp(i-stabilityRange + 1);
        varRange = max(insideTemp(i - stabilityRange + 1:i)) - min(insideTemp(i - stabilityRange + 1:i));
        if varRange < tol
            nt = counter;
            terminate = 1;
        end
    end
    counter = counter + 1;
end
