function [xmax,x_out,t_out,u_out,iT] = ThicknessFinder2mat(tempMax,tmax,tol, graphData, tileData)

%% ThicknessFinder function description
% This function computes the required thickness for a required transition layer maxiumum
%temperature

% Inputs
% tempMax = Maximum temperature transition layer temperature
% tmax = The maxium time the model computes over
% tol = tolerance of the model ( +-tol unit temperature)
% graphData = The time  and temperature data
% TileData = All of the tile properties needed to run the model

% Outputs
% xmax = the thickness required to have a transition layer temperature
% lower than tempMax
% x_out = the nx used for the last guess
% t_out = the nt used for the last guess
% iT = the transition point in nx

% Initialise the counter for the shooting method
n = 2;
tolerance_thick = 0.001; %thickness accuracy of +-1mm
Al_thickness_ratio = 0.1; %10% of total thickness

% A first guess of the thickness
xmax_guess_nminus1 = 0.01;


% An appropriate time step and spacial step is found for the parameters of
% the first guess
[nt,nx,iT] = StabilityTileAl(tmax, xmax_guess_nminus1, tol, Al_thickness_ratio, graphData, tileData);

% The maximum outside temperature is found for the first thickness guess
[~, ~, u] = shuttleNeumann2mat(tmax, nt,iT, xmax_guess_nminus1, nx,'crank nicolson', graphData, tileData);
Inner_Temperature_nminus1 = max(u(:,iT));

% The error in temperature for the first guess is found
error_nminus1 = Inner_Temperature_nminus1 - tempMax;

% A second guess of the thickness
xmax_guess_n = 0.06;

% The maximum outside temperature is found for the second thickness guess
[~, ~, u] = shuttleNeumann2mat(tmax, nt,iT, xmax_guess_n, nx,'crank nicolson', graphData, tileData);
Inner_Temperature_n = max(u(:,iT));


% The error in the temperature for the second guess is found
error_n = Inner_Temperature_n - tempMax;


% Sets the error margin that terminates the shooting method
error_margin = tolerance_thick;
error_nplus1 = error_margin + error_margin; %initialise variable larger than error_margin

while abs(error_nplus1) > error_margin

    xmax_guess_nplus1 = abs(xmax_guess_n - error_n * ((xmax_guess_n - xmax_guess_nminus1)/(error_n - error_nminus1)));
    [x_out, t_out, u_out] = shuttleNeumann2mat(tmax, nt,iT,xmax_guess_nplus1, nx,'crank nicolson', graphData, tileData);
    Inner_Temperature_nplus1 = max(u_out(:,iT));

    %new error
    error_nplus1 = Inner_Temperature_nplus1 - tempMax;

    %reset variables to find the next guess
    xmax_guess_nminus1 = xmax_guess_n;
    xmax_guess_n = xmax_guess_nplus1;
    error_nminus1 = error_n;
    error_n = error_nplus1;
    n = n + 1;
    
    % if the program runs for too long terminate the while loop as the
    % solution is not possible
    if n > 2200
        xmax = 0;
        return
    end

end

% The final thickness is the final guess
xmax = xmax_guess_nplus1;

end
