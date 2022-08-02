function [x, t, u, L] = shuttleNeumann2mat(tmax, nt, iT, xmax, nx, method, graphData, tileData)
% Function for modelling temperature in a space shuttle tile
% Modelling Group 12
%
% Input arguments:
% tmax   - maximum time (s)
% nt     - number of timesteps
% xmax   - total thickness (m)
% nx     - number of spatial steps
% method - solution method ('forward', 'backward' etc)
% timeData - time vector for surface temperatures (s)
% tempData - surface temperature vector (C or K)
%
% Return arguments:
% x      - distance vector (m)
% t      - time vector (s)
% u      - temperature matrix (C or K or F)
% L      - outside temperature data


% Set tile properties of silica fiber foam
thermConTile = tileData (1); % W/(m K)densityTile  = 144;   % kg/m3 ~ 9 lb/ft^3
specHeatTile = tileData (2);  % J/kgK ~ 0.3 Btu/lb/F at 500F
density = tileData (3);

%Set aluminium NASA 398 properties
thermConAl = 128; % W/(m K)
densityAl =  2760; %kg/m3
specHeatAl = 915; % J/kgK

%Initialise everything.
dt = tmax / (nt-1);
t = (0:nt-1) * dt;
dx = xmax / (nx-1);
x = (0:nx-1) * dx;
u = zeros(nt, nx);
alphaTile = thermConTile /(density * specHeatTile);
alphaAl = thermConAl / (densityAl * specHeatAl);
pTile = alphaTile * dt / dx^2;
pAl = alphaAl * dt / dx^2;

% vector indexing
ivec1 = 2:(iT-1);
ivec2 = (iT+1):nx-1;

% For the transition point, use the average density and specific heat
pT1 = thermConTile * 2 / (density*specHeatTile+densityAl*specHeatAl) * dt /(dx^2);
pT2 = thermConAl * 2 / (density*specHeatTile+densityAl*specHeatAl) * dt /(dx^2);

% Use interpolation to get outside temperature at time vector t
% and store it as left-hand boundary vector L.
L = interp1(graphData(2,:), graphData(1,:), t, "linear", "extrap");

% set initial conditions equal to boundary temperature at t=0.
u(1, :) = L(1);

% boundary conditions

u(:,1) = L;

% Select method and run simulation.
switch method
    case 'forward'
        for n=1:nt-1
            % left hand half: material 1
            u(n+1,ivec1) = (1 - 2*pTile) * u(n,ivec1) + pTile * (u(n,ivec1-1) + u(n,ivec1+1));

            % right-hand half: material 2
            u(n+1,ivec2) = (1 - 2*pAl) * u(n,ivec2) + pAl * (u(n,ivec2-1) + u(n,ivec2+1));

            % transition point
            u(n+1,iT) = u(n,iT) + pT2*(u(n,iT+1) - u(n,iT)) - pT1*(u(n,iT) - u(n,iT-1));
        end

    case 'dufort'
        for n=1:nt-1
            % set index for 'old' point
            if n == 1
                nm = 1; % at first timestep, old point doesn't exist as n-1 = 0.
                % Use value at timestep 1 instead.
            else
                nm = n-1; % after first timestep, proceed normally.
            end

            % left hand half: material 1
            u(n+1,ivec1) = ((1 - 2*pTile) * u(nm,ivec1) + ...
                2 * pTile * (u(n,ivec1-1) + u(n,ivec1+1))) / (1 + 2*pTile);

            % transition point
            u(n+1,iT) = ((1 - pT1 - pT2) * u(nm,iT) + ...
                2 * pT1 * u(n,iT-1) + 2 * pT2 * u(n,iT+1)) / (1 + pT1 + pT2);

            % right-hand half: material 2
            u(n+1,ivec2) = ((1 - 2*pAl) * u(nm,ivec2) + ...
                2 * pAl * (u(n,ivec2-1) + u(n,ivec2+1))) / (1 + 2*pAl);
        end

    case 'backward'
        for n=1:nt-1
            %left Boundary
            b(1) = 1;
            c(1) = 0;
            d(1) = L(n+1);

            % left hand half: material 1
            a(ivec1) = -pTile;
            b(ivec1) = 1 + 2*pTile;
            c(ivec1) = -pTile;
            d(ivec1) = u(n,ivec1);

            % transition point
            a(iT) = - pT1;
            b(iT) = 1 + pT1 + pT2;
            c(iT) = - pT2;
            d(iT) = u(n,iT);

            % right-hand half: material 2
            a(ivec2) = -pAl;
            b(ivec2) = 1 + 2*pAl;
            c(ivec2) = -pAl;
            d(ivec2) = u(n,ivec2);

            %Right Boundary
            a(nx) = -2*pAl;
            b(nx) = 1 + 2*pAl;
            d(nx) = u(n,nx);

            u(n+1,:) = tdm(a,b,c,d);
        end
    case 'crank nicolson'
        for n=1:nt-1
            %left boundary consition
            b(1) = 1;
            c(1) = 0;
            d(1) = L(n+1);

            % left hand half: material 1
            a(ivec1) = -pTile/2;
            b(ivec1) = 1 + pTile;
            c(ivec1) = -pTile/2;
            d(ivec1) = u(n,ivec1);

            % transition point
            a(iT) = - pT1/2;
            b(iT) = 1 + pT1/2 + pT2/2;
            c(iT) = - pT2/2;
            d(iT) = (1 - pT1/2 - pT2/2)*u(n,iT) + (pT1/2)*u(n,iT-1) + (pT2/2)*u(n,iT+1);

            % right-hand half: material 2
            a(ivec2) = -pAl/2;
            b(ivec2) = 1 + pAl;
            c(ivec2) = -pAl/2;
            d(ivec2) = u(n,ivec2);

            %right boundary points
            a(nx) = -2*pAl;
            b(nx) = 1 + 2*pAl;
            d(nx) = u(n,nx);

            u(n+1,:) = tdm(a,b,c,d);
        end
    otherwise
        error (['Undefined method: ' method])

end
% Tri-diagonal matrix solution
function x = tdm(a,b,c,d)
n = length(b);

% Eliminate a terms
for i = 2:n
    factor = a(i) / b(i-1);
    b(i) = b(i) - factor * c(i-1);
    d(i) = d(i) - factor * d(i-1);
end

x(n) = d(n) / b(n);

% Loop backwards to find other x values by back-substitution
for i = n-1:-1:1
    x(i) = (d(i) - c(i) * x(i+1)) / b(i);
end

