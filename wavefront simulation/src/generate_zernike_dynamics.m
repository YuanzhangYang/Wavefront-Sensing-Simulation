function generate_zernike_dynamics()
% Script file: generate_zernike_dynamics.m
%
% Purpose:
%  This function is to imulate dynamic interferograms for selected Zernike aberrations
% (defocus, hcoma, vcoma, spherical, astigmatism) as their strength varies,
% and saves each sequence as a GIF animation.

% Record of revisions:
%     Date        Programmer          Description of change
%     ====        ==========          =====================
%   07/09/25      Yuanzhang Yang         Original code 
%
% Define variables:
%   wavelength    -- Wavelength of monochromatic light (meters)
%   k             -- Wave number (2π/λ)
%   D             -- Diameter of the aperture (meters)
%   N             -- Size of the simulation grid (NxN)
%   Delta         -- Fixed optical path difference (meters)
%   linear_coeff  -- Linear phase coefficient for wedge angle
%   x, y          -- 1D coordinate arrays for grid (meters)
%   X, Y          -- 2D meshgrid coordinates (meters)
%   rho           -- Normalized radial coordinate (0-1)
%   theta         -- Azimuthal angle (radians)
%   aperture      -- Circular aperture mask (logical)
%   phase_types   -- Cell array of aberration types
%   Z             -- Zernike polynomial value
%   delta_phi     -- Phase difference (radians)
%   I             -- Interference intensity

%% ----- Simulation Parameters -----
    wavelength = 632.8e-9;         % Wavelength of light (meters)
    k = 2 * pi / wavelength;       % Wave number
    D = 20e-3;                     % Diameter of the circular aperture (meters)
    N = 500;                       % Size of the simulation grid (NxN)
    Delta = 1e-6;                  % Fixed optical path difference (meters)

    % Generate spatial grid
    x = linspace(-D/2, D/2, N);    % 1D coordinate (x-axis)
    y = linspace(-D/2, D/2, N);    % 1D coordinate (y-axis)
    [X, Y] = meshgrid(x, y);      % 2D coordinate grid
    rho = sqrt(X.^2 + Y.^2) / (D/2);   % Normalized radial coordinate
    theta = atan2(Y, X);              % Azimuthal angle in polar coordinates
    aperture = rho <= 1;              % Circular aperture mask

    % Define aberration types to simulate
    aberrations = {'defocus', 'hcoma', 'vcoma', 'spherical', 'astig0', 'astig45'};

    % Coefficient values to simulate dynamically
    coef_range = linspace(-2, 2, 40); % Sweep range for Zernike coefficients

    %% ----- Loop Over Each Aberration Type -----
    for i = 1:length(aberrations)
        type = aberrations{i};                          % Current aberration type
        filename = sprintf('dynamic_%s.gif', type);     % Output GIF filename

        figure('Position', [100, 100, 600, 600]);

        for idx = 1:length(coef_range)
            coef = coef_range(idx);                     % Current Zernike coefficient

            % Generate Zernike phase map for current aberration and coefficient
            Z = getZernike(type, rho, theta, coef);

            % Compute phase difference and resulting interferogram intensity
            delta_phi = k * Delta * Z;
            I = 2 * (1 + cos(delta_phi));               % Interference intensity
            I(~aperture) = 0;                           % Apply aperture mask

            % Display interferogram
            imagesc(x*1e3, y*1e3, I);                   % Convert to mm for display
            axis image off;
            colormap(parula);
            clim([0 2]);                               % Set consistent color limits
            title(sprintf('%s  Coef: %.2f', upper(type), coef), ...
                  'FontSize', 14, 'FontWeight', 'bold');
            drawnow;

            % Convert frame to image for GIF
            frame = getframe(gcf);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im, 256);

            % Write frame to GIF
            if idx == 1
                imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.15);
            else
                imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.15);
            end
        end

        fprintf('Saved: %s\n\n', filename);
        close; % Close figure after saving
    end
end

%% ----- Zernike Polynomial Generator -----
function Z = getZernike(type, rho, theta, coef)
% getZernike
% ----------
% Returns the Zernike phase term for a given aberration type and coefficient.
%
% Inputs:
%   - type  : string specifying aberration type
%   - rho   : normalized radial coordinate
%   - theta : azimuthal angle
%   - coef  : scalar Zernike coefficient(dimensionless)
%
% Output:
%   - Z : 2D matrix of Zernike phase values

    switch lower(type)
        case 'defocus'
            Z = coef * (sqrt(3)*rho.^2 - 1);
        case 'hcoma'
            Z = coef * sqrt(8)*(3*rho.^3 - 2*rho).*cos(theta);
        case 'vcoma'
            Z = coef * sqrt(8)*(3*rho.^3 - 2*rho).*sin(theta);
        case 'spherical'
            Z = coef * sqrt(5)*(6*rho.^4 - 6*rho.^2 + 1);
        case 'astig0'
            Z = coef * sqrt(6)*rho.^2 .* cos(2*theta);
        case 'astig45'
            Z = coef * sqrt(6)*rho.^2 .* sin(2*theta);
        otherwise
            error('Unsupported aberration type: %s', type);
    end
end