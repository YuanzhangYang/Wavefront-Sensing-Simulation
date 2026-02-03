% Script file: interferogram_simulation.m
%
% Purpose:
%   This script simulates and visualizes interferograms for various optical aberrations
%   and plane wavefronts, with and without a slight wedge angle, using a circular
%   aperture and Zernike polynomials to model phase differences.
%
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
%

% Parameter settings
wavelength = 632.8e-9;  % Wavelength of monochromatic light (632.8 nm)
k = 2*pi/wavelength;    % Wave number (k = 2π/λ)
D = 20e-3;              % Diameter of the aperture (20 mm)
N = 500;                % Size of the simulation grid (500x500)
Delta = 1e-6;           % Fixed optical path difference (1 µm)
linear_coeff = 1e-3;    % Linear phase coefficient to introduce a small wedge angle

% Define the grid range
x = linspace(-D/2, D/2, N);  % Create linear array for x-axis from -D/2 to D/2
y = linspace(-D/2, D/2, N);  % Create linear array for y-axis from -D/2 to D/2
[X, Y] = meshgrid(x, y);     % Generate 2D grid coordinates

% Calculate polar coordinates
rho = sqrt(X.^2 + Y.^2) / (D/2);  % Normalized radius (ρ = sqrt(x² + y²) / (D/2))
theta = atan2(Y, X);              % Angular coordinate (θ = atan2(y, x))

% Apply circular aperture mask
aperture = rho <= 1;  % Create circular aperture mask (1 inside radius, 0 outside)

% Array of phase aberration types
phase_types = {'defocus', 'hcoma', 'vcoma', 'spherical', 'astig0', 'astig45', 'plane_wave', 'plane_with_linear'};

% Create figure layout for plotting
figure('Position', [100, 100, 1600, 600]); % Set figure window size and position
t = tiledlayout(2, 8, 'TileSpacing', 'compact', 'Padding', 'compact'); % Create 2x8 tiled layout

% Plot interferograms for perfect plane wavefront with eight phase types
for i = 1:8
    nexttile; % Move to next tile in layout
    type = phase_types{i}; % Select current phase type
    switch lower(type) % Define Zernike polynomial or phase for each type
        case 'defocus'
            Z = sqrt(3)*rho.^2 - 1; % Defocus aberration (Zernike polynomial)
        case 'hcoma'
            Z = sqrt(8) * (3 * rho.^3 - 2 * rho) .* cos(theta); % Horizontal coma
        case 'vcoma'
            Z = sqrt(8) * (3 * rho.^3 - 2 * rho) .* sin(theta); % Vertical coma
        case 'spherical'
            Z = sqrt(5)*(6*rho.^4 - 6*rho.^2 + 1); % Spherical aberration
        case 'astig0'
            Z = sqrt(6)*rho.^2 .* cos(2*theta); % Astigmatism at 0 degrees
        case 'astig45'
            Z = sqrt(6)*rho.^2 .* sin(2*theta); % Astigmatism at 45 degrees
        case 'plane_wave'
            Z = 0; % Perfect plane wave (no aberration)
        case 'plane_with_linear'
            Z = 0; % Plane wave with linear phase (wedge angle)
    end
    % Calculate phase difference
    if strcmpi(type, 'plane_wave')
        delta_phi = 0; % No phase difference for perfect plane wave
    elseif strcmpi(type, 'plane_with_linear')
        delta_phi = k * linear_coeff * X; % Linear phase for wedge angle
    else
        delta_phi = k * Delta * Z; % Phase difference for aberrations
    end
    I = 2*(1 + cos(delta_phi)); % Interference intensity (I = 2(1 + cos(Δφ)))
    I = I .* aperture; % Apply circular aperture mask
    % Plot 2D contour or surface plot
    if i <= 6
        [C, h] = contourf(x*1e3, y*1e3, I, 20, 'LineColor', 'none'); % Contour plot for aberrations
        set(h, 'EdgeColor', 'none'); % Remove contour edges
        colormap(parula); % Apply parula colormap
        colorbar('FontSize', 10, 'FontName', 'Times'); % Add colorbar
        title(['Perfect Plane Wavefront - ', upper(type)], 'FontSize', 12, 'FontWeight', 'bold'); % Set title
    else
        h = surf(x*1e3, y*1e3, I, 'EdgeColor', 'none'); % Surface plot for plane wave cases
        colormap(jet); % Apply jet colormap
        colorbar('FontSize', 10, 'FontName', 'Times'); % Add colorbar
        if i == 7
            title('Perfect Plane WavefrontFleet Interferogram', 'FontSize', 12, 'FontWeight', 'bold'); % Title for plane wave
        elseif i == 8
            title('Slight Wedge Angle Interferogram', 'FontSize', 12, 'FontWeight', 'bold'); % Title for wedge angle
        end
        zlabel('Intensity (a.u.)', 'FontSize', 12, 'FontWeight', 'bold'); % Label z-axis
        view(2); % Set top-down view
    end
    set(gca, 'FontSize', 10, 'FontName', 'Times', 'LineWidth', 1); % Set axes properties
    xlabel('x (mm)', 'FontSize', 12, 'FontWeight', 'bold'); % Label x-axis
    ylabel('y (mm)', 'FontSize', 12, 'FontWeight', 'bold'); % Label y-axis
    axis equal; % Set equal aspect ratio
    clim([0 2]); % Set color limits for intensity % Set color limits for intensity
    if i == 4 || i == 7
        colorbar('FontSize', 10, 'FontName', 'Times'); % Add colorbar for specific tiles
    end
end

% Plot interferograms with linear phase (wedge angle) for six phase types
for i = 1:6
    nexttile; % Move to next tile in layout
    type = phase_types{i}; % Select current phase type
    switch lower(type) % Define Zernike polynomial for each type
        case 'defocus'
            Z = sqrt(3)*rho.^2 - 1; % Defocus aberration
        case 'hcoma'
            Z = sqrt(8) * (3 * rho.^3 - 2 * rho) .* cos(theta); % Horizontal coma
        case 'vcoma'
            Z = sqrt(8) * (3 * rho.^3 - 2 * rho) .* sin(theta); % Vertical coma
        case 'spherical'
            Z = sqrt(5)*(6*rho.^4 - 6*rho.^2 + 1); % Spherical aberration
        case 'astig0'
            Z = sqrt(6)*rho.^2 .* cos(2*theta); % Astigmatism at 0 degrees
        case 'astig45'
            Z = sqrt(6)*rho.^2 .* sin(2*theta); % Astigmatism at 45 degrees
    end
    linear_phase = k * linear_coeff * X; % Linear phase for wedge angle
    delta_phi = k * Delta * Z + linear_phase; % Total phase difference
    I = 2*(1 + cos(delta_phi)); % Interference intensity
    I = I .* aperture; % Apply circular aperture mask
    % Plot 2D contour plot
    [C, h] = contourf(x*1e3, y*1e3, I, 20, 'LineColor', 'none'); % Contour plot
    set(h, 'EdgeColor', 'none'); % Remove contour edges
    colormap(parula); % Apply parula colormap
    colorbar('FontSize', 10, 'FontName', 'Times'); % Add colorbar
    set(gca, 'FontSize', 10, 'FontName', 'Times', 'LineWidth', 1); % Set axes properties
    xlabel('x (mm)', 'FontSize', 12, 'FontWeight', 'bold'); % Label x-axis
    ylabel('y (mm)', 'FontSize', 12, 'FontWeight', 'bold'); % Label y-axis
    title(['Slight Wedge Angle - ', upper(type)], 'FontSize', 12, 'FontWeight', 'bold'); % Set title
    axis equal; % Set equal aspect ratio
    clim([0 2]); % Set color limits for intensity
    if i == 4
        colorbar('FontSize', 10, 'FontName', 'Times'); % Add colorbar for specific tile
    end
end
