function [betasw, beta90sw, bsw] = betasw(lambda, T, theta, S, delta)
% betasw_ZHH2009: Calculates volume scattering parameters of seawater.
% 
% Inputs:
%   lambda: Wavelength (nm) (e.g., 700 for FLBB)
%   T: Temperature in degree Celsius (can be a scalar or an array)
%   theta: Scattering angle in degrees (a single value)
%   S: Salinity (can be a scalar or an array of the same size as T)
%   delta: Depolarization ratio (optional, default = 0.039)
%
% Outputs:
%   betasw: Volume scattering at the angle defined by theta
%   beta90sw: Volume scattering at 90 degrees
%   bsw: Total scattering coefficient

    if nargin < 5
        delta = 0.039; % Default value if not provided
    end

    % Ensure T and S are column vectors and have the same size
    T = T(:);
    S = S(:);
    if length(T) ~= length(S)
        error('Temperature (T) and Salinity (S) must have the same number of elements.');
    end

    % Initialize outputs to NaN
    betasw = NaN(size(T));
    beta90sw = NaN(size(T));
    bsw = NaN(size(T));

    % Filter out NaN values in T and S
    validIdx = ~isnan(T) & ~isnan(S);
    T_valid = T(validIdx);
    S_valid = S(validIdx);
    
    % Constants
    Na = 6.0221417930e23;   % Avogadro's constant
    Kbz = 1.3806503e-23;    % Boltzmann constant
    M0 = 18e-3;             % Molecular weight of water in kg/mol
    
    % Convert temperature to Kelvin
    Tk = T_valid + 273.15;
    
    % Calculate refractive indices and derivatives
    [nsw, dnds] = RInw(lambda, T_valid, S_valid);
    
    % Calculate isothermal compressibility
    IsoComp = BetaT(T_valid, S_valid);
    
    % Calculate density of seawater
    density_sw = rhou_sw(T_valid, S_valid);
    
    % Calculate partial derivative of natural logarithm of water activity w.r.t. salinity
    dlnawds = dlnasw_ds(T_valid, S_valid);
    
    % Calculate density derivative of refractive index from PMH model
    DFRI = PMH(nsw);
    
    % Pre-calculate common terms
    lambda_m = lambda * 1e-9;
    theta_rad = theta * pi / 180;
    cosTheta2 = cos(theta_rad).^2;
    
    % Volume scattering at 90 degree due to the density fluctuation
    beta_df = (pi^2 / 2) * (lambda_m.^(-4)) .* Kbz .* Tk .* IsoComp .* DFRI.^2 .* (6 + 6 * delta) ./ (6 - 7 * delta);
    
    % Volume scattering at 90 degree due to the concentration fluctuation
    flu_con = S_valid .* M0 .* dnds.^2 ./ density_sw ./ (-dlnawds) ./ Na;
    beta_cf = 2 * pi^2 * (lambda_m.^(-4)) .* nsw.^2 .* flu_con .* (6 + 6 * delta) ./ (6 - 7 * delta);
    
    % Total volume scattering at 90 degree
    beta90sw_valid = beta_df + beta_cf;
    
    % Total scattering coefficient
    bsw_valid = 8 * pi / 3 * beta90sw_valid * (2 + delta) / (1 + delta);
    
    % Volume scattering at angles defined by theta
    betasw_valid = beta90sw_valid .* (1 + cosTheta2 .* (1 - delta) ./ (1 + delta));

    % Assign valid results back to the corresponding positions
    betasw(validIdx) = betasw_valid;
    beta90sw(validIdx) = beta90sw_valid;
    bsw(validIdx) = bsw_valid;

end

function [nsw, dnswds]= RInw(lambda, Tc, S)
    % Ensure inputs are column vectors
    Tc = Tc(:);
    S = S(:);

    % refractive index of air is from Ciddor (1996, Applied Optics)
    n_air = 1.0 + (5792105.0./(238.0185 - 1./(lambda/1e3).^2) + 167917.0./(57.362 - 1./(lambda/1e3).^2)) / 1e8;
    
    % refractive index of seawater is from Quan and Fry (1994, Applied Optics)
    n0 = 1.31405; n1 = 1.779e-4 ; n2 = -1.05e-6 ; n3 = 1.6e-8 ; n4 = -2.02e-6 ;
    n5 = 15.868; n6 = 0.01155;  n7 = -0.00423;  n8 = -4382 ; n9 = 1.1455e6;
    
    nsw = n0 + (n1 + n2.*Tc + n3.*Tc.^2).*S + n4.*Tc.^2 + (n5 + n6.*S + n7.*Tc)./lambda + n8./lambda.^2 + n9./lambda.^3; % pure seawater
    nsw = nsw.*n_air;
    dnswds = (n1 + n2.*Tc + n3.*Tc.^2 + n6./lambda).*n_air;
end

function IsoComp = BetaT(Tc, S)
    % Ensure inputs are column vectors
    Tc = Tc(:);
    S = S(:);

    % pure water secant bulk Millero (1980, Deep-sea Research)
    kw = 19652.21 + 148.4206.*Tc - 2.327105.*Tc.^2 + 1.360477e-2.*Tc.^3 - 5.155288e-5.*Tc.^4;
    Btw_cal = 1./kw;
    
    % seawater secant bulk
    a0 = 54.6746 - 0.603459.*Tc + 1.09987e-2.*Tc.^2 - 6.167e-5.*Tc.^3;
    b0 = 7.944e-2 + 1.6483e-2.*Tc - 5.3009e-4.*Tc.^2;
    
    Ks = kw + a0.*S + b0.*S.^1.5;
    
    % calculate seawater isothermal compressibility from the secant bulk
    IsoComp = 1./(Ks*1e5); % unit is pa
end

function density_sw = rhou_sw(Tc, S)
    % Ensure inputs are column vectors
    Tc = Tc(:);
    S = S(:);

    % density of water and seawater,unit is Kg/m^3, from UNESCO,38,1981
    a0 = 8.24493e-1;  a1 = -4.0899e-3; a2 = 7.6438e-5; a3 = -8.2467e-7; a4 = 5.3875e-9;
    a5 = -5.72466e-3; a6 = 1.0227e-4;  a7 = -1.6546e-6; a8 = 4.8314e-4;
    b0 = 999.842594; b1 = 6.793952e-2; b2 = -9.09529e-3; b3 = 1.001685e-4;
    b4 = -1.120083e-6; b5 = 6.536332e-9;
     
    % density for pure water 
    density_w = b0 + b1.*Tc + b2.*Tc.^2 + b3.*Tc.^3 + b4.*Tc.^4 + b5.*Tc.^5;
    
    % density for pure seawater
    density_sw = density_w + ((a0 + a1.*Tc + a2.*Tc.^2 + a3.*Tc.^3 + a4.*Tc.^4).*S + ...
                 (a5 + a6.*Tc + a7.*Tc.^2).*S.^1.5 + a8.*S.^2);
end

function dlnawds = dlnasw_ds(Tc, S)
    % Ensure inputs are column vectors
    Tc = Tc(:);
    S = S(:);

    dlnawds = (-5.58651e-4 + 2.40452e-7.*Tc - 3.12165e-9.*Tc.^2 + 2.40808e-11.*Tc.^3) +...
               1.5.*(1.79613e-5 - 9.9422e-8.*Tc + 2.08919e-9.*Tc.^2 - 1.39872e-11.*Tc.^3).*S.^0.5 +...
               2.*(-2.31065e-6 - 1.37674e-9.*Tc - 1.93316e-11.*Tc.^2).*S;
end

function n_density_derivative = PMH(n_wat)
    n_wat2 = n_wat.^2;
    n_density_derivative = (n_wat2 - 1).*(1 + 2/3.*(n_wat2 + 2).*(n_wat./3 - 1./3./n_wat).^2);
end