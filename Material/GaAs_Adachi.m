function epsilon = GaAs_Adachi(lambda)
    
    %% Function calculates the dielectric constant of GaAs from the wavelength in meters 
    % using the Adachi model in 
    % S. Ozaki & S. Adachi Spectroscopic ellipsometry and thermoreflectance of GaAs
    % J. Appl. Phys. 78, 3380–3386 (1995) https://doi.org/10.1063/1.359966
    
    % convert wavelength in meters to eV
    EPhoton = 6.62607015e-34*299792458/(lambda*1.602176634e-19); % in eV

    % Fitted Parameters
    E0 = 1.42;
    DELTA0 = 0.33;
    E1 = 2.91;
    DELTA1 = 0.23;

    A = 7;
    % Gamma_E0, Gamma_E1
    GAMMA = [0.03,0.12];

    % B1 B2 B1x B2x
    B = [3.5, 1.75, 1.2, 0.6];
    
    % E'_0, E2(1), E2(2), E'_1
    E = [4.45,4.77,5,6.6];
    C = [0.8,1.35,0.35,0.7];
    GAMMA_E = [0.67,0.62,0.4,0.6];

    % Derived Parameters
    % epsilon E0 and E0+Delta0 Contributions
    chi0 = (EPhoton + 1i*GAMMA(1))/E0;
    chiso = (EPhoton + 1i*GAMMA(1))/(E0 + DELTA0);
    fchi0 = chi0^-2*(2 - (1+chi0)^0.5 - (1-chi0)^0.5);
    fchiso = chiso^-2*(2 - (1+chiso)^0.5 - (1-chiso)^0.5);
    epsilon = A*E0^-1.5*(fchi0 + 0.5*fchiso*(E0/(E0+DELTA0))^1.5);

    % E1 and E1+Delta1 Contributions
    chi1 = (EPhoton + 1i*GAMMA(2))/E1;
    chi1s = (EPhoton + 1i*GAMMA(2))/(E1 + DELTA1); 
    epsilon = epsilon + -B(1)*chi1^-2*log(1-chi1^2) - B(2)*chi1s^-2*log(1-chi1s^2);

    epsilon = epsilon + (B(3)/(E1 - EPhoton - 1i*GAMMA(2)) + ...
        B(4)/(E1 + DELTA1 - EPhoton - 1i*GAMMA(2)))...
        *1.051799790264625; % 1.05...  is the limit of the sum, sum(1./((2.*n-1).^3))
    
    % E'_0, E2(1), E2(2) and E'_1 Contributions
    epsilon = epsilon + sum(C.*E.^2./(E.^2 - EPhoton^2 - 1i*EPhoton.*GAMMA_E));
end