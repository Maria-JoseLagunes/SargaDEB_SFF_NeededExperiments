%Estimate Irradiance as a function of light intensity and t based on
%Fourier equations


function Irradiance = get_irradiance(t, lightIntensity)

Irradiance = lightIntensity/2 + lightIntensity/2  * sin(2 * pi * (t + 12)/ 24);


% Irradiance = lightIntensity / 2 * (1 + sin(2 * pi * (t - 6) / 24));
% % shift by -6 so that max is at t = 12 and zero at t = 0 and 24
% Irradiance(Irradiance < 0) = 0;  % set negative values to 0 (night)

end
