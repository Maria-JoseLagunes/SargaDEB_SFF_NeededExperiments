function pastelJet = pastel_jet(n)
    if nargin < 1
        n = 256; % Default number of colors
    end
    
    % Get the standard jet colormap
    jetMap = jet(n);
    
    % Convert to HSV to adjust saturation and brightness
    hsvMap = rgb2hsv(jetMap);
    
    % Reduce saturation and increase brightness for a pastel effect
    hsvMap(:,2) = hsvMap(:,2) * 0.8;  % Reduce saturation (pastel effect)
    hsvMap(:,3) = hsvMap(:,3) * 0.9;  % Increase brightness slightly
    
    % Ensure values stay within valid range
    hsvMap(hsvMap > 1) = 1;  
    
    % Convert back to RGB
    pastelJet = hsv2rgb(hsvMap);
end
