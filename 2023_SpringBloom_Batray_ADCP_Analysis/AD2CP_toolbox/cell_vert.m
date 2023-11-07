function z = cell_vert(pitch, roll, velocity_range, beam_number)
    % CELL_VERT Calculate a vertical displacement below instrument for
    % each ADCP bin adjusting for pitch and roll (in degrees)
    % Positive roll: Port wing up
    % Positive pitch: Pitch up
    %
    % Beam 1: Forward   (47.5 degrees off horizontal)
    % Beam 2: Port      (25 degrees off horizontal)
    % Beam 3: Aft       (47.5 degrees off horizontal)
    % Beam 4: Starboard (25 degrees off horizontal)
    %
    % Beam angle is only incorporated in pitch for Beams 1 & 3 and
    % in roll for Beams 2 & 4
    
    % Ensure beam_number is valid
    if beam_number < 1 || beam_number > 4
        error('Must specify beam number as an integer between 1 and 4.');
    end
    
    % Define beam angles
    if beam_number == 1 || beam_number == 3
        beam_angle = 47.5;
    else
        beam_angle = 25;
    end
    
    % Calculate z based on beam_number
    switch beam_number
        case 1
            pitch_adjusted = velocity_range .* sin(deg2rad(90 + beam_angle + pitch));
            z = (pitch_adjusted .* sin(deg2rad(90 - roll)));
        case 2
            pitch_adjusted = velocity_range .* sin(deg2rad(90 - pitch));
            z = (pitch_adjusted .* sin(deg2rad(90 + roll + beam_angle)));
        case 3
            pitch_adjusted = velocity_range .* sin(deg2rad(90 - beam_angle + pitch));
            z = (pitch_adjusted .* sin(deg2rad(90 - roll)));
        case 4
            pitch_adjusted = velocity_range .* sin(deg2rad(90 - pitch));
            z = (pitch_adjusted .* sin(deg2rad(90 - roll + beam_angle)));
    end
    
    % Transpose z to match Python function output
    z = z';
end