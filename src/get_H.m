function H = get_H(default_H)
% get_H - function to ask the user values for the tail LE y-coordinate
% if he does not want to use the default ones (taken from the paper).
% 
% INPUTS: 
% - default_H, float: default H value if not specified
%
% OUTPUTS:
% - H, float: output value for H, either specified or default
%
% CALLED FUNCTIONS: -
%
% REVISIONS:
% - #v0 28/11/24, Boscariol Jacopo
%               Changes: release.

    % asking the user if they want to specify a value
    answer = input(['Do you want to input a value for H which is different ' ...
        'from the default one (', num2str(default_H), ' m)? (y/n): '], 's');
    
    % checking if the answer is 'y' or 'Y'
    if strcmpi(answer, 'y')
        % defining max_H to get a reasonable valid H
        max_H = 0.2;

        % init H for the control loop
        H = max_H + 1;
    
        % control loop
        while abs(H) > max_H
            % requesting the value of H
            H = input(['Please enter the value for H in m (H must be ' ...
                'between +-', num2str(max_H), ' m): ']);
    
            % checking H validity
            if abs(H) > max_H
                disp(['H must be between +-', num2str(max_H), ' m. ' ...
                    'Please try again.']);
            end
        end
    
        % display valid H
        disp(['You entered a valid value for H: ', num2str(H), ' m']);

    else
        disp(['The default value for H (', num2str(default_H), ' m) will be ' ...
            'used.']);

        % assigning output
        H = default_H;
    end

end
