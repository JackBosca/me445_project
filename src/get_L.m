function L = get_L(cW, default_L)
% get_L - function to ask the user values for the tail LE x-coordinate
% if he does not want to use the default ones (taken from the paper).
% 
% INPUTS: 
% - cW, float: wing chord to test validity of requested L
% - default_L, float: default L value if not specified
%
% OUTPUTS:
% - L, float: output value for L, either specified or default
%
% CALLED FUNCTIONS: -
%
% REVISIONS:
% - #v0 28/11/24, Boscariol Jacopo
%               Changes: release.

    % asking the user if they want to specify a value
    answer = input(['Do you want to input a value for L which is different ' ...
        'from the default one (', num2str(default_L), ' m)? (y/n): '], 's');
    
    % checking if the answer is 'y' or 'Y'
    if strcmpi(answer, 'y')
        % init L for the control loop
        L = 0;

        % defining delta, max_L to get a reasonable valid L
        delta = cW/2;
        max_L = 0.3;
    
        % control loop
        while L < (cW + delta) || L > max_L
            % requesting the value of L
            L = input(['Please enter the value for L in m (' ...
                , num2str(cW + delta), ' m <= L <= ' ...
                , num2str(max_L), ' m): ']);
    
            % checking L validity
            if L < (cW + delta) || L > max_L
                disp(['L must be greater than ', num2str(cW + delta), ' m ' ...
                    'and smaller than ', num2str(max_L), ' m. ' ...
                    'Please try again.']);
            end
        end
    
        % display valid L
        disp(['You entered a valid value for L: ', num2str(L), ' m']);

    else
        disp(['The default value for L (', num2str(default_L), ' m) will be ' ...
            'used.']);

        % assigning output
        L = default_L;
    end

end
