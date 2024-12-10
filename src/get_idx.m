function idx = get_idx(alpha_vec)
% get_idx - function to ask the user which undisturbed angle of attack
% should be used in plots for fig10 scripts.
% 
% INPUTS: 
% - alpha_vec, float: vector of possible alphas (rad)
%
% OUTPUTS:
% - idx, int: index corresponding to the desiderd alpha
%
% CALLED FUNCTIONS: -
%
% REVISIONS:
% - #v0 10/12/24, Boscariol Jacopo
%               Changes: release.

    % requesting the value of alpha
    alpha = input(['Please enter the undisturbed integer alpha value ' ...
        'you would like to be adopted in plots in deg: ' ...
        '(-20° <= alpha <= 30°): ']);
    
    % control loop
    while alpha < -20 || alpha > 30 || mod(alpha, 1)
        % not valid alpha message
        disp(['Alpha must be an integer greater than -20° and smaller ' ...
                'than 30°.']);

        % requesting the value of alpha
        alpha = input(['Please enter the undisturbed integer alpha value ' ...
            'you would like to be adopted in plots in deg: ' ...
            '(-20° <= alpha <= 30°): ']);
    end

    % display valid alpha
    disp(['You entered a valid value for alpha: ', num2str(alpha), '°']);

    % converting alpha to corresponding index
    differences = abs(alpha_vec - alpha);
    [~, idx] = min(differences);

end
