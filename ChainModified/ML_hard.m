function x_out = ML_hard(p, H_matrix, y_vector)
% ML_hard: Maximum Likelihood (ML) symbol detection with hard decision.
%
%   x_out = ML_hard(p, H_matrix, y_vector) performs ML detection
%   to estimate the transmitted symbols using hard decision.
%
%   Inputs:
%       - p: Structure containing system parameters.
%       - H_matrix: Channel matrix for each tone.
%       - y_vector: Received signal vector for each tone.
%
%   Output:
%       - x_out: Estimated transmitted symbols vector.

    % Initialize output matrix
    x_out = complex(zeros(p.Mr * p.Nr, p.frmLen / p.sq), zeros(p.Mr * p.Nr, p.frmLen / p.sq));
 
    % Loop over tones to detect symbols
    for n = 1:p.frmLen / p.sq

        % Initialize maximum distance for ML hypothesis
        d_ml = inf;
        num_lattice_pts = prod(p.Q); % Number of lattice points

        % Generate constellation vector B
        B = [zeros(1, p.Mt * p.Nt - 1), 1];
        for i = p.Mt * p.Nt - 1:-1:1
            B(i) = B(i + 1) * p.Q(i + 1);
        end

        % Initialize iteration count and symbol vectors
        iterat  = 0;
        svx = complex(zeros(p.Nt * p.Mt, 1), zeros(p.Nt * p.Mt, 1));
        svx_ml = complex(zeros(p.Nt * p.Mt, 1), zeros(p.Nt * p.Mt, 1)); 

        % Extract channel matrix and received signal vector for the current tone
        H = H_matrix(((n - 1) * p.Mr * p.Nr + 1):(n * p.Mr * p.Nr), 1:p.Nt * p.Mt);
        y = y_vector(((n - 1) * p.Mr * p.Nr + 1):(n * p.Mr * p.Nr), 1);

        % Iterate over lattice points
        while iterat  < num_lattice_pts
               
            % Generate symbol vector (in decimal)
            svd = mod(floor(iterat  ./ B), p.Q);
            iterat  = iterat  + 1;

            % Get corresponding complex modulation symbol vector
            for i = 1:p.Nt * p.Mt
                svx(i, 1) = p.syms(i, svd(i) + 1);
            end  

            % Compute Euclidean distance
            d = norm(y - H * svx);
            
            % Update ML hypothesis if a smaller distance is found
            if d < d_ml
                svx_ml = svx;
                d_ml = d;
            end
        
        end

        % Store ML estimate for current tone in output matrix
        x_out(:, n) = svx_ml;
    end
end
