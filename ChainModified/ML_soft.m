function LLRout = ML_soft(p, H_matrix, y_vector, N0_vector, LLR_in)
% ML_soft: Maximum Likelihood (ML) symbol detection with soft decision.
%
%   LLRout = ML_soft(p, H_matrix, y_vector, N0_vector, LLR_in)
%   performs ML detection with soft decision to estimate the transmitted
%   symbols.
%
%   Inputs:
%       - p: Structure containing system parameters.
%       - H_matrix: Channel matrix.
%       - y_vector: Received signal vector .
%       - N0_vector: Noise power for each tone.
%       - LLR_in: Input Log-Likelihood Ratio (LLR) information.
%
%   Output:
%       - LLRout: Output LLR information.

    % Initialize LLRout matrix
    LLRout = complex(zeros(p.frmLen / p.sq, p.sq), zeros(p.frmLen / p.sq, p.sq));  

    N0_vector2 = complex(size(N0_vector), size(N0_vector)); 
    N0_vector2 = N0_vector;

    % Initialize vector to store ML estimates
    svx_out_ML = complex(zeros(p.Mr * p.Nr * p.frmLen / p.sq, 1), zeros(p.Mr * p.Nr * p.frmLen / p.sq, 1));

    % Loop over tones to detect symbols
    for n = 1:p.frmLen / p.sq

        q = log2(p.Q);
        sq = sum(p.q);
        svd_ml = zeros(1, p.Nt * p.Mt);  % ML symbol vector (decimal)
        svb_ml = zeros(1, sq);
        d_ml = complex(inf, inf);
        dc_ml = complex(inf * ones(1, sq), inf * ones(1, sq));
        num_lattice_pts = prod(p.Q);
        
        % Generate constellation vector B
        B = [zeros(1, p.Mt * p.Nt - 1), 1];
        for i = p.Mt * p.Nt - 1:-1:1
            B(i) = B(i + 1) * p.Q(i + 1);
        end

        iterat  = 0;
        svx = complex(zeros(p.Nt * p.Mt, 1), zeros(p.Nt * p.Mt, 1));
        svx_ml = complex(zeros(p.Nt * p.Mt, 1), zeros(p.Nt * p.Mt, 1)); 

        % Extract channel matrix and received signal vector for the current tone
        H = H_matrix(((n - 1) * p.Mr * p.Nr + 1):(n * p.Mr * p.Nr), 1:p.Nt * p.Mt); 
        y = y_vector(((n - 1) * p.Mr * p.Nr + 1):(n * p.Mr * p.Nr), 1);
        [Qq, R] = qrd1(H);
        yth = Qq' * y;
        W = Qq'; 
        N0_vector2(((n - 1) * p.Mr * p.Nr + 1):(n * p.Mr * p.Nr), 1) = N0_vector .* diag(W * W');
        
        % Iterate over lattice points
        while iterat  < num_lattice_pts
            % Generate symbol vector (in decimal)
            svd = mod(floor(iterat  ./ B), p.Q);
            iterat  = iterat  + 1;

            % Get corresponding complex modulation symbol vector
            for i = 1:p.Nt * p.Mt
                svx(i, 1) = p.syms(i, svd(i) + 1);
            end  

            % Rotate by R and compute distance
            d = (yth - R * svx)' * ((yth - R * svx) .* (1 ./ N0_vector));
            
            % Convert decimal to binary representation
            svb_tem = cell(length(svd), 1);
            for i = 1:length(svd)
                svb_tem{i} = de2bi(svd(i), q(i), 'left-msb');
            end
            svb = [svb_tem{:}];
            
            % Update ML hypothesis if smaller distance is found
            s = 1;
            for i = 1:sq
                if i > sum(q(1:s))
                    s = s + 1;
                end
                if svb(i) + svb_ml(i) == 1
                    if d < d_ml
                        dc_ml(i) = d_ml;                
                    elseif d < dc_ml(i)
                        dc_ml(i) = d;
                    end            
                end    
            end
            if d < d_ml
                svd_ml = svd;
                svx_ml = svx;
                d_ml = d;
                cell_array_svd_ml = num2cell(svd_ml);
                for i = 1:length(svd_ml)
                    cell_array_svd_ml{i} = de2bi(svd_ml(i), q(i), 'left-msb');
                end
                svb_ml = [cell_array_svd_ml{:}];
            end    
        end

        % Compute LLRs from d_ml and dc_ml
        LLRout(n, :) = (1 - 2 * svb_ml) .* (d_ml - dc_ml); % + LLR_in when you have llrin we should add it 
        
        % Store ML estimate for current tone
        svx_out_ML(((n - 1) * p.Mr * p.Nr + 1):(n * p.Mr * p.Nr), 1) = svx_ml;
    end
end
