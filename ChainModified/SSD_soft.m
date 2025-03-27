function LLRout = SSD_soft(p, H_matrix, y_vector, N0_vector)   
    
    LLRout = complex(zeros(p.frmLen / p.sq, p.sq), zeros(p.frmLen / p.sq, p.sq));  
    svx_out_ML = complex(zeros(p.Mr * p.Nr * p.frmLen / p.sq, 1), zeros(p.Mr * p.Nr * p.frmLen / p.sq, 1));

    for n = 1:p.frmLen / p.sq
        
        q = log2(p.Q);
        sq = sum(p.q);
        svd_ml = zeros(1, p.Nt * p.Mt);  % ML symbol vector (decimal)
        svb_ml = zeros(1, sq);
        d_ml = complex(inf, inf);
        dc_ml = complex(inf * ones(1, sq), inf * ones(1, sq));
  
        num_lattice_pts = prod(p.Q);
        B = [zeros(1, p.Mt * p.Nt - 1), 1]; 
        for i = p.Mt * p.Nt - 1:-1:1
            B(i) = B(i + 1) * p.Q(i + 1); 
        end
        % Initialize iteration count and symbol vectors
        iterat  = 0;  % Iterator for lattice points
        svx = complex(zeros(p.Nt * p.Mt, 1), zeros(p.Nt * p.Mt, 1));  
        svx_ml = complex(zeros(p.Nt * p.Mt, 1), zeros(p.Nt * p.Mt, 1)); 

        % Extract the channel matrix and received signal vector for the current tone
        H = H_matrix(((n - 1) * p.Mr * p.Nr + 1):(n * p.Mr * p.Nr), 1:p.Nt * p.Mt);
        [W, R] = compute_W_R(H);  % Perform WR decomposition on H to get W and R
        y = y_vector(((n - 1) * p.Mr * p.Nr + 1):(n * p.Mr * p.Nr), 1);  % Received signal for the tone
        % Iterate over all lattice points for symbol detection
        while iterat < num_lattice_pts
               
            % Generate the symbol vector (in decimal format)
            svd = mod(floor(iterat ./ B), p.Q);  % Symbol vector in decimal format
            iterat = iterat + 1;

            % Convert decimal symbols to modulation symbols (complex symbols)
            for i = 1:p.Nt * p.Mt
                svx(i, 1) = p.syms(i, svd(i) + 1);  % Map decimal values to modulation symbols
            end  
            
            % Rotate by R and compute distance
            d = (W'*y - R * svx)' * ((W'*y - R * svx) .* (1 ./ N0_vector));
            % Convert decimal to binary representation
            svb_tem = cell(length(svd), 1);
            %:
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


function [Q, R] = compute_W_R(H)
% compute_WR  decomposition of the channel matrix H with updates.

    [Q, R] = qr(H);
    
    N = size(H, 1);  % Number of rows in H (assuming square or rectangular)

    % Iterative procedure to update Q and R
    for m = N-1:-1:1
        for n = m+1:N-1
            % rho_mn = rmn / rnn
            rho_mn = R(m, n) / R(n, n);            
            % q_m = q_m - q_n * rho_mn*
            Q(:, m) = Q(:, m) - Q(:, n) * conj(rho_mn);            
            %  rmn = rmn - rnn * rho_mn
            R(m, n) = R(m, n) - R(n, n) * rho_mn;
            %  rmN = rmN - rnN * rho_mn
            R(m, N) = R(m, N) - R(n, N) * rho_mn;
        end
        
        % Normalize q_m 
        R(m, m) = R(m, m) / norm(Q(:, m));  
        R(m, N) = R(m, N) / norm(Q(:, m));  
        Q(:, m) = Q(:, m) / norm(Q(:, m));  
    end
    
end
