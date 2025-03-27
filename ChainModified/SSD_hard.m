function x_out = SSD_hard(p, H_matrix, y_vector)

x_out = complex(zeros(p.Mr * p.Nr, p.frmLen / p.sq), zeros(p.Mr * p.Nr, p.frmLen / p.sq));
 
    for n = 1:p.frmLen / p.sq

        d_ml = inf;
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
            
            % Euclidean distance in the SSD ||W*y-Rx)||
            d = norm((W'*y - R * svx));  
            
            if d < d_ml
                svx_ml = svx;  
                d_ml = d;     
            end
        end

        x_out(:, n) = svx_ml;  % Store the estimated symbols for this tone
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
