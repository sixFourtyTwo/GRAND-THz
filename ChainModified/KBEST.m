%% K-Best detector
function x_out = KBEST(p, Hmat, yvec, Kpar)

    % Initialize output
    x_out = complex(zeros(p.Mr * p.Nr * p.frmLen / p.sq, 1), zeros(p.Mr * p.Nr * p.frmLen / p.sq, 1));
    
    % Define xer as matrix (row vector)
    xer = zeros(1, p.Mt*p.Nt);  
    
    for n = 1 : p.frmLen / p.sq  % Loop over tones to detect
        % Extract relevant portions of H and y for the current tone
        H = Hmat(((n-1) * p.Mr * p.Nr + 1) : (n * p.Mr * p.Nr), 1 : p.Nt * p.Mt); 
        y = yvec(((n-1) * p.Mr * p.Nr + 1) : (n * p.Mr * p.Nr), 1);
        
        % Assign symbols to xer
        xer = p.syms(1, :);
        
        % Preprocessing
        [Q, R] = qr(H);
        y_hat = Q' * y;
        
        % Initialize Partial Euclidean Distance (PED) with last TX symbol
        PED_list = abs(xer * R(p.Mt * p.Nt, p.Mt * p.Nt) - y_hat(p.Mt * p.Nt)).^2;
        [PED_list, idx] = sort(PED_list);
        s = xer(:, idx);
        
        % Take the K-best
        s = s(:, 1 : min(Kpar, length(PED_list)));
        Kbest_PED_list = PED_list(1 : min(Kpar, length(PED_list)));
        
        % For each TX symbol
        for Level = (p.Mt * p.Nt) - 1 : -1 : 1
            PED_list = [];
            
            % Obtain the cumulative Euclidean distance considering the K-best previous nodes
            for k = 1 : length(Kbest_PED_list)
                tmp = Kbest_PED_list(k) + abs(xer * R(Level, Level) - y_hat(Level) + R(Level, Level+1 : (p.Mt * p.Nt)) * s(:, k)).^2;
                PED_list = [PED_list, tmp];
            end
            
            % Sort in ascending order
            s = [kron(ones(1, length(Kbest_PED_list)), xer); kron(s, ones(1, length(xer)))];
            [PED_list, idx] = sort(PED_list);
            s = s(:, idx);
            
            % Take the K-best
            s = s(:, 1 : min(Kpar, length(PED_list)));
            Kbest_PED_list = PED_list(1 : min(Kpar, length(PED_list)));
        end
        
        % Take the best
        x_out(((n-1) * p.Mr * p.Nr + 1) : (n * p.Mr * p.Nr), 1) = s(:, 1);
    end
end
