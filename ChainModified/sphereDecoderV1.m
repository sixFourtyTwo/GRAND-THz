function decodedSymbols = sphereDecoder(yvec, p)
    % yvec: received signal vector
    % p: structure containing system parameters, including p.syms and other necessary parameters
    
    latticeMatrix = p.syms; % Get lattice basis matrix from system parameters
    radius = p.radius; % Get decoding radius from system parameters
    
    numSymbols = size(latticeMatrix, 2); % Number of symbols in the lattice
    
    % Initialize decoded symbols and minimum distance
    decodedSymbols = zeros(size(yvec));
    minDistanceSquared = Inf;
    
    % Iterate through all lattice points
    for i = 1:numSymbols
        latticePoint = latticeMatrix(:, i); % Get a lattice point
        
        % Calculate distance between received signal and lattice point
        distanceSquared = norm(yvec - latticePoint)^2;
        
        % If the distance is within the decoding sphere
        if distanceSquared <= radius^2
            % Check if it's the closest point
            if distanceSquared < minDistanceSquared
                % Update closest point and minimum distance
                minDistanceSquared = distanceSquared;
                decodedSymbols = latticePoint;
            end
        end
    end
end
