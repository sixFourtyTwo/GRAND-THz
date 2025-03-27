function x_out = sphereDecoderHardType(p, Hmat, yvec)
    % sphereDecoderHardType performs sphere decoding with hard decision for a given set of parameters
    % Inputs:
    %   p - Parameters structure
    %   Hmat - Channel matrix
    %   yvec - Received symbols
    % Output:
    %   x_out - Decoded symbol sequence

    % Initialize output symbol sequence
    x_out = complex(zeros(p.Mr * p.Nr * p.frmLen / p.sq, 1), zeros(p.Mr * p.Nr * p.frmLen / p.sq, 1));
    Modtemp=lower(p.qam{1});
    for n = 1:p.frmLen / p.sq  % Loop over tones to detect
        % Extract channel matrix and received symbols for the current tone
        H = Hmat(((n - 1) * p.Mr * p.Nr + 1):(n * p.Mr * p.Nr), 1:p.Nt * p.Mt);
        rxSymbs = yvec(((n - 1) * p.Mr * p.Nr + 1):(n * p.Mr * p.Nr), 1);

        % Determine modulation types for each transmit antenna
        Nt = p.Nt * p.Mt;
        moduTypes = cell(1, Nt);
        for i = 1:ceil(Nt / 2)
            moduTypes{1, i} = Modtemp ;
        end
        for i = ceil(Nt / 2) + 1: Nt
            moduTypes{1, i} = Modtemp;
        end

        % Perform sphere decoding with hard decision
        [out, svx, nVisitedNodes] = nrSphereDecoderHard(H, rxSymbs, moduTypes, p);

        % Update the decoded symbol sequence
        x_out(((n - 1) * p.Mr * p.Nr + 1):(n * p.Mr * p.Nr), 1) = svx;
    end
end
