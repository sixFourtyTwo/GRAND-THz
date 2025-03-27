function [LLRout, x_out, nVisitedNodes] = sphereDecoderSoftType(p, Hmat, yvec, N0)
    % sphereDecoderSoftType performs sphere decoding with soft decision for a given set of parameters
    % Inputs:
    %   p - Parameters structure
    %   Hmat - Channel matrix
    %   yvec - Received symbols
    %   N0 - Noise variance
    % Outputs:
    %   LLRout - Output LLRs
    %   x_out - Decoded symbol sequence
    %   nVisitedNodes - Number of visited nodes in the decoding process

    % Initialize output variables
    x_out = complex(zeros(p.Mr * p.Nr * p.frmLen / p.sq, 1), zeros(p.Mr * p.Nr * p.frmLen / p.sq, 1));
    LLRout = zeros(p.frmLen / p.sq, p.sq);
    Modtemp=lower(p.qam{1});

    for n = 1:p.frmLen / p.sq  % Loop over tones to detect
        % Extract channel matrix and received symbols for the current tone
        H = Hmat(((n - 1) * p.Mr * p.Nr + 1):(n * p.Mr * p.Nr), 1:p.Nt * p.Mt);
        rxSymbs = yvec(((n - 1) * p.Mr * p.Nr + 1):(n * p.Mr * p.Nr), 1);
        
        % Determine modulation types for each transmit antenna
        Nt = p.Nt * p.Mt;
        moduTypes = cell(1, Nt);
        for i = 1:ceil(Nt / 2)
            moduTypes{1, i} =Modtemp;
        end
        for i = ceil(Nt / 2) + 1: Nt
            moduTypes{1, i} = Modtemp;
        end


        % Perform sphere decoding with soft decision
        [outSoft, svx, nVisitedNodes] = nrSphereDecoderSoft(H, rxSymbs, moduTypes, p, N0);

        % Store LLRs and decoded symbols
        LLRout(n, :) = outSoft;
        x_out(((n - 1) * p.Mr * p.Nr + 1):(n * p.Mr * p.Nr), 1) = svx;

    end
end
