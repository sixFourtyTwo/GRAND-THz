function [out, svx_hard, nVisitedNodes] = nrSphereDecoderHard(H, rxSymbs, moduTypes, p)
            % nrSphereDecoderHard performs sphere decoding for seeking the maximum-likelihood solution
            % for a set of symbols transmitted over the MIMO channel.
            % Inputs:
            %   H - Channel matrix Nr x Nt
            %   rxSymbs - Received symbols: Nr x 1
            %   moduTypes - Modulation for each TX antenna in a cell array of 1 x Nt
            %   p - Parameters structure
            % Outputs:
            %   out - Output bits
            %   svx_hard - Decoded symbols
            %   nVisitedNodes - Number of leaf nodes traversed
        
            [Nr, Nt] = size(H);
            assert(Nr == length(rxSymbs), 'Number of rows in H must be equal to the length of rxSymbs');
            assert(Nr >= Nt, 'Nr must be greater or equal to Nt');
        
            % Generate constellation
            constls = genConstls(moduTypes);
        
            % Initialize variables
            nOutBits = 0;
            offsets = zeros(1, Nt);
            for i = 1:Nt
                offsets(i) = nOutBits;
                nOutBits = nOutBits + constls{1, i};
            end
        
            bitsML = false(1, nOutBits);
            for i = 1:Nt
                bitsML(offsets(i) + (1:constls{1, i})) = constls{3, i}(1, :);
            end
        
            lambdaML = inf;
            out = false(1, nOutBits);
        
            % QR decomposition of H
            [Q, R] = qr(H, 0);
            Qy = Q' * rxSymbs;
        
            % Initialize partial symbol vector of the visiting node & PD (start from the root)
            psv = constls{2, Nt}(1);
            psvIdx = 1;
            psvBitLabels = constls{3, Nt}(1, :);
            traverseDone = false;
            nVisitedNodes = 0;
        
            while (~traverseDone)
                curLevel = Nt - length(psv) + 1;
                di = abs(Qy(curLevel) - sum(R(curLevel, curLevel:Nt) .* psv)).^2;
        
                if curLevel == Nt
                    ped(curLevel) = di;
                else
                    ped(curLevel) = ped(curLevel + 1) + di;
                end
        
                if curLevel == 1 % Leaf node
                    nVisitedNodes = nVisitedNodes + 1;
                    if ped(curLevel) < lambdaML % Smaller Euclidean distance found
                        lambdaML = ped(1);
                        bitsML = psvBitLabels;
                    end 
                    % Go to next node (right or up)
                    [psvIdx, psv, psvBitLabels, traverseDone] = moveToNextNode(psvIdx, psv, psvBitLabels, constls);
                else % Non-leaf node
                    if ped(curLevel) >= lambdaML % Ped > search radius, prune subtree
                        % Go to next node (right or up)
                        [psvIdx, psv, psvBitLabels, traverseDone] = moveToNextNode(psvIdx, psv, psvBitLabels, constls);
                    else % Go down a level
                        psvIdx = [1, psvIdx];
                        psv = [constls{2, curLevel - 1}(1), psv]; % Add a symbol to the left
                        psvBitLabels = [constls{3, curLevel - 1}(1, :), psvBitLabels]; % Add the corresponding bit to the left
                    end
                end
            end
        
            out = bitsML';
            q = log2(p.Q);
            svb_ml = zeros(1, p.Nt * p.Mt); 
            svx_hard = complex(zeros(p.Nt * p.Mt, 1), zeros(p.Nt * p.Mt, 1));
        
            for i = 1:numel(svb_ml)
                svb_ml(i) = bi2de(bitsML((i - 1) * q(1) + 1 : i * q(1)), 'left-msb');
            end
        
            for i = 1:p.Nt * p.Mt
                svx_hard(i, 1) = p.syms(i, svb_ml(i) + 1);
            end
end

function [psvIdx, psv, psvBitLabels, traverseDone] = moveToNextNode(psvIdx, psv, psvBitLabels, constls)
        % Move to the right node of the same level if symbol idx < M. Else go up.
    
        Nt = size(constls, 2);
        curLevel = Nt - length(psvIdx) + 1;
    
        while true
            if psvIdx(1) < 2^constls{1, curLevel}
                psvIdx(1) = psvIdx(1) + 1;
                psv(1) = constls{2, curLevel}(psvIdx(1));
                psvBitLabels(1:constls{1, curLevel}) = constls{3, curLevel}(psvIdx(1), :);
                traverseDone = false;
                return;
            else
                if curLevel == Nt
                    traverseDone = true;
                    return;
                else % Up a level
                    psvIdx(1) = [];
                    psv(1) = [];
                    psvBitLabels(1:constls{1, curLevel}) = [];
                    curLevel = curLevel + 1;
                end
            end
        end
end

function constls = genConstls(moduTypes)
        % Generate constellation info for a cell
    
        constls = cell(3, length(moduTypes));
    
        for i = 1:length(moduTypes)
            moduType = moduTypes{i};
            switch lower(moduType)
                case 'bpsk'
                    M = 2;
                case 'qpsk'
                    M = 4;
                case '16qam'
                    M = 16;
                case '64qam'
                    M = 64;
                case '256qam'
                    M = 256;
            end
    
            K = log2(M);
            constls{1, i} = K;
    
            bitLabels = de2bi(0:M-1, 'left-msb');
            constls{3, i} = bitLabels;
    
            symbBitsIn = bitLabels.';
            constSymbs = nrModuMapper(symbBitsIn(:), lower(moduType));
    
            constls{2, i} = constSymbs;
        end
end
