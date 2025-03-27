function [outtemp,svx_soft,nVistedNodes] = nrSphereDecoderSoft(H,rxSymbs,moduTypes,p,N0)
% [out,nVistedNodes] = nrSphereDecoder(H,rxSymbs,moduTypes,outType) uses sphere decoding algorithms
% (single tree seach)for seeking the maximum-likelihood solution 
%  for a set of symbols transmitted over the MIMO channel. 
% Input:
%   H - channel matrix Nr X Nt
%   rxSymbs - received symbols: Nr X 1
%   moduTypes - nodulation for each tx antenna in a cell array of 1 X Nt
%   outType - 'soft' or 'hard'
% Output:
%   out - soft bits or hard bits output, soft bits (not sacled by 1/N0)
%   nVistedNodes - number of leaf nodes traversed

    [Nr,Nt] = size(H);
    assert(Nr == length(rxSymbs));
    assert(Nr >= Nt,'Nr must be greater or equal to Nt!')

    constls = genConstls(moduTypes);
    KModu = zeros(1,Nt);
    for i= 1 : Nt
        KModu(i) = constls{1,i};
    end
    
    % total number of bits 
    nOutBits = 0;
    offsets = zeros(1,Nt); % offset bit index for each symbol in symbol vector
    for i = 1: Nt
        offsets(i) = nOutBits;
        nOutBits = nOutBits+constls{1,i};
    end
    
    
    
    % hypothesis ML bits
    bitsML = false(1,nOutBits);
    for i = 1: Nt
        bitsML(offsets(i)+ (1:constls{1,i})) = constls{3,i}(1,:);
    end
    
    
    % hypothesis minimum distances for ML bits and its counterpart
    lambdaML = inf;% sacler
    lambdaMLBar = inf*ones(1,nOutBits);
    

    out = zeros(1,nOutBits);

    
    
    % QR decomposition of H
    [Q, R] = qr(H,0);
    Qy = Q'*rxSymbs;%Q^H*y
    
    
    % intinitalize patial symbol vector of the visiting node & PD
    % start form the root
    psv = constls{2,Nt}(1);
    psvIdx = 1;
    psvBitLabels = constls{3,Nt}(1,:);
    ped = zeros(1,Nt);
    traverseDone = false;
    nVistedNodes = 0;
    while (~traverseDone)
    
        
        curLevl = Nt - length(psv) +1;
        %disp(psvIdx);
        % compute distance increment
        di = (abs(Qy(curLevl) -  sum(R(curLevl,curLevl:Nt).*psv)).^2) *(1 ./ N0);
    
        % compute partial Eucidean distance
        if curLevl == Nt
            ped(curLevl) = di;
        else
            ped(curLevl) = ped(curLevl +1) + di;
        end
    
    
            if curLevl == 1 % leaf node has a complete symbol vector
                nVistedNodes = nVistedNodes + 1;
                if ped(curLevl) < lambdaML %smaller euclidean is found
                    % update the counter hypotheses
                    lambdaMLBar(psvBitLabels ~= bitsML) = lambdaML;
                    % update the hypotheses
                    lambdaML = ped(1);
                    bitsML = psvBitLabels;
                else
                    % update the counter hypotheses
                    lambdaMLBar(psvBitLabels ~= bitsML & lambdaMLBar > ped(1)) = ped(1);  
                end
                % move right or up
                [psvIdx,psv,psvBitLabels,traverseDone] = moveToNextNode(psvIdx,psv,psvBitLabels,constls);
                    
            else % non-leaf nodes
                if ped(curLevl) < lambdaML % go down a level
                    psvIdx = [1, psvIdx];
                    psv = [constls{2,curLevl-1}(1),psv]; % add a symbol to left
                    psvBitLabels = [constls{3,curLevl-1}(1,:),psvBitLabels]; % add the corresponding bit to the left
                else  % Check if there is a smaller Euclidean distance for the couther hypothesis in the sub-tree
                    
                    % for those bits not yet traversed
                    nBitsNotVisited = sum(KModu(1:curLevl-1));
                    lamdaMax1 = max(lambdaMLBar(1:nBitsNotVisited));
    
                    % for those bits already traversed
                    lambdaMLBarVisted =  lambdaMLBar(nBitsNotVisited+1:end);
                    lamdaMax2 = max(lambdaMLBarVisted(psvBitLabels ~= bitsML(nBitsNotVisited+1:end)));
    
                    if ped(curLevl) < max(lamdaMax1,lamdaMax2)
                        % Continue searching counter-hypotheses, down a level
                        psvIdx = [1, psvIdx];
                        psv = [constls{2,curLevl-1}(1),psv]; % add a symbol to left
                        psvBitLabels = [constls{3,curLevl-1}(1,:),psvBitLabels]; % add the corresponding bit to the left
                    else % go right or up
                        [psvIdx,psv,psvBitLabels,traverseDone] = moveToNextNode(psvIdx,psv,psvBitLabels,constls);
                    end
    
                end
            end
    
        
        
    end

        bitsML = (bitsML == 1);
        % out(bitsML) = (lambdaML - lambdaMLBar(bitsML));
        % out(~bitsML) = (lambdaMLBar(~bitsML) - lambdaML);
        outtemp=(1 - 2 * bitsML) .* (lambdaML - lambdaMLBar);% u can add LLrin
        



        q = log2(p.Q);
        svb_ml = zeros(1, p.Nt * p.Mt); 
        svx_soft = complex(zeros(p.Nt * p.Mt, 1), zeros(p.Nt * p.Mt, 1));
        
        % Convert the chunk of bits to decimal
        for i = 1:numel(svb_ml) ,  svb_ml(i) = bi2de(bitsML((i - 1) * q(1) + 1 : i * q(1)), 'left-msb');end
        for i = 1:p.Nt * p.Mt , svx_soft(i, 1) = p.syms(i, svb_ml(i) + 1);  end
        
      

end

function [psvIdx,psv,psvBitLabels,traverseDone] = moveToNextNode(psvIdx,psv,psvBitLabels,constls)
    % Move to the right node of the same level if symbol idx < M. Else
    % go up.

    Nt = size(constls,2);
    curLevl = Nt - length(psvIdx) +1;
    while 1
        if psvIdx(1) < 2^constls{1,curLevl}
            psvIdx(1) = psvIdx(1) + 1;
            psv(1) = constls{2,curLevl}(psvIdx(1));
            psvBitLabels(1:constls{1,curLevl}) = constls{3,curLevl}(psvIdx(1),:);
            traverseDone = false;
            return;
        else
            if curLevl == Nt
                traverseDone =  true;
                return;
            else % up a level
                psvIdx(1) = [];psv(1) = [];
                psvBitLabels(1:constls{1,curLevl})= [];
                curLevl = curLevl+1;
            end

        end
    end

end

function constls = genConstls(moduTypes)
% constls = genConstls(moduTypes) output the constellation info of a cell
% list of moduTypes.

    constls = cell(3,length(moduTypes));
    for i = 1: length(moduTypes)
        moduType = moduTypes{i};
        switch lower(moduType)
        case 'bpsk'
            M =2;
        case 'qpsk'
            M = 4;
        case '16qam'
            M = 16;
        case  '64qam'
            M = 64;
        case '256qam'
            M = 256;
        end
    
        K = log2(M);
        constls{1,i} = K;
    
        bitLabels = de2bi(0:M-1,'left-msb');
        constls{3,i} = bitLabels;
    
        symbBitsIn = bitLabels.';
        constSymbs = nrModuMapper(symbBitsIn(:),lower(moduType));
    
        constls{2,i} = constSymbs;
    
    end
end
