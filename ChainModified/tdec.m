% turbo decoder: generates extrinsic llrs for encoded bits
function [y, llr, crcpass, iter] = tdec(obj,x)
interlvrIndices = obj.InterleaverIndices;
blkLen = length(interlvrIndices);    
pK = log2(obj.TrellisStructure.numInputSymbols);
assert(pK==1, 'tdec function:','Only rate 1/n trellises are supported.');
pN = log2(obj.TrellisStructure.numOutputSymbols);
pMLen = log2(obj.TrellisStructure.numStates);
pNumTails = pMLen*(pN);

% Bit order
dIdx   = (2*pN-1)*blkLen;

yD     = reshape(x((1:dIdx).', 1), 2*pN-1, blkLen);
lc1D   = yD(1:pN, :);
y1T    = x(dIdx + (1:pNumTails).', 1);
Lc1_in = [lc1D(:); y1T];

Lu1_in = zeros(blkLen+pMLen, 1);

lc2D1  = zeros(1, blkLen);
lc2D2  = yD(pN+1:2*pN-1, :);
lc2D   = [lc2D1; lc2D2];
y2T    = x(dIdx + pNumTails + (1:pNumTails).', 1);
Lc2_in = [lc2D(:); y2T];

% Turbo Decode: Lu used in turbo-decoding; Lc's used in iterative detection
out1 = zeros(blkLen, 1);  
for iterIdx = 1:obj.NumIterations
    release(obj.app1); % Release the object
Lc1_in = real(Lc1_in);
    [Lu1_out,Lc1_out] = step(obj.app1, Lu1_in, Lc1_in);
    tmp = Lu1_out((1:blkLen).', 1);
    tmp2 = tmp(:);
Lc2_in = real(Lc2_in);

    [Lu2_out,Lc2_out] = step(obj.app2, [tmp2(interlvrIndices(:)); zeros(pMLen,1)], Lc2_in);
    
    out1(interlvrIndices(:), 1) = Lu2_out((1:blkLen).', 1);
    Lu1_in = [out1; zeros(pMLen,1)];
    % crc check
    llr1 = out1 + tmp2; y = (llr1>=0);
    [~, crcfail] = step(obj.hCRCDet,y); if crcfail==0, break; end  
end
out2(interlvrIndices(:), 1) = Lc2_out((1:2:2*blkLen).', 1);
llr1 = Lc1_out((1:2:2*blkLen).', 1) + out2; % combine Lc1 from both encoders
llr2 = Lc1_out((2:2:2*blkLen).', 1); 
llr3 = Lc2_out((2:2:2*blkLen).', 1); 
llr = reshape([llr1 llr2 llr3]',3*blkLen,1);
crcpass = 1 - crcfail;
iter = iterIdx;
end