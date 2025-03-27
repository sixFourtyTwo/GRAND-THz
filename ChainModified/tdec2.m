% turbo decoder: generates extrinsic llrs for encoded bits
function [list,y, llr, crcpass, iter] = tdec2(obj,x,p)
interlvrIndices = obj.InterleaverIndices;

blkLen = length(interlvrIndices);  
% blkLen=256;

pK = log2(obj.TrellisStructure.numInputSymbols);
% pK=1;

% assert(pK==1, 'tdec function:','Only rate 1/n trellises are supported.');
pN = log2(obj.TrellisStructure.numOutputSymbols);
% pN=2;

pMLen = log2(obj.TrellisStructure.numStates);
% pMLen=3;

% number of tailing bits
pNumTails = pMLen*(pN);
% pNumTails=6;

list = zeros(p.dataLen*(2*length(p.codegen)-1),1+obj.NumIterations);
% Bit order

% size of the encoded bits
dIdx   = (2*pN-1)*blkLen;
% dIdx=768;

%yD is the input shrinked to the encoding size and reshaped into blocks each of input
% size
yD     = reshape(x((1:dIdx).', 1), 2*pN-1, blkLen);
% yD of size (3,256)

lc1D   = yD(1:pN, :);
% lc1D of size (2,256)

y1T    = x(dIdx + (1:pNumTails).', 1);
% The first 6 tail bits

Lc1_in = [lc1D(:); y1T];
% combining the first decoder input: sys bits and tailing bits

Lu1_in = zeros(blkLen+pMLen, 1);

lc2D1  = zeros(1, blkLen);

lc2D2  = yD(pN+1:2*pN-1, :);
% Last row of yD

lc2D   = [lc2D1; lc2D2];
% its first line 0, the systematic for the second decoder.
y2T    = x(dIdx + pNumTails + (1:pNumTails).', 1);
% second tail bits corresponding to second encoder

Lc2_in = [lc2D(:); y2T];
% second decoder input combined

outt2(interlvrIndices(:), 1) = Lc2_in((1:2:2*blkLen).', 1);
lllr1 = Lc1_in((1:2:2*blkLen).', 1) + outt2; % combine Lc1 from both encoders
lllr2 = Lc1_in((2:2:2*blkLen).', 1);
lllr3 = Lc2_in((2:2:2*blkLen).', 1);
lllr = reshape([lllr1 lllr2 lllr3]',3*blkLen,1);

list(:,1) = lllr;


% Turbo Decode: Lu used in turbo-decoding; Lc's used in iterative detection
out1 = zeros(blkLen, 1); 
indx=2;
for iterIdx = 1:obj.NumIterations
    [Lu1_out,Lc1_out] = step(obj.app1, Lu1_in, Lc1_in);
    % Lu1 is the extrinsic equivalent to E in the diagram
    % Lc1 is the equivalent to D1, the result
    tmp = Lu1_out((1:blkLen).', 1);
    tmp2 = tmp(:);
    % size (256,1)
    [Lu2_out,Lc2_out] = step(obj.app2, [tmp2(interlvrIndices(:)); zeros(pMLen,1)], Lc2_in);
    out1(interlvrIndices(:), 1) = Lu2_out((1:blkLen).', 1);
    Lu1_in = [out1; zeros(pMLen,1)];
    % crc check
    llr1 = out1 + tmp2; y = (llr1>=0);
    [~, crcfail] = step(obj.hCRCDet,y); %if crcfail==0, break; end 
    out2(interlvrIndices(:), 1) = Lc2_out((1:2:2*blkLen).', 1);
    % double check the jump
    llr1 = Lc1_out((1:2:2*blkLen).', 1) + out2; % combine Lc1 from both encoders
    llr2 = Lc1_out((2:2:2*blkLen).', 1); 
    llr3 = Lc2_out((2:2:2*blkLen).', 1); 
    llr = reshape([llr1 llr2 llr3]',3*blkLen,1);
    list(:,indx) = llr;
    indx = indx + 1;
end

crcpass = 1 - crcfail;
iter = iterIdx;

end