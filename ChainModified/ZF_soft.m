function LLRout = ZF_soft(p,Hmat,yvec,N0)

q1 = p.q(1); % number of bits in the modulation type
Q1 = p.Q(1); % number of symbols
bits = de2bi(0:p.Q(1)-1,p.q(1),'left-msb'); % the vector of the bits (00/01..)
% disp("the bits vector is");
% disp(bits);
mask_0 = zeros(Q1/q1,q1); mask_1 = zeros(Q1/q1,q1);
% disp("the mask_0 vector is")
% disp(mask_0);
% disp("the mask_1 vector is")
% disp(mask_1);
for i=1:q1
    mask_0(:,i)=find(bits(:,i)==0); 
    mask_1(:,i)=find(bits(:,i)==1); 
end % masks


LLRout = complex(zeros(p.frmLen/p.sq,p.sq),zeros(p.frmLen/p.sq,p.sq));  
for n=1:p.frmLen/(p.sq) % number of tones to detect
    H = Hmat(((n-1)*p.Mr*p.Nr+1):(n*p.Mr*p.Mr),1:p.Nt*p.Mt); y = yvec(((n-1)*p.Mr*p.Nr+1):(n*p.Mr*p.Nr),1);
    xhat=((H')*H)\(H'*y);
    % xhat is the hard decision of z transformed and quantized
    invN0scaled =  diag(1./(N0));  
    % create a matrix (ones(p.Nt*p.Mt,1)*p.syms(1,1:p.Q(1))
    dist1 = invN0scaled*abs(xhat*ones(1,p.Q(1))-ones(p.Nt*p.Mt,1)*p.syms(1,1:p.Q(1))).^2;
    ofst=0;  % only for uniform constellations; modify p.masks0/1 to generalize
    for i=1:p.Nt*p.Mt
        dist = dist1(i,:);%- ((-1+2*bits)*LLRin(n,(1+ofst):(p.q(i)+ofst))'/2)';
        % this (-1+2*bits) is not clear and needs to be removed or replaced
        for j=1:p.q(i), LLRout(n,ofst+j) = min(dist(mask_0(:,j))) - min(dist(mask_1(:,j))); end
        ofst = ofst + p.q(i);
    end
end
end

 