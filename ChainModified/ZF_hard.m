function x_out=ZF_hard(p,Hmat,yvec) %TYPE
    x_out=complex(zeros(p.Mr*p.Nr*p.frmLen/p.sq,1),zeros(p.Mr*p.Nr*p.frmLen/p.sq,1));
    for n=1:p.frmLen/(p.sq) % number of tones to detect
        H = Hmat(((n-1)*p.Mr*p.Nr+1):(n*p.Mr*p.Nr),1:p.Nt*p.Mt); y = yvec(((n-1)*p.Mr*p.Nr+1):(n*p.Mr*p.Nr),1);
        %HH=Hmat(((n-1)*p.Mr*p.Nr+1):(n*p.Mr*p.Nr),1:p.Nt*p.Mt);
        x_hat=((H')*H)\(H'*y);
        for i=1:length(x_hat)
             x_hat(i)=slice2(x_hat(i),2*(1/sqrt(2*(p.Q(1)-1)/3)),p.Q(1));
        end
        x_out(((n-1)*p.Mr*p.Nr+1):(n*p.Mr*p.Nr),1)=x_hat;
    end
end