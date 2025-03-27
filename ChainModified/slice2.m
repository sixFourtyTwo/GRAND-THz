% -------------------------------------------------------------------------
% (2^j)-QAM slicer, j=0,2,4,etc (modify to support other constellations)
% -------------------------------------------------------------------------
% r: complex point in space
% d: minimum distance between two adjacent constellation points (=2)
% Q: order of QAM constellation
% s: nearest QAM point to r in lattice
% -------------------------------------------------------------------------
% Round each dimension to nearest (infinite) lattice point, then saturate

% Round each dimension to nearest (infinite) lattice point, then saturate

function s = slice2(r,d,Q)
if(Q==2)
    const = (1/sqrt(2*(2-1)/3))*[ +1+1i -1-1i ];
    dist1=(real(r)-real(const(1)))^2+(imag(r)-imag(const(1)))^2;
    dist2=(real(r)-real(const(2)))^2+(imag(r)-imag(const(2)))^2;

    if(dist1<=dist2)
        s=const(1);
    else
        s=const(2);
    end

else
    maxL = d*(sqrt(Q)-1)/2; % max dist of a lattice point from origin along each dimension
    sr = real(r); srq = d*(floor(sr/d)+1/2); srq(srq>maxL) = maxL; srq(srq <-maxL) = -maxL;
    si = imag(r); siq = d*(floor(si/d)+1/2); siq(siq>maxL) = maxL; siq(siq <-maxL) = -maxL;
    s = srq + 1i*siq;
end


