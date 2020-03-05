function [ Gs,Gms ] = srl( numroots,denroots,gain )
%   Draw the symmetric root locus for the transfer function G(s)
%   solve: 1+rG(s)G(-s), for r going from 0 to inf.
Gs=zpk([numroots],[denroots],gain);
[a,b,c,d]=ssdata(Gs);
Gms=zpk(ss(-a,-b,c,d));
rlocus(Gs*Gms);
axis equal;


end

