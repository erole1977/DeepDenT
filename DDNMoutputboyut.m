function [boyut]=DDNMoutputboyut(pbest,p,h,m,q)
% p girdi sayýsý
% h gizli tabaka birim sayýsý
% m dendrit sayýsý
% q gizli tabaka sayýsý
Wx{1,1}=pbest(1:(p*m));
Dx{1,1}=pbest((p*m)+1:2*(p*m));
Wh{1,1}=pbest(2*(p*m)+1:2*(p*m)+m*h);
Dh{1,1}=pbest(2*p*m+m*h+1:2*p*m+2*m*h);
ksoma{1,1}=pbest(2*p*m+2*m*h+1:2*p*m+2*m*h+h);
dsoma{1,1}=pbest(2*p*m+2*m*h+h+1:2*p*m+2*m*h+h+h);
ks{1,1}=pbest(2*p*m+2*m*h+h+h+1);
ds{1,1}=pbest(2*p*m+2*m*h+h+h+2);
nhh=2*p*m+2*m*h+h+h+2;
if q>1
for i=2:q
    Wx{i,1}=pbest(nhh+1:nhh+(p*m));
    Dx{i,1}=pbest(nhh+(p*m)+1:nhh+2*(p*m));
    Wh{i,1}=pbest(nhh+2*(p*m)+1:nhh+2*(p*m)+m*h);   
    Dh{i,1}=pbest(nhh+2*p*m+m*h+1:nhh+2*p*m+2*m*h);
    ksoma{i,1}=pbest(nhh+2*p*m+2*m*h+1:nhh+2*p*m+2*m*h+h);
    dsoma{i,1}=pbest(nhh+2*p*m+2*m*h+h+1:nhh+2*p*m+2*m*h+h+h);
    ks{i,1}=pbest(nhh+2*p*m+2*m*h+h+h+1);
    ds{i,1}=pbest(nhh+2*p*m+2*m*h+h+h+2);
    nhh=nhh+2*p*m+2*m*h+h+h+2;
end
end
boyut=nhh+h+1;
end