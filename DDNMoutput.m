function [yhat]=DDNMoutput(pbest,xt,p,h,m,q)
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
V=pbest(nhh+1:nhh+h);
V=reshape(V,[h,1]);
by=pbest(nhh+h+1);
Wx{1,1}=reshape(Wx{1,1},[p,m]);
Dx{1,1}=reshape(Dx{1,1},[p,m]);
Wh{1,1}=reshape(Wh{1,1},[h,m]);
Dh{1,1}=reshape(Dh{1,1},[h,m]);
ksoma{1,1}=reshape(ksoma{1,1},[1,h]);
dsoma{1,1}=reshape(dsoma{1,1},[1,h]);
if q>1
for i=2:q
    Wx{i,1}=reshape(Wx{i,1},[p,m]);
    Dx{i,1}=reshape(Dx{i,1},[p,m]);
    Wh{i,1}=reshape(Wh{i,1},[h,m]);
    Dh{i,1}=reshape(Dh{i,1},[h,m]);
    ksoma{i,1}=reshape(ksoma{i,1},[1,h]);
    dsoma{i,1}=reshape(dsoma{i,1},[1,h]);
end
end
X=lagmatrix(xt,[1:p,0]);
Xgirdi=X(p+1:end,1:p);
f=@(x) 1./(1+exp(-x));
%nogr: Öðrenme örneði sayýsý
nogr=size(Xgirdi,1);
for t=1:nogr-h
    if t==1
    for i1=1:q
        if i1==1
           for j=1:h
                if j==1
                    Y1=f(Wx{i1,1}.*Xgirdi(t+j-1,:)'+Dx{i1,1});
                    V=sum(prod(Y1));
                    H{t,j}=f(-ks{i1,1}*V-ds{i1,1});
                else
                    Y1=f(Wx{i1,1}.*Xgirdi(t+j-1,:)'+Dx{i1,1});
                    Y2=f(Wh{i1,1}.*H{t,j-1}'+Dh{i1,1});
                    V=sum((prod(Y1)).*(prod(Y2)));
                    H{t,j}=f(-ks{i1,1}*V-ds{i1,1});
                end
           end
        else
            for j=1:h
                if j==1
                    Y1=f(Wx{i1,1}.*H{t,j}'+Dx{i1,1});
                    V=sum(prod(Y1));
                    H{t,j}=f(-ks{i1,1}*V-ds{i1,1});
                else
                    Y1=f(Wx{i1,1}.*H{t,j}'+Dx{i1,1});
                    Y2=f(Wh{i1,1}.*H{t,j-1}'+Dh{i1,1});
                    V=sum((prod(Y1)).*(prod(Y2)));
                    H{t,j}=f(-ks{i1,1}*V-ds{i1,1});
                end
           end
        end
    end
    else
    for i1=1:q
        if i1==1
           for j=1:h
                if j==1
                    Y1=f(Wx{i1,1}.*Xgirdi(t+j-1,:)'+Dx{i1,1});
                    Y2=f(Wh{i1,1}.*H{t-1,j}+Dh{i1,1});
                    V=sum((prod(Y1)).*(prod(Y2)));
                    H{t,j}=f(-ks{i1,1}*V-ds{i1,1});
                else
                    Y1=f(Wx{i1,1}.*Xgirdi(t+j-1,:)'+Dx{i1,1});
                    Y2=f(Wh{i1,1}.*H{t,j-1}'+Dh{i1,1});
                    V=sum((prod(Y1)).*(prod(Y2)));
                    H{t,j}=f(-ks{i1,1}*V-ds{i1,1});
                end
           end
        else
            for j=1:h
                if j==1
                    Y1=f(Wx{i1,1}.*H{t,j}'+Dx{i1,1});
                    Y2=f(Wh{i1,1}.*H{t-1,j}'+Dh{i1,1});
                    V=sum((prod(Y1)).*(prod(Y2)));
                    H{t,j}=f(-ks{i1,1}*V-ds{i1,1});
                else
                    Y1=f(Wx{i1,1}.*H{t,j}'+Dx{i1,1});
                    Y2=f(Wh{i1,1}.*H{t-1,j}'+Dh{i1,1});
                    V=sum((prod(Y1)).*(prod(Y2)));
                    H{t,j}=f(-ks{i1,1}*V-ds{i1,1});
                end
           end
        end
    end
    end
    yhat(t)=f(H{t,h}*V+by);
end
end