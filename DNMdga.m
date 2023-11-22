function [rmsetest,mapetest,msebest,bestkr,rmse,yhatg,yhattestg,rn,mape,itn,p]=DNMdga(xt,pg,h,m,q,ntest)
%xt: time series,pg- girdi sayısı, h- gizli tabaka birim sayısı ,m dendrit sayısı, q-gizli tabaka sayısı 
itr=1000;ps=30;
msebestEski=10^6;kkr=0;cor=0.3;
x0=xt;
n0=length(xt);
xtest=x0((n0-ntest+1):n0);
xt=x0(1:(n0-ntest));
x=(xt-min(xt))/(max(xt)-min(xt));
saat=clock;
rand('seed',saat(6)*10000000);
rn=rand('seed');
p=DDNMoutputboyut(unifrnd(0,1,1000,1),pg,h,m,q);
%p parametre sayısı
n=length(x);
A=unifrnd(0,1,ps,p);
%ps kromozom sayısı
%p gen sayısı
%Pso paramterleri program hızlı çalışması için başta oluşturuluyor;
for k=1:ps
    yhat=DDNMoutput(A(k,:),x,pg,h,m,q);
    nh=length(yhat);
    for i=1:nh
        hata(i)=(yhat(i)-x(n-nh+i))^2;
    end
    mse(k)=mean(hata);
end
MSEegt=min(mse);
for i=1:ps
    if MSEegt==mse(i)
        dd=i;
        break
    end
end
bestkr=A(dd,:);
msebest=mse(dd);
i22=0;
for i1=1:itr
    if i22>=30
        A=unifrnd(0,1,ps,p);
        i22=0;
     end
    for i2=1:ps
        dizi=randperm(ps);
        i3=0;
        while i3<=10
            i3=i3+1;
            if dizi(i3)==i2
               kk(i3)=dizi(i3+1);
               i3=i3+1;
            else
                kk(i3)=dizi(i3);
            end
        end
        i4=0;
        for i5=1:10
            if (kk(i5)<1)||(kk(i5)>30)
            else
                i4=i4+1;
                kk2(i4)=kk(i5);
            end
            if i4==3
                break
            end
        end
        %mutasyon
        Ayeni=0.8*(A(kk2(1),:)-A(kk2(2),:))+A(kk2(3),:);
        %Çaprazlama
        for i6=1:p
            if rand<cor
                Ayeni2(i6)=Ayeni(i6);
            else
                Ayeni2(i6)=A(i2,i6);
            end
        end
        yhat=DDNMoutput(Ayeni2,x,pg,h,m,q);
        nh=length(yhat);
        for i=1:nh
            hata(i)=(yhat(i)-x(n-nh+i))^2;
        end
        mseyeni=mean(hata);
        if mseyeni<mse(i2)
            A(i2,:)=Ayeni2;
            mse(i2)=mseyeni;
        end   
    end    
    MSEegt=min(mse);
    for i=1:ps
        if MSEegt==mse(i)
         dd=i;
         break
        end
    end
    msebest=mse(dd);
    bestkr=A(dd,:);
    if abs((msebestEski-msebest)/msebest)<10^-3
        kkr=kkr+1;       
    end
    if kkr>50
        break
    end
    msebestEski=msebest;
end
itn=i1;
        yhat=DDNMoutput(bestkr,x,pg,h,m,q);
        nh=length(yhat);
        yhatg=(yhat)*(max(xt)-min(xt))+min(xt);
        for i=1:nh
            hata(i)=(yhatg(i)-xt(n-nh+i))^2;
            hata2(i)=abs((yhatg(i)-xt(n-nh+i))/xt(n-nh+i));
        end
    rmse=mean(hata)^0.5;
    mape=mean(hata2);
    x=(x0-min(x0))/(max(x0)-min(x0));
    yhattum=DDNMoutput(bestkr,x,pg,h,m,q);
    nt=length(yhattum);
    yhattest=yhattum((nt-ntest+1):nt);
    yhattestg=(yhattest)*(max(xt)-min(xt))+min(xt);
    for i=1:ntest
        hata3(i)=(xtest(i)-yhattestg(i))^2;
        hata4(i)=abs((xtest(i)-yhattestg(i))/xtest(i));
    end
    rmsetest=(mean(hata3))^0.5;
    mapetest=mean(hata4);
end