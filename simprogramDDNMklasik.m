function [Stats1,Stats2,ONG,B1,B,RMSETEST,MAPETEST,bestkr,pgbest,hhbest,mmbest,qqbest]=simprogramDDNMklasik(x,nval,ntest)
xval=x(1:end-ntest);
pg=1:8;hh=1:4;mm=1:2;qq=1:4;
kk=0;
B=[1,2,3,4,5,6];
for i1=1:length(pg)
    for i2=1:length(hh)
        for i3=1:length(mm)
            for i4=1:length(qq)
               for i5=1:5
                   kk=kk+1
                   [rmsetest(kk),~,~,~,~,~,~,rn(kk)]=DNMdga(xval,pg(i1),hh(i2),mm(i3),qq(i4),nval);
                   B=[B;[pg(i1),hh(i2),mm(i3),qq(i4),rmsetest(kk),rn(kk)]];
               end
            end
        end
    end
end
B=B(2:end,:);
minB=min(B(:,5));
for i=1:length(B(:,5))
    if minB(1)==B(i,5)
        pgbest=B(i,1);hhbest=B(i,2);mmbest=B(i,3);qqbest=B(i,4);
        break
    end
end
B1=[1,2,3];
for i4=1:30
        [rmsetest(i4),mapetest(i4),~,bestkr,~,~,yhattestg{i4,1},rn(i4)]=DNMdga(x,pgbest,hhbest,mmbest,qqbest,ntest);
        B1=[B1;[rmsetest(i4),mapetest(i4),rn(i4)]];
end  
B1=B1(2:end,:);
Stats1=[mean(B1(:,1)),median(B1(:,1)),std(B1(:,1)),iqr(B1(:,1)),min(B1(:,1)),max(B1(:,1))];
Stats2=[mean(B1(:,2)),median(B1(:,2)),std(B1(:,2)),iqr(B1(:,2)),min(B1(:,2)),max(B1(:,2))];
minB1=min(B1(:,1));
for i=1:length(B1(:,1))
    if minB1(1)==B1(i,1)
       ONG=yhattestg{i,1};
       RMSETEST=B1(i,1);MAPETEST=B1(i,2);
        break
    end
end
ONG=ONG';
end