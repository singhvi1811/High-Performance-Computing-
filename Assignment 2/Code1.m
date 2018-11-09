clear 
a=[1 0 0 0;1 2 0 0; 1 2 3 0; 1 2 3 4];
b=[1; 5; 14; 30];
n=length(b);
x=zeros(n,1);
b(1)=b(1)/a(1,1);
for i=1:n-1
   b(i+1:n)=b(i+1:n)-a(i+1:n,i)*b(i);
   b(i+1)=b(i+1)/a(i+1,i+1);
end

a1=[2 2 3; 0 1 2; 0 0 1];
b1=[14; 5; 1];
b11=b1;
n1=length(b1);
b1(n1)=b1(n1)/a1(n1,n1);
for j=n1:-1:2
    b1(j) = b1(j)/a1(j,j);
    b1(1:j-1) = b1(1:j-1)-a1(1:j-1,j)*b1(j);
    
end
sol=inv(a1)*b11;