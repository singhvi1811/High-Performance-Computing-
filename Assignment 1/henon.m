% The henon map function is applied to one entire column vector which contains n
% values per iteration.The basin matrix is created and updated after each
% each updated value of the n by 2 vector.
clear 
close all
n=1000;
xmin = -3;
xmax = 3;
ymin = -3;
ymax = 3;
x0=linspace(xmin,xmax,n); %n initial values of x1 
y0=linspace(ymin,ymax,n);%n initial values of x2
param=struct('a',2.12,'b',-0.3,'K',50,'L',50); %struct that defines all the constant parameters required.
tic
basin=zeros(1000,1000); %initalizing the basin matrix 
for i=1 :n
    ytemp=y0(i)*ones(n,1); %generates one column at a time 
    X=[x0',ytemp]; %initial value vector (nx2)
        %this loop find out the next vector and updates the basin matrix whenever the elemnts go above the lockout value
        for j=1:param.K
            Y(:,2)=X(:,1);
            Y(:,1)=param.a-X(:,1).^2+param.b*X(:,2);
            X=Y; %updates(overwrites) the current value of X vector    
            ix=find(sqrt(Y(:,1).^2+Y(:,2).^2)>param.L); %check for the norm of the vectors and returns the indices for norms greter than lockout value
            basin(i,ix)=1; %updates the basin matrix values to 1 wherever the norm croses the lockout value
        end            
end
surf(x0,y0,basin,'edgecolor','none');
view(0,90);
grid on
colormap jet
toc
