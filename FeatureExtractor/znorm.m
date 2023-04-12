function [x1]=znorm(x)
[ss1 ss2]=size(x);

x1=reshape(x,[1,ss1*ss2]);
x1=(x1-mean(x1))/std(x1);

x1=reshape(x1,[ss1,ss2]);


end