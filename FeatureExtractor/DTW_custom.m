function [DTW_val]= DTW_custom(a,bb)


step=15;
step2=10;

ll=length(a);
dat=[];
i_trc=[];
for i=1:step:length(bb)-ll+1
    b=bb(i:i+ll-1);

    dat=[dat dtw(a,b)];
    i_trc=[i_trc i];
end

[val,ind]=min(dat);
index=i_trc(ind);


DTW_val=val;

end
