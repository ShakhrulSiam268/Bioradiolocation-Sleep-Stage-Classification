function[dtw_val dfw_val] = Dynamic_Wrapping(signal)


%% Downsampling


down_ratio=5;
xin=signal;

[rw cl]=size(signal);
nx=cl;
xin_d=[];   %downsampled signal

for i=1:nx
    tmp=[];
    for j=1:down_ratio:rw
        tmp=[tmp;xin(j,i)];
    end
    xin_d(:,i)=tmp;
end

%% DTW

epc=10;
feat=[]

for i=1:nx
    a=xin_d(:,i)';
    if i<=epc+1
        b=xin_d(:,1:2*epc+1);
        b(:,i)=[];
    else if i>(nx-epc)
            b=xin_d(:,nx-2*epc:nx);
            b(:,2*epc-nx+i)=[];
        else b=xin_d(:,i-epc:i+epc);
            b(:,epc+1)=[];
        end
    end
    
    [rw clm]=size(b);
    b=reshape(b,[1,rw*clm]);
    
    feat(i)=DTW_custom(a,b);
    
    i
end

%% DFW
DFW_mat=[];
DFW_mat2=[];

for i=1:nx
    
    if i==1 || nx
        a=xin(:,i);
        a2=a;
    else
        a1=xin(751:1500,i-1);
        a2=xin(:,i);
        a3=xin(1:750,i+1);
        
        a=[a1;a2;a3];
    end
    
    FT=fft(hanning(length(a)).*a);
    FTA=abs(FT);
    FW=FTA(1:40);
    DFW_mat(i,:)=FW;
    
    
    FT2=fft(hanning(length(a2)).*a2);
    FTA2=abs(FT2);
    FW2=FTA2(1:40);
    DFW_mat2(i,:)=FW2;
    
    
    
    i
end
%%

bbf=reshape(DFW_mat2,[1,nx*40]);

for i=1:nx
    
    a=DFW_mat(i,:);
    
    
    feat2(i)=DFW_custom(a,bbf);
    i
end

dtw_val=feat;
dfw_val=feat2;

end