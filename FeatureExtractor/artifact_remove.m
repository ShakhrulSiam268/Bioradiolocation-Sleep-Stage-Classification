function [xin_art]=artifact_remove(SB,t_hold)

SB2=SB;
win_size=5;
step=2;
total_entropy=0;
count=0;

n=length(SB)/1500;


for i=1:step:length(SB)-10
    startpoint=i;
    endpoint=i+win_size;
    entropy=0;
    
    for j=startpoint:endpoint
        entropy=entropy+SB(j)*SB(j)*log10(SB(j)*SB(j));
    end
    total_entropy=total_entropy+entropy;
    count=count+1;
    ent(count)=entropy;
end

threshold=(t_hold*total_entropy)/count;


for i=1:step:length(SB)-10
    startpoint=i;
    endpoint=i+win_size;
    entropy=0;
    
    for j=startpoint:endpoint
        entropy=entropy+SB(j)*SB(j)*log10(SB(j)*SB(j));
    end
    
    if (entropy)>threshold
        for jj=startpoint:endpoint
            SB2(jj)=0;
        end
    end
end

%% replace IAI elements + z-normalization + flipping
disp('replace IAI elements + z-normalization + flipping')


IAI(1,1)=1;
diffSB=SB-SB2;
count=1;
for i=1:n*1500-1
    
    if diffSB(i)==0 && diffSB(i+1)~=0
        IAI(count,2)=i;
        count=count+1;
        
    else if diffSB(i)~=0 && diffSB(i+1)==0
            IAI(count,1)=i;
            
        end
    end
end

[a b]=size(IAI);
IAI(a,2)=n*1500;
SB3=SB2;

for i=1:a
    
    for jj=1:16
        
        tempdat(jj,:)=SB(1,IAI(i,1):IAI(i,2));
        temp_mean(jj)=mean(tempdat(jj,:));
    end
    
    
    [val index]=max(temp_mean);
    
    
    zz=tempdat(index,:);
    z_norm=(zz-mean(zz))/(sqrt(var(zz)));
    
    
    [pky,pkix] = findpeaks( z_norm);
    [vly,vlix] = findpeaks(-z_norm);
    if length(pkix)<length(vlix)
        z_norm=-z_norm;
    end
    
    
    SB3(IAI(i,1):IAI(i,2))=z_norm;
    
    clear tempdat temp_mean
    
    xin_art=SB3;
end



end