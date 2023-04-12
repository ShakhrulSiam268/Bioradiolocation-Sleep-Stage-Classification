function[median_amp median_wid] = Resp(s)

[pky,pkixx] = findpeaks( s);
[vly,vlixx] = findpeaks(-s);


if pkixx(1)>vlixx(1)
    vlixx(1)=[];
    vly(1)=[];
    
end


vly=-vly;
plot(s,'b')
br_len=min(length(pkixx),length(vlixx));

% for i=1:br_len
%     hold on
%     plot(pkixx(i),pky(i),'.r','Markersize',15)
%     hold on
%     plot(vlixx(i),vly(i),'.g','Markersize',15)
%     
% end
amp=[];
wid=[];
for i=1:br_len-1
    
    amp=[amp (pky(i)-vly(i)) (pky(i+1)-vly(i)) ];
    wid=[wid (vlixx(i)-pkixx(i)) (vlixx(i+1)-pkixx(i))];
    
end
med_amp=median(amp);
med_wid=median(wid);

factor=0.15;
amp2=[];
wid2=[];

for i=1:length(amp)-1
    
    if amp(i)>med_amp*factor
        amp2=[amp2 amp(i)];
    end
end

for i=1:length(wid)-1
    
    if wid(i)>med_wid*factor
        wid2=[wid2 wid(i)];
    end
end

median_amp=median(amp2);
median_wid=median(wid2);

end





