function[p_iqr_med t_iqr_med Pse Tse PTdiff Vbr Vin Vex Fbr Fin Fex RTfr]=depth_feat(s)

factor=2.5;


%win=12;
%s=xin(:,index-win:index+win);
%s=reshape(s,[1,1500*(2*win+1)])';

s=highpass(s,0.15,50);
%s=s/std(s);

[pky,pkixx] = findpeaks( s);
[vly,vlixx] = findpeaks(-s);
s=s(vlixx(2)-1:vlixx(end)+1);

[pky,pkixx] = findpeaks( s);
[vly,vlixx] = findpeaks(-s);
vly=-vly;

x=[1:length(s)];

%% Rectify signal
clear s_rec
for i=1:length(s)
    if s(i)>=0
        s_rec(i)=s(i);
    else s_rec(i)=-s(i);
    end
end
s_rec=s_rec';

%%

pv_mat=[];

br_len=min(length(pkixx),length(vlixx));

for i=1:br_len-1
    
    pv_mat(i,1)=vlixx(i);
    pv_mat(i,2)=pkixx(i);
    pv_mat(i,3)=vlixx(i+1);
    
    pv_mat(i,4)=vly(i);
    pv_mat(i,5)=pky(i);
    pv_mat(i,6)=vly(i+1);
    
end
[rw clm]=size(pv_mat);

pv_length=[];

for i=1:rw
    pv_length=[pv_length pv_mat(i,5)-pv_mat(i,4)];
    pv_length=[pv_length pv_mat(i,5)-pv_mat(i,6)];
    
end
th=median(pv_length)/factor;

pv_track=[];

pvc=1;
for i=1:2:2*br_len-1
    
    pv_track(1,i)=vlixx(pvc);
    pv_track(2,i)=vly(pvc);
    pv_track(3,i)=0;
    
    pv_track(1,i+1)=pkixx(pvc);
    pv_track(2,i+1)=pky(pvc);
    pv_track(3,i+1)=1;
    pvc=pvc+1;
    
end

[rw clm]=size(pv_track);
pv_t_valid=pv_track(:,1);
parity=0;

for i=2:clm-1
    
    if (abs(pv_track(2,i)-pv_track(2,i-1))>th || abs(pv_track(2,i)-pv_track(2,i+1))>th ) && (parity~=pv_track(3,i))
        pv_t_valid=[pv_t_valid pv_track(:,i)];
        parity=pv_track(3,i);
    end
    
end

% val=skewness(sort(pv_length));

%% Area find

[rw clm]=size(pv_t_valid);
val_index=[];
val_amp=[];
pk_index=[];
pk_amp=[];


for i=1:clm-1
    
    if pv_t_valid(3,i)==0
        val_index=[val_index pv_t_valid(1,i)];
        val_amp=[val_amp pv_t_valid(2,i)];
    else pk_index=[pk_index pv_t_valid(1,i)];
        pk_amp=[pk_amp pv_t_valid(2,i)];
        
    end
end

eff_len=min(length(val_amp),length(pk_amp));

brth_mat_valid=[];


%%
jj=1;
for i=1:eff_len-1
    
    brth_mat_valid(i,1)=val_index(i);
    brth_mat_valid(i,2)=pk_index(i);
    brth_mat_valid(i,3)=val_index(i+1);
    
    brth_mat_valid(i,4)=val_amp(i);
    brth_mat_valid(i,5)=pk_amp(i);
    brth_mat_valid(i,6)=val_amp(i+1);
    
    width=brth_mat_valid(i,3)-brth_mat_valid(i,1);
    
    w_in=brth_mat_valid(i,2)-brth_mat_valid(i,1);
    w_ex=brth_mat_valid(i,3)-brth_mat_valid(i,2);
    
    
    %if (brth_mat(i,3)-brth_mat(i,2))>(brth_mat(i,2)-brth_mat(i,1))
    amp=brth_mat_valid(i,5)-brth_mat_valid(i,4);
    
    
    brth_mat_valid(i,7)=amp;
    brth_mat_valid(i,8)=width;
    brth_mat_valid(i,11)=w_in;
    brth_mat_valid(i,12)=w_ex;
    
    jj=jj+1;
end


%% Area

[rw clm]=size(brth_mat_valid);
modex=[];
valid_pk=1;

for k1 = 1:rw
    if brth_mat_valid(k1,4)<0 && brth_mat_valid(k1,5)>0 && brth_mat_valid(k1,6)<0
        
        ixrng1 = brth_mat_valid(k1,1):brth_mat_valid(k1,2);
        ixrng2 = brth_mat_valid(k1,2):brth_mat_valid(k1,3);
        int_s1(k1) = trapz(x(ixrng1), s_rec(ixrng1));
        int_s2(k1) = trapz(x(ixrng2), s_rec(ixrng2));
        %         brth_mat_valid(k1,9)=int_s1(k1);
        %         brth_mat_valid(k1,10)=int_s2(k1);
        modex(valid_pk,1)=brth_mat_valid(k1,1);
        modex(valid_pk,2)=brth_mat_valid(k1,2);
        modex(valid_pk,3)=brth_mat_valid(k1,3);
        modex(valid_pk,4)=int_s1(k1);
        modex(valid_pk,5)=int_s2(k1);
        modex(valid_pk,6)=brth_mat_valid(k1,11);
        modex(valid_pk,7)=brth_mat_valid(k1,12);
        modex(valid_pk,8)=brth_mat_valid(k1,4);
        modex(valid_pk,9)=brth_mat_valid(k1,5);
        modex(valid_pk,10)=brth_mat_valid(k1,6);
        
        valid_pk=valid_pk+1;
        
    end
    brth_mat_valid(k1,9)=0;
    brth_mat_valid(k1,10)=0;
    
    
end
%% graph area
% cm=[1 0 0;0 1 1];
%
% figure(1)
% plot(x, s)
% hold on
% for k1 = 1:rw
%
%     ixrng1 = [brth_mat_valid(k1,1):brth_mat_valid(k1,2)]';
%     ixrng2 = [brth_mat_valid(k1,2):brth_mat_valid(k1,3)]';
%     hp{k1} = patch(x([ixrng1; flipud(ixrng1)])', [s(ixrng1); zeros(size(ixrng1))]', cm(1,:));
%     hp{k1} = patch(x([ixrng2; flipud(ixrng2)])', [s(ixrng2); zeros(size(ixrng2))]', cm(2,:));
%
% end
% hold off
% grid on


%% breathing features

if valid_pk~=1
    
    p_iqr_med=median(modex(:,9))/(prctile(modex(:,9),75)-prctile(modex(:,9),25));
    t_iqr_med=-median(modex(:,8))/(prctile(modex(:,8),75)-prctile(modex(:,8),25));
    
    
    %     Se=sampen2(s,100)
    
    Pse=sampen2(modex(:,9),2);
    Tse=sampen2(modex(:,8),2);
    
    p_isn=isnan(Pse);
    t_isn=isnan(Tse);
    if p_isn==1
        Pse=0;
    end
    if t_isn==1
        Tse=0;
    end
    
    PT=[modex(:,9)-modex(:,8) ; modex(:,9)-modex(:,10)];
    PTdiff=median(PT);
    [rw clm]=size(modex);
    
    depth_mat=[];
    c1=zeros(1,rw);
    t1=zeros(1,rw);
    for kk=1:rw
        c1(kk)=c1(kk)+sum([s(modex(kk,1):modex(kk,3))]);
        t1(kk)=length(s(modex(kk,1):modex(kk,3)));
        %         figure
        %         plot([s(modex(kk,1):modex(kk,3))])
        %         ylabel(c1(kk))
    end
    
    c2=zeros(1,rw);
    t2=zeros(1,rw);
    for kk=1:rw
        c2(kk)=c2(kk)+sum([s(modex(kk,1):modex(kk,2))]);
        t2(kk)=length(s(modex(kk,1):modex(kk,2)));
    end
    
    c3=zeros(1,rw);
    t3=zeros(1,rw);
    for kk=1:rw
        c3(kk)=c3(kk)+sum([s(modex(kk,2):modex(kk,3))]);
        t3(kk)=length(s(modex(kk,2):modex(kk,3)));
    end
    
    Vbr=median(c1);
    Vin=median(c2);
    Vex=median(c3);
    
    Fbr=median(c1./t1);
    Fin=median(c2./t2);
    Fex=median(c3./t3);
    
    if Fex~=0
        RTfr=Fin/Fex;
    else RTfr=0;
    end
else
    %     Se=0;
    Pse=0;
    Tse=0;
    PTdiff=0;
    Vbr=0;
    Vin=0;
    Vex=0;
    Fbr=0;
    Fin=0;
    Fex=0;
    
    RTfr=0;
    p_iqr_med=0;
    t_iqr_med=0;
    
    
    
    
    
end

end