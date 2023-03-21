function [angle]=pk(Vab,Vrc,T,tstep)
pk1=findpeaks(Vab);
pk2=findpeaks(Vrc);
if isempty(pk1)==1
    pk1=max(Vab);
    idx_pk1=(find(Vab==pk1(end)));
else
    idx_pk1=(find(Vab==pk1(end)));
end
if isempty(pk2)==1
    pk2=max(Vrc);
    idx_pk2=(find(Vrc==pk2(end)));
else
    idx_pk2=(find(Vrc==pk2(end)));
end
time_diff=abs(idx_pk1-idx_pk2)*tstep;
angle=time_diff*360/T;
end

