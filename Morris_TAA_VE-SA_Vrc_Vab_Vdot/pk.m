function [angle]=pk(Vab,Vrc,T,tstep)
pk1=findpeaks(Vab);
idx_pk1=(find(Vab==pk1(end)));
pk2=findpeaks(Vrc);
idx_pk2=(find(Vrc==pk2(end)));
time_diff=abs(idx_pk1-idx_pk2)*tstep;
angle=time_diff*360/T;
% lag=idx_pk1-idx_pk2
%disp(angle);
end

