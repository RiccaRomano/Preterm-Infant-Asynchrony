for n=1:11
    currentFile = sprintf('biomarkers%d.mat',n-1);
    load(currentFile)
    VAmaxes(n)=VAmax*1000;
    VAmins(n)=VAmin*1000;
    Pelmaxes(n)=maxPel;
    Pelmins(n)=minPel;
    Pplmaxes(n)=Pplmax;
    Pplmins(n)=Pplmin;
    Vdotmaxes(n)=Vdotmax*1000;
    Vdotmins(n)=Vdotmin*1000;
end
close all
Pao=[0.01 1 2 3 4 5 6 7 8 9 10]
figure
plot(Pao,VAmaxes,'*-',Pao,VAmins,'*-',Pao,VAmaxes-VAmins,'*-')
xlabel('Pao'); ylabel('VA')
legend('max','min','difference')
figure
plot(Pao,Pelmaxes,'*-',Pao,Pelmins,'*-',Pao,Pelmaxes-Pelmins,'*-')
xlabel('Pao'); ylabel('Pel')
legend('max','min','difference')
figure
plot(Pao,Pplmaxes,'*-',Pao,Pplmins,'*-',Pao,Pplmaxes-Pplmins,'*-')
xlabel('Pao'); ylabel('Ppl')
legend('max','min','difference')
figure
plot(Pao,Vdotmaxes,'*-',Pao,Vdotmins,'*-',Pao,Vdotmaxes-Vdotmins,'*-')
xlabel('Pao'); ylabel('Vdot')
legend('max','min','difference')