global sim

switch sim
    case 0
        disp('baseline')
    case 1
        disp('high Cw')
        stateCw='high';
    case 2
        disp('mid Rum')
        Rum=80;
    case 3
        disp('high Rum')
        Rum=230;
    case 4
        disp('mid Rsm')
        Rsm=60;
    case 5
        disp('high Rsm')
        Rsm=120;
    case 6
        disp('0 Apic')
        Apic=0.01;
    case 7
        disp('low Pfrac')
        Pfrac=0.28;
    case 8
        disp('low RR')
        RR=40;
    case 9
        disp('high RR')
        RR=70;
    case 10
        disp('Preterm')
        stateCw='high';
        Rum=230;
        Rsm=120;
        Apic=0.01;
        Pfrac=0.28;
        RR=70;
    case 11
        disp('Unhealthy lung, low Cw')
        Phalf=15;
        Ptau=3;
        beta=0.1;
    case 12
        disp('Unhealthy lung, high Cw')
        stateCw='high';
        Phalf=15;
        Ptau=3;
        beta=0.1;
    case 13
        disp('Obstructed airway')
        Rum=2000;
    otherwise
        disp('inappropriate output')
end
