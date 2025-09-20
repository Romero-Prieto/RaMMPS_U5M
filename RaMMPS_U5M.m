clear
pATh                  = "/Users/lshjr3/Documents/RaMMPS/under-five mortality/";
lISt                  = {'RaMMPS','RaMMPSdst','RaMMPSHH','MICSmalawi','DHSmalawi','DHSmalawidst','IGMEmalawi','RaMMPScalls'};
for i = 1:numel(lISt)
    options    = detectImportOptions(char(pATh + string(lISt{i}) + ".csv"));
    for j = 1:numel(options.VariableTypes)
        if isequal(options.VariableTypes{j},'char')
            options.VariableTypes{j} = 'categorical';
        end
        if isequal(options.VariableTypes{j},'datetime') & ~isequal(lISt{i},'RaMMPScalls')
            options.VariableOptions(1,j).InputFormat = 'dd/MM/yyyy';
        end
    end
    sET        = readtable(char(pATh + string(lISt{i}) + ".csv"),options);
    assignin('base',lISt{i},sET);
    clear options j sET
end
clear lISt i
RESolUTioN            = 300;
save(char(pATh + "Results/RaMMPSf.mat"),'RaMMPS','RESolUTioN');
RaMMPS                = RaMMPS(RaMMPS.caseid ~= "0A0" & isundefined(RaMMPS.missing),:);
RaMMPSHH              = RaMMPSHH(RaMMPSHH.caseid ~= "0A0" & isundefined(RaMMPSHH.missing),:);
RaMMPSdst             = RaMMPSdst(isundefined(RaMMPSdst.missing),:);
save(char(pATh + "Results/RaMMPS.mat"));

clear
pATh                  = "/Users/lshjr3/Documents/RaMMPS/under-five mortality/";
load(char(pATh + "Results/RaMMPS.mat"),'RaMMPS','RESolUTioN');
x{1}                  = [(0:7:28)/365.25,[(2:1:12),(15:3:24),(36:12:60)]/12]';
x{2}                  = [string(0);string([(7:7:28)';(2:1:11)';(12:3:24)';(36:12:60)']) + char([kron('d',ones(4,1));kron('m',ones(18,1))])];
x{3}                  = round(x{1}*365.25*24);
n                     = diff(x{1},1);
instrument            = {'TPH','FPH'};
date                  = max(RaMMPS.interview);
Ts                    = {datetime([2014 year(date)]',[1 month(date)]',[1 day(date)]'),datetime([2014 2016]',1,1),datetime([2016 2018]',1,1),datetime([2018 2020]',1,1),datetime([2020 2022]',1,1),datetime([2022 year(date)]',[1 month(date)]',[1 day(date)]')};
sET                   = [RaMMPS.WR(RaMMPS.k == 1),ones(sum(RaMMPS.k == 1),1)];
models                = {'post-strat.','selected'};
R                     = 50;
aGEs                  = [5 16 23];


bLOKs                 = {'SBH','TPH','FPH'};
data                  = RaMMPS(RaMMPS.k == 1,{'sample','caseid','K'});
for i = 1:numel(models)
    rng(0);
    S          = [ones(size(data,1),1),rand(size(data,1),R)];
    for j = 1:numel(bLOKs)
        sEL             = find(data.sample == bLOKs{j});
        W               = sET(sEL,i);
        W               = [0;cumsum(W(1:end - 1))]/sum(W);
        for r = 2:R + 1
            temp     = tabulate([sum(W < S(sEL,r)',1)';numel(sEL) + 1]);
            S(sEL,r) = temp(1:end - 1,2);
            clear temp
            clc;
            r/(R + 1)
        end
    clear W sEL    
    end
    WRaMMPS{i} = S(repelem((1:size(data,1))',data.K),:);
    clear S
end

for i = 1:numel(instrument)
    for j = 1:numel(models)
        h          = j + (i - 1)*numel(models);
        dATa{h}    = {RaMMPS.sample == instrument{i},WRaMMPS{j}};
        for k = 1:numel(Ts)
            date      = eXAcTTime(Ts{k});
            date      = string(instrument{i}) + ", " + string(sprintf('%0.1f',date(1))) + "-" + string(sprintf('%0.1f',date(2)));
            pOPs{k,h} = {char(date);''};
            clear date
        end
    end
end

sex              = [2 1];
for i = 1:numel(dATa)
    pOPs{numel(Ts) + 1,i}    = pOPs{1,i};
    pOPs{numel(Ts) + 1,i}{1} = char(pOPs{1,i}{1} + " (female)");
    pOPs{numel(Ts) + 2,i}    = pOPs{1,i};
    pOPs{numel(Ts) + 2,i}{1} = char(pOPs{1,i}{1} + " (male)");
end

p                     = rand(size(RaMMPS,1),R + 1);
d                     = (p.*RaMMPS.D_min + (1 - p).*RaMMPS.D_max)/365.25;
p                     = rand(size(RaMMPS,1),R + 1);
bR                    = (RaMMPS.birth ~= 'livebirth');
B                     = p.*eXAcTTime(RaMMPS.B_min) + (1 - p).*eXAcTTime(RaMMPS.B_max);
B(bR,:)               = NaN;
bR                    = (RaMMPS.birth ~= 'stillbirth');
sB                    = p.*eXAcTTime(RaMMPS.B_min) + (1 - p).*eXAcTTime(RaMMPS.B_max);
sB(bR,:)              = NaN;

interview             = eXAcTTime(RaMMPS.interview);
D                     = B + d;
O                     = min(D,interview);
p                     = rand(size(RaMMPS,1),R + 1);
ageS                  = p.*RaMMPS.age + (1 - p).*(RaMMPS.age + 1);
ageB                  = ageS - (interview - B);
ageSB                 = ageS - (interview - sB);
sEL                   = (RaMMPS.survival == 'dead' | RaMMPS.survival == 'alive');
        
for i = 1:numel(dATa)
    s = dATa{i}{1};
    w = dATa{i}{2};
    for h = 1:size(pOPs,1)
        if h <= numel(Ts)
            TS            = Ts{h};
            T             = eXAcTTime(TS);
            sEX           = ones(size(RaMMPS,1),1) == 1;
        else
            TS            = Ts{1};
            T             = eXAcTTime(TS);
            sEX           = (RaMMPS.sex == sex(h - numel(Ts)));
        end
        sH              = dATa{i}{1} & sEX;
        
        for j = 1:numel(x{1}) - 1
            a             = max(min(B(sH,:) + x{1}(j),min(O(sH,:),T(2))),T(1));
            o             = max(min(B(sH,:) + x{1}(j + 1),min(O(sH,:),T(2))),T(1)); 
            exposure(j,:) = sum((o - a).*(w(sH,:).*sEL(sH)));
            events(j,:)   = sum((d(sH,:) >= x{1}(j) & d(sH,:) < x{1}(j + 1) & D(sH,:) >= T(1) & D(sH,:) < T(2)).*(w(sH,:).*sEL(sH,:)));
            clear a o
            clc;
            j/(numel(x{1}) - 1)
        end
        
        stillbirths     = sum((RaMMPS.FlagBI(s) ~= 1 & RaMMPS.gestation(s) >= 28 & sB(s,:) >= T(1) & sB(s,:) < T(2)).*w(s,:));
        births          = sum((B(s,:) >= T(1) & B(s,:) < T(2)).*w(s,:));
        stillbirths     = stillbirths./(births + stillbirths);        
        SRB             = sum((B(s,:) >= T(1) & B(s,:) < T(2) & RaMMPS.sex(s) == 1).*w(s,:))./sum((B(s,:) >= T(1) & B(s,:) < T(2) & RaMMPS.sex(s) == 2).*w(s,:));
        m               = events./exposure;
        q               = 1 - [ones(1,R + 1);exp(-cumsum(m.*n))];
        
        TaBle.q{h,i}    = q;
        TaBle.s{h,i}    = 1 - (1 - q([1 2 5],:)).*(1 - stillbirths);
        TaBle.SNR{h,i}  = [TaBle.s{h,i}(1,:);TaBle.q{h,i}(2,:);TaBle.q{h,i}(5,:);TaBle.s{h,i}(1,:)./TaBle.q{h,i}(5,:)/1000;TaBle.q{h,i}(2,:)./TaBle.q{h,i}(5,:)/1000];
        TaBle.SRB{h,i}  = SRB;
        TaBle.tAU{h,i}  = TS;
        TaBle.BTHW{h,i} = births./sum(w(s,:));
        clear T TS sEX events exposure q m births stillbirths SRB sH
    end
    
    A                    = [max(min(RaMMPS.age(dATa{i}{1})),15),min(max(RaMMPS.age(dATa{i}{1})),49)];
    TFR                  = 0;
    for j = A(1):A(2)
        a               = max(min(ageS(s,:) - 3,j + 1),j);
        o               = max(min(ageS(s,:),j + 1),j);
        exposure        = sum((RaMMPS.k(s) == 1).*(o - a).*w(s,:));
        events          = sum((interview(s) - B(s,:) <= 3 & interview(s) - B(s,:) > 0).*(ageB(s,:) >= j & ageB(s,:) < j + 1).*w(s,:));
        TFR             = TFR + events./max(exposure,eps);
        clear events exposure a o
        clc;
        (j - A(1) + 1)/(A(2) - A(1) + 1)
    end
    
    TSR                  = 0;
    for j = A(1):A(2)
        a               = max(min(ageS(s,:) - 3,j + 1),j);
        o               = max(min(ageS(s,:),j + 1),j);
        exposure        = sum((RaMMPS.k(s) == 1).*(o - a).*w(s,:));
        events          = sum((interview(s) - sB(s,:) <= 3 & interview(s) - sB(s,:) > 0).*(ageSB(s,:) >= j & ageSB(s,:) < j + 1).*w(s,:));
        TSR             = TSR + events./max(exposure,eps);
        clear events exposure a o
        clc;
        (j - A(1) + 1)/(A(2) - A(1) + 1)
    end
        
    TaBle.TFR{1,i}       = TFR;
    TaBle.TSR{1,i}       = TSR;
    TaBle.sRB{1,i}       = sum((RaMMPS.birth(s) == 'livebirth' & RaMMPS.sex(s) == 1).*w(s,:))./sum((RaMMPS.birth(s) == 'livebirth' & RaMMPS.sex(s) == 2).*w(s,:)); 
    TaBle.parity{1,i}    = sum((RaMMPS.birth(s) == 'livebirth').*w(s,:))./sum((RaMMPS.k(s) == 1).*w(s,:));
    TaBle.childless{1,i} = sum((RaMMPS.mother(s) ~= 1 & RaMMPS.k(s) == 1).*w(s,:))./sum((RaMMPS.k(s) == 1).*w(s,:))*100;
    clear TFR TSR w s A
end
clear b ageB ageSB ageS bR h i j k r

for i = 1:numel(dATa)
    TaBle.sEX{1,i} = log(1 - TaBle.q{numel(Ts) + 2,i})./log(1 - TaBle.q{numel(Ts) + 1,i});
end

lABelS                = {'Sex','Place of residence','Region','Education','Electricity'};
lABelSd               = {{'female' 'male'} {'urban' 'rural'} {'North' 'Central' 'South'} {'less than complete secondary' 'complete secondary or more'} {'access' 'no access'}}; 
pOPh{1}               = {'All women'};
for i = 1:numel(lABelS)
    pOPh{end + 1,1} = {char(lABelS{i} + ": " + lABelSd{i}{1})};
    for j = 2:numel(lABelSd{i})
        pOPh{end + 1,1} = {char(lABelSd{i}{j})};
    end
end


load(char(pATh + "Results/RaMMPS.mat"),'RaMMPS');
RaMMPS.s              = RaMMPS.sample;
s                     = find(RaMMPS.s == 'FPH');
RaMMPS.s(s)           = 'SBH';
instrument            = {'SBH'};
BrassTimeD            = datetime([kron((2008:2017)',ones(12,1)),kron(ones(2017 - 2008 + 1,1),(1:12)'),15*ones((2017 - 2008 + 1)*12,1)]);
BrassTime             = (2008.0 + 1/24:1/12:2018.0 - 1/24)';
BrassP                = [BrassTime < 2013.0,BrassTime >= 2013.0];
BrassP                = BrassP./sum(BrassP);
CD                    = {'North','South','East','West'};
family                = 1;

for i = 1:numel(instrument)
    for j = 1:size(models,2)
        h             = j + (i - 1)*numel(models);
        H             = h + numel(dATa);
        dATaB{h}      = {RaMMPS.s == instrument{i} & RaMMPS.k == 1,WRaMMPS{j}};
        
        for k = 1:size(BrassP,2)
            date       = string(sprintf('%0.1f',BrassP(:,k)'*BrassTime - 2.5)) + "-" + string(sprintf('%0.1f',BrassP(:,k)'*BrassTime + 2.5));
            pOPs{k,H}  = {char("SBH (Brass: " + string(CD{family}) + "), " + date);''};
            clear date
        end
        pOPs{k + 1,H} = {char("SBH (Brass: " + string(CD{family}) + ")");''};
        clear h H
    end
end

for i = 1:numel(dATaB)
    s                       = dATaB{i}{1};
    w                       = dATaB{i}{2}(s,:);
        
    T                       = eXAcTTime(RaMMPS.interview(s));
    T                       = (T'*w)./sum(w);
    xx                      = RaMMPS.age(s) >= 15:5:45 & RaMMPS.age(s) < 20:5:50;
    W                       = xx'*w;
    B                       = xx'*((RaMMPS.sons(s) + RaMMPS.daughters(s)).*w);
    D                       = xx'*((RaMMPS.sonsD(s) + RaMMPS.daughtersD(s)).*w);
    
    h                       = i + numel(dATa);
    TaBle.sRB{1,h}          = (RaMMPS.sons(s)'*w)./(RaMMPS.daughters(s)'*w);
    children                = RaMMPS.sons(s) + RaMMPS.daughters(s);
    TaBle.parity{1,h}       = sum(children.*w)./sum(w);        
    TaBle.childless{1,h}    = ((children == 0)'*w)./((children >= 0)'*w)*100;
    
    temp                    = BraSsTruSseLl(W,B,D,T);
    for k = 1:numel(CD)
        for j = 1:2
            BraSs_s{k,i}{j} = LinInterPol(flip(temp{k,end}),flip(temp{k,j + 1}),BrassTime);
            BraSs_r{k,i}{j} = flip(temp{k,j + 1});
        end
        BraSs_r{k,i}{j + 1} = flip(temp{k,end});
    end
    
    temp                    = temp(family,:);
    for j = 1:numel(temp) - 1
        q                   = LinInterPol(flip(temp{end}),flip(temp{j}),BrassTime);
        for k = 1:size(BrassP,2)
            if isequal(j,1)
                TaBle.q{k,h} = NaN(max(aGEs),R + 1);
            end
            TaBle.q{k,h}(aGEs(j),:) = BrassP(:,k)'*q;
            TaBle.tAU{k,h}          = BrassP(:,k)'*BrassTime;
        end
        TaBle.q{k + 1,h}{j} = q; 
    end
    TaBle.q{k + 1,h}{j + 1} = BrassTimeD;
    TaBle.tAU{k + 1,h}      = BrassTimeD;
    clear s w W T xx W B D h children temp q
end


load(char(pATh + "Results/RaMMPS.mat"),'DHSmalawi');
mAx                   = median(DHSmalawi.interview);
mIn                   = datetime([year(mAx) - 5,month(mAx),day(mAx)],'Format','dd/MM/yyyy');
dATe{1}               = [mIn mAx];
dATe{2}               = eXAcTTime(dATe{1});
Tdhs                  = dATe{2};
models                = {'Women 15-49','Women 18-49','Women 18-49, mobile owners','Women 18-49, mobile owners'};
dATaDHS{1}            = {DHSmalawi.age >= 15,DHSmalawi.W};
dATaDHS{2}            = {DHSmalawi.age >= 18,DHSmalawi.W};
dATaDHS{3}            = {DHSmalawi.age >= 18 & DHSmalawi.mobile == 1 & DHSmalawi.jure == 1,DHSmalawi.WR};
dATaDHS{4}            = {dATaDHS{3}{1},DHSmalawi.W};

for i = 1:numel(dATaDHS)
    data      = "DHS VII, " + string(sprintf('%0.1f',dATe{2}(1))) + "-" + string(sprintf('%0.1f',dATe{2}(2)));
    h         = i + numel(dATa) + numel(dATaB);
    pOPs{1,h} = {data;models{i}};
    pOPs{2,h} = {data + " (female)";models{i}};
    pOPs{3,h} = {data + " (male)";models{i}};
end
clear mIn mAx models data

rng(0);
s                     = DHSmalawi(DHSmalawi.k == 1, {'cluster','K','woman','Women','iNDeX'});
w                     = rand(size(s,1),R);
w                     = [(1:size(s,1))',ceil(s.Women.*w) + s.iNDeX];
Wdhs                  = NaN(size(DHSmalawi,1),R + 1);
for r = 1:R + 1
    S         = tabulate([w(:,1);w(:,r)]);
    Wdhs(:,r) = repelem(S(:,2) - 1,s.K);
    clc;
    r/(R + 1)
end
clear r s S w

for i = 1:numel(dATaDHS)
    dATaDHS{i}{2} = Wdhs.*dATaDHS{i}{2};
end
clear Wdhs
save(char(pATh + "Results/ReSAmPLiNG.mat"),'dATa','dATaDHS','Ts','-v7.3','-nocompression');

DHSmalawi.date        = eXAcTTime(DHSmalawi.interview);
DHSmalawi.births      = zeros(size(DHSmalawi,1),1);  
DHSmalawi.stillbirths = zeros(size(DHSmalawi.births,1),1);
sET                   = find(~isnan(DHSmalawi.row));
for j = 1:numel(sET)
    cal = char(DHSmalawi.CAL(sET(j)));
    if ismember('T',cal) || ismember('B',cal)
        cal = cal(DHSmalawi.row(sET(j)) + 1:end);
        sT  = find(ismember(cal,'T'));
        for h = 1:numel(sT)
            dT = datetime(year(DHSmalawi.interview(sET(j))),month((DHSmalawi.interview(sET(j)))) - sT(h) + 1,day(DHSmalawi.interview(sET(j))),'Format','dd/MM/yyyy');
            gA = cal(sT(h):min(sT(h) + 6,end));
            if isequal(gA,'TPPPPPP') && dT >= dATe{1}(1) && dT  < dATe{1}(2)
               DHSmalawi.stillbirths(sET(j)) = DHSmalawi.stillbirths(sET(j)) + 1;
            end
            clear dT gA
        end
        sB  = find(ismember(cal,'B'));
        for h = 1:numel(sB)
            dB = datetime(year(DHSmalawi.interview(sET(j))),month((DHSmalawi.interview(sET(j)))) - sB(h) + 1,day(DHSmalawi.interview(sET(j))),'Format','dd/MM/yyyy');
            gA = cal(sB(h):min(sB(h) + 6,end));
            if dB >= dATe{1}(1) && dB  < dATe{1}(2)
               DHSmalawi.births(sET(j))      = DHSmalawi.births(sET(j)) + 1;
            end
            clear dB gA
        end
        clear sT sB
    end
    clear cal
    j/numel(sET)
end

p                     = rand(size(DHSmalawi,1),R + 1);
d                     = (p.*DHSmalawi.D_min + (1 - p).*(DHSmalawi.D_max - eps))/365.25;
B                     = eXAcTTime(DHSmalawi.Birth);
B                     = B*ones(1,R + 1);
D                     = B + d;
O                     = min(D,Tdhs(2));

p                     = rand(size(DHSmalawi,1),R + 1);
dob                   = p.*eXAcTTime(DHSmalawi.DOB) + (1 - p).*eXAcTTime(datetime(year(DHSmalawi.DOB),month(DHSmalawi.DOB) + 1,day(DHSmalawi.DOB)));
ageS                  = dATe{2}(2) - dob;
ageB                  = B - dob;

for i = 1:numel(dATaDHS)
    h                    = i + numel(dATa) + numel(dATaB);
    s                    = dATaDHS{i}{1};
    w                    = dATaDHS{i}{2};
    for k = 1:3
        if k == 1
            sEX           = ones(size(DHSmalawi,1),1);
        else
            sEX           = (DHSmalawi.sex == sex(k - 1));
        end
        sH                   = dATaDHS{i}{1} & sEX;
        
        for j = 1:numel(x{1}) - 1
            a             = max(min(B(sH,:) + x{1}(j),min(O(sH,:),Tdhs(2))),Tdhs(1));
            o             = max(min(B(sH,:) + x{1}(j + 1),min(O(sH,:),Tdhs(2))),Tdhs(1));
            exposure(j,:) = sum((o - a).*w(sH,:));
            events(j,:)   = sum((d(sH,:) >= x{1}(j) & d(sH,:) < x{1}(j + 1) & D(sH,:) >= Tdhs(1) & D(sH,:) < Tdhs(2)).*w(sH,:));
            clear a o
            clc;
            [21 j/(numel(x{1}) - 1)]
        end

        births               = sum(DHSmalawi.births(s).*w(s,:));
        stillbirths          = sum(DHSmalawi.stillbirths(s).*w(s,:));
        stillbirths          = stillbirths./(stillbirths + births);
        SRB                  = ((DHSmalawi.sex(s) == 1 & B(s,:) >= Tdhs(1) & B(s,:) < Tdhs(2))'*w(s,:))./((DHSmalawi.sex(s) == 2 & B(s,:) >= Tdhs(1) & B(s,:) < Tdhs(2))'*w(s,:));
        m                    = events./exposure;
        q                    = 1 - [ones(1,R + 1);exp(-cumsum(m.*n))];

        TaBle.q{k,h}         = q;    
        TaBle.s{k,h}         = 1 - (1 - q([1 2 5],:)).*(1 - stillbirths);
        TaBle.SNR{k,h}       = [TaBle.s{k,h}(1,:);TaBle.q{k,h}(2,:);TaBle.q{k,h}(5,:);TaBle.s{k,h}(1,:)./TaBle.q{k,h}(5,:)/1000;TaBle.q{k,h}(2,:)./TaBle.q{k,h}(5,:)/1000];
        TaBle.SRB{k,h}       = SRB;
        TaBle.tAU{k,h}       = dATe{1}';
        clear sEX events exposure births stillbirths m q SRB
    end
    
    A                    = [max(min(DHSmalawi.age(dATaDHS{i}{1})),15),min(max(DHSmalawi.age(dATaDHS{i}{1})),49)];
    TFR                  = 0;
    for j = A(1):A(2)
        a        = max(min(ageS(s,:) - 3,j + 1),j);
        o        = max(min(ageS(s,:),j + 1),j);
        exposure = sum((DHSmalawi.k(s) == 1).*(o - a).*w(s,:));
        events   = sum((dATe{2}(2) - B(s,:) <= 3 & dATe{2}(2) - B(s,:) > 0).*(ageB(s,:) >= j & ageB(s,:) < j + 1).*w(s,:));
        TFR      = TFR + events./exposure;
        clc;
        [22 (j - A(1) + 1)/(A(2) - A(1) + 1)]
    end
    
    TaBle.TFR{1,h}       = TFR;
    TaBle.sRB{1,h}       = sum((DHSmalawi.sex(s) == 1).*w(s,:))./sum((DHSmalawi.sex(s) == 2).*w(s,:));
    TaBle.parity{1,h}    = sum(~isnan(DHSmalawi.bidx(s)).*w(s,:))./sum((DHSmalawi.k(s) == 1).*w(s,:));
    TaBle.childless{1,h} = sum((DHSmalawi.mother(s) ~= 1 & DHSmalawi.k(s) == 1).*w(s,:))./sum((DHSmalawi.k(s) == 1).*w(s,:))*100;
    clear h j TFR w A
end

for i = 1:numel(dATaDHS)
    h              = i + numel(dATa) + numel(dATaB);
    TaBle.sEX{1,h} = log(1 - TaBle.q{3,h})./log(1 - TaBle.q{2,h});
end


load(char(pATh + "Results/RaMMPS.mat"),'MICSmalawi');
mAx                   = median(MICSmalawi.interview);
mIn                   = datetime([year(mAx) - 5,month(mAx),day(mAx)],'Format','dd/MM/yyyy');
dATe{1}               = [mIn mAx];
dATe{2}               = eXAcTTime(dATe{1});
Tmics                 = dATe{2};
models                = {'Women 15-49','Women 18-49','Women 18-49, mobile owners','Women 18-49, mobile owners'};
dATaMICS{1}           = {MICSmalawi.age >= 15,MICSmalawi.W};
dATaMICS{2}           = {MICSmalawi.age >= 18,MICSmalawi.W};
dATaMICS{3}           = {MICSmalawi.age >= 18 & MICSmalawi.mobile == 1,MICSmalawi.WR};
dATaMICS{4}           = {dATaMICS{3}{1},MICSmalawi.W};

for i = 1:numel(dATaMICS)
    data      = "MICS 6, " + string(sprintf('%0.1f',dATe{2}(1))) + "-" + string(sprintf('%0.1f',dATe{2}(2)));
    h         = i + numel(dATa) + numel(dATaB) + numel(dATaDHS);
    pOPs{1,h} = {data;models{i}};
    pOPs{2,h} = {data + " (female)";models{i}};
    pOPs{3,h} = {data + " (male)";models{i}};
end
clear mIn mAx models data

rng(0);
s                     = MICSmalawi(MICSmalawi.k == 1, {'cluster','K','woman','Women','iNDeX'});
w                     = rand(size(s,1),R);
w                     = [(1:size(s,1))',ceil(s.Women.*w) + s.iNDeX];
Wmics                 = NaN(size(MICSmalawi,1),R + 1);
for r = 1:R + 1
    S          = tabulate([w(:,1);w(:,r)]);
    Wmics(:,r) = repelem(S(:,2) - 1,s.K);
    clc;
    r/(R + 1)
end
clear r s S w

for i = 1:numel(dATaMICS)
    dATaMICS{i}{2} = Wmics.*dATaMICS{i}{2};
end
clear Wdhs
load(char(pATh + "Results/ReSAmPLiNG.mat"),'dATa','dATaDHS','Ts'); 
save(char(pATh + "Results/ReSAmPLiNG.mat"),'dATa','dATaDHS','dATaMICS','Ts','-v7.3','-nocompression');


p                     = rand(size(MICSmalawi,1),R + 1);
d                     = (p.*MICSmalawi.D_min + (1 - p).*(MICSmalawi.D_max - eps))/365.25;
p                     = rand(size(MICSmalawi,1),R + 1);
B                     = datetime(datenum(MICSmalawi.B_min) + datenum(MICSmalawi.B_max - MICSmalawi.B_min).*p,'ConvertFrom','datenum');
B                     = eXAcTTime(B);
D                     = B + d;
O                     = min(D,Tmics(2));

p                     = rand(size(MICSmalawi,1),R + 1);
dob                   = p.*eXAcTTime(MICSmalawi.DOB_min) + (1 - p).*eXAcTTime(MICSmalawi.DOB_max);
ageS                  = dATe{2}(2) - dob;
ageB                  = B - dob;

for i = 1:numel(dATaMICS)
    h                    = i + numel(dATa) + numel(dATaB) + numel(dATaDHS);
    s                    = dATaMICS{i}{1};
    w                    = dATaMICS{i}{2};
    for k = 1:3
        if k == 1
            sEX           = ones(size(MICSmalawi,1),1);
        else
            sEX           = (MICSmalawi.sex == sex(k - 1));
        end
        sH                   = dATaMICS{i}{1} & sEX;
        
        for j = 1:numel(x{1}) - 1
            a             = max(min(B(sH,:) + x{1}(j),min(O(sH,:),Tmics(2))),Tmics(1));
            o             = max(min(B(sH,:) + x{1}(j + 1),min(O(sH,:),Tmics(2))),Tmics(1));
            exposure(j,:) = sum((o - a).*w(sH,:));
            events(j,:)   = sum((d(sH,:) >= x{1}(j) & d(sH,:) < x{1}(j + 1) & D(sH,:) >= Tmics(1) & D(sH,:) < Tmics(2)).*w(sH,:));
            clear a o
            clc;
            [31 j/(numel(x{1}) - 1)]
        end

        SRB                  = ((MICSmalawi.sex(s) == 1 & B(s,:) >= Tmics(1) & B(s,:) < Tmics(2))'*w(s,:))./((MICSmalawi.sex(s) == 2 & B(s,:) >= Tmics(1) & B(s,:) < Tmics(2))'*w(s,:));
        m                    = events./exposure;
        q                    = 1 - [ones(1,R + 1);exp(-cumsum(m.*n))];

        TaBle.q{k,h}         = q;    
        TaBle.SRB{k,h}       = SRB;
        TaBle.tAU{k,h}       = dATe{1}';
        clear sEX events exposure m q SRB
    end
    
    A                    = [max(min(MICSmalawi.age(dATaMICS{i}{1})),15),min(max(MICSmalawi.age(dATaMICS{i}{1})),49)];
    TFR                  = 0;
    for j = A(1):A(2)
        a        = max(min(ageS(s,:) - 3,j + 1),j);
        o        = max(min(ageS(s,:),j + 1),j);
        exposure = sum((MICSmalawi.k(s) == 1).*(o - a).*w(s,:));
        events   = sum((dATe{2}(2) - B(s,:) <= 3 & dATe{2}(2) - B(s,:) > 0).*(ageB(s,:) >= j & ageB(s,:) < j + 1).*w(s,:));
        TFR      = TFR + events./exposure;
        clc;
        [32 (j - A(1) + 1)/(A(2) - A(1) + 1)]
    end
    
    TaBle.TFR{1,h}       = TFR;
    TaBle.sRB{1,h}       = sum((MICSmalawi.sex(s) == 1).*w(s,:))./sum((MICSmalawi.sex(s) == 2).*w(s,:));
    TaBle.parity{1,h}    = sum(~isnan(MICSmalawi.bidx(s)).*w(s,:))./sum((MICSmalawi.k(s) == 1).*w(s,:));
    TaBle.childless{1,h} = sum((MICSmalawi.mother(s) ~= 1 & MICSmalawi.k(s) == 1).*w(s,:))./sum((MICSmalawi.k(s) == 1).*w(s,:))*100;
    clear h j TFR w A
end

for i = 1:numel(dATaMICS)
    h              = i + numel(dATa) + numel(dATaB) + numel(dATaDHS);
    TaBle.sEX{1,h} = log(1 - TaBle.q{3,h})./log(1 - TaBle.q{2,h});
end



load(char(pATh + "Results/RaMMPS.mat"),'IGMEmalawi');
IGME             = mat2cell(table2array(IGMEmalawi),size(IGMEmalawi,1),[1 ones(1,(size(IGMEmalawi,2) - 1)/3)*3]);
for i = 2:numel(IGME)
    IGME{i} = IGME{i}/1000;
end

A{1}             = log(IGME{2}(1:22,:));
A{2}             = log(IGME{5}(1:22,:));
for i = 1:2
    for j = 1:size(A{i},1)
        delta = 0.005;
        mu    = (A{i}(j,2) + A{i}(j,3))/2;
        sigma = 1;
        E     = 10;
        while E^2 > eps
            dE    = (normcdf(A{i}(j,2),mu,sigma*exp(delta)) - normcdf(A{i}(j,2),mu,sigma*exp(-0.005)))/(2*delta);
            sigma = sigma*exp(-dE*E);
            E     = normcdf(A{i}(j,2),mu,sigma) - 0.05;
        end
        A{2 + i}(j,:) = exp(normrnd(mu,sigma,1,10000));
        clear E sigma mu delta
    end
end
A{4}             = A{4}./(1 + A{4});
A{3}             = A{4} + (1 - A{4}).*A{3};
IGME{7}          = prctile(A{3}',[50 5 95])';
IGME{5}          = IGME{5}(1:22,:)./(1 + IGME{5}(1:22,:));
IGME{6}          = NaN(size(IGME{5}));

h                = 1 + numel(dATa) + numel(dATaB) + numel(dATaDHS) + numel(dATaMICS);
k                = 4;
for i = 1:k
    q               = NaN(numel(x{1}),3);
    q(aGEs,:)       = [IGME{2}(k - i + 1,:);IGME{3}(k - i + 1,:);IGME{4}(k - i + 1,:)];
    TaBle.q{i,h}    = q;
    s               = [IGME{5}(k - i + 1,:);IGME{6}(k - i + 1,:);IGME{7}(k - i + 1,:)];
    TaBle.s{i,h}    = s;
    pOPs{i,h}       = {char("UN IGME, " + string(IGME{1}(k - i + 1)));''};
    TaBle.tAU{i,h}  = datetime(IGME{1}(k - i + 1) - .5,7,1);
    clear q s 
end
TaBle.q{i + 1,h}   = {IGME{2},IGME{3},IGME{4},datetime(IGME{1} - .5,7,1)};
TaBle.s{i + 1,h}   = {IGME{5},IGME{6},IGME{7},datetime(IGME{1} - .5,7,1)};
pOPs{i + 1,h}      = {'UN IGME';''};
TaBle.tAU{i + 1,h} = datetime(IGME{1} - .5,7,1);

[qS,qP]            = RaMMPS_Bayes(000,0); %Only run this line once, as RaMMPS_Bayes(5000,250)%
TaBleEs.q          = TaBle.q;
pOPsex             = pOPs(1,:);
for i = 1:numel(qS)
    TaBleEs.q(:,i) = [qP{i}';qS{i}{1};qS{i}{2}];
    TaBle.sEX{2,i} = qS{i}{2}./qS{i}{1};
    TaBle.sEX{3,i} = qS{i}{3};
    pOPsex{1,i}    = {pOPsex{1,i}{1};'Direct estimation'};
    pOPsex{2,i}    = {pOPsex{1,i}{1};'Predicted: Poisson Regression'};
    pOPsex{3,i}    = {pOPsex{1,i}{1};'Predicted: Poisson, via Equation 4'};
end


TaBle.childless{1} = NaN(1,R + 1);
TaBle.childless{2} = NaN(1,R + 1);
TaBle.parity{1}    = NaN(1,R + 1);
TaBle.parity{2}    = NaN(1,R + 1);

mAX                = 2;
tHtiles            = [50 2.5 97.5];
selection          = {[1 1 2],[1 3 2],[1 5 2],[1 7 1],[1 8 1],[1 9 2],[1 11 1],[1 12 1],[1 13 2]};
bASe               = {TaBle.childless,TaBle.parity,TaBle.TFR,TaBle.sRB};
for i = 1:numel(selection)
    temp       = pOPs{selection{i}(1),selection{i}(2)};
    temp{1}    = char(temp{1});
    if numel(strfind(temp{1},"(")) > 0
        temp{1}    = temp{1}(1:strfind(temp{1},"(") - 2);
    end
    if numel(strfind(temp{1},",")) > 0
        temp{1}    = temp{1}(1:strfind(temp{1},",") - 1);
    end
    label{i,1} = temp;
    clear temp
    for j = 1:numel(bASe)
        for k = 1:mAX
            h          = k + (j - 1)*mAX;
            table{i,h} = NaN(1,3);
        end

        for k = 1:selection{i}(3)
            h          = k + (j - 1)*mAX;
            temp       = bASe{j}{selection{i}(1),selection{i}(2) + k - 1};
            if size(temp,2) == R + 1 
                table{i,h} = prctile(temp(:,2:end)',tHtiles);
            elseif numel(temp) > 0
                table{i,h} = temp;
            end
        end
    end
end

sEt          = {'$\textrm{Childless Women \%}$','$\textrm{Average Parity}$','$\textrm{TFR}$','$\textrm{SRB}$'};
models       = {'$\textit{post-strat.}$','$\textit{selected}$'};
vARs         = {models models models models};
foRMaT       = {'%0.2f','%0.2f','%0.2f'};
nOTe         = {'$\textrm{Instrument}$','$\textit{Bootstrapping}$ $\mathrm{p50}$/$\mathit{[p2.5,p97.5]}$'};
lABs         = {{1} {2} {3} {4 5 6} {7 8 9}};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,label,cell2mat(table),0.140,0.055,[]);
exportgraphics(gcf,char(pATh + "Results/Table_3.png"),'Resolution',RESolUTioN);
s            = [1 2 3 4 5 6 7 8 9];
bOx          = table(s,:);
pOPsd        = label(s);
sets         = sEt;
vars         = vARs;
save(char(pATh + "Results/Tables/Table_3.mat"),'bOx','pOPsd','vars','sets');
clear table labels temp bOx pOPsd vars sets s



selection  = {[2 3 2],[3 3 2],[4 3 2],[5 3 2],[6 3 2],[1 15 1],[2 15 1],[3 15 1]};
mAX        = 2;
for i = 1:numel(selection)
    label{i,1} = pOPs{selection{i}(1),selection{i}(2)};
    for j = 1:numel(aGEs)
        for k = 1:mAX
            h          = k + (j - 1)*mAX;
            table{i,h} = NaN(1,3);
        end
        
        for k = 1:selection{i}(3)
            h          = k + (j - 1)*mAX;
            temp       = TaBleEs.q{selection{i}(1),selection{i}(2) + k - 1}(aGEs(j),:);
            if size(temp,2) > 3
                table{i,h} = prctile(temp(:,2:end)',tHtiles)*1000;
            else
                table{i,h} = temp*1000;
            end
        end
    end
end

sEt          = {'NMR - $\mathit{q}\mathrm{(28}\mathit{d}\mathrm{)}$','IMR - $\mathit{q}\mathrm{(12}\mathit{m}\mathrm{)}$','U5MR - $\mathit{q}\mathrm{(60}\mathit{m}\mathrm{)}$'};
vARs         = {models models models};
nOTePredict  = {'$\textrm{Instrument}$/$\textit{method}$','$\mathrm{RaMMPS\,predictions\,are\,Poisson-Hamiltonian\,MCMC\,estimation.}$ $\mathrm{p50}$/$\mathit{[p2.5,p97.5]}$ $\mathrm{posterior}$ $\mathrm{distribution.}$'};
lABs         = {{1 2 3 4 5} {6 7 8}};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTePredict,label,cell2mat(table),0.175,0.060,[]);
exportgraphics(gcf,char(pATh + "Results/Table_A2.png"),'Resolution',RESolUTioN);
s            = [1 2 3 4 5 6 7 8];
bOx          = table(s,:);
pOPsd        = label(s);
sets         = sEt;
vars         = vARs;
save(char(pATh + "Results/Tables/Table_A2.mat"),'bOx','pOPsd','vars','sets');
clear table labels temp bOx pOPsd vars sets s


selection  = {[1 1 2],[1 3 2],[2 1 2],[2 3 2],[3 1 2],[3 3 2],[1 8 1],[1 9 2],[1 12 1],[1 13 2]};
for i = 1:numel(selection)
    label{i,1} = pOPsex{selection{i}(1),selection{i}(2)};
    for j = 1:numel(aGEs)
        for k = 1:mAX
            h          = k + (j - 1)*mAX;
            table{i,h} = NaN(1,3);
        end
        
        for k = 1:selection{i}(3)
            h          = k + (j - 1)*mAX;
            temp       = TaBle.sEX{selection{i}(1),selection{i}(2) + k - 1}(aGEs(j),:);
            if size(temp,2) > 3
                table{i,h} = prctile(temp(:,2:end)',tHtiles);
            else
                table{i,h} = temp;
            end
        end
    end
end

nOTePredict2 = {'$\textrm{Instrument}$/$\textit{method}$','$\mathrm{if\,direct\,estimation,}$ $\mathit{Bootstrapping\,CI.}$ $\mathrm{If\,predicted,\,Poisson-Hamiltonian\,MCMC\,estimation.}$ $\mathrm{p50}$/$\mathit{[p2.5,p97.5]}$ $\mathrm{posterior\,distribution.}$'};
lABs         = {{1 2} {3 4} {5 6} {7 8} {9 10}};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTePredict2,label,cell2mat(table),0.175,0.060,[]);
exportgraphics(gcf,char(pATh + "Results/Table_5.png"),'Resolution',RESolUTioN);
s            = [1 2 3 4 5 6 7 8 9 10];
bOx          = table(s,:);
pOPsd        = label(s);
sets         = sEt;
vars         = vARs;
save(char(pATh + "Results/Tables/Table_5.mat"),'bOx','pOPsd','vars','sets');
clear table labels temp bOx pOPsd vars sets s



selection  = {[1 1 2],[1 3 2],[1 5 2],[2 5 2],[1 7 1],[1 8 1],[1 9 2],[1 11 1],[1 12 1],[1 13 2],[1 15 1],[2 15 1],[3 15 1],[4 15 1]};
for i = 1:numel(selection)
    label{i,1} = pOPs{selection{i}(1),selection{i}(2)};
    for j = 1:numel(aGEs)
        for k = 1:mAX
            h          = k + (j - 1)*mAX;
            table{i,h} = NaN(1,3);
        end
        
        for k = 1:selection{i}(3)
            h          = k + (j - 1)*mAX;
            temp       = TaBle.q{selection{i}(1),selection{i}(2) + k - 1}(aGEs(j),:);
            if size(temp,2) == R + 1
                table{i,h} = prctile(temp(:,2:end)',tHtiles)*1000;
            else
                table{i,h} = temp*1000;
            end
        end
    end
end

sEt          = {'NMR - $\mathit{q}\mathrm{(28}\mathit{d}\mathrm{)}$','IMR - $\mathit{q}\mathrm{(12}\mathit{m}\mathrm{)}$','U5MR - $\mathit{q}\mathrm{(60}\mathit{m}\mathrm{)}$'};;
vARs         = {models models models};
lABs         = {{1} {2} {3 4} {5 6 7} {8 9 10} {11 12 13 14}};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,label,cell2mat(table),0.175,0.060,[]);
exportgraphics(gcf,char(pATh + "Results/Table_4.png"),'Resolution',RESolUTioN);
s            = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
bOx          = table(s,:);
pOPsd        = label(s);
sets         = sEt;
vars         = vARs;
save(char(pATh + "Results/Tables/Table_4.mat"),'bOx','pOPsd','vars','sets');
clear table labels temp bOx pOPsd vars sets s

selection  = {[7 1 2],[8 1 2],[7 3 2],[8 3 2],[2 8 1],[3 8 1],[2 9 2],[3 9 2],[2 12 1],[3 12 1],[2 13 2],[3 13 2]};
for i = 1:numel(selection)
    label{i,1} = pOPs{selection{i}(1),selection{i}(2)};
    for j = 1:numel(aGEs)
        for k = 1:mAX
            h          = k + (j - 1)*mAX;
            table{i,h} = NaN(1,3);
        end
        
        for k = 1:selection{i}(3)
            h          = k + (j - 1)*mAX;
            temp       = TaBle.q{selection{i}(1),selection{i}(2) + k - 1}(aGEs(j),:);
            if size(temp,2) == R + 1
                table{i,h} = prctile(temp(:,2:end)',tHtiles)*1000;
            else
                table{i,h} = temp*1000;
            end
        end
    end
end

sEt          = {'NMR - $\mathit{q}\mathrm{(28}\mathit{d}\mathrm{)}$','IMR - $\mathit{q}\mathrm{(12}\mathit{m}\mathrm{)}$','U5MR - $\mathit{q}\mathrm{(60}\mathit{m}\mathrm{)}$'};
vARs         = {models models models};
lABs         = {{1 2} {3 4} {5 6} {7 8} {9 10} {11 12}};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,label,cell2mat(table),0.175,0.060,[]);
exportgraphics(gcf,char(pATh + "Results/Table_O2.png"),'Resolution',RESolUTioN);
s            = [1 2 3 4 5 6 7 8 9 10 11 12];
bOx          = table(s,:);
pOPsd        = label(s);
sets         = sEt;
vars         = vARs;
save(char(pATh + "Results/Tables/Table_O2.mat"),'bOx','pOPsd','vars','sets');
clear table labels temp bOx pOPsd vars sets s


selection  = {[2 3 2],[3 3 2],[4 3 2],[5 3 2],[6 3 2],[1 15 1],[2 15 1],[3 15 1],[4 15 1]};
mAX        = 2;
for i = 1:numel(selection)
    label{i,1} = pOPs{selection{i}(1),selection{i}(2)};
    for j = 1:numel(aGEs)
        for k = 1:mAX
            h          = k + (j - 1)*mAX;
            table{i,h} = NaN(1,3);
        end
        
        for k = 1:selection{i}(3)
            h          = k + (j - 1)*mAX;
            temp       = TaBle.q{selection{i}(1),selection{i}(2) + k - 1}(aGEs(j),:);
            if size(temp,2) == R + 1
                table{i,h} = prctile(temp(:,2:end)',tHtiles)*1000;
            else
                table{i,h} = temp*1000;
            end
        end
    end
end

lABs         = {{1 2 3 4 5} {6 7 8 9}};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,label,cell2mat(table),0.175,0.060,[]);
exportgraphics(gcf,char(pATh + "Results/Table_O1.png"),'Resolution',RESolUTioN);
s            = [1 2 3 4 5 6 7 8 9];
bOx          = table(s,:);
pOPsd        = label(s);
sets         = sEt;
vars         = vARs;
save(char(pATh + "Results/Tables/Table_O1.mat"),'bOx','pOPsd','vars','sets');
clear table labels temp bOx pOPsd vars sets s





P                        = cell(0);
selection                = {[1 1],[1 3],[2 5],[1 8],[1 12],[1 15]};
selectionT               = {[1 1],[1 3],[3 5],[1 8],[1 12],[5 15]};
LAB                      = {'Neonatal mortality rate $\mathbf{\it{q}}\mathbf{(28}\mathbf{\it{d}}\mathbf{)}$','Infant mortality rate $\mathbf{\it{q}}\mathbf{(12}\mathbf{\it{m}}\mathbf{)}$','Under-five mortality rate $\mathbf{\it{q}}\mathbf{(60}\mathbf{\it{m}}\mathbf{)}$'};
load(char(pATh + "Results/paleTTe.mat"),'paleTTe');
coloR                    = paleTTe([1 2 3 7 4 8]);

mPIX                     = 538756;
pix                      = 1/37.7952755906;
z                        = min(sqrt(mPIX/((10*3)*(10*2)/pix^2)),1);
fi                       = figure('Color',[1 1 1],'Position',z*10*[0 0 3 2]/pix,'Theme','light');
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
for i = 1:6
    nexttile(i)
    if i > 3
        ax{i}                       = gca;
        ax{i}.FontSize              = 10*z;
        ax{i}.XScale                = 'log';
        ax{i}.XAxis.TickLabelFormat = '%.1f';
        ax{i}.YAxis.TickLabelFormat = '%.2f';
        ax{i}.XTick                 = .1./(2.^(9:-1:-2))*1000;
        ax{i}.XAxis.MinorTickValues = 10:10:200;
        xlim([6.25 100])

        t                           = ax{i}.XAxis.TickValues;
        for j = 1:numel(t)
            T{j} = char("$\mathbf{" + string(sprintf(ax{i}.XAxis.TickLabelFormat,t(j))) + "}$");
        end
        ax{i}.XTickLabelRotation         = 0;
        ax{i}.XAxis.TickLabels           = T;
        ax{i}.XAxis.TickLabelInterpreter = 'latex';
        clear t T

        xlabel('$\textbf{Deaths per 1000 births (log scale)}$','Interpreter','latex','FontSize',11*z);
        if isequal(mod(i,3),1)
            ylabel('$\textbf{Kernel density}$','Interpreter','latex','FontSize',11*z);
        end
    else
        d                           = scatter(datetime('01-Jan-2021'),0);
        delete(d)
        ax{i}                       = gca;
        ax{i}.FontSize              = 10*z;
        ax{i}.YAxis.TickLabelFormat = '%.1f';
        ax{i}.YMinorGrid            = 'on';
        ax{i}.XMinorGrid            = 'on';
        ax{i}.XAxis.TickLabelFormat = 'yyyy';
        ax{i}.XTick                 = datetime(2010:2:2024,1,1);
        ax{i}.XAxis.MinorTickValues = datetime(2010:2024,1,1);
        ax{i}.YScale                = 'log';
        ax{i}.YAxis.TickLabelFormat = '%.1f';
        ax{i}.YTick                 = .1./(2.^(9:-1:-2))*1000;
        ax{i}.YAxis.MinorTickValues = 10:10:200;
        ylim([6.25 100])
        xlim([datetime(2010,1,1) datetime(2024,1,1)]);
        
        T                           = ax{i}.XAxis.TickLabels;
        for j = 1:numel(T)
            T{j} = char("$\mathbf{" + T{j} + "}$");
        end
        ax{i}.XTickLabelRotation         = 0;
        ax{i}.XAxis.TickLabels           = T;
        ax{i}.XAxis.TickLabelInterpreter = 'latex';
        clear T

        t                           = ax{i}.YAxis.TickValues;
        for j = 1:numel(t)
            T{j} = char("$\mathbf{" + string(sprintf(ax{i}.YAxis.TickLabelFormat,t(j))) + "}$");
        end
        ax{i}.YTickLabelRotation         = 0;
        ax{i}.YAxis.TickLabels           = T;
        ax{i}.YAxis.TickLabelInterpreter = 'latex';
        clear t T

        xlabel('$\textbf{Year}$','Interpreter','latex','FontSize',11*z);
        if isequal(mod(i,3),1)
            ylabel('$\textbf{Deaths per 1000 births (log scale)}$','Interpreter','latex','FontSize',11*z);
        end        
    end
    title(char("$\textbf{" + string(char(96 + i)) + ". " + string(LAB{mod(i - 1,3) + 1}) + "}$"),'Interpreter','latex','FontSize',12*z);
    grid on;
    box on;
    hold on;
end

for i = 1:6
    nexttile(i)
    if i > 3
        mf  = 4.0;
        ylim([0 mf])
        for j = 1:numel(selection)
            q{j} = TaBle.q{selection{j}(1),selection{j}(2)}(aGEs(i - 3),:);
            q{j} = recode(q{j},NaN,eps);
            if numel(q{j}) > 3
                [f{j},xi{j}] = ksdensity(log(max(q{j}(2:end),eps)));
                xi{j}        = exp(xi{j})*1000;
                pc{j}        = prctile(q{j}(2:end),tHtiles)*1000;
                [g{j},~]     = ksdensity(log(q{j}(2:end)),log(pc{j}/1000));
            else
                xi{j}        = kron(q{j}(2:3),ones(1,2))*1000;
                f{j}         = [0,mf,mf,0];
                pc{j}        = q{j}*1000;
                g{j}         = ones(1,numel(tHtiles))*mf;
            end
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j},'LineWidth',1.15*z);
        end
        for j = 1:numel(selection)
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.25*z,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.75*z);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.75*z);
        end
        t                                = ax{i}.YAxis.TickValues;
        for j = 1:numel(t)
            T{j} = char("$\mathbf{" + string(sprintf(ax{i}.XAxis.TickLabelFormat,t(j))) + "}$");
        end
        ax{i}.YTickLabelRotation         = 0;
        ax{i}.YAxis.TickLabels           = T;
        ax{i}.YAxis.TickLabelInterpreter = 'latex';
        clear t T

        if i == 5
            for j = 1:numel(selection)
                leGend{j} = char("$\textbf{" + pOPs{selection{j}(1),selection{j}(2)}{1} + "}$");
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',9*z,'Location','southoutside','NumColumns',3,'Box','off');
        end
    else
        for j = 1:numel(selectionT)
            tAU{j} = TaBle.tAU{selectionT{j}(1),selectionT{j}(2)};
            if numel(tAU{j}) > 2
                q{j} = TaBle.q{selectionT{j}(1),selectionT{j}(2)}{1 + mod(i + 2,3)};
            else
                q{j} = kron(ones(2,1),TaBle.q{selectionT{j}(1),selectionT{j}(2)}(aGEs(1 + mod(i + 2,3)),:));
            end
            if size(q{j},2) > 3
                q{j} = prctile(q{j}(:,2:end)',tHtiles)'*1000;
            else
                q{j} = q{j}*1000;
            end
            P{end + 1} = fill([tAU{j};flip(tAU{j})],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j},'LineWidth',1.15*z);
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j},q{j}(:,1),'color',coloR{j},'LineWidth',1.25*z,'LineStyle',':');
        end
        P{end + 1} = plot(datetime(ones(2,1)*[2020 1 1]),[3.125 200],'color',[0.1 0.4 0.5],'LineStyle','-.','LineWidth',1.00*z);
        if i == 2
            for j = 1:numel(selectionT)
                leGend{j} = char("$\textbf{" + pOPs{selectionT{j}(1),selectionT{j}(2)}{1} + "}$");
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',9*z,'Location','southoutside','NumColumns',3,'Box','off');
        end
    end
end
exportgraphics(gcf,char(pATh + "Results/Figure_2.png"),'Resolution',RESolUTioN);
for i = 1:numel(P)
    delete(P{i})
end
clear leGend P

coloR                    = paleTTe([1 2 3 7 4 8]);
P                        = cell(0);
selection                = {[7 3],[8 3],[2 8],[3 8],[2 12],[3 12]};
selectionT               = {[7 3],[8 3],[2 8],[3 8],[2 12],[3 12]};
for i = 1:6
    nexttile(i)
    if i > 3
        mf  = 4.0;
        ylim([0 mf])
        for j = 1:numel(selection)
            q{j} = TaBle.q{selection{j}(1),selection{j}(2)}(aGEs(i - 3),:);
            q{j} = recode(q{j},NaN,eps);
            if numel(q{j}) > 3
                [f{j},xi{j}] = ksdensity(log(max(q{j}(2:end),eps)));
                xi{j}        = exp(xi{j})*1000;
                pc{j}        = prctile(q{j}(2:end),tHtiles)*1000;
                [g{j},~]     = ksdensity(log(q{j}(2:end)),log(pc{j}/1000));
            else
                xi{j}        = kron(q{j}(2:3),ones(1,2))*1000;
                f{j}         = [0,mf,mf,0];
                pc{j}        = q{j}*1000;
                g{j}         = ones(1,numel(tHtiles))*mf;
            end
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j},'LineWidth',1.15*z);
        end
        for j = 1:numel(selection)
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.0*z,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.50*z);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.50*z);
        end
        t                                = ax{i}.YAxis.TickValues;
        for j = 1:numel(t)
            T{j} = char("$\mathbf{" + string(sprintf(ax{i}.XAxis.TickLabelFormat,t(j))) + "}$");
        end
        ax{i}.YTickLabelRotation         = 0;
        ax{i}.YAxis.TickLabels           = T;
        ax{i}.YAxis.TickLabelInterpreter = 'latex';
        clear t T

        if i == 5
            for j = 1:numel(selection)
                leGend{j} = char("$\textbf{" + pOPs{selection{j}(1),selection{j}(2)}{1} + "}$");
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',9*z,'Location','southoutside','NumColumns',4,'Box','off');
        end
    else
        for j = 1:numel(selectionT)
            tAU{j} = TaBle.tAU{selectionT{j}(1),selectionT{j}(2)};
            if numel(tAU{j}) > 2
                q{j} = TaBle.q{selectionT{j}(1),selectionT{j}(2)}{1 + mod(i + 2,3)};
            else
                q{j} = kron(ones(2,1),TaBle.q{selectionT{j}(1),selectionT{j}(2)}(aGEs(1 + mod(i + 2,3)),:));
            end
            if size(q{j},2) > 3
                q{j} = prctile(q{j}(:,2:end)',tHtiles)'*1000;
            else
                q{j} = q{j}*1000;
            end
            P{end + 1} = fill([tAU{j};flip(tAU{j})],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j},'LineWidth',1.15*z);
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j},q{j}(:,1),'color',coloR{j},'LineWidth',1.0,'LineStyle',':');
        end
        P{end + 1} = plot(datetime(ones(2,1)*[2020 1 1]),[3.125 200],'color',[0.1 0.4 0.5],'LineStyle','-.','LineWidth',0.75);
        if i == 2
            for j = 1:numel(selectionT)
                leGend{j} = char("$\textbf{" + pOPs{selectionT{j}(1),selectionT{j}(2)}{1} + "}$");
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',9*z,'Location','southoutside','NumColumns',3,'Box','off');
        end
    end
end
exportgraphics(gcf,char(pATh + "Results/Figure_3.png"),'Resolution',RESolUTioN);
for i = 1:numel(P)
    delete(P{i})
end
clear leGend P

coloR                    = paleTTe([1 2 3 5 6 7 4 8]);
P                        = cell(0);
selection                = {[2 3],[3 3],[4 3],[5 3],[6 3],[1 8],[1 12],[1 15]};
selectionT               = {[2 3],[3 3],[4 3],[5 3],[6 3],[1 8],[1 12],[5 15]};
for i = 1:6
    nexttile(i)
    if i > 3
        mf  = 2.5;
        ylim([0 mf])
        for j = 1:numel(selection)
            q{j} = TaBle.q{selection{j}(1),selection{j}(2)}(aGEs(i - 3),:);
            q{j} = recode(q{j},NaN,eps);
            if numel(q{j}) > 3
                [f{j},xi{j}] = ksdensity(log(max(q{j}(2:end),eps)));
                xi{j}        = exp(xi{j})*1000;
                pc{j}        = prctile(q{j}(2:end),tHtiles)*1000;
                [g{j},~]     = ksdensity(log(q{j}(2:end)),log(pc{j}/1000));
            else
                xi{j}        = kron(q{j}(2:3),ones(1,2))*1000;
                f{j}         = [0,mf,mf,0];
                pc{j}        = q{j}*1000;
                g{j}         = ones(1,numel(tHtiles))*mf;
            end
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j},'LineWidth',1.15*z);
        end
        for j = 1:numel(selection)
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.0*z,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.50*z);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.50*z);
        end
        t                                = ax{i}.YAxis.TickValues;
        for j = 1:numel(t)
            T{j} = char("$\mathbf{" + string(sprintf(ax{i}.XAxis.TickLabelFormat,t(j))) + "}$");
        end
        ax{i}.YTickLabelRotation         = 0;
        ax{i}.YAxis.TickLabels           = T;
        ax{i}.YAxis.TickLabelInterpreter = 'latex';
        clear t T

        if i == 5
            for j = 1:numel(selection)
                leGend{j} = char("$\textbf{" + pOPs{selection{j}(1),selection{j}(2)}{1} + "}$");
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',9*z,'Location','southoutside','NumColumns',4,'Box','off');
        end
    else
        for j = 1:numel(selectionT)
            tAU{j} = TaBle.tAU{selectionT{j}(1),selectionT{j}(2)};
            if numel(tAU{j}) > 2
                q{j} = TaBle.q{selectionT{j}(1),selectionT{j}(2)}{1 + mod(i + 2,3)};
            else
                q{j} = kron(ones(2,1),TaBle.q{selectionT{j}(1),selectionT{j}(2)}(aGEs(1 + mod(i + 2,3)),:));
            end
            if size(q{j},2) > 3
                q{j} = prctile(q{j}(:,2:end)',tHtiles)'*1000;
            else
                q{j} = q{j}*1000;
            end
            P{end + 1} = fill([tAU{j};flip(tAU{j})],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j},'LineWidth',1.15*z);
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j},q{j}(:,1),'color',coloR{j},'LineWidth',1.0,'LineStyle',':');
        end
        P{end + 1} = plot(datetime(ones(2,1)*[2020 1 1]),[3.125 200],'color',[0.1 0.4 0.5],'LineStyle','-.','LineWidth',0.75);
        if i == 2
            for j = 1:numel(selectionT)
                leGend{j} = char("$\textbf{" + pOPs{selectionT{j}(1),selectionT{j}(2)}{1} + "}$");
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',9*z,'Location','southoutside','NumColumns',4,'Box','off');
        end
    end
end
exportgraphics(gcf,char(pATh + "Results/Figure_4.png"),'Resolution',RESolUTioN);
for i = 1:numel(P)
    delete(P{i})
end
clear leGend P

coloR                    = paleTTe([1 2 3 7 4 8]);
P                        = cell(0);
selection                = {[7 3],[8 3],[2 8],[3 8],[2 12],[3 12]};
selectionT               = {[7 3],[8 3],[2 8],[3 8],[2 12],[3 12]};
for i = 1:6
    nexttile(i)
    if i > 3
        mf  = 4.0;
        ylim([0 mf])
        for j = 1:numel(selection)
            q{j} = TaBleEs.q{selection{j}(1),selection{j}(2)}(aGEs(i - 3),:);
            q{j} = recode(q{j},NaN,eps);
            if numel(q{j}) > 3
                [f{j},xi{j}] = ksdensity(log(max(q{j}(2:end),eps)));
                xi{j}        = exp(xi{j})*1000;
                pc{j}        = prctile(q{j}(2:end),tHtiles)*1000;
                [g{j},~]     = ksdensity(log(q{j}(2:end)),log(pc{j}/1000));
            else
                xi{j}        = kron(q{j}(2:3),ones(1,2))*1000;
                f{j}         = [0,mf,mf,0];
                pc{j}        = q{j}*1000;
                g{j}         = ones(1,numel(tHtiles))*mf;
            end
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j},'LineWidth',1.15*z);
        end
        for j = 1:numel(selection)
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.0*z,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.50*z);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.50*z);
        end
        t                                = ax{i}.YAxis.TickValues;
        for j = 1:numel(t)
            T{j} = char("$\mathbf{" + string(sprintf(ax{i}.XAxis.TickLabelFormat,t(j))) + "}$");
        end
        ax{i}.YTickLabelRotation         = 0;
        ax{i}.YAxis.TickLabels           = T;
        ax{i}.YAxis.TickLabelInterpreter = 'latex';
        clear t T

        if i == 5
            for j = 1:numel(selection)
                leGend{j} = char("$\textbf{" + pOPs{selection{j}(1),selection{j}(2)}{1} + "}$");
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',9*z,'Location','southoutside','NumColumns',3,'Box','off');
        end
    else
        for j = 1:numel(selectionT)
            tAU{j} = TaBle.tAU{selectionT{j}(1),selectionT{j}(2)};
            if numel(tAU{j}) > 2
                q{j} = TaBleEs.q{selectionT{j}(1),selectionT{j}(2)}{1 + mod(i + 2,3)};
            else
                q{j} = kron(ones(2,1),TaBleEs.q{selectionT{j}(1),selectionT{j}(2)}(aGEs(1 + mod(i + 2,3)),:));
            end
            if size(q{j},2) > 3
                q{j} = prctile(q{j}(:,2:end)',tHtiles)'*1000;
            else
                q{j} = q{j}*1000;
            end
            P{end + 1} = fill([tAU{j};flip(tAU{j})],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j},'LineWidth',1.15*z);
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j},q{j}(:,1),'color',coloR{j},'LineWidth',1.0,'LineStyle',':');
        end
        P{end + 1} = plot(datetime(ones(2,1)*[2020 1 1]),[3.125 200],'color',[0.1 0.4 0.5],'LineStyle','-.','LineWidth',0.75);
        if i == 2
            for j = 1:numel(selectionT)
                leGend{j} = char("$\textbf{" + pOPs{selectionT{j}(1),selectionT{j}(2)}{1} + "}$");
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',9*z,'Location','southoutside','NumColumns',3,'Box','off');
        end
    end
end
exportgraphics(gcf,char(pATh + "Results/Figure_3Reg.png"),'Resolution',RESolUTioN);
for i = 1:numel(P)
    delete(P{i})
end
clear leGend P

coloR                    = paleTTe([1 2 3 5 6 7 4 8]);
P                        = cell(0);
selection                = {[2 3],[3 3],[4 3],[5 3],[6 3],[1 8],[1 12],[1 15]};
selectionT               = {[2 3],[3 3],[4 3],[5 3],[6 3],[1 8],[1 12],[5 15]};
for i = 1:6
    nexttile(i)
    if i > 3
        mf  = 2.5;
        ylim([0 mf])
        for j = 1:numel(selection)
            q{j} = TaBleEs.q{selection{j}(1),selection{j}(2)}(aGEs(i - 3),:);
            q{j} = recode(q{j},NaN,eps);
            if numel(q{j}) > 3
                [f{j},xi{j}] = ksdensity(log(max(q{j}(2:end),eps)));
                xi{j}        = exp(xi{j})*1000;
                pc{j}        = prctile(q{j}(2:end),tHtiles)*1000;
                [g{j},~]     = ksdensity(log(q{j}(2:end)),log(pc{j}/1000));
            else
                xi{j}        = kron(q{j}(2:3),ones(1,2))*1000;
                f{j}         = [0,mf,mf,0];
                pc{j}        = q{j}*1000;
                g{j}         = ones(1,numel(tHtiles))*mf;
            end
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j},'LineWidth',1.15*z);
        end
        for j = 1:numel(selection)
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.0*z,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.50*z);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.50*z);
        end
        t                                = ax{i}.YAxis.TickValues;
        for j = 1:numel(t)
            T{j} = char("$\mathbf{" + string(sprintf(ax{i}.XAxis.TickLabelFormat,t(j))) + "}$");
        end
        ax{i}.YTickLabelRotation         = 0;
        ax{i}.YAxis.TickLabels           = T;
        ax{i}.YAxis.TickLabelInterpreter = 'latex';
        clear t T

        if i == 5
            for j = 1:numel(selection)
                leGend{j} = char("$\textbf{" + pOPs{selection{j}(1),selection{j}(2)}{1} + "}$");
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',9*z,'Location','southoutside','NumColumns',4,'Box','off');
        end
    else
        for j = 1:numel(selectionT)
            tAU{j} = TaBle.tAU{selectionT{j}(1),selectionT{j}(2)};
            if numel(tAU{j}) > 2
                q{j} = TaBleEs.q{selectionT{j}(1),selectionT{j}(2)}{1 + mod(i + 2,3)};
            else
                q{j} = kron(ones(2,1),TaBleEs.q{selectionT{j}(1),selectionT{j}(2)}(aGEs(1 + mod(i + 2,3)),:));
            end
            if size(q{j},2) > 3
                q{j} = prctile(q{j}(:,2:end)',tHtiles)'*1000;
            else
                q{j} = q{j}*1000;
            end
            P{end + 1} = fill([tAU{j};flip(tAU{j})],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j},'LineWidth',1.15*z);
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j},q{j}(:,1),'color',coloR{j},'LineWidth',1.0,'LineStyle',':');
        end
        P{end + 1} = plot(datetime(ones(2,1)*[2020 1 1]),[3.125 200],'color',[0.1 0.4 0.5],'LineStyle','-.','LineWidth',0.75);
        if i == 2
            for j = 1:numel(selectionT)
                leGend{j} = char("$\textbf{" + pOPs{selectionT{j}(1),selectionT{j}(2)}{1} + "}$");
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontSize',9*z,'Location','southoutside','NumColumns',4,'Box','off');
        end
    end
end
exportgraphics(gcf,char(pATh + "Results/Figure_4Reg.png"),'Resolution',RESolUTioN);
for i = 1:numel(P)
    delete(P{i})
end
clear leGend P


LAB                      = {'RaMMPS','DHS VII','MICS 6'};
z                        = min(sqrt(mPIX/((10*3)*(10*1)/pix^2)),1);
fi                       = figure('Color',[1 1 1],'Position',z*10*[0 0 3 1]/pix,'Theme','light');
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
wEIGhTS                  = {RaMMPS.WR(RaMMPS.k == 1),DHSmalawi.WR(DHSmalawi.k == 1 & ~isnan(DHSmalawi.WR)),MICSmalawi.WR(MICSmalawi.k == 1 & ~isnan(MICSmalawi.WR))}; 
for i = 1:3
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontSize              = 10*z;
    ax{i}.XAxis.TickLabelFormat = '%.1f';
    ax{i}.YAxis.TickLabelFormat = '%.2f';
    xlabel('$\textbf{Post-stratified weights}$','Interpreter','latex','FontSize',11*z);
    ylabel('$\textbf{Probability density function}$','Interpreter','latex','FontSize',11*z);
    xlim([0 15])    
    grid on;
    box on;
    hold on;
    
    y      = wEIGhTS{i};
    x      = min(y):.25:max(y);
    histogram(y,x,'Normalization','pdf','FaceColor',coloR{i},'EdgeColor',[0 0 0],'FaceAlpha',0.25,'EdgeAlpha',0.45);
    title(char("$\textbf{" + string(char(96 + i)) + ". " + string(LAB{i}) + ", N = " + string(numel(y)) + "}$"),'Interpreter','latex','FontSize',12*z);
    f      = ksdensity(y,x,'Function','pdf','Bandwidth',.1);
    plot(x,f,'LineWidth',0.75,'Color',[coloR{i} 0.75]);
    clear x y f

    T                           = ax{i}.XAxis.TickLabels;
    for j = 1:numel(T)
        T{j} = char("$\mathbf{" + T{j} + "}$");
    end
    ax{i}.XAxis.TickLabels           = T;
    ax{i}.XAxis.TickLabelInterpreter = 'latex';
    clear T

    T                           = ax{i}.YAxis.TickLabels;
    for j = 1:numel(T)
        T{j} = char("$\mathbf{" + T{j} + "}$");
    end
    ax{i}.YAxis.TickLabels           = T;
    ax{i}.YAxis.TickLabelInterpreter = 'latex';
    clear T
end
exportgraphics(gcf,char(pATh + "Results/Figure_A2.png"),'Resolution',RESolUTioN);



tAU                      = (2006:.25:2018)';
LAB                      = {'Infant mortality rate $\mathbf{\it{q}}\mathbf{(12}\mathbf{\it{m}}\mathbf{)}$','Under-five mortality rate $\mathbf{\it{q}}\mathbf{(60}\mathbf{\it{m}}\mathbf{)}$'};
coloR                    = {[0.00 0.00 0.75],[0.95 0.00 0.95],[0.85 0.35 0.01],[0.45 0.65 0.20],[0.65 0.10 0.20],[0.00 0.55 0.65],[0.05 0.05 0.05]};
z                        = min(sqrt(mPIX/((10*2)*(10*2)/pix^2)),1);
fi                       = figure('Color',[1 1 1],'Position',z*10*[0 0 2 2]/pix,'Theme','light');
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
for i = 1:4
    nexttile(i)
    delete(d)
    ax{i}                       = gca;
    ax{i}.FontSize              = 10*z;
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    ax{i}.YMinorGrid            = 'on';
    ax{i}.XMinorGrid            = 'on';
    ax{i}.XAxis.TickLabelFormat = '%.0f';
    ax{i}.XTick                 = 2006:2:2020;
    ax{i}.XAxis.MinorTickValues = 2006:2020;
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    xlim([2006 2020]);

    if i > 2
        xlabel('$\textbf{Year}$','Interpreter','latex','FontSize',11*z);
    end

    if 1 + mod(i + 1,2) == 1
        ax{i}.YTick                 = 0:20:200;
        ax{i}.YAxis.MinorTickValues = 0:20:200;
        ylim([0 120])
        ylabel('$\textbf{Deaths per 1000 births}$','Interpreter','latex','FontSize',11*z);
    else
        ax{i}.YScale                = 'log';
        ax{i}.YTick                 = .1./(2.^(9:-1:-2))*1000;
        ax{i}.YAxis.MinorTickValues = 10:10:100;
        ylim([25 150])
        ylabel('$\textbf{Deaths per 1000 births (log scale)}$','Interpreter','latex','FontSize',11*z);
    end

    t                           = ax{i}.XAxis.TickValues;
    for j = 1:numel(t)
        T{j} = char("$\mathbf{" + string(sprintf(ax{i}.XAxis.TickLabelFormat,t(j))) + "}$");
    end
    ax{i}.XTickLabelRotation         = 0;
    ax{i}.XAxis.TickLabels           = T;
    ax{i}.XAxis.TickLabelInterpreter = 'latex';
    clear T

    t                           = ax{i}.YAxis.TickValues;
    for j = 1:numel(t)
        T{j} = char("$\mathbf{" + string(sprintf(ax{i}.YAxis.TickLabelFormat,t(j))) + "}$");
    end
    ax{i}.YTickLabelRotation         = 0;
    ax{i}.YAxis.TickLabels           = T;
    ax{i}.YAxis.TickLabelInterpreter = 'latex';
    clear t T

    title(char("$\textbf{" + string(char(96 + i)) + ". " + string(LAB{ceil(i/2)}) + "}$"),'Interpreter','latex','FontSize',12*z);
    grid on;
    box on;
    hold on;
end

for i = 1:4
    nexttile(i)
    h = ceil(i/2);
    if 1 + mod(i + 1,2) == 1
        for j = 1:numel(CD)
            plot(prctile(BraSs_r{j,1}{3}(:,2:end)',50)',prctile(BraSs_r{j,1}{h}(:,2:end)',50)'*1000,'color',coloR{j},'LineWidth',1.00*z);
        end
        for j = 1:numel(CD)
            plot(BraSs_r{j,1}{3}(:,2:50),BraSs_r{j,1}{h}(:,2:50)*1000,'color',[coloR{j} 0.1],'LineWidth',0.75*z);
        end
    else
        for j = 1:numel(CD)
            temp{j} = max(LinInterPol(BraSs_r{j,1}{3},BraSs_r{j,1}{h},tAU),0);
            temp{j} = prctile(temp{j}(:,2:end)',tHtiles)'*1000;
            fill([tAU;flip(tAU)],max([temp{j}(:,2);flip(temp{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j});
        end
        for j = 1:numel(CD)
            plot(tAU,temp{j}(:,1),'color',coloR{j},'lineWidth',1.00*z,'LineStyle',':');
        end
    end
    if h == 1
        legend(CD,'Interpreter','latex','FontSize',9*z,'Location','southoutside','NumColumns',3,'Box','off');
    end
end
exportgraphics(gcf,char(pATh + "Results/Figure_A3.png"),'Resolution',RESolUTioN);
clear leGend tAU




clear
pATh                     = "/Users/lshjr3/Documents/RaMMPS/under-five mortality/";
load(char(pATh + "Results/RaMMPS.mat"),'RaMMPS','RESolUTioN');
sET                      = {RaMMPS.TPH_duration(RaMMPS.k == 1),RaMMPS.FPH_duration(RaMMPS.k == 1)};
LAB                      = {'Truncated Pregnancy History','Full Pregnancy History','Household Questionnaire','All RaMMPS'};
coloR                    = {[0.00 0.00 0.75],[0.95 0.00 0.95],[0.85 0.35 0.01],[0.45 0.65 0.20],[0.65 0.10 0.20],[0.00 0.55 0.65],[0.05 0.05 0.05]};

mPIX                     = 538756;
pix                      = 1/37.7952755906;
z                        = min(sqrt(mPIX/((10*2)*(10*1)/pix^2)),1);
fi                       = figure('Color',[1 1 1],'Position',z*10*[0 0 2 1]/pix,'Theme','light');
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
for i = 1:2
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontSize              = 10*z;
    ax{i}.XAxis.TickLabelFormat = '%.1f';
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    ax{i}.XTick                 = 0:2.50:40;
    ax{i}.XAxis.MinorTickValues = 0:1.25:40;
    ax{i}.YTick                 = 0:0.10:1;
    ax{i}.YAxis.MinorTickValues = 0:0.05:1;

    xlabel('$\textbf{Duration (in minutes)}$','Interpreter','latex','FontSize',11*z);
    ylabel('$\textbf{Probability density function}$','Interpreter','latex','FontSize',11*z);
    grid on;
    grid minor;
    box on;
    hold on;
    
    xi     = 0:0.25:40;
    y      = kron(sET{i}(~isnan(sET{i}))/60,ones(1,1));
    histogram(y,xi,'Normalization','pdf','FaceColor',coloR{i},'EdgeColor',[0 0 0],'FaceAlpha',0.1,'EdgeAlpha',0.25);
    title("$\textbf{" + string(char(96 + i)) + ". " + char(LAB{i} + ", N = " + string(numel(y))) + "}$",'Interpreter','latex','FontSize',12*z);
    f      = ksdensity(y,xi,'Function','pdf','Bandwidth',.1);
    plot(xi,f,'LineWidth',0.50,'Color',[coloR{i} 0.5]);
    W      = mean(y);
    round(W,2)
    Z(i,:) = prctile(y,[50 25 75]);
    [f,W]  = ksdensity(y,W,'Function','pdf','Bandwidth',.1);
    plot(W(1)*ones(1,2),[0 f(1)],'LineWidth',1.0,'Color',coloR{i});
    xlim([-.5 15])
    ylim([0 .80]);

    t                           = ax{i}.XAxis.TickValues;
    for j = 1:numel(t)
        T{j} = char("$\mathbf{" + string(sprintf(ax{i}.XAxis.TickLabelFormat,t(j))) + "}$");
    end
    ax{i}.XAxis.TickLabels           = T;
    ax{i}.XAxis.TickLabelInterpreter = 'latex';
    clear t T

    t                           = ax{i}.YAxis.TickValues;
    for j = 1:numel(t)
        T{j} = char("$\mathbf{" + string(sprintf(ax{i}.YAxis.TickLabelFormat,t(j))) + "}$");
    end
    ax{i}.YAxis.TickLabels           = T;
    ax{i}.YAxis.TickLabelInterpreter = 'latex';
    clear t T
end
exportgraphics(gcf,char(pATh + "Results/Figure_A1.png"),'Resolution',RESolUTioN);
Z(:,4)       = Z(:,3) - Z(:,2);


clear
pATh         = "/Users/lshjr3/Documents/RaMMPS/under-five mortality/";
load(char(pATh + "Results/RaMMPSf.mat"),'RaMMPS','RESolUTioN');
RaMMPS       = RaMMPS(RaMMPS.caseid ~= "0A0" & RaMMPS.k == 1,:);
sample       = {'SBH','TPH','FPH'};

date         = [min(RaMMPS.interview) - 2;datetime(2022,5,25);datetime(2022,9,13);datetime(2023,3,20);max(RaMMPS.interview)];
for i = 1:numel(date) - 1
    t      = char("$\textrm{" + string(i) + ". from " + string(datetime(date(i) + 1,'Format','MMM d, yyyy')) + " to " + string(datetime(date(i + 1),'Format','MMM d, yyyy')) + "}$");
    pOP{i} = {t;''};
end
pOP{end + 1} = {'$\textrm{Total}$';''};

for i = 1:numel(pOP) - 1
    for j = 1:numel(sample)
        N                = sum(RaMMPS.sample == sample{j} & RaMMPS.group == i);
        Nk               = sum(RaMMPS.sample == sample{j} & RaMMPS.group == i & isundefined(RaMMPS.missing));
        taBLe(i,j*2 - 1) = N;
        taBLe(i,j*2)     = N - Nk;
    end
end
taBLe        = [taBLe;sum(taBLe,1)];
taBLe        = [taBLe,sum(taBLe(:,1:2:end),2),sum(taBLe(:,2:2:end),2)];
taBLe        = kron(taBLe,[1 NaN(1,2)]);
sEt          = {{'$\textrm{SBH}$'},{'$\textrm{TPH}$'},{'$\textrm{FPH/SBH}$'},{'$\textrm{Total}$'}};
numbers      = {'$\textit{interviews}$','$\textit{excluded}$'};

lABs         = {{1 2 3 4} {5}};
vARs         = {numbers numbers numbers numbers};
foRMaT       = {'%0.0f','%0.0f','%0.0f'};

nOTe         = {'$\textrm{Round of data collection}$/$\textrm{Instrument}$','$\textrm{For consistency, FPHs also collect SBHs (i.e., children alive and ever born). Due to missing values, excluded interviews cannot be used to estimate post-stratification weights.}$'};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,pOP,taBLe,0.190,0.065,[]);
exportgraphics(gcf,char(pATh + "Results/Table_1.png"),'Resolution',RESolUTioN);
s            = [1 2 3 4 5];
bOx          = taBLe(s,:);
pOPsd        = pOP(s);
sets         = sEt;
vars         = vARs;
save(char(pATh + "Results/Tables/Table_1.mat"),'bOx','pOPsd','vars','sets');


clear
pATh         = "/Users/lshjr3/Documents/RaMMPS/under-five mortality/";
load(char(pATh + "Results/RaMMPS.mat"),'RaMMPS','DHSmalawi','MICSmalawi','RESolUTioN');
load(char(pATh + "Results/ReSAmPLiNG.mat"));
dATaDHS      = dATaDHS(end - 2:end);
dATaMICS     = dATaMICS(end - 2:end);

mobile{1}    = ((DHSmalawi.mobile == 1 & DHSmalawi.k == 1 & dATaDHS{1}{1})'*dATaDHS{1}{2})./((DHSmalawi.k == 1 & dATaDHS{1}{1})'*dATaDHS{1}{2})*100;
mobile{2}    = ((MICSmalawi.mobile == 1 & MICSmalawi.k == 1 & dATaMICS{1}{1})'*dATaMICS{1}{2})./((MICSmalawi.k == 1 & dATaMICS{1}{1})'*dATaMICS{1}{2})*100;
round(prctile(mobile{1}(2:end),[50 2.5 97.5]),2)
round(prctile(mobile{2}(2:end),[50 2.5 97.5]),2)

list         = {'$\textrm{SBH}$','$\textrm{TPH}$','$\textrm{FPH}$'};
models       = {'$\textit{post-strat.}$','$\textit{selected}$'};
s            = RaMMPS.sample == 'SBH' | RaMMPS.sample == 'FPH';
dATa         = [{{s,dATa{1}{2}},{s,dATa{2}{2}}} dATa];
listDHS      = {'$\textrm{DHS VII, women 18-49}$'};
listMICS     = {'$\textrm{MICS 6, women 18-49}$'};
modelsDHS    = {'$\textit{all}$','$\textit{m.\,owner*}$','$\textit{m.\,owner}$'};
modelsMICS   = modelsDHS;

for i = 1:numel(list)
    for j = 1:numel(models)
        I          = numel(models)*(i - 1) + j;
        date       = max(RaMMPS.interview);
        ts{I,1}    = datetime([2014 year(date)]',[1 month(date)]',[1 day(date)]');
        date       = eXAcTTime(ts{I,1});
        date       = string(list{i}) + ", " + string(sprintf('%0.1f',date(1))) + "-" + string(sprintf('%0.1f',date(2)));
        pOPsD{I,1} = {char("RaMMPS " + date);models{j}};
        pOPsS{I,1} = list{i};
        clear data I
    end
end

for i = 1:numel(modelsDHS)
    mIn              = min(DHSmalawi.interview);
    mAx              = max(DHSmalawi.interview);
    mAx              = mean([mIn mAx]);
    mIn              = datetime([year(mAx) - 5,month(mAx),day(mAx)],'Format','dd/MM/yyyy');
    date             = [mIn mAx]';
    ts{end + 1,1}    = date;
    date             = eXAcTTime(date);    
    date             = ", " + string(sprintf('%0.1f',date(1))) + "-" + string(sprintf('%0.1f',date(2)));
    date             = char("DHS VII" + date);
    pOPsD{end + 1,1} = {date;modelsDHS{i}};
    pOPsS{end + 1,1} = listDHS{1};
    clear deta mIn mAx
end

for i = 1:numel(modelsMICS)
    mIn              = min(MICSmalawi.interview);
    mAx              = max(MICSmalawi.interview);
    mAx              = mean([mIn mAx]);
    mIn              = datetime([year(mAx) - 5,month(mAx),day(mAx)],'Format','dd/MM/yyyy');
    date             = [mIn mAx]';
    ts{end + 1,1}    = date;
    date             = eXAcTTime(date);    
    date             = ", " + string(sprintf('%0.1f',date(1))) + "-" + string(sprintf('%0.1f',date(2)));
    date             = char("MICS 6" + date);
    pOPsD{end + 1,1} = {date;modelsMICS{i}};
    pOPsS{end + 1,1} = listMICS{1};
    clear deta mIn mAx
end


pACk         = {{RaMMPS,dATa},{DHSmalawi,dATaDHS},{MICSmalawi,dATaMICS}};
lABelS       = {'Place of residence','Region','Age','Education','Household size','Electricity','Drinking water','Roofing'};
lABelSd      = {{'urban' 'rural'} {'North' 'Central' 'South'} {'18-29' '30-39' '40-49'} {'less than complete primary' 'incomplete secondary' 'complete secondary or more'} {'1-4' '5-8' '9+'} {'access' 'no access'} {'safe source' 'other source'} {'durable material' 'other material'}}; 
outcomes     = {[1 2],[1 2 3],[1 2 3],[1 2 3],[1 2 3],[2 1],[2 1],[2 1]};

H            = 0;
for h = 1:numel(pACk)
    d    = pACk{h}{1};
    data = [d.UR,d.Region,d.GO,d.Education,d.household,d.Electricity,d.Water,d.Roofing];
    for i = 1:numel(pACk{h}{2})
        H  = H + 1;
        s  = pACk{h}{2}{i}{1} & (d.k == 1);
        sW = pACk{h}{2}{i}{2}(s,:);
        I  = 0;
        for j = 1:numel(outcomes)
            for k = 1:numel(outcomes{j})
                I        = I + 1;
                BST      = ((data(s,j) == outcomes{j}(k))'*sW)./sum(sW,1);
                bOx{I,H} = prctile(100*BST(2:end),[50 2.5 97.5]);
                
                if isequal(h,1) & isequal(i,1)
                    if k == 1
                        pOPsd{I,1} = {char(string(lABelS{j}) + ": " + string(lABelSd{j}{k}))};
                    else
                        pOPsd{I,1} = {lABelSd{j}{k}};
                    end
                end
            end
        end
        
        if isequal(h,1) & isequal(i,1)
            bOx{end + 1,H}   = [sum(s) NaN NaN];
            pOPsd{end + 1,1} = {'Observations'};            
        else
            bOx{end,H}       = [sum(s) NaN NaN];
        end
    end
end

selection    = 1 + [0 cumsum([ones(1,2)*numel(models) numel(modelsDHS) numel(modelsMICS)])];
for i = 1:numel(selection)
    sEt{i} = {pOPsS{selection(i)}};
end

lABs         = {{1} {3 4 5} {6 7 8} {9 10 11} {12 13 14} {15} {17} {19} {21}};
vARs         = {models models models modelsMICS([1 3]) modelsDHS([1 3])};
foRMaT       = {'%0.2f','%0.2f','%0.2f'};

nOTe         = {'$\textrm{Attributes}$/$\textrm{Instrument}$','$\mathit{Bootstrapping}$ $\mathrm{p50}$/$\mathit{[p2.5,p97.5]}$'};
tABleBAyEs(sEt([1 2 3 5 4]),vARs,foRMaT,lABs,nOTe,pOPsd,cell2mat(bOx(:,[1 2 3 4 5 6 10 12 7 9])),0.190,0.065,[]);
exportgraphics(gcf,char(pATh + "Results/Table_2.png"),'Resolution',RESolUTioN);
s            = [1 3 4 5 6 7 8 9 10 11 12 13 14 15 17 19 21];
bOx          = bOx(s,:);
pOPsd        = pOPsd(s);
sets         = sEt;
vars         = vARs;
save(char(pATh + "Results/Tables/Table_2.mat"),'bOx','pOPsd','vars','sets');


clear
pATh       = "/Users/lshjr3/Documents/RaMMPS/under-five mortality/";
load(char(pATh + "Results/RaMMPS.mat"),'RaMMPS','RESolUTioN');
list       = {'TPH','FPH'};
RaMMPS     = RaMMPS(RaMMPS.birth == 'livebirth' & (RaMMPS.survival == 'dead' | RaMMPS.survival == 'alive'),:);
RaMMPS     = RaMMPS(~isnat(RaMMPS.B_min) & RaMMPS.B_min <= RaMMPS.interview,:);
leGend     = {'$\textbf{Date of death}$','$\textbf{Date of interview}$'};

for i = 1:numel(list)
    alpha      = min(RaMMPS.B_min + (RaMMPS.B_max - RaMMPS.B_min)/2,RaMMPS.interview);
    omega      = min(alpha + (RaMMPS.D_max + RaMMPS.D_min)/2,RaMMPS.interview);
    dATa{i}    = RaMMPS(RaMMPS.sample == list{i} & alpha >= datetime(2009,1,1) & omega >= datetime(2014,1,1) & alpha < RaMMPS.interview,:);
    m          = tabulate(dATa{i}.caseid);
    mOTheRs{i} = sum(cell2mat(m(:,2)) > 0);
    bIRtHs{i}  = size(dATa{i},1);
end

mPIX                     = 538756;
pix                      = 1/37.7952755906;
z                        = min(sqrt(mPIX/((10*1)*(10*numel(dATa))/pix^2)),1);
fi                       = figure('Color',[1 1 1],'Position',z*10*[0 0 numel(dATa) 1]/pix,'Theme','light');
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(1,numel(dATa),'Padding','compact','TileSpacing','compact');
for i = 1:numel(dATa)
    nexttile(i)
    d                           = scatter(datetime('01-Jan-2021'),0);
    delete(d)
    ax{i}                       = gca;
    ax{i}.FontSize              = 10*z;
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    ax{i}.YMinorGrid            = 'on';
    ax{i}.XMinorGrid            = 'on';
    ax{i}.XAxis.TickLabelFormat = 'yyyy';
    ax{i}.XTick                 = datetime(2014:2:2024,1,1);
    ax{i}.XAxis.MinorTickValues = datetime(2014:2024,1,1);
    ax{i}.YTick                 = 0:1:5;
    ax{i}.YAxis.MinorTickValues = 0:.5:5;
    if i == 1
        ylabel('$\textbf{Age (in years)}$','Interpreter','latex','FontSize',11*z);
    end
    xlabel('$\textbf{Year}$','Interpreter','latex','FontSize',11*z);
    xlim([datetime(2014,1,1) datetime(2024,1,1)]);
    ylim([0 8])

    T                           = ax{i}.XAxis.TickLabels;
    for j = 1:numel(T)
        T{j} = char("$\mathbf{" + T{j} + "}$");
    end
    ax{i}.XAxis.TickLabels           = T;
    ax{i}.XAxis.TickLabelInterpreter = 'latex';
    clear T
        
    t                           = ax{i}.YAxis.TickValues;
    for j = 1:numel(t)
        T{j} = char("$\mathbf{" + string(sprintf(ax{i}.YAxis.TickLabelFormat,t(j))) + "}$");
    end
    ax{i}.YAxis.TickLabels           = T;
    ax{i}.YAxis.TickLabelInterpreter = 'latex';
    clear t T

    temp                        = max(dATa{i}.interview);
    eXACt{i}                    = [datetime(2014,1,1) temp];
    temp                        = round(year(temp) + (day(temp,'dayofyear') - 1)./days(datetime(year(temp) + 1,1,1) - datetime(year(temp),1,1)),1);
    temp                        = char(string(list{i}) + ", " + string(sprintf('%0.1f',2014)) + "-" +string(sprintf('%0.1f',temp)));
    pOPs{i}                     = char("RaMMPS " + string(temp));
    temp                        = char("$\textbf{" + string(char(96 + i)) + ". " + string(temp) + "}$");
    temp2                       = char("$\textbf{" + string(mOTheRs{i}) + " mothers and " + string(bIRtHs{i}) + " births}$");
    title({temp;temp2},'Interpreter','latex','FontSize',11*z);
    clear d temp temp2
    grid on;
    box on;
    hold on;
end
for i = 1:numel(dATa)
    nexttile(i);
    scatter(NaT,NaN,10,'filled','MarkerFaceColor','r','MarkerFaceAlpha',1,'MarkerEdgeColor','m','MarkerEdgeAlpha',.5);
    scatter(NaT,NaN,2.5,'filled','MarkerFaceColor','b','MarkerFaceAlpha',1,'MarkerEdgeColor','k','MarkerEdgeAlpha',.0);
    s = dATa{i};
    s = s(s.survival ~= 'unknown',:);
    for j = 1:size(s,1)
        if isequal(s.survival(j),'dead')
            r  = rand(2);
            r  = ones(size(r))/2;
            B  = r(1)*datenum(s.B_min(j)) + (1 - r(1))*datenum(s.B_max(j));
            xd = r(2)*s.D_min(j) + (1 - r(2))*s.D_max(j);
            D  = datetime(datevec(B + xd));
            B  = datetime(datevec(B));
            I  = s.interview(j);
            D  = min(D,I - 30*rand);
            B  = min(B,D - xd);
            xd = xd/365.25;
            xi = (datenum(I) - datenum(B))/365.25;
            plot([B D],[0 xd],'color',[1.0 0.0 0.0 1.00],'LineWidth',0.75);
            plot([D I],[xd xi],'color',[1.0 0.0 0.0 0.25],'LineStyle','-.','LineWidth',0.5);
            scatter(D,xd,10,'filled','MarkerFaceColor','r','MarkerFaceAlpha',1,'MarkerEdgeColor','m','MarkerEdgeAlpha',.5);
            scatter(I,xi,2.5,'filled','MarkerFaceColor','b','MarkerFaceAlpha',1,'MarkerEdgeColor','k','MarkerEdgeAlpha',.0);
        else
            r  = rand(1);
            r  = ones(size(r))/2;
            B  = r*datenum(s.B_min(j)) + (1 - r)*datenum(s.B_max(j));
            B  = datetime(datevec(B));
            I  = s.interview(j);
            xi = (datenum(I) - datenum(B))/365.25;
            plot([B I],[0 xi],'color',[0.0 0.0 0.0 0.025],'LineWidth',0.5);
            scatter(I,xi,2.5,'filled','MarkerFaceColor','b','MarkerFaceAlpha',.25,'MarkerEdgeColor','k','MarkerEdgeAlpha',.0);
        end
    end
    plot(datetime(ones(2,1)*[2020 1 1]),[0 8],'color',[0.1 0.4 0.25],'LineStyle','--','LineWidth',0.75);
    if i == 2
        legend(leGend(1:2),'Interpreter','latex','FontSize',9*z,'Location','southoutside','NumColumns',2,'Box','off');
    end
    clear B D I xi xd r
end
exportgraphics(gcf,char(pATh + "Results/Figure_1.png"),'Resolution',RESolUTioN);