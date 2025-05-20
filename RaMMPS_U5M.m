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
save(char(pATh + "Results/RaMMPSf.mat"),'RaMMPS');
RaMMPS                = RaMMPS(RaMMPS.caseid ~= "0A0" & isundefined(RaMMPS.missing),:);
RaMMPSHH              = RaMMPSHH(RaMMPSHH.caseid ~= "0A0" & isundefined(RaMMPSHH.missing),:);
RaMMPSdst             = RaMMPSdst(isundefined(RaMMPSdst.missing),:);
save(char(pATh + "Results/RaMMPS.mat"));

clear
pATh                  = "/Users/lshjr3/Documents/RaMMPS/under-five mortality/";
load(char(pATh + "Results/RaMMPS.mat"),'RaMMPS');
x{1}                  = [(0:7:28)/365.25,[(2:1:12),(15:3:24),(36:12:60)]/12]';
x{2}                  = [string(0);string([(7:7:28)';(2:1:11)';(12:3:24)';(36:12:60)']) + char([kron('d',ones(4,1));kron('m',ones(18,1))])];
n                     = x{1}(2:end) - x{1}(1:end - 1);
instrument            = {'TPH','FPH'};
date                  = max(RaMMPS.interview);
Ts                    = {datetime([2014 year(date)]',[1 month(date)]',[1 day(date)]'),datetime([2014 2016]',1,1),datetime([2016 2018]',1,1),datetime([2018 2020]',1,1),datetime([2020 2022]',1,1),datetime([2022 year(date)]',[1 month(date)]',[1 day(date)]')};
sET                   = [RaMMPS.WR(RaMMPS.k == 1),ones(sum(RaMMPS.k == 1),1)];
models                = {'post-strat.','selected'};
R                     = 1000;
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
    for j = 1:size(models,2)
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


setH                  = ones(size(RaMMPS.WR));
setH(:,end + 1)       = (RaMMPS.sex == 2);
setH(:,end + 1)       = (RaMMPS.sex == 1);
setH(:,end + 1)       = (RaMMPS.UR == 1);
setH(:,end + 1)       = (RaMMPS.UR == 2);
setH(:,end + 1)       = (RaMMPS.Region == 1);
setH(:,end + 1)       = (RaMMPS.Region == 2);
setH(:,end + 1)       = (RaMMPS.Region == 3);
setH(:,end + 1)       = (RaMMPS.Education == 1 | RaMMPS.Education == 2);
setH(:,end + 1)       = (RaMMPS.Education == 3);
setH(:,end + 1)       = (RaMMPS.Electricity == 2);
setH(:,end + 1)       = (RaMMPS.Electricity == 1);
T                     = eXAcTTime(Ts{1});

for i = 1:numel(dATa)
    for h = 1:size(setH,2)
        s               = dATa{i}{1} & setH(:,h); 
        wH              = dATa{i}{2}(s,:);
        for j = 1:numel(x{1}) - 1
            a             = max(min(B(s,:) + x{1}(j),min(O(s,:),T(2))),T(1));
            o             = max(min(B(s,:) + x{1}(j + 1),min(O(s,:),T(2))),T(1)); 
            exposure(j,:) = sum((o - a).*wH.*sEL(s));
            events(j,:)   = sum((d(s,:) >= x{1}(j) & d(s,:) < x{1}(j + 1) & D(s,:) >= T(1) & D(s,:) < T(2)).*wH.*sEL(s));
            clear a o
            clc;
            j/(numel(x{1}) - 1)
        end
        
        stillbirths     = sum((RaMMPS.FlagBI(s) ~= 1 & RaMMPS.gestation(s) >= 28 & sB(s,:) >= T(1) & sB(s,:) < T(2)).*wH);
        births          = sum((B(s,:) >= T(1) & B(s,:) < T(2)).*wH);
        stillbirths     = stillbirths./(births + stillbirths);        
        SRB             = sum((B(s,:) >= T(1) & B(s,:) < T(2) & RaMMPS.sex(s) == 1).*wH)./sum((B(s,:) >= T(1) & B(s,:) < T(2) & RaMMPS.sex(s) == 2).*wH);
        m               = events./exposure;
        q               = 1 - [ones(1,R + 1);exp(-cumsum(m.*n))];
        
        TaBlE.q{h,i}    = q;
        TaBlE.s{h,i}    = 1 - (1 - q([1 2 5],:)).*(1 - stillbirths);
        TaBlE.tAU{h,i}  = T;
        clear events exposure q m births stillbirths SRB wH s
    end
end
clear T i sH interview sH sEL p d D B O j h sB


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
mIn                   = min(DHSmalawi.interview);
mAx                   = max(DHSmalawi.interview);
mAx                   = mean([mIn mAx]);
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
    clc;
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
            j/(numel(x{1}) - 1)
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
        (j - A(1) + 1)/(A(2) - A(1) + 1)
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

setH             = ones(size(DHSmalawi.W));
setH(:,end + 1)  = (DHSmalawi.sex == 2);
setH(:,end + 1)  = (DHSmalawi.sex == 1);
setH(:,end + 1)  = (DHSmalawi.UR == 1);
setH(:,end + 1)  = (DHSmalawi.UR == 2);
setH(:,end + 1)  = (DHSmalawi.Region == 1);
setH(:,end + 1)  = (DHSmalawi.Region == 2);
setH(:,end + 1)  = (DHSmalawi.Region == 3);
setH(:,end + 1)  = (DHSmalawi.Education == 1 | DHSmalawi.Education == 2);
setH(:,end + 1)  = (DHSmalawi.Education == 3);
setH(:,end + 1)  = (DHSmalawi.Electricity == 2);
setH(:,end + 1)  = (DHSmalawi.Electricity == 1);

for i = 1:numel(dATaDHS)
    k                = i + numel(dATa) + numel(dATaB);
    for h = 1:size(setH,2)
        s               = dATaDHS{i}{1} & setH(:,h); 
        wH              = dATaDHS{i}{2}(s,:);
        for j = 1:numel(x{1}) - 1
            a             = max(min(B(s,:) + x{1}(j),min(O(s,:),Tdhs(2))),Tdhs(1));
            o             = max(min(B(s,:) + x{1}(j + 1),min(O(s,:),Tdhs(2))),Tdhs(1));
            exposure(j,:) = sum((o - a).*wH);
            events(j,:)   = sum((d(s,:) >= x{1}(j) & d(s,:) < x{1}(j + 1) & D(s,:) >= Tdhs(1) & D(s,:) < Tdhs(2)).*wH);
            clear a o
            clc;
            j/(numel(x{1}) - 1)
        end
        
        births          = sum(DHSmalawi.births(s).*wH);
        stillbirths     = sum(DHSmalawi.stillbirths(s).*wH);
        stillbirths     = stillbirths./(births + stillbirths);        
        SRB             = ((DHSmalawi.sex(s) == 1 & B(s,:) >= Tdhs(1) & B(s,:) < Tdhs(2))'*wH)./((DHSmalawi.sex(s) == 2 & B(s,:) >= Tdhs(1) & B(s,:) < Tdhs(2))'*wH);
        m               = events./exposure;
        q               = 1 - [ones(1,R + 1);exp(-cumsum(m.*n))];
        
        TaBlE.q{h,k}    = q;
        TaBlE.s{h,k}    = 1 - (1 - q([1 2 5],:)).*(1 - stillbirths);
        TaBlE.tAU{h,k}  = Tdhs;
        clear events exposure q m births stillbirths SRB s wH
    end
    clear h j k
end
clear p d b B D O dob ageS ageB sET i



load(char(pATh + "Results/RaMMPS.mat"),'MICSmalawi');
mIn                   = min(MICSmalawi.interview);
mAx                   = max(MICSmalawi.interview);
mAx                   = mean([mIn mAx]);
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
            j/(numel(x{1}) - 1)
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
        (j - A(1) + 1)/(A(2) - A(1) + 1)
    end
    
    TaBle.TFR{1,h}       = TFR;
    TaBle.sRB{1,h}       = sum((MICSmalawi.sex(s) == 1).*w(s,:))./sum((MICSmalawi.sex(s) == 2).*w(s,:));
    TaBle.parity{1,h}    = sum(~isnan(MICSmalawi.bidx(s)).*w(s,:))./sum((MICSmalawi.k(s) == 1).*w(s,:));
    TaBle.childless{1,h} = sum((MICSmalawi.mother(s) ~= 1 & MICSmalawi.k(s) == 1).*w(s,:))./sum((MICSmalawi.k(s) == 1).*w(s,:))*100;
    clear h j TFR w A
end

for i = 1:numel(dATaDHS)
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

[qS,qP]            = RaMMPS_Bayes(000,0);
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
saveas(gcf,char(pATh + "Results/RaMMPS-U5M Table 2.png"));
clear table label temp

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
saveas(gcf,char(pATh + "Results/RaMMPS-U5M Table 3.png"));
clear table label temp


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
saveas(gcf,char(pATh + "Results/RaMMPS-U5M Table 4.png"));
clear table label temp



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
saveas(gcf,char(pATh + "Results/RaMMPS-U5M Table A1.png"));
clear table label temp

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
saveas(gcf,char(pATh + "Results/RaMMPS-U5M Table A2.png"));
clear table label temp


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
saveas(gcf,char(pATh + "Results/RaMMPS-U5M Table A3.png"));
clear table label temp


selection  = {[1 1 2],[1 3 2],[1 8 1],[1 9 2],[1 15 1],[2 15 1],[3 15 1]};
for i = 1:numel(selection)
    label{i,1} = pOPs{selection{i}(1),selection{i}(2)};
    for j = 1:3
        for k = 1:mAX
            h          = k + (j - 1)*mAX;
            table{i,h} = NaN(1,3);
        end
        
        for k = 1:selection{i}(3)
            h          = k + (j - 1)*mAX;
            temp       = TaBle.s{selection{i}(1),selection{i}(2) + k - 1}(j,:);
            if size(temp,2) == R + 1
                table{i,h} = prctile(temp(:,2:end)',tHtiles)*1000;
            else
                table{i,h} = temp*1000;
            end
        end
    end
end

nOTe         = {'$\textrm{Instrument}$/$\textit{method}$',''};
sEt          = {'Stillbirth $\mathrm{[28}\mathit{w}\mathrm{,}\mathrm{B}\mathrm{)}$','Perinatal $\mathrm{[28}\mathit{w}\mathrm{,}\mathrm{7}\mathit{d}\mathrm{)}$','Perinatal $\mathrm{[28}\mathit{w}\mathrm{,}\mathrm{28}\mathit{d}\mathrm{)}$'};
vARs         = {models models models};
lABs         = {{1} {2} {3 4} {5 6 7}};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,label,cell2mat(table),0.175,0.060,[]);
saveas(gcf,char(pATh + "Results/RaMMPS-perinatal Table 2.png"));
s            = [1 2 3 4 5 6 7];
bOx          = table(s,:);
pOPsd        = label(s);
sets         = sEt;
vars         = vARs;
save(char(pATh + "Results/Tables/perinatal_Table_2.mat"),'bOx','pOPsd','vars','sets');
clear table labels temp bOx pOPsd vars sets s


selection  = {[2 1 2],[3 1 2],[4 1 2],[5 1 2],[6 1 2],[2 3 2],[3 3 2],[4 3 2],[5 3 2],[6 3 2],[1 15 1],[2 15 1],[3 15 1]};
mAX        = 2;
for i = 1:numel(selection)
    label{i,1} = pOPs{selection{i}(1),selection{i}(2)};
    for j = 1:3
        for k = 1:mAX
            h          = k + (j - 1)*mAX;
            table{i,h} = NaN(1,3);
        end
        
        for k = 1:selection{i}(3)
            h          = k + (j - 1)*mAX;
            temp       = TaBle.s{selection{i}(1),selection{i}(2) + k - 1}(j,:);
            if size(temp,2) == R + 1
                table{i,h} = prctile(temp(:,2:end)',tHtiles)*1000;
            else
                table{i,h} = temp*1000;
            end
        end
    end
end

lABs         = {{1 2 3 4 5} {6 7 8 9 10} {11 12 13}};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,label,cell2mat(table),0.175,0.060,[]);
saveas(gcf,char(pATh + "Results/RaMMPS-perinatal Table S5.png"));
s            = [1 2 3 4 5 6 7 8 9 10 11 12 13];
bOx          = table(s,:);
pOPsd        = label(s);
sets         = sEt;
vars         = vARs;
save(char(pATh + "Results/Tables/perinatal_Table_S5.mat"),'bOx','pOPsd','vars','sets');
clear table temp bOx pOPsd vars sets s



selection  = {[1 1 2],[1 3 2],[1 8 1],[1 9 2]};
for i = 1:numel(selection)
    label{i,1} = pOPs{selection{i}(1),selection{i}(2)};
    for j = 1:5
        for k = 1:mAX
            h          = k + (j - 1)*mAX;
            table{i,h} = NaN(1,3);
        end
        
        for k = 1:selection{i}(3)
            h          = k + (j - 1)*mAX;
            temp       = TaBle.SNR{selection{i}(1),selection{i}(2) + k - 1}(j,:);
            if size(temp,2) == R + 1
                table{i,h} = prctile(temp(:,2:end)',tHtiles)*1000;
            else
                table{i,h} = temp*1000;
            end
        end
    end
end

sEt          = {'Stillbirth $\mathrm{[28}\mathit{w}\mathrm{,}\mathrm{B}\mathrm{)}$','ENMR - $\mathit{q}\mathrm{(7}\mathit{d}\mathrm{)}$','NMR - $\mathit{q}\mathrm{(28}\mathit{d}\mathrm{)}$','Stillbirth/NMR','ENMR/NMR'};
vARs         = {models models models models models};
lABs         = {{1} {2} {3 4}};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,label,cell2mat(table),0.135,0.060,[]);
saveas(gcf,char(pATh + "Results/RaMMPS-perinatal Table 3.png"));
s            = [1 2 3 4];
bOx          = table(s,:);
pOPsd        = label(s);
sets         = sEt;
vars         = vARs;
save(char(pATh + "Results/Tables/perinatal_Table_3.mat"),'bOx','pOPsd','vars','sets');
clear table labels temp bOx pOPsd vars sets s


tAU{1}       = ", " + string(sprintf('%0.1f',TaBlE.tAU{1,1}(1))) + "-" + string(sprintf('%0.1f',TaBlE.tAU{1,1}(2)));
tAU{2}       = ", " + string(sprintf('%0.1f',TaBlE.tAU{1,3}(1))) + "-" + string(sprintf('%0.1f',TaBlE.tAU{1,3}(2)));
tAU{3}       = ", " + string(sprintf('%0.1f',TaBlE.tAU{1,7}(1))) + "-" + string(sprintf('%0.1f',TaBlE.tAU{1,7}(2)));
sEt          = {"$\textrm{TPH" + tAU{1} + "}$","$\textrm{FPH" + tAU{2} + "}$","$\textrm{DHS-VII" + tAU{3} + "}$"};
DHS          = {'$\textit{All women}$','$\textit{MP owners*}$','$\textit{MP owners}$'};
vARs         = {models models DHS};
selection    = [1 2 3 4 7 8 9];
for i = 1:size(TaBlE.s,1)
    for j = 1:numel(selection)
        bOx{i,j}  = prctile(TaBlE.s{i,selection(j)}(1,2:end),[50 2.5 97.5])*1000;
        bOx2{i,j} = prctile(TaBlE.s{i,selection(j)}(2,2:end),[50 2.5 97.5])*1000;
    end
end
lABs         = {{1} {4 5} {6 7 8} {9 10} {11 12}};
nOTe         = {'$\textrm{Attributes}$/$\textrm{Instrument\&}\textit{method}$',''};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,pOPh,cell2mat(bOx),0.175,0.070,[]);
saveas(gcf,char(pATh + "Results/RaMMPS-perinatal Table S4a.png"));
s            = [1 4 5 6 7 8 9 10 11 12];
bOx          = bOx(s,:);
pOPsd        = pOPh(s);
sets         = sEt;
vars         = vARs;
save(char(pATh + "Results/Tables/perinatal_Table_S4a.mat"),'bOx','pOPsd','vars','sets');

tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,pOPh,cell2mat(bOx2),0.175,0.070,[]);
saveas(gcf,char(pATh + "Results/RaMMPS-perinatal Table S4b.png"));
bOx          = bOx2(s,:);
pOPsd        = pOPh(s);
sets         = sEt;
vars         = vARs;
save(char(pATh + "Results/Tables/perinatal_Table_S4b.mat"),'bOx','pOPsd','vars','sets');
clear bOx pOPsd vars sets s


P                        = cell(0);
selection                = {[1 1],[1 3],[2 5],[1 8],[1 12],[1 15]};
selectionT               = {[1 1],[1 3],[3 5],[1 8],[1 12],[5 15]};
LAB                      = {'Neonatal mortality rate $\mathit{q}\mathrm{(28}\mathit{d}\mathrm{)}$','Infant mortality rate $\mathit{q}\mathrm{(12}\mathit{m}\mathrm{)}$','Under-five mortality rate $\mathit{q}\mathrm{(60}\mathit{m}\mathrm{)}$'};
load(char(pATh + "Results/paleTTe.mat"),'paleTTe');
coloR                    = paleTTe([1 2 3 7 4 8]);
pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 21 14]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
for i = 1:6
    nexttile(i)
    if i > 3
        ax{i}                       = gca;
        ax{i}.FontName              = 'Times New Roman';
        ax{i}.FontSize              = 10;
        ax{i}.XScale                = 'log';
        ax{i}.XAxis.TickLabelFormat = '%.1f';
        ax{i}.YAxis.TickLabelFormat = '%.2f';
        ax{i}.XTick                 = .1./(2.^(9:-1:-2))*1000;
        ax{i}.XAxis.MinorTickValues = 10:10:200;
        xlim([6.25 100])
        xlabel('$\mathit{deaths}$ $\mathit{per}$ $\mathit{1000}$ $\mathit{births}$ (log scale)','Interpreter','latex','FontName','Times New Roman','FontSize',11);
        ylabel('$\mathit{kernel}$ $\mathit{density}$','Interpreter','latex','FontName','Times New Roman','FontSize',11);
        title(char(string(char(96 + i)) + ". " + string(LAB{i - 3})),'Interpreter','latex');
    else
        d                           = scatter(datetime('01-Jan-2021'),0);
        delete(d)
        ax{i}                       = gca;
        ax{i}.FontName              = 'Times New Roman';
        ax{i}.FontSize              = 10;
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
        xlabel('$\mathit{year}$','Interpreter','latex','FontName','Times New Roman','FontSize',11);
        ylabel('$\mathit{deaths}$ $\mathit{per}$ $\mathit{1000}$ $\mathit{births}$ (log scale)','Interpreter','latex','FontName','Times New Roman','FontSize',11);
        xlim([datetime(2010,1,1) datetime(2024,1,1)]);
        title(char(string(char(96 + i)) + ". " + string(LAB{i})),'Interpreter','latex','FontName','Times New Roman','FontSize',12);
    end
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
            if numel(q{j}) == R + 1
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
            P{end + 1} = plot(max(xi{j},eps),f{j},'color',coloR{j},'LineWidth',1.00);
        end
        for j = 1:numel(selection)
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.75,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.50);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.50);
        end
        if i == 5
            for j = 1:numel(selection)
                leGend{j} = pOPs{selection{j}(1),selection{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    else
        for j = 1:numel(selectionT)
            tAU{j} = TaBle.tAU{selectionT{j}(1),selectionT{j}(2)};
            if numel(tAU{j}) > 2
                q{j} = TaBle.q{selectionT{j}(1),selectionT{j}(2)}{1 + mod(i + 2,3)};
            else
                q{j} = kron(ones(2,1),TaBle.q{selectionT{j}(1),selectionT{j}(2)}(aGEs(1 + mod(i + 2,3)),:));
            end
            if size(q{j},2) == R + 1
                q{j} = prctile(q{j}(:,2:end)',tHtiles)'*1000;
            else
                q{j} = q{j}*1000;
            end
            P{end + 1} = fill([tAU{j};flip(tAU{j})],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j});
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j},q{j}(:,1),'color',coloR{j},'LineWidth',1.0,'LineStyle',':');
        end
        P{end + 1} = plot(datetime(ones(2,1)*[2020 1 1]),[3.125 200],'color',[0.1 0.4 0.5],'LineStyle','--','LineWidth',0.75);
        if i == 2
            for j = 1:numel(selectionT)
                leGend{j} = pOPs{selectionT{j}(1),selectionT{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    end
end
saveas(gcf,char(pATh + "Results/RaMMPS-U5M Figure 2.png"));
for i = 1:numel(P)
    delete(P{i})
end
clear leGend P


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
            if numel(q{j}) == R + 1
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
            P{end + 1} = plot(max(xi{j},eps),f{j},'color',coloR{j},'LineWidth',1.00);
        end
        for j = 1:numel(selection)
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.75,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.50);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.50);
        end
        if i == 5
            for j = 1:numel(selection)
                leGend{j} = pOPs{selection{j}(1),selection{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    else
        for j = 1:numel(selectionT)
            tAU{j} = TaBle.tAU{selectionT{j}(1),selectionT{j}(2)};
            if numel(tAU{j}) > 2
                q{j} = TaBle.q{selectionT{j}(1),selectionT{j}(2)}{1 + mod(i + 2,3)};
            else
                q{j} = kron(ones(2,1),TaBle.q{selectionT{j}(1),selectionT{j}(2)}(aGEs(1 + mod(i + 2,3)),:));
            end
            if size(q{j},2) == R + 1
                q{j} = prctile(q{j}(:,2:end)',tHtiles)'*1000;
            else
                q{j} = q{j}*1000;
            end
            P{end + 1} = fill([tAU{j};flip(tAU{j})],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j});
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j},q{j}(:,1),'color',coloR{j},'LineWidth',1.0,'LineStyle',':');
        end
        P{end + 1} = plot(datetime(ones(2,1)*[2020 1 1]),[3.125 200],'color',[0.1 0.4 0.5],'LineStyle','--','LineWidth',0.75);
        if i == 2
            for j = 1:numel(selectionT)
                leGend{j} = pOPs{selectionT{j}(1),selectionT{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    end
end
saveas(gcf,char(pATh + "Results/RaMMPS-U5M Figure 3.png"));
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
            if numel(q{j}) == R + 1
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
            P{end + 1} = plot(max(xi{j},eps),f{j},'color',coloR{j},'LineWidth',1.00);
        end
        for j = 1:numel(selection)
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.75,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.50);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.50);
        end
        if i == 5
            for j = 1:numel(selection)
                leGend{j} = pOPs{selection{j}(1),selection{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',4,'Box','off');
        end
    else
        for j = 1:numel(selectionT)
            tAU{j} = TaBle.tAU{selectionT{j}(1),selectionT{j}(2)};
            if numel(tAU{j}) > 2
                q{j} = TaBle.q{selectionT{j}(1),selectionT{j}(2)}{1 + mod(i + 2,3)};
            else
                q{j} = kron(ones(2,1),TaBle.q{selectionT{j}(1),selectionT{j}(2)}(aGEs(1 + mod(i + 2,3)),:));
            end
            if size(q{j},2) == R + 1
                q{j} = prctile(q{j}(:,2:end)',tHtiles)'*1000;
            else
                q{j} = q{j}*1000;
            end
            P{end + 1} = fill([tAU{j};flip(tAU{j})],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j});
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j},q{j}(:,1),'color',coloR{j},'LineWidth',1.0,'LineStyle',':');
        end
        P{end + 1} = plot(datetime(ones(2,1)*[2020 1 1]),[3.125 200],'color',[0.1 0.4 0.5],'LineStyle','--','LineWidth',0.75);
        if i == 2
            for j = 1:numel(selectionT)
                leGend{j} = pOPs{selectionT{j}(1),selectionT{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',4,'Box','off');
        end
    end
end
saveas(gcf,char(pATh + "Results/RaMMPS-U5M Figure 4.png"));
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
            P{end + 1} = plot(max(xi{j},eps),f{j},'color',coloR{j},'LineWidth',1.00);
        end
        for j = 1:numel(selection)
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.75,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.50);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.50);
        end
        if i == 5
            for j = 1:numel(selection)
                leGend{j} = pOPs{selection{j}(1),selection{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
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
            P{end + 1} = fill([tAU{j};flip(tAU{j})],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j});
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j},q{j}(:,1),'color',coloR{j},'LineWidth',1.0,'LineStyle',':');
        end
        P{end + 1} = plot(datetime(ones(2,1)*[2020 1 1]),[3.125 200],'color',[0.1 0.4 0.5],'LineStyle','--','LineWidth',0.75);
        if i == 2
            for j = 1:numel(selectionT)
                leGend{j} = pOPs{selectionT{j}(1),selectionT{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    end
end
saveas(gcf,char(pATh + "Results/RaMMPS-U5M Figure 3-Reg.png"));
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
            P{end + 1} = plot(max(xi{j},eps),f{j},'color',coloR{j},'LineWidth',1.00);
        end
        for j = 1:numel(selection)
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.75,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.50);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.50);
        end
        if i == 5
            for j = 1:numel(selection)
                leGend{j} = pOPs{selection{j}(1),selection{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',4,'Box','off');
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
            P{end + 1} = fill([tAU{j};flip(tAU{j})],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j});
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j},q{j}(:,1),'color',coloR{j},'LineWidth',1.0,'LineStyle',':');
        end
        P{end + 1} = plot(datetime(ones(2,1)*[2020 1 1]),[3.125 200],'color',[0.1 0.4 0.5],'LineStyle','--','LineWidth',0.75);
        if i == 2
            for j = 1:numel(selectionT)
                leGend{j} = pOPs{selectionT{j}(1),selectionT{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',4,'Box','off');
        end
    end
end
saveas(gcf,char(pATh + "Results/RaMMPS-U5M Figure 4-Reg.png"));
for i = 1:numel(P)
    delete(P{i})
end
clear leGend P


LAB                      = {'$\textrm{RaMMPS}$','$\textrm{DHS VII}$','$\textrm{MICS 6}$'};
pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 21 7]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
wEIGhTS                  = {RaMMPS.WR(RaMMPS.k == 1),DHSmalawi.WR(DHSmalawi.k == 1 & ~isnan(DHSmalawi.WR)),MICSmalawi.WR(MICSmalawi.k == 1 & ~isnan(MICSmalawi.WR))}; 
for i = 1:3
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 10;
    ax{i}.XAxis.TickLabelFormat = '%.1f';
    ax{i}.YAxis.TickLabelFormat = '%.2f';
    xlabel('$\textit{post-stratified weights}$','Interpreter','latex','FontName','Times New Roman','FontSize',11);
    ylabel('$\textit{probability density function}$','Interpreter','latex','FontName','Times New Roman','FontSize',11);
    xlim([0 15])    
    grid on;
    box on;
    hold on;
    
    y      = wEIGhTS{i};
    x      = min(y):.25:max(y);
    histogram(y,x,'Normalization','pdf','FaceColor',coloR{i},'EdgeColor',[0 0 0],'FaceAlpha',0.25,'EdgeAlpha',0.45);
    title(char(string(char(96 + i)) + ". " + string(LAB{i}) + ", N = " + string(numel(y))),'Interpreter','latex','FontName','Times New Roman','FontSize',12);
    f      = ksdensity(y,x,'Function','pdf','Bandwidth',.1);
    plot(x,f,'LineWidth',0.75,'Color',[coloR{i} 0.75]);
    clear x y f
end
saveas(gcf,char(pATh + "Results/RaMMPS-U5M Figure 5.png"));









P                        = cell(0);
selection                = {[1 1],[1 3],[1 8],[1 15]};
selectionT               = {[1 1],[1 3],[1 8],[5 15]};
LAB                      = {'Stillbirth $\mathrm{[28}\mathit{w}\mathrm{,}\mathrm{B}\mathrm{)}$','Perinatal $\mathrm{[28}\mathit{w}\mathrm{,}\mathrm{7}\mathit{d}\mathrm{)}$','Perinatal $\mathrm{[28}\mathit{w}\mathrm{,}\mathrm{28}\mathit{d}\mathrm{)}$'};
coloR                    = {[0.00 0.00 0.75],[0.95 0.00 0.95],[0.00 0.55 0.65],[0.05 0.05 0.05]};
pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 21 14]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
for i = 1:6
    nexttile(i)
    if i <= 3
        ax{i}                       = gca;
        ax{i}.FontName              = 'Times New Roman';
        ax{i}.FontSize              = 10;
        ax{i}.XScale                = 'log';
        ax{i}.XAxis.TickLabelFormat = '%.1f';
        ax{i}.YAxis.TickLabelFormat = '%.2f';
        ax{i}.XTick                 = .1./(2.^(9:-1:-2))*1000;
        ax{i}.XAxis.MinorTickValues = 10:10:200;
        xlim([6.25 100])
        xlabel('$\mathit{deaths}$ $\mathit{per}$ $\mathit{1000}$ $\mathit{births}$ (log scale)','Interpreter','latex','FontName','Times New Roman','FontSize',11);
        ylabel('$\mathit{kernel}$ $\mathit{density}$','Interpreter','latex','FontName','Times New Roman','FontSize',11);
        title(char(string(char(64 + i)) + ". " + string(LAB{i})),'Interpreter','latex');
    else
        d                           = scatter(datetime('01-Jan-2021'),0);
        delete(d)
        ax{i}                       = gca;
        ax{i}.FontName              = 'Times New Roman';
        ax{i}.FontSize              = 10;
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
        xlabel('$\mathit{year}$','Interpreter','latex','FontName','Times New Roman','FontSize',11);
        ylabel('$\mathit{deaths}$ $\mathit{per}$ $\mathit{1000}$ $\mathit{births}$ (log scale)','Interpreter','latex','FontName','Times New Roman','FontSize',11);
        xlim([datetime(2010,1,1) datetime(2024,1,1)]);
        title(char(string(char(64 + i)) + ". " + string(LAB{i - 3})),'Interpreter','latex','FontName','Times New Roman','FontSize',12);
    end
    grid on;
    box on;
    hold on;
end

for i = 1:6
    nexttile(i)
    if i <= 3
        mf  = 4.0;
        ylim([0 mf])
        for j = 1:numel(selection)
            q{j} = TaBle.s{selection{j}(1),selection{j}(2)}(i,:);
            q{j} = recode(q{j},NaN,eps);
            if numel(q{j}) == R + 1
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
            P{end + 1} = plot(max(xi{j},eps),f{j},'color',coloR{j},'LineWidth',1.00);
        end
        for j = 1:numel(selection)
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.75,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.50);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.50);
        end
        if i == 2
            for j = 1:numel(selection)
                leGend{j} = pOPs{selection{j}(1),selection{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    else
        for j = 1:numel(selectionT)
            tAU{j} = TaBle.tAU{selectionT{j}(1),selectionT{j}(2)};
            if numel(tAU{j}) > 2
                q{j} = TaBle.s{selectionT{j}(1),selectionT{j}(2)}{1 + mod(i + 2,3)};
            else
                q{j} = kron(ones(2,1),TaBle.s{selectionT{j}(1),selectionT{j}(2)}(1 + mod(i + 2,3),:));
            end
            if size(q{j},2) == R + 1
                q{j} = prctile(q{j}(:,2:end)',tHtiles)'*1000;
            else
                q{j} = q{j}*1000;
            end
            P{end + 1} = fill([tAU{j}(1:numel(q{j}(:,2)));flip(tAU{j}(1:numel(q{j}(:,3))))],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j});
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j}(1:numel(q{j}(:,1))),q{j}(:,1),'color',coloR{j},'LineWidth',1.0,'LineStyle',':');
        end
        P{end + 1} = plot(datetime(ones(2,1)*[2020 1 1]),[3.125 200],'color',[0.1 0.4 0.5],'LineStyle','--','LineWidth',0.75);
        if i == 5
            for j = 1:numel(selectionT)
                leGend{j} = pOPs{selectionT{j}(1),selectionT{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    end
end
saveas(gcf,char(pATh + "Results/RaMMPS-perinatal Figure 1.png"));
for i = 1:numel(P)
    delete(P{i})
end
clear leGend P


P                        = cell(0);
selection                = {[2 3],[3 3],[4 3],[5 3],[6 3],[1 8],[1 15]};
selectionT               = {[2 3],[3 3],[4 3],[5 3],[6 3],[1 8],[5 15]};
coloR                    = {[0.00 0.00 0.75],[0.95 0.00 0.95],[0.85 0.35 0.01],[0.45 0.65 0.20],[0.65 0.10 0.20],[0.00 0.55 0.65],[0.05 0.05 0.05]};
for i = 1:6
    nexttile(i)
    if i <= 3
        mf  = 2.0;
        ylim([0 mf])
        for j = 1:numel(selection)
            q{j} = TaBle.s{selection{j}(1),selection{j}(2)}(i,:);
            q{j} = recode(q{j},NaN,eps);
            if numel(q{j}) == R + 1
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
            P{end + 1} = plot(max(xi{j},eps),f{j},'color',coloR{j},'LineWidth',1.00);
        end
        for j = 1:numel(selection)
            P{end + 1} = fill(max(xi{j},eps),f{j},coloR{j},'FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',coloR{j});
            P{end + 1} = plot(ones(2,1)*max(pc{j}(1),eps),[0 g{j}(1)],'color',coloR{j},'LineWidth',1.75,'LineStyle',':');
            P{end + 1} = plot(ones(2,1)*max(pc{j}(2),eps),[0 g{j}(2)],'color',coloR{j},'LineWidth',0.50);
            P{end + 1} = plot(ones(2,1)*max(pc{j}(3),eps),[0 g{j}(3)],'color',coloR{j},'LineWidth',0.50);
        end
        if i == 2
            for j = 1:numel(selection)
                leGend{j} = pOPs{selection{j}(1),selection{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    else
        for j = 1:numel(selectionT)
            tAU{j} = TaBle.tAU{selectionT{j}(1),selectionT{j}(2)};
            if numel(tAU{j}) > 2
                q{j} = TaBle.s{selectionT{j}(1),selectionT{j}(2)}{1 + mod(i + 2,3)};
            else
                q{j} = kron(ones(2,1),TaBle.s{selectionT{j}(1),selectionT{j}(2)}(1 + mod(i + 2,3),:));
            end
            if size(q{j},2) == R + 1
                q{j} = prctile(q{j}(:,2:end)',tHtiles)'*1000;
            else
                q{j} = q{j}*1000;
            end
            P{end + 1} = fill([tAU{j}(1:numel(q{j}(:,2)));flip(tAU{j}(1:numel(q{j}(:,3))))],max([q{j}(:,2);flip(q{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j});
        end
        
        for j = 1:numel(selectionT)
            P{end + 1} = plot(tAU{j}(1:numel(q{j}(:,1))),q{j}(:,1),'color',coloR{j},'LineWidth',1.0,'LineStyle',':');
        end
        P{end + 1} = plot(datetime(ones(2,1)*[2020 1 1]),[3.125/2 200],'color',[0.1 0.4 0.5],'LineStyle','--','LineWidth',0.75);
        if i == 5
            for j = 1:numel(selectionT)
                leGend{j} = pOPs{selectionT{j}(1),selectionT{j}(2)}{1};
            end
            P{end + 1} = legend(leGend,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
        end
    end
end
saveas(gcf,char(pATh + "Results/RaMMPS-perinatal Figure 2.png"));
for i = 1:numel(P)
    delete(P{i})
end
clear leGend P



tAU                      = (2006:.25:2018)';
LAB                      = {'Infant mortality rate $\mathit{q}\mathrm{(12}\mathit{m}\mathrm{)}$','Under-five mortality rate $\mathit{q}\mathrm{(60}\mathit{m}\mathrm{)}$'};
coloR                    = {[0.00 0.00 0.75],[0.95 0.00 0.95],[0.85 0.35 0.01],[0.45 0.65 0.20],[0.65 0.10 0.20],[0.00 0.55 0.65],[0.05 0.05 0.05]};
pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 14 14]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
for i = 1:4
    nexttile(i)
    d                           = scatter(datetime('01-Jan-2021'),0);
    delete(d)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 10;
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    ax{i}.YMinorGrid            = 'on';
    ax{i}.XMinorGrid            = 'on';
    ax{i}.XAxis.TickLabelFormat = 'yyyy';
    ax{i}.XTick                 = datetime(2006:2:2020,1,1);
    ax{i}.XAxis.MinorTickValues = datetime(2006:2020,1,1);
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    if 1 + mod(i + 1,2) == 1
        ax{i}.YTick                 = 0:20:200;
        ax{i}.YAxis.MinorTickValues = 0:20:200;
        %ylim([0 100])
        ylabel('$\mathit{deaths}$ $\mathit{per}$ $\mathit{1000}$ $\mathit{births}$','Interpreter','latex','FontName','Times New Roman','FontSize',11);    
    else
        ax{i}.YScale                = 'log';
        ax{i}.YTick                 = .1./(2.^(9:-1:-2))*1000;
        ax{i}.YAxis.MinorTickValues = 10:10:100;
        ylim([25 100])
        ylabel('$\mathit{deaths}$ $\mathit{per}$ $\mathit{1000}$ $\mathit{births}$ (log scale)','Interpreter','latex','FontName','Times New Roman','FontSize',11);
    end
    xlabel('$\mathit{year}$','Interpreter','latex','FontName','Times New Roman','FontSize',11);
    xlim([datetime(2006,1,1) datetime(2020,1,1)]);
    title(char(string(char(96 + i)) + ". " + string(LAB{ceil(i/2)})),'Interpreter','latex','FontName','Times New Roman','FontSize',12);
    grid on;
    box on;
    hold on;
end

for i = 1:4
    nexttile(i)
    h = ceil(i/2);
    if 1 + mod(i + 1,2) == 1
        for j = 1:numel(CD)
            plot(prctile(BraSs_r{j,1}{3}(:,2:end)',50)',prctile(BraSs_r{j,1}{h}(:,2:end)',50)'*1000,'color',coloR{j},'LineWidth',1.00);
        end
        for j = 1:numel(CD)
            plot(BraSs_r{j,1}{3}(:,2:50),BraSs_r{j,1}{h}(:,2:50)*1000,'color',[coloR{j} 0.1],'LineWidth',0.75);
        end
        xlim([min(tAU) 2020])
    else
        for j = 1:numel(CD)
            temp{j} = max(LinInterPol(BraSs_r{j,1}{3},BraSs_r{j,1}{h},tAU),0);
            temp{j} = prctile(temp{j}(:,2:end)',tHtiles)'*1000;
            fill([tAU;flip(tAU)],max([temp{j}(:,2);flip(temp{j}(:,3))],eps),coloR{j},'FaceAlpha',.05,'EdgeAlpha',1.00,'LineStyle','-','EdgeColor',coloR{j});
        end
        for j = 1:numel(CD)
            plot(tAU,temp{j}(:,1),'color',coloR{j},'lineWidth',1.00,'LineStyle',':');
        end
        if h == 1
            ylim([25 150])
        else
            ylim([25 200])
        end
        xlim([min(tAU) 2020])
    end
    if h == 1
        legend(CD,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
    end
end
saveas(gcf,char(pATh + "Results/RaMMPS-U5M Figure A1.png"));
clear leGend tAU




clear
pATh         = "/Users/lshjr3/Documents/RaMMPS/under-five mortality/";
load(char(pATh + "Results/RaMMPS.mat"),'RaMMPS');
%sET         = {RaMMPS.TPH_duration(RaMMPS.k == 1),RaMMPS.FPH_duration(RaMMPS.k == 1),RaMMPS.HHD_duration(RaMMPS.k == 1),RaMMPS.SVY_duration(RaMMPS.k == 1)};
sET          = {RaMMPS.TPH_duration(RaMMPS.k == 1),RaMMPS.FPH_duration(RaMMPS.k == 1)};
LAB          = {'$\textrm{a. Truncated Pregnancy History}$','$\textrm{b. Full Pregnancy History}$','$\textrm{c. Household Questionnaire}$','$\textrm{d. All RaMMPS}$'};

coloR                    = {[0.00 0.00 0.75],[0.95 0.00 0.95],[0.85 0.35 0.01],[0.45 0.65 0.20],[0.65 0.10 0.20],[0.00 0.55 0.65],[0.05 0.05 0.05]};
pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 14 7]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
for i = 1:2
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 10;
    ax{i}.XAxis.TickLabelFormat = '%.1f';
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    ax{i}.XTick                 = 0:2.50:40;
    ax{i}.XAxis.MinorTickValues = 0:1.25:40;
    ax{i}.YTick                 = 0:0.10:1;
    ax{i}.YAxis.MinorTickValues = 0:0.05:1;


    xlabel('$\textit{duration (in minutes)}$','Interpreter','latex','FontName','Times New Roman','FontSize',11);
    ylabel('$\textit{probability density function}$','Interpreter','latex','FontName','Times New Roman','FontSize',11);
    grid on;
    grid minor;
    box on;
    hold on;
    
    xi     = 0:0.25:40;
    y      = kron(sET{i}(~isnan(sET{i}))/60,ones(1,1));
    histogram(y,xi,'Normalization','pdf','FaceColor',coloR{i},'EdgeColor',[0 0 0],'FaceAlpha',0.1,'EdgeAlpha',0.25);
    title(char(LAB{i} + ", N = " + string(numel(y))),'Interpreter','latex','FontName','Times New Roman','FontSize',12);
    f      = ksdensity(y,xi,'Function','pdf','Bandwidth',.1);
    plot(xi,f,'LineWidth',0.50,'Color',[coloR{i} 0.5]);
    z      = mean(y);
    round(z,2)
    Z(i,:) = prctile(y,[50 25 75]);
    [f,z]  = ksdensity(y,z,'Function','pdf','Bandwidth',.1);
    plot(z(1)*ones(1,2),[0 f(1)],'LineWidth',1.0,'Color',coloR{i});
    xlim([-.5 15])
    ylim([0 .80]);
end
saveas(gcf,char(pATh + "Results/RaMMPS-duration.png"));
Z(:,4)       = Z(:,3) - Z(:,2);

lABs         = {{1} {2}};
vARs         = {{'median','Q1','Q3','IQR'}};
foRMaT       = {'%0.2f','%0.2f','%0.2f'};
nOTe         = {'$\textit{Survey Instrument}$',''};
tABleBAyEs({'$\textrm{Duration (in minutes)}$'},vARs,foRMaT,lABs,nOTe,{{'TPH',''};{'FPH',''}},kron(Z,[1 NaN NaN]),0.120,0.055,[]);
saveas(gcf,char(pATh + "Results/RaMMPS-duration Table.png"));



clear
pATh         = "/Users/lshjr3/Documents/RaMMPS/under-five mortality/";
load(char(pATh + "Results/RaMMPSf.mat"),'RaMMPS');
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
saveas(gcf,char(pATh + "Results/RaMMPS-U5M Table 0.png"));


clear
pATh         = "/Users/lshjr3/Documents/RaMMPS/under-five mortality/";
load(char(pATh + "Results/RaMMPS.mat"),'RaMMPS','DHSmalawi','MICSmalawi');
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
saveas(gcf,char(pATh + "Results/RaMMPS-U5M Table 1.png"));

lABs         = {{1} {3 4 5} {6 7 8} {9 10 11} {12 13 14} {15} {21}};
vARs         = {models models models modelsDHS([1 3])};
nOTe         = {'$\textrm{Attributes}$/$\textrm{Instrument\&}\textit{method}$',''};
tABleBAyEs(sEt([2 3 4]),vARs(2:end),foRMaT,lABs,nOTe,pOPsd,cell2mat(bOx(:,[3 4 5 6 7 9])),0.190,0.065,[]);
saveas(gcf,char(pATh + "Results/RaMMPS-perinatal Table 1.png"));
s            = [1 3 4 5 6 7 8 9 10 11 12 13 14 15 21];
bOx          = bOx(s,[3 4 5 6 7 9]);
pOPsd        = pOPsd(s);
sets         = sEt([2 3 5]);
vars         = vARs(2:end);
save(char(pATh + "Results/Tables/perinatal_Table_1.mat"),'bOx','pOPsd','vars','sets');
clear bOx pACk lABelS lABelSd pOPsD pOPsS vars sets s



load(char(pATh + "Results/RaMMPS.mat"),'RaMMPSdst');
list         = {'$\textrm{RaMMPS, all 18-64}$'};
models       = {'$\textit{unweighted}$','$\textit{weighted}$'};

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

sET                   = [ones(size(RaMMPSdst.WT)),RaMMPSdst.WT];
R                     = size(dATa{1}{2},2) - 1;
bLOKs                 = {1,2,3,4};
data                  = RaMMPSdst;
for i = 1:numel(models)
    rng(0);
    S          = [ones(size(data,1),1),rand(size(data,1),R)];
    for j = 1:numel(bLOKs)
        sEL             = find(data.group == bLOKs{j});
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
    dATadst{i} = {data.group  > 0,S};
    clear S
end

load(char(pATh + "Results/RaMMPS.mat"),'DHSmalawidst');
listDHS               = {'$\textrm{DHS VII, all 18-64}$'};
modelsDHS             = {'$\textit{all}$','$\textit{m.\,owner}$','$\textit{re-weighted}$'};
dATaDHSdst{1}         = {DHSmalawidst.age >= 18 & DHSmalawidst.age < 65 & DHSmalawidst.jure == 1,DHSmalawidst.Wh};
dATaDHSdst{2}         = {DHSmalawidst.age >= 18 & DHSmalawidst.age < 65 & DHSmalawidst.mobile == 1 & DHSmalawidst.jure == 1,DHSmalawidst.Wh};
dATaDHSdst{3}         = {dATaDHSdst{2}{1},DHSmalawidst.WT};

rng(0);
s                     = DHSmalawidst(DHSmalawidst.k == 1, {'cluster','k','K','h','H','iNDeX'});
w                     = rand(size(s,1),R);
w                     = [(1:size(s,1))',ceil(s.H.*w) + s.iNDeX];
Wdhs                  = NaN(size(DHSmalawidst,1),R + 1);
for r = 1:R + 1
    S                  = tabulate([w(:,1);w(:,r)]);
    Wdhs(:,r)          = repelem(S(:,2) - 1,s.K);
    clc;
    r/(R + 1)
end
clear r s S w

for i = 1:numel(dATaDHSdst)
    dATaDHSdst{i}{1,2} = Wdhs.*dATaDHSdst{i}{1,2};
end

for i = 1:numel(modelsDHS)
    mIn              = min(DHSmalawidst.interview);
    mAx              = max(DHSmalawidst.interview);
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


pACk         = {{RaMMPSdst,dATadst},{DHSmalawidst,dATaDHSdst}};
lABelS       = {'Place of residence','Region','Sex','Age','Education','Household size','Electricity','Drinking water','Roofing'};
lABelSd      = {{'urban' 'rural'} {'North' 'Central' 'South'} {'Female' 'Male'} {'18-29' '30-39' '40-49' '50-64'} {'less than complete primary' 'incomplete secondary' 'complete secondary or more'} {'1-4' '5-8' '9+'} {'access' 'no access'} {'safe source' 'other source'} {'durable material' 'other material'}}; 
outcomes     = {[1 2],[1 2 3],[2 1],[1 2 3 4],[1 2 3],[1 2 3],[2 1],[2 1],[2 1]};

H            = 0;
for h = 1:numel(pACk)
    d    = pACk{h}{1};
    data = [d.UR,d.Region,d.sex,d.GO,d.Education,d.household,d.Electricity,d.Water,d.Roofing];
    for i = 1:numel(pACk{h}{2})
        H  = H + 1;
        s  = pACk{h}{2}{i}{1};
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

selection    = 1 + [0 cumsum([numel(models)])];
for i = 1:numel(selection)
    sEt{i} = {pOPsS{selection(i)}};
end

lABs       = {{1} {3 4 5} {6} {8 9 10 11} {12 13 14} {15 16 17} {18} {24}};
vARs       = {models modelsDHS};
foRMaT     = {'%0.2f','%0.2f','%0.2f'};

nOTe       = {'$\textrm{Attributes}$/$\textit{method}$',''};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,pOPsd,cell2mat(bOx(:,[1 2 3 4 5])),0.190,0.065,[]);
saveas(gcf,char(pATh + "Results/RaMMPS-Data Table 3.png"));
s            = [1 3 4 5 6 8 9 10 11 12 13 14 15 16 17 18 24];
bOx          = bOx(s,[1 2 3 4 5]);
pOPsd        = pOPsd(s);
sets         = sEt([1 2 3 4 5]);
vars         = vARs;
save(char(pATh + "Results/Tables/data_Table_3.mat"),'bOx','pOPsd','vars','sets');
clear bOx  pOPsd vars sets s



clear
pATh             = "/Users/lshjr3/Documents/RaMMPS/under-five mortality/";
load(char(pATh + "Results/RaMMPS.mat"),'RaMMPS');
load(char(pATh + "Results/ReSAmPLiNG.mat"),'dATa','Ts');
instrument       = {'TPH','FPH'};
models           = {'weighted','unweighted'};
Ts               = Ts(1);

for i = 1:numel(instrument)
    for j = 1:size(models,2)
        h          = j + (i - 1)*size(models,2);
        for k = 1:numel(Ts)
            date      = eXAcTTime(Ts{k});
            date      = string(instrument{i}) + ", " + string(sprintf('%0.1f',date(1))) + "-" + string(sprintf('%0.1f',date(2)));
            pOPs{k,h} = {char(date);''};
            clear date
        end
        vARs{h}    = "$\mathit{" + string(models{j}) + "}$";
    end
    sEt{i} = "$\mathrm{" + string(instrument{i}) + "}$";
end

rng(0);
R                = size(dATa{1}{2},2) - 1;
p                = rand(size(RaMMPS,1),R + 1);
Event            = p.*datenum(RaMMPS.B_min) + (1 - p).*datenum(RaMMPS.B_max);
Event            = reshape(datetime(datevec(reshape(Event,numel(Event),1)),'Format','dd/MM/yyyy'),size(RaMMPS,1),R + 1);

sEL              = (Event >= datetime(year(Ts{1}(1)),1,1));
Gestation        = (RaMMPS.birth == 'stillbirth') & (RaMMPS.gestation >= 28);
NGestation       = (RaMMPS.birth == 'stillbirth') & (RaMMPS.gestation  < 28);
ClearBI          = (RaMMPS.birth == 'stillbirth') & (RaMMPS.FlagBI ~= 1);
NClearBI         = (RaMMPS.birth == 'stillbirth') & (RaMMPS.FlagBI == 1);

fRoMA            = char("$\textrm{reporting events from 01/01/" + year(Ts{1}(1)) + " onward}$");
fRoMB            = char("$\textrm{from 01/01/" + year(Ts{1}(1)) + "}$");
for i = 1:numel(dATa)
    s                = dATa{i}{1};
    wS               = dATa{i}{2}(s,:);    
    
    TaBle.WT{1,i}    = sum((RaMMPS.k(s) == 1).*wS,1);
    TaBle.W{1,i}     = sum((RaMMPS.k(s) == RaMMPS.K(s)).*sEL(s,:).*wS,1);
    TaBle.B{1,i}     = sum((RaMMPS.birth(s) == 'livebirth' & RaMMPS.flagPQ(s) ~= 're-classified').*sEL(s,:).*wS,1);
    TaBle.F{1,i}     = sum((RaMMPS.birth(s) == 'stillbirth' | RaMMPS.flagPQ(s) == 're-classified').*sEL(s,:).*wS,1);
    TaBle.R{1,i}     = sum((RaMMPS.flagPQ(s) == 're-classified').*sEL(s,:).*wS,1);
    TaBle.Miss{1,i}  = sum((RaMMPS.birth(s) == 'stillbirth').*isnan(RaMMPS.gestation(s)).*sEL(s,:).*wS,1);
    TaBle.M{1,i}     = sum(NGestation(s).*sEL(s,:).*wS,1);
    TaBle.StilN{1,i} = sum(Gestation(s).*sEL(s,:).*wS,1);
    TaBle.Over{1,i}  = sum(NClearBI(s).*Gestation(s).*sEL(s,:).*wS,1);    
    TaBle.Still{1,i} = sum(ClearBI(s).*Gestation(s).*sEL(s,:).*wS,1);   
    TaBle.SR{1,i}    = TaBle.Still{1,i}./(TaBle.B{1,i} + TaBle.R{1,i} + TaBle.Still{1,i})*1000;
    clear wS s
end

nAMeS{1,1}   = {'$\textrm{Women}$';'$\textrm{total number}$'};
nAMeS{2,1}   = {'$\textrm{Women}$';fRoMA};
nAMeS{3,1}   = {'$\textrm{Livebirths}$';fRoMB};
nAMeS{4,1}   = {'$\textrm{Fetal deaths}$';fRoMB};
nAMeS{5,1}   = {'$\textrm{(\textendash) Actual livebirths: breathing, moving, or crying}$';fRoMB};
nAMeS{6,1}   = {'$\textrm{(\textendash) Fetal deaths with missing gestation}$';fRoMB};
nAMeS{7,1}   = {'$\textrm{(\textendash) Fetal deaths with less than 28 weeks (miscarriages)}$';fRoMB};
nAMeS{8,1}   = {'$\textrm{(=) Stillbirths (28 weeks or more)}$';'$\textrm{from 01/01/2014}$'};
nAMeS{9,1}   = {'$\textrm{(\textendash) Stillbirths (overlapping a birth interval)}$';fRoMB};
nAMeS{10,1}  = {'$\textrm{(=) Stillbirths (28 weeks or more \& clear birth interval)}$';fRoMB};
nAMeS{11,1}  = {'$\textrm{Stillbirth Rate (per 1,000 pregnancies)}$';fRoMB};

for i = 1:numel(dATa)
    bOx{1,i}  = [mean(TaBle.WT{1,i}(2:end)) prctile(TaBle.WT{1,i}(2:end),[ 2.5 97.5])];
    bOx{2,i}  = [mean(TaBle.W{1,i}(2:end)) prctile(TaBle.W{1,i}(2:end),[ 2.5 97.5])];
    bOx{3,i}  = [mean(TaBle.B{1,i}(2:end)) prctile(TaBle.B{1,i}(2:end),[ 2.5 97.5])];
    bOx{4,i}  = [mean(TaBle.F{1,i}(2:end)) prctile(TaBle.F{1,i}(2:end),[ 2.5 97.5])];
    bOx{5,i}  = [mean(TaBle.R{1,i}(2:end)) prctile(TaBle.R{1,i}(2:end),[ 2.5 97.5])];
    bOx{6,i}  = [mean(TaBle.Miss{1,i}(2:end)) prctile(TaBle.Miss{1,i}(2:end),[ 2.5 97.5])];
    bOx{7,i}  = [mean(TaBle.M{1,i}(2:end)) prctile(TaBle.M{1,i}(2:end),[ 2.5 97.5])];
    bOx{8,i}  = [mean(TaBle.StilN{1,i}(2:end)) prctile(TaBle.StilN{1,i}(2:end),[ 2.5 97.5])];
    bOx{9,i}  = [mean(TaBle.Over{1,i}(2:end)) prctile(TaBle.Over{1,i}(2:end),[ 2.5 97.5])];    
    bOx{10,i} = [mean(TaBle.Still{1,i}(2:end)) prctile(TaBle.Still{1,i}(2:end),[ 2.5 97.5])];
    bOx{11,i} = [prctile(TaBle.SR{1,i}(2:end),[50 2.5 97.5])];
end
    
lABs         = {{1 2 3},{4 5 6 7 8},{9 10},{11}};
vARs         = {vARs(1:2) vARs(1:2)};
foRMaT       = {'%0.2f','%0.2f','%0.2f','%0.2f'};
nOTe         = {'$\textrm{Selection Group}$/$\textrm{Instrument\&}\textit{method}$','$\textit{}$'};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,nAMeS,cell2mat(bOx),0.250,0.070,[]);
saveas(gcf,char(pATh + "Results/RaMMPS-Perinatal Table S2.png"));
s            = [1 2 3 4 5 6 7 8 9 10 11];
bOx          = bOx(s,:);
pOPsd        = nAMeS(s);
sets         = sEt;
vars         = vARs;
save(char(pATh + "Results/Tables/perinatal_Table_S2.mat"),'bOx','pOPsd','vars','sets');
clear bOx  pOPsd vars sets s








clear
pATh                  = "/Users/lshjr3/Documents/RaMMPS/under-five mortality/";
load(char(pATh + "Results/RaMMPS.mat"),'RaMMPScalls');
R                     = 0;
bLOKs                 = {1,2,3,4};
data                  = RaMMPScalls(RaMMPScalls.k == RaMMPScalls.K,{'group','caseid','K','catioutcome','catiOUTcome'});

rng(0);
WS                    = [ones(size(data,1),1),rand(size(data,1),R)];
for j = 1:numel(bLOKs)
    sEL             = find(data.group == bLOKs{j});
    W               = ones(size(sEL));
    W               = [0;cumsum(W(1:end - 1))]/sum(W);
    for r = 2:R + 1
        temp      = tabulate([sum(W < WS(sEL,r)',1)';numel(sEL) + 1]);
        WS(sEL,r) = temp(1:end - 1,2);
        clear temp
        clc;
        r/(R + 1)
    end
clear W sEL    
end
W                     = WS(repelem((1:size(data,1))',data.K),:);

data.complete         = (data.catiOUTcome ==  1);
data.refusal          = (data.catiOUTcome ==  3);
data.eligible         = (data.catiOUTcome ==  1 | data.catiOUTcome ==  2 | data.catiOUTcome ==  3 | data.catiOUTcome == 4 | data.catiOUTcome == 7);

CATI                  = data.complete'*WS;
cases                 = size(data,1);
calls                 = size(RaMMPScalls,1);
response              = (data.complete'*WS)./(data.eligible'*WS)*100;
refusal               = (data.refusal'*WS)./(data.eligible'*WS)*100;
callperCATI           = calls./CATI;
numbersCATI           = cases./CATI;

nAMeS{1,1}            = {'$\textrm{Complete CATI}$';'$\textrm{}$'};
nAMeS{2,1}            = {'$\textrm{Cases}$';'$\textrm{}$'};
nAMeS{3,1}            = {'$\textrm{Calls placed}$';'$\textrm{}$'};
nAMeS{4,1}            = {'$\textrm{Response rate (\%)}$';'$\textit{excludes ineligible and numbers not in use}$'};
nAMeS{5,1}            = {'$\textrm{Refusal rate (\%)}$';'$\textit{excludes ineligible and numbers not in use}$'};
nAMeS{6,1}            = {'$\textrm{Calls per complete CATI}$';'$\textrm{}$'};
nAMeS{7,1}            = {'$\textrm{Cases per complete CATI}$';'$\textrm{}$'};

bOx{1,1}              = prctile(CATI,[50 2.5 97.5]);
bOx{2,1}              = [cases,NaN(1,2)];
bOx{3,1}              = [calls,NaN(1,2)];
bOx{4,1}              = prctile(response,[50 2.5 97.5]);
bOx{5,1}              = prctile(refusal,[50 2.5 97.5]);
bOx{6,1}              = prctile(callperCATI,[50 2.5 97.5]);
bOx{7,1}              = prctile(numbersCATI,[50 2.5 97.5]);

sEt                   = {'$\textrm{}$'};
lABs                  = {{1 2 3} {4 5} {6 7}};
models                = {'$\textit{}$'};
vARs                  = {models};
foRMaT                = {'%0.2f'};
nOTe                  = {'$\textrm{Summary}$',''};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,nAMeS,cell2mat(bOx),0.250,0.070,[]);
saveas(gcf,char(pATh + "Results/RaMMPS-Data Table 1.png"));
s                     = [1 2 3 4 5 6 7];
bOx                   = bOx(s,:);
pOPsd                 = nAMeS(s);
sets                  = sEt;
vars                  = vARs;
save(char(pATh + "Results/Tables/data_Table_1.mat"),'bOx','pOPsd','vars','sets');
clear nAMeS bOx pOPsd vars sets s


TAB                   = tabulate(data.catiOUTcome);
TAB                   = TAB(:,1);
for i = 1:numel(TAB)
    sEL        = (data.catiOUTcome == TAB(i));    
    bOx{i,1}   = [sum(sEL),NaN(1,2)];
    bOx{i,2}   = prctile((sEL'*WS)./(sum(WS,1))*100,[50 2.5 97.5]);
    temp       = data.catioutcome(sEL);
    temp       = char(temp(1));
    j          = find(temp == '(');
    if numel(j) > 0
        nAMeS{i,1} = {char("$\textrm{" + temp(1:j - 2) + "}$");char("$\textit{" + temp(j:end) + "}$")};
    else
        nAMeS{i,1} = {char("$\textrm{" + temp + "}$");'$\textit{}$'};
    end
    clear sEL temp j
end

bOx{i + 1,1}        = [size(data,1),NaN(1,2)];
bOx{i + 1,2}        = [100,NaN(1,2)];
nAMeS{i + 1,1}      = {'$\textrm{Total cases}$';'$\textrm{}$'};

lABs                  = {{1 2 3 4},{5 6},{7},{8}};
models                = {'$\textit{Number}$','$\textit{\%}$'};
vARs                  = {models};
foRMaT                = {'%0.0f','%0.2f'};
nOTe                  = {'$\textrm{Call Outcome}\textit{ (final dispositions)}$',''};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,nAMeS,cell2mat(bOx),0.250,0.070,[]);
saveas(gcf,char(pATh + "Results/RaMMPS-Data Table 2.png"));
s                     = [1 2 3 4 5 6 7 8];
bOx                   = bOx(s,:);
pOPsd                 = nAMeS(s);
sets                  = sEt;
vars                  = vARs;
save(char(pATh + "Results/Tables/data_Table_2.mat"),'bOx','pOPsd','vars','sets');
clear nAMeS bOx pOPsd vars sets s






clear
pATh       = "/Users/lshjr3/Documents/RaMMPS/under-five mortality/";
load(char(pATh + "Results/RaMMPS.mat"),'RaMMPS');
x{1}       = [(0:7:28)/365.25,[(2:1:12),(15:3:24),(36:12:60)]/12]';
x{2}       = [string(0);string([(7:7:28)';(2:1:11)';(12:3:24)';(36:12:60)']) + char([kron('d',ones(4,1));kron('m',ones(18,1))])];
n          = x{1}(2:end) - x{1}(1:end - 1);
list       = {'TPH','FPH'};
RaMMPS.W   = RaMMPS.WR;
RaMMPS     = RaMMPS(RaMMPS.birth == 'livebirth' & (RaMMPS.survival == 'dead' | RaMMPS.survival == 'alive'),:);
RaMMPS     = RaMMPS(~isnat(RaMMPS.B_min),:);
R          = 1000;
date       = max(RaMMPS.interview);
Ts         = {datetime([2014 year(date)]',[1 month(date)]',[1 day(date)]'),datetime([2014 2017]',1,1),datetime([2017 2020]',1,1),datetime([2020 year(date)]',[1 month(date)]',[1 day(date)]')};
leGend     = {'Date of death','Date of interview','$\mathit{q}\mathrm{(}\mathit{x}\mathrm{)}$','Bootstrapping [p2.5, p97.5]'};

for i = 1:numel(list)
    alpha   = datetime((datenum(RaMMPS.B_min) + datenum(RaMMPS.B_max))/2,'ConvertFrom','datenum');
    omega   = datetime(min(datenum(alpha) + (datenum(RaMMPS.D_min) + datenum(RaMMPS.D_max))/2,datenum(RaMMPS.interview)),'ConvertFrom','datenum');
    dATa{i} = RaMMPS(RaMMPS.sample == list{i} & alpha >= datetime(2009,1,1) & omega >= datetime(2014,1,1) & alpha < RaMMPS.interview,:);
end
for i = 1:numel(dATa)
    s       = dATa{i};
    k       = 1;
    K       = 1;
    for j = 2:numel(s.caseid)
        if isequal(s.caseid(j),s.caseid(j - 1))
            k(end,1)     = k(end) + 1;
        else
            k(end + 1,1) = 1;
            K(end + 1,1) = j;
        end
    end
    rng(0);
    W       = s.W(K);
    W       = [0;cumsum(W(1:end - 1))]/sum(W);
    S       = [ones(numel(K),1),rand(numel(K),R)];
    w       = ones(size(s.W));
    for r = 2:R + 1
        temp   = tabulate([sum(W < S(:,r)',1)';numel(K) + 1]);
        S(:,r) = temp(1:end - 1,2);
        w(:,r) = repelem(S(:,r),k);
        clear temp
        clc;
        r/(R + 1)
    end
    p       = rand(size(s,1),R + 1);
    d       = p.*datenum(s.D_min) + (1 - p).*datenum(s.D_max);
    p       = rand(size(s,1),R + 1);
    b       = p.*datenum(s.B_min) + (1 - p).*datenum(s.B_max);
    B       = reshape(datetime(datevec(reshape(b,numel(b),1))),size(s,1),R + 1);
    D       = B + d;
    O       = reshape(datetime(datevec(reshape(min(datenum(D),datenum(s.interview)),numel(b),1))),size(s,1),R + 1);
    for h = 1:numel(Ts)
        T           = Ts{h};
        for j = 1:numel(x{1}) - 1
            a             = max(min(B + x{1}(j)*365.25,min(O,T(2))),T(1));
            o             = max(min(B + x{1}(j + 1)*365.25,min(O,T(2))),T(1));
            exposure(j,:) = sum((datenum(o) - datenum(a)).*w/365.25);
            events(j,:)   = sum((d/365.25 >= x{1}(j) & d/365.25 < x{1}(j + 1) & D >= T(1) & D < T(2)).*w);
            clear a o
        end
        births      = sum(k);
        mothers     = numel(k);
        m           = events./exposure;
        q           = 1 - [ones(1,R + 1);exp(-cumsum(m.*n))];
        OuT{h}{1,i} = q;
        OuT{h}{2,i} = m;
        OuT{h}{3,i} = mothers;
        OuT{h}{4,i} = births;
        OuT{h}{5,i} = events;
        OuT{h}{6,i} = exposure;
        clear T m q
    end
    clear B D O events exposure births mothers b d
end
for i = 1:numel(dATa)
    q        = OuT{1}{1,i};
    m        = OuT{1}{2,i};
    events   = OuT{1}{5,i};
    exposure = OuT{1}{6,i};
    d        = q(2:end,:) - q(1:end - 1,:);
    L        = min(d./m,(1 - q(1:end - 1,:)).*n);
    a        = (L - (1 - q(2:end,:)).*n)./d;
    t        = NaN(1,R + 1);
    X        = {[x{1},NaN(size(x{1},1),2)],[[n;NaN],NaN(size(x{1},1),2)],[exposure;t],[events;t],[m;t],q,[d;t],[a;t],[L;t]};
    for j = 3:numel(X)
        X{j} = [X{j}(:,1),prctile(X{j}(:,3:end)',[2.5 97.5])'];
    end
    for j = 1:numel(x{2})
       leGendA{j,1} = {char(x{2}(j));''};
    end
    temp     = max(dATa{i}.interview);
    temp     = round(year(temp) + (day(temp,'dayofyear') - 1)./days(datetime(year(temp) + 1,1,1) - datetime(year(temp),1,1)),1);
    sEt      = {char("RaMMPS Under-five mortality, Malawi 2014.0-" + string(temp) + " (" + string(list{i}) + ": " + string(OuT{1}{3,i}(1,1)) + " mothers and " + string(OuT{1}{4,i}(1,1)) + " births)")};
    vARs     = {'$\mathit{x}$ (in years)','$\mathit{n}$','$\mathit{exposure}$','$\mathit{events}$','$\mathit{_n}\mathit{m}\mathit{_x}$','$\mathit{q}\mathrm{(}\mathit{x}\mathrm{)}$','$\mathit{_n}\mathit{d}\mathit{_x}$','$\mathit{_n}\mathit{a}\mathit{_x}$','$\mathit{_n}\mathit{L}\mathit{_x}$'};
    foRMaT   = {'%0.4f','%0.4f','%0.1f','%0.1f','%0.4f','%0.4f','%0.4f','%0.4f','%0.4f'};
    lABs     = {num2cell(1:23)};
    nOTe     = {'','$\mathit{Bootstrapping}$ $\mathrm{p50}$/$\mathit{[p2.5,p97.5]}$'};
    tABleLT(sEt,{vARs},foRMaT,lABs,nOTe,leGendA,cell2mat(X),[],0.0125,0.0750);
    saveas(gcf,char(pATh + "Results/RaMMPS LT " + string(list{i}) + ".png"));
    clear events exposure m q d L a l t X
end


pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 numel(dATa)*7 7]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(1,numel(dATa),'Padding','compact','TileSpacing','compact');
for i = 1:numel(dATa)
    nexttile(i)
    d                           = scatter(datetime('01-Jan-2021'),0);
    delete(d)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 10;
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    ax{i}.YMinorGrid            = 'on';
    ax{i}.XMinorGrid            = 'on';
    ax{i}.XAxis.TickLabelFormat = 'yyyy';
    ax{i}.XTick                 = datetime(2014:2:2024,1,1);
    ax{i}.XAxis.MinorTickValues = datetime(2014:2024,1,1);
    ax{i}.YTick                 = 0:1:5;
    ax{i}.YAxis.MinorTickValues = 0:.5:5;
    ylabel('$\mathit{age}$ (in years)','Interpreter','latex','FontName','Times New Roman','FontSize',11);
    xlabel('$\mathit{year}$','Interpreter','latex','FontName','Times New Roman','FontSize',11);
    xlim([datetime(2014,1,1) datetime(2024,1,1)]);
    ylim([0 8])
    temp                        = max(dATa{i}.interview);
    eXACt{i}                    = [datetime(2014,1,1) temp];
    temp                        = round(year(temp) + (day(temp,'dayofyear') - 1)./days(datetime(year(temp) + 1,1,1) - datetime(year(temp),1,1)),1);
    temp                        = char(string(list{i}) + ", " + string(sprintf('%0.1f',2014)) + "-" +string(sprintf('%0.1f',temp)));
    pOPs{i}                     = char("RaMMPS " + string(temp));
    temp                        = char(string(char(96 + i)) + ". " + string(temp));
    temp2                       = char(string(OuT{1}{3,i}(1,1)) + " mothers and " + string(OuT{1}{4,i}(1,1)) + " births");
    title({temp;temp2},'Interpreter','latex','FontName','Times New Roman','FontSize',10);
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
        legend(leGend(1:2),'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',2,'Box','off');
    end
end
clear B D I xi xd r
saveas(gcf,char(pATh + "Results/RaMMPS-U5M Figure 1.png"));



clear
pATh       = "/Users/lshjr3/Documents/RaMMPS/under-five mortality/";
load(char(pATh + "Results/RaMMPS.mat"),'RaMMPSHH');
R          = 1000;
mIn        = min(RaMMPSHH.interview);
mAx        = max(RaMMPSHH.interview);

T          = datetime(year(mIn),month(mIn) + 1,1);
while datetime(year(T(end,1)),month(T(end,1)) + 3,1) < mAx
    T(end + 1,1) = datetime(year(T(end,1)),month(T(end,1)) + 3,1);
end
T(:,2)     = datetime(year(T),month(T) + 3,1);
T(:,3)     = T(:,1) + (T(:,2) - T(:,1))/2;


data       = RaMMPSHH(RaMMPSHH.WT > 0,:);
K          = find(data.k == 1);
k          = [K(2:end);size(data,1) + 1] - K;

rng(0);
W          = data.WT(K);
W          = [0;cumsum(W(1:end - 1))]/sum(W);
S          = [ones(numel(K),1),rand(numel(K),R)];
w          = ones(size(data,1),1);
for r = 2:R + 1
    temp   = tabulate([sum(W < S(:,r)',1)';numel(K) + 1]);
    S(:,r) = temp(1:end - 1,2);
    w(:,r) = repelem(S(:,r),k);
    clear temp
    clc;
    r/(R + 1)
end

data.Omega = data.O;
dead       = (data.status == 'dead' | data.status == 'dead migrant');
HH         = (data.k == 1);
for i = 0:3
    if isequal(i,0)
        region = 1;
    else
        region = (data.Region == i);
    end
        
    for j = 1:size(T,1)
        set             = region.*(data.interview >= T(j,1) & data.interview < T(j,2));
        setD            = region.*(data.Omega >= T(j,1) & data.Omega < T(j,2));
        ReP{i + 1}(j,:) = (setD'*(dead.*w))./(set'*(HH.*w));
    end
end



freddy                   = datetime(2023,[2 3],[5 14]);
blocks                   = [mIn datetime([2022 2022 2023],[5 9 3],[25 13 20])];
LAB                      = {'$\mathit{Malawi}$','$\mathit{North}$','$\mathit{Central}$','$\mathit{South}$'};
coloR                    = {[0.00 0.00 0.75],[0.95 0.00 0.95],[0.85 0.35 0.01],[0.45 0.65 0.20],[0.65 0.10 0.20],[0.00 0.55 0.65],[0.05 0.05 0.05]};
pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 14 14]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
for i = 1:4
    nexttile(i)
    d                           = scatter(datetime('01-Jan-2021'),0);
    delete(d)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 10;
    ax{i}.YAxis.TickLabelFormat = '%.1f';
    ax{i}.YMinorGrid            = 'on';
    ax{i}.XMinorGrid            = 'off';

    xlabel('$\mathit{date}$','Interpreter','latex','FontName','Times New Roman','FontSize',11);
    ylabel('$\mathit{quarterly\,deaths\,per\,housuhold\,(\%)}$','Interpreter','latex','FontName','Times New Roman','FontSize',11);
    title(char(string(i) + ". " + string(LAB{i})),'Interpreter','latex','FontName','Times New Roman','FontSize',12);    
    xlim([mIn mAx]);
    
    grid on;
    box on;
    hold on;
end


for i = 1:4
    nexttile(i)
    CB  = prctile(ReP{i}',[50 2.5 97.5])'*100;
    plot(T(:,3),CB(:,1),'color','b','LineWidth',1.0);
    fill([T(:,3);flip(T(:,3))],[CB(:,2);flip(CB(:,3))],'b','FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor','b');
    plot(T(:,3),CB(:,2:3),':','color','b','LineWidth',0.50);
    
    ylim([0 max(CB(:,3)*1.25)]);
    ylim([0 10]);

    for j = 2021:2023
        season = datetime([j j + 1],[11 5],1);
        fill([season([1 1]) season([2 2])],[100 0 0 100],'k','FaceAlpha',.075,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor','k');
    end
    fill([freddy([1 1]) freddy([2 2])],[100 0 0 100],'r','FaceAlpha',.15,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor','r');
    
    for j = 1:numel(blocks)
        plot(blocks([j j]),[0 100],'-.','color','k','LineWidth',0.50);    
    end
end
saveas(gcf,char(pATh + "Results/HHd.png"));