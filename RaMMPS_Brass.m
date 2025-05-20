clear
load('RaMMPS.mat','RaMMPSHH');

R          = 10;
K          = find(RaMMPSHH.k == 1);
k          = [K(2:end);size(RaMMPSHH,1) + 1] - K;
data       = RaMMPSHH;
data.w     = data.w/sum(data.w);
rng(0);
W          = data.w(K);
W          = [0;cumsum(W(1:end - 1))]/sum(W);
S          = [ones(numel(K),1),rand(numel(K),R)];
w          = ones(size(data,1),1);
for r = 2:R + 1
    temp   = tabulate([(1:numel(K))';sum(W < S(:,r)',1)']);
    S(:,r) = temp(:,2) - 1;
    w(:,r) = repelem(S(:,r),k);
    clear temp
    clc;
    r/(R + 1)
end
mig        = find(data.status == "migrant");
birth      = find(data.agegroup == "A. birth");
P          = ones(sum(k),R + 1);
P(mig,:)   = rand(numel(mig),R + 1);
Q          = ones(sum(k),R + 1);
Q(birth,:) = rand(numel(birth),R + 1);
P          = P.*Q;
clear Q mig birth
miss       = find(isnan(data.A));
ex         = data.O - data.A;
ex(miss)   = max(ex);
ex         = P.*ex;
O          = data.O*ones(1,R + 1);
A          = O - ex;
U5         = (data.agegroup == "A. under5" | data.agegroup == "A. birth");
Sample     = (data.status ~= "dead migrant" & data.status ~= "migrant");
temp       = datetime(kron([2021 2022]',ones(12,1)),kron(ones(2,1),(1:12)'),ones(24,1));
Tw{1}{1,1} = temp(12:22);
Tw{1}{2,1} = Tw{1}{1,1}([1 end]);
for j = 1:size(Tw{1},1)
    Tw{1}{j,2}      = Tw{1}{j,1}(1:end - 1) + days(Tw{1}{j,1}(2:end) - Tw{1}{j,1}(1:end - 1))/2;
    Tw{2}{j,1}      = year(Tw{1}{j,1}) + (day(Tw{1}{j,1},'dayofyear') - 1)./days(datetime(year(Tw{1}{j,1}) + 1,1,1) - datetime(year(Tw{1}{j,1}),1,1));
    Tw{2}{j,2}      = (Tw{2}{j,1}(2:end) + Tw{2}{j,1}(1:end - 1))/2;
    for i = 2:numel(Tw{2}{j,1})
        exposure     = sparse(max(min(O,Tw{2}{j,1}(i)) - max(A,Tw{2}{j,1}(i - 1)),0)).*Sample;
        events       = sparse(O >= Tw{2}{j,1}(i - 1)).*sparse(O < Tw{2}{j,1}(i)).*sparse(data.status == "dead" | data.status == "dead migrant").*Sample;
        M(i - 1,:)   = full(sum(events.*w.*(1 - U5))./sum(exposure.*w.*(1 - U5)));
        U5D(i - 1,:) = full(sum(events.*w.*U5)./sum(exposure.*w.*U5));
        clear exposure events
    end
    mortalityR{j,1} = M;
    mortalityR{j,2} = U5D;
    clear M U5D
end

D              = [mortalityR{2,2}' mortalityR{2,1}']*1000;
clear A data ex i j k K miss mortalityR O P r R RaMMPSHH S temp Tw U5 w W ans

load('prospects.mat')
country        = "Malawi";
population     = population(population.Name == country,8:end);
mortality      = mortality(mortality.Name == country,8:end);
sET            = [];
for i = 1950:2021
    temp1 = reshape(table2array(population(population.year == i,2:end))',101,2);
    temp2 = table2array(mortality(mortality.year == i,3:end));
    sET   = [sET;[i*ones(101,1),(0:100)',temp1,temp2]];
    clear temp1 temp2
end

% Newton-Raphson %
sET            = sET(:,[1:4 9:10]);
n              = size(D,1);
set            = sET(sET(:,1) == 2021,2:end);
pop            = sum(set(:,2:3),2);
population     = set(:,2:3);
m              = set(:,end - 1);
q              = set(:,end);
A              = 1 - 1./q + 1./m;
X              = cumprod(1 - q(1:end - 1));
X              = 0.5*log((1 - X)./X);
B              = [zeros(n,1) ones(n,1)];
set            = (1:n);
delta          = 0.001;
for i = 1:10
    for j = 1:3
        S    = zeros(1,3);
        S(j) = delta;
        y    = (B(set,1) + S(1,1))' + X*(B(set,2).*exp(S(1,2)))';
        s    = [ones(1,numel(set));1./(1 + exp(2*y));zeros(1,numel(set))];
        s    = 1 - s(2:end,:)./s(1:end - 1,:);
        M    = s./(1 - (1 - A).*s);
        d    = pop.*M;
        d    = [sum(d(1:5,:),1)'/sum(pop(1:5,:),1) sum(d(6:end,:),1)'/sum(pop(6:end,:),1)]*1000;
        F{j} = log(D(set,:)) - log(d);
        clear S
        clear y s M d S
    end
    
    for k = 1:2
        F{k} = (F{k} - F{3})/delta;
    end
    for k = 1:numel(set)
        S(k,:) = (pinv([F{1}(k,:)' F{2}(k,:)'])*F{3}(k,:)')';
        clear temp
    end
    
    s        = [B(:,1) log(B(:,2))];
    B(set,1) = B(set,1) - S(:,1);
    B(set,2) = B(set,2).*exp(-S(:,2));
    s        = [B(:,1) log(B(:,2))] - s;
    s        = max(abs(s),[],2);
    set      = find(s > 10^-10);
    s
    clear s S F 
    if numel(set) == 0
        break
    end
 end
clear set m q i j k delta n pop 

% Figures %
n              = 100;
a              = (-1:3/(n - 1):2)';
b              = (0:2/(n - 1):2)';
B              = [[kron(a,ones(n,1)),kron(ones(n,1),b)];[0,1];B];
for i = 1950:2021
    set           = sET(sET(:,1) == i,2:end);
    m             = set(:,end - 1);
    q             = set(:,end);
    A             = 1 - 1./q + 1./m;
    X             = cumprod(1 - q(1:end - 1));
    X             = 0.5*log((1 - X)./X);
    y             = B(:,1)' + X*B(:,2)';
    X             = [ones(1,size(y,2));1./(1 + exp(2*y));zeros(1,size(y,2))];
    Q             = 1 - X(:,n^2 + 1:end);
    X             = 1 - X(2:end,:)./X(1:end - 1,:);    
    M             = X./(1 - (1 - A).*X);
    X             = sum(set(:,2:3),2);
    d             = X.*M;
    d             = [sum(d(1:5,:),1)'/sum(X(1:5,:),1) sum(d(6:end,:),1)'/sum(X(6:end,:),1)]*1000;
    for j = 1:size(d,2)
        CO{i - 1949,j} = reshape(d(1:n^2,j),n,n);
    end
    Z(i - 1949,:) = d(n^2 + 1,:);
    Mortality     = M(:,n^2 + 1:end);
    clear d M X y m q set    
end
B              = B(n^2 + 1:end,:);

pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 14 21]/pix;
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
sources                  = {'$\mathit{UN - WPP2021}$','$\mathit{RaMMPS - HH}$'};
coloR                    = {[0.00 0.00 0.75],[0.95 0.00 0.95],[0.85 0.35 0.01],[0.05 0.05 0.05],[0.00 0.55 0.65],[0.85 0.30 0.10]};
LAB                      = {'A. Population pyramid','B. Matching parameters','C. Crude death rate - $\mathit{_5}\mathit{m}\mathit{_0}$ (under-5)','D. Crude death rate - $\mathit{_9}\mathit{_5}\mathit{m}\mathit{_5}$ (above-5)','E. Survival function - $\mathit{S}\mathrm{(}\mathit{x}\mathrm{)}$','F. Mortality rates - $\mathit{_1}\mathit{m}\mathit{_x}$'};
xlab                     = {'$\mathit{male}$ (000)                            $\mathit{female}$ (000)','$\mathit{\alpha}$ (level)','$\mathit{deaths}$ $\mathit{per}$ $\mathit{1000}$ $\mathit{person-years}$','$\mathit{deaths}$ $\mathit{per}$ $\mathit{1000}$ $\mathit{person-years}$','$\mathit{x}$ (age in years)','$\mathit{x}$ (age in years)'};
ylab                     = {'$\mathit{age}$','$\mathit{\beta}$ (shape)','$\mathit{kernel}$ $\mathit{density}$','$\mathit{kernel}$ $\mathit{density}$','',''};
for i = 1:6
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontName              = 'Times New Roman';
    ax{i}.FontSize              = 10;
    ax{i}.XAxis.TickLabelFormat = '%.1f';
    ax{i}.YAxis.TickLabelFormat = '%.2f';
    if i == 6
        ax{i}.YScale                = 'log';
    end
    if i  > 4
        xlim([0 100]);
        ax{i}.XAxis.TickLabelFormat = '%.0f';
    end
    xlabel(xlab{i},'Interpreter','latex','FontName','Times New Roman','FontSize',11);
    ylabel(ylab{i},'Interpreter','latex','FontName','Times New Roman','FontSize',11);
    title(LAB{i},'Interpreter','latex');
    grid on;
    box on;
    hold on;
end

nexttile(1)
ax{1}.YAxis.TickLabelFormat = '%.0f';
ax{1}.XAxis.TickLabelFormat = '%.0f';
set(gca,'FontName','Times New Roman','XTick',(-400:200:400),'XTickLabel',[(400:-200:0),(200:200:400)]);
barh((0:100) + .5,max(population(:,1),population(:,2)),1,'FaceColor','k','FaceAlpha',.10,'EdgeAlpha',.00)
barh((0:100) + .5,population(:,1),1,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',.20,'EdgeAlpha',.20)
barh((0:100) + .5,-max(population(:,1),population(:,2)),1,'FaceColor','k','FaceAlpha',.10,'EdgeAlpha',.00)
barh((0:100) + .5,-population(:,2),1,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',.20,'EdgeAlpha',.20)
ylim([0 100])

nexttile(2)
ax{2}.YAxis.TickLabelFormat = '%.1f';
[u{1},u{2}] = contour(a,b,CO{end,1},round([Z(end,1) D(1,1)],1),'color',coloR{3},'LineStyle','-.','ShowText','on','LabelSpacing',150);
clabel(u{1},u{2},'FontSize',8,'FontName','Times New Roman','color',coloR{3},'Interpreter','latex');
[u{1},u{2}] = contour(a,b,CO{end,2},round([Z(end,2) D(1,2)],1),'color',coloR{5},'LineStyle','-.','ShowText','on','LabelSpacing',150);
clabel(u{1},u{2},'FontSize',8,'FontName','Times New Roman','color',coloR{5},'Interpreter','latex');

plot([B(1,1) -2],B(1,2)*ones(2,1),'color','m','Linewidth',1.0);
plot(B(1,1)*ones(2,1),[B(1,2) 0],'color','m','Linewidth',1.0);
plot([B(2,1) -2],B(2,2)*ones(2,1),'color','k','Linewidth',1.0);
plot(B(2,1)*ones(2,1),[B(2,2) 0],'color','k','Linewidth',1.0);
for i = 3:size(D,1)
    plot([B(i,1) -2],B(i,2)*ones(2,1),'color',[0.1 0.1 0.1 0.025],'Linewidth',0.5);
    plot(B(i,1)*ones(2,1),[B(i,2) 0],'color',[0.1 0.1 0.1 0.025],'Linewidth',0.5);
end
legend({'$\mathit{_5}\mathit{m}\mathit{_0}$ (per 1000)','$\mathit{_9}\mathit{_5}\mathit{m}\mathit{_5}$ (per 1000)'},'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
xlim([-1 2])
ylim([0 2])


for i = 3:4
    nexttile(i)
    [f{i},xi{i}]  = ksdensity(D(2:end,i - 2));
    plot(Z(end,i - 2)*ones(2,1),[max(f{i})*1.05 0],'color','m','LineWidth',1.00);
    plot(xi{i},f{i},'color','k','LineWidth',1.00);
    fill(xi{i},f{i},[0.1 0.1 0.1],'FaceAlpha',.05,'EdgeAlpha',0.25,'LineStyle',':','EdgeColor',[0.1 0.1 0.1]);
    
    [f{i},xi{i}]  = ksdensity(D(2:end,i - 2),[D(1,i - 2) prctile(D(2:end,i - 2),[50 2.5 97.5])]);
    %plot(xi{i}(1)*ones(2,1),[f{i}(1) 0],'color','k','LineWidth',0.75,'LineStyle',':')
    for j = 3:4
        plot(xi{i}(j)*ones(2,1),[f{i}(j) 0],'color',[0.1 0.1 0.1],'LineWidth',0.75,'LineStyle','-')
    end
    ylim([0 max(f{i})*1.05])
    legend(sources,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
end

d              = Q(2:end,:) - Q(1:end - 1,:);
l              = 1 - Q;
L              = l(2:end,:) + d.*A;
T              = cumsum(L,'reverse');
e              = T./l(1:end - 1,:);

nexttile(5)
plot((0:101),1 - Q(:,1),'color','m','Linewidth',1.0);
plot((0:101),1 - Q(:,2),'color','k','Linewidth',1.0);
plot((0:101),1 - Q(:,3:end),'color',[0.1 0.1 0.1 0.0125],'Linewidth',0.5);
legend({char("$\mathit{e}\mathit{_0 = " + string(round(e(1,1),1)) + "}$"),char("$\mathit{e}\mathit{_0 = " + string(round(e(1,2),1)) + "}$")},'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');

nexttile(6)
plot((0:100) + .5,Mortality(:,1),'color','m','Linewidth',1.0);
plot((0:100) + .5,Mortality(:,2),'color','k','Linewidth',1.0);
plot((0:100) + .5,Mortality(:,3:end),'color',[0.1 0.1 0.1 0.025],'Linewidth',0.5);
legend(sources,'Interpreter','latex','FontName','Times New Roman','FontSize',9,'FontAngle','oblique','Location','southoutside','NumColumns',3,'Box','off');
saveas(gcf,'RaMMPS/Brass relational.png');
