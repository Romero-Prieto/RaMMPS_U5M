function [qS,qP,XB] = RaMMPS_Bayes(N,warmUp)

if N > 0
    pATh       = "/Users/lshjr3/Documents/RaMMPS/under-five mortality/";
    load(char(pATh + "Results/RaMMPS.mat"),'RaMMPS','RESolUTioN');

    x{1}       = [(0:7:28)/365.25,[(2:1:12),(15:3:24),(36:12:60)]/12]';
    x{2}       = [string(0);string([(7:7:28)';(2:1:11)';(12:3:24)';(36:12:60)']) + char([kron('d',ones(4,1));kron('m',ones(18,1))])];
    n          = diff(x{1},1);
    list       = {'TPH','FPH'};
    date       = max(RaMMPS.interview);
    Ts         = {datetime([2014 year(date)]',[1 month(date)]',[1 day(date)]'),datetime([2014 2016]',1,1),datetime([2016 2018]',1,1),datetime([2018 2020]',1,1),datetime([2020 2022]',1,1),datetime([2022 year(date)]',[1 month(date)]',[1 day(date)]')};

    for i = 1:5
        Tx{i} = min(datetime([2014 2014]',1 + 24*[i - 1 i]',[1 1]'),Ts{1}(end));
        tx(i) = datetime(2014,1 + 24*(i - 1/2),1);
    end

    sET        = [RaMMPS.WR,ones(size(RaMMPS.WR))];
    models     = {'post-strat.','selected'};
    epsilon    = [0.125 0.125 0.125 0.125];

    for i = 1:numel(list)
        for j = 1:size(sET,2)
            h          = j + (i - 1)*size(sET,2);
            RaMMPS.W   = sET(:,j);
            dATa{h}    = RaMMPS(RaMMPS.sample == list{i},:);
            for k = 1:numel(Ts)
                date      = year(Ts{k}) + (day(Ts{k},'dayofyear') - 1)./days(datetime(year(Ts{k}) + 1,1,1) - datetime(year(Ts{k}),1,1));
                date      = string(list{i}) + ", " + string(sprintf('%0.1f',date(1))) + "-" + string(sprintf('%0.1f',date(2)));
                date      = char("RaMMPS " + date);
                pOPs{h,k} = {date;models{j}};
                clear date
            end
            vARs{h}    = "$\textit{" + string(models{j}) + "}$";
        end
        sEt{i} = "$\textrm{" + string(list{i}) + "}$";
    end

    for i = 1:numel(dATa)
        s            = dATa{i}(dATa{i}.W > 0,:);
        s.index      = (1:size(s,1))';
        K            = s.index(s.k == 1);
        k            = s.K(s.k == 1);

        rng(0);
        p            = 0.5;   %%%%%%
        p            = rand(size(s,1),1);
        d            = p.*s.D_min + (1 - p).*(datenum(s.D_max) + 1/2);
        
        p            = 0.5;   %%%%%%
        p            = rand(size(s,1),1);
        bR           = (s.birth ~= 'livebirth');
        B            = p.*datenum(s.B_min) + (1 - p).*datenum(s.B_max);    
        B(bR,:)      = NaN;

        bR           = (s.birth ~= 'stillbirth');
        sB           = p.*datenum(s.B_min) + (1 - p).*datenum(s.B_max);
        sB(bR,:)     = NaN;

        sB           = datetime(datevec(sB),'Format','dd/MM/yyyy');
        B            = datetime(datevec(B),'Format','dd/MM/yyyy');
        D            = B + d;
        O            = min(D,s.interview);
        p            = rand(size(s,1),1);
        p            = 0.5;   %%%%%%
        ageS         = p.*s.age + (1 - p).*(s.age + 1);
        ageB         = ageS - years(s.interview - B);

        exposure     = zeros(size(s,1),numel(x{1}) - 1);
        events       = zeros(size(s,1),numel(x{1}) - 1);
        T            = Ts{1};
        for j = 1:numel(x{1}) - 1
            a               = max(min(B + x{1}(j)*365.25,min(O,T(2))),T(1));
            o               = max(min(B + x{1}(j + 1)*365.25,min(O,T(2))),T(1));
            exposure(:,j)   = exposure(:,j) + (datenum(o) - datenum(a))/365.25;
            events(:,j)     = events(:,j) + (d/365.25 >= x{1}(j) & d/365.25 < x{1}(j + 1) & D >= T(1) & D < T(2));
            clear a o
        end

        w            = (exposure > 0).*s.W;
        w            = (w./sum(w,1)).*sum(exposure > 0);
        sEL          = find(sum(exposure,2) > 0 & (s.survival == 'dead' | s.survival == 'alive'));
        sex          = recode(s.sex(sEL),[2 NaN],[-1 0]);
        UR           = recode(s.UR(sEL),2,-1);
        Region       = s.Region(sEL) == [1 3];
        Age          = 2*(s.age(sEL) < 30) - 1;
        Education    = s.Education(sEL) > [1 2];    

        Xn           = [sex UR Region Age Education];

        nX           = [1 2 3 5 6 8 9 9 + (1:numel(Tx))];
        X            = kron(ones(1,numel(x{1}) - 1),[Xn zeros(numel(sEL),numel(Tx))]);
        Xn           = size(Xn,2); 
        X            = mat2cell(X,numel(sEL),ones(1,numel(x{1}) - 1)*(Xn + numel(Tx)));

        for h = 1:numel(Tx)
            T           = Tx{h};
            for j = 1:numel(x{1}) - 1
                a              = max(min(B(sEL) + x{1}(j)*365.25,min(O(sEL),T(2))),T(1));
                o              = max(min(B(sEL) + x{1}(j + 1)*365.25,min(O(sEL),T(2))),T(1));
                X{j}(:,h + Xn) = (datenum(o) - datenum(a))/365.25 > 0;
                clear a o
            end
        end
        
        mo.type      = 'Poisson';
        mo.events    = events(sEL,:);
        mo.exposure  = exposure(sEL,:);
        mo.w         = w(sEL,:);
        mo.X         = X;
        mo.Beta      = [log(sum(mo.events.*mo.w)./sum(mo.exposure.*mo.w))';NaN(size(X{1},2),1)];
        mo.sEL       = find(isnan(mo.Beta));
        [~,mo.B,~,e] = CG_PoissonNB(mo.events,mo.exposure,mo.X,mo.w,1000,mo.type,mo.Beta);

        rEG          = zeros(max(nX) + 1,3);
        rEG(nX,:)    = [mo.B(mo.sEL),mo.B(mo.sEL) + 1.96*e(mo.sEL)*[-1 1]];
        rEG(end,:)   = [size(mo.events,1) NaN(1,2)];
        MLE{i}       = rEG;
        mle{i}       = [mo.B(mo.sEL),e(mo.sEL)];

        [mo,~,a,sam] = HMC(mo,warmUp + N,epsilon(i));
        sam          = sam(warmUp + 1:end,:);
        sam          = sam(~isnan(sam(:,1)),:);
        sample{i}    = sam;
        
        rEG          = zeros(max(nX) + 1,3);
        rEG(nX,:)    = prctile(sam,[50 2.5 97.5])';
        rEG(end,:)   = [size(mo.events,1) NaN(1,2)];
        bOx{i}       = rEG;        
        simulated    = 0;
        
        if isequal(simulated,1)
            B          = (datetime([year(Ts{1}(1)) - 5,month(Ts{1}(1)),day(Ts{1}(1))]):Ts{end}(2) - 1)';
            B          = [B;B];
            O          = datetime(kron(datevec(max(s.interview)),ones(size(B))),'Format','dd/MM/yyyy');
            sEL        = (1:size(B,1));
            w          = ones(size(B,1),numel(x{1}) - 1);

            exposure   = zeros(size(B,1),numel(x{1}) - 1);
            events     = zeros(size(B,1),numel(x{1}) - 1);
            T          = [Tx{1}(1) Tx{end}(2)]';
            for j = 1:numel(x{1}) - 1
                a               = max(min(B + x{1}(j)*365.25,min(O,T(2))),T(1));
                o               = max(min(B + x{1}(j + 1)*365.25,min(O,T(2))),T(1));
                exposure(:,j)   = exposure(:,j) + (datenum(o) - datenum(a))/365.25;
                clear a o
            end

            period     = zeros(size(B,1),numel(Tx));
            for h = 1:numel(Tx) 
                T           = Tx{h};
                for j = 1:numel(x{1}) - 1
                    a               = max(min(B + x{1}(j)*365.25,min(O,T(2))),T(1));
                    o               = max(min(B + x{1}(j + 1)*365.25,min(O,T(2))),T(1));
                    period(:,h)     = period(:,h) + (datenum(o) - datenum(a))/365.25;
                    clear a o
                end
            end
            period     = period./max(sum(exposure,2),eps);
            X          = [[-ones(numel(B)/2);ones(numel(B)/2)] period];
            mo.X       = X(:,nX);
            mo.w       = w;
        end

        pE         = {};
        for h = 1:numel(Ts) 
            T           = Ts{h};
            for j = 1:numel(x{1}) - 1
                a             = max(min(B + x{1}(j)*365.25,min(O,T(2))),T(1));
                o             = max(min(B + x{1}(j + 1)*365.25,min(O,T(2))),T(1));
                pE{h}(:,j)    = (datenum(o(sEL)) - datenum(a(sEL)))/365.25;
                clear a o
            end
        end

        sex        = [-1 1];
        for r = 1:size(sam,1)
            for j = 1:numel(pE)
                events_hat    = exp(mo.B(1:size(mo.events,2))' + cell2mat(mo.X)*kron(eye(size(mo.events,2)),sam(r,:)')).*pE{j};
                mP{i}{j}(:,r) = sum(events_hat.*mo.w)'./sum(pE{j}.*mo.w)';
                qP{i}{j}(:,r) = [0;1 - exp(-cumsum(mP{i}{j}(:,r).*n))];
            end

            for j = 1:numel(sex)
                events_hat    = exp(mo.B(1:size(mo.events,2))' + cell2mat(mo.X)*kron(eye(size(mo.events,2)),sam(r,:)')).*pE{1};
                sEX           = mo.X{1}(:,1) == sex(j);
                mS{i}{j}(:,r) = sum(events_hat.*mo.w.*sEX)'./sum(pE{1}.*mo.w.*sEX)';
                qS{i}{j}(:,r) = [0;1 - exp(-cumsum(mS{i}{j}(:,r).*n))];
            end

            sams          = sam(r,:);
            sams(1)       = 0;
            events_hat    = exp(mo.B(1:size(mo.events,2))' + cell2mat(mo.X)*kron(eye(size(mo.events,2)),sams')).*pE{1};
            mS{i}{3}(:,r) = sum(events_hat.*mo.w)'./sum(pE{1}.*mo.w)';        
            qS{i}{3}(:,r) = [0;1 - exp(-cumsum(mS{i}{j}(:,r).*n)*exp(sam(r,1)))]./[0;1 - exp(-cumsum(mS{i}{j}(:,r).*n)*exp(-sam(r,1)))];
        end
        save(char(pATh + "Results/RaMMPS_Bayes.mat"));
    end
else
    pATh       = "/Users/lshjr3/Documents/RaMMPS/under-five mortality/";
    load(char(pATh + "Results/RaMMPS_Bayes.mat"));
end

XB{1,1}       = "$\textrm{Sex: male = 1; female = -1}$";
XB{2,1}       = "$\textrm{Place of Residence: urban = 1; rural = -1}$";
XB{3,1}       = "$\textrm{Region: North = 1}$";
XB{4,1}       = "$\textrm{Central = 1}$";
XB{5,1}       = "$\textrm{South = 1}$";
XB{6,1}       = "$\textrm{Mother's Age: less than 30 = 1; 30 or more = -1}$";
XB{7,1}       = "$\textrm{Mother's Education: incomplete primary = 1}$";
XB{8,1}       = "$\textrm{incomplete secondary = 1}$";
XB{9,1}       = "$\textrm{complete secondary or more = 1}$";

lABs{1}{1}    = 1;
lABs{2}{1}    = 2;
lABs{3}{1}    = 3;
lABs{3}{2}    = 4;
lABs{3}{3}    = 5;
lABs{4}{1}    = 6;
lABs{5}{1}    = 7;
lABs{5}{2}    = 8;
lABs{5}{3}    = 9;
for i = 1:numel(Tx)
    date          = year(Tx{i}) + (day(Tx{i},'dayofyear') - 1)./days(datetime(year(Tx{i}) + 1,1,1) - datetime(year(Tx{i}),1,1));
    date          = string(sprintf('%0.1f',date(1))) + "-" + string(sprintf('%0.1f',date(2)));    
    XB{end + 1,1} = "$\textrm{period [\it{" + date + "}\rm{)}}$";
    lABs{6}{i}    = max(cell2mat(lABs{5})) + i;
    clear date
end
XB{end + 1,1} = "Reported Children";
lABs{7}{1}    = max(cell2mat(lABs{end})) + 1;

vARs          = {vARs(1:2) vARs(1:2)};
foRMaT        = {'%0.4f','%0.4f','%0.4f'};
nOTe          = {'$\textrm{ln}[\lambda_{j}(x)] = \textrm{ln}[m(x)] + \sum_{}^{}\mathop{}_{\mkern-5mu k} z^{k}_{j}\cdot\beta_{k}$ + $\sum_{}^{}\mathop{}_{\mkern-5mu p} y^{p}_{j}\cdot\gamma_{p} + t_{j}(x)$','$\textrm{Poisson-(Bayesian) MCMC Hamiltonian estimation.}$ $\mathrm{p50}$/$\mathit{[p2.5,p97.5]}$ $\textrm{posterior distribution.}$'};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,XB,cell2mat(bOx),0.230,0.080,[]);
exportgraphics(gcf,char(pATh + "Results/Table_A1.png"),'Resolution',RESolUTioN);


nOTe          = {'$\textrm{ln}[\lambda_{j}(x)] = \textrm{ln}[m(x)] + \sum_{}^{}\mathop{}_{\mkern-5mu k} z^{k}_{j}\cdot\beta_{k}$ + $\sum_{}^{}\mathop{}_{\mkern-5mu p} y^{p}_{j}\cdot\gamma_{p} + t_{j}(x)$','$\textrm{Maximum Likelihood estimation (conjugate gradient methods). 95\% CIs, assuming a normal distribution of the coefficients.}$'};
tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,XB,cell2mat(MLE),0.230,0.080,[]);
exportgraphics(gcf,char(pATh + "Results/Table_A3.png"),'Resolution',RESolUTioN);



models                   = {'Post-stratified','Selected'};
LAB                      = XB(nX);
for i = 1:numel(LAB)
    temp   = char(LAB{i});
    A      = find(ismember(temp,'{}'));
    LAB{i} = temp(A(1) + 1:A(end) - 1);
end
coloR                    = {[0.00 0.00 0.75],[0.95 0.00 0.95],[0.85 0.35 0.01],[0.45 0.65 0.20],[0.65 0.10 0.20],[0.00 0.55 0.65],[0.00 0.00 0.75],[0.95 0.00 0.95],[0.85 0.35 0.01],[0.45 0.65 0.20],[0.65 0.10 0.20],[0.00 0.55 0.65]};

mPIX                     = 538756;
pix                      = 1/37.7952755906;
z                        = min(sqrt(mPIX/((10*4)*(10*3)/pix^2)),1);
fi                       = figure('Color',[1 1 1],'Position',z*10*[0 0 4 3]/pix,'Theme','light');
axes1                    = axes('Parent',fi,'Position',[0.025 0.025 0.975 0.975]);
hold(axes1,'on');
TL                       = tiledlayout(3,4,'Padding','compact','TileSpacing','compact');
for i = 1:size(sample{1},2)
    nexttile(i)
    ax{i}                       = gca;
    ax{i}.FontSize              = 10*z;
    ax{i}.XAxis.TickLabelFormat = '%.2f';
    ax{i}.YAxis.TickLabelFormat = '%.2f';
    
    if i >= 9
        xlabel('$\textbf{Estimated coefficient}$','Interpreter','latex','FontSize',11*z);
    end
    if isequal(mod(i,4),1)
        ylabel('$\textbf{Probability density function}$','Interpreter','latex','FontSize',11*z);
    end
    title(char(string("$\textbf{" + char(96 + i)) + ". " + LAB{i}) + "}$",'Interpreter','latex');
    grid on;
    box on;
    hold on;
end

for j = 1:numel(sample)
    for i = 1:size(sample{1},2)
        nexttile(i)
        ax{i}.XTickMode  = 'auto';
        ax{i}.YTickMode  = 'auto';
        ax{i}.XLabelMode = 'auto';
        ax{i}.YLabelMode = 'auto';

        H{i}             = histogram(sample{j}(:,i),25,'Normalization','pdf','FaceColor',coloR{i},'EdgeColor',[0 0 0],'FaceAlpha',0.1,'EdgeAlpha',0.25);
        [f,xi]           = ksdensity(sample{j}(:,i));
        P{i}             = plot(xi,f,'LineWidth',1.25*z,'Color',coloR{i});        
        w                = prctile(sample{j}(:,i),[0.1 99.9]);
        xlim([w(1) w(2)])
        ylim([0 max(f)*1.15])
        
        y                = prctile(sample{j}(:,i),[50 2.5 97.5]);
        [f,y]            = ksdensity(sample{j}(:,i),y);
        for k = 1:numel(y)
             G{i,k} = plot(y([k k]),[0 f(k)],'LineWidth',1.00*z,'Color',coloR{i},'LineStyle','-.');
        end
        G{i,end + 1}     = plot([0 0],[0 max(f)*1.50],'LineWidth',1.15*z,'Color','k','LineStyle','-');

        t                = ax{i}.XAxis.TickValues;
        for k = 1:numel(t)
            Tk{k} = char("$\mathbf{" + string(sprintf(ax{i}.XAxis.TickLabelFormat,t(k))) + "}$");
        end
        ax{i}.XTickLabelRotation         = 0;
        ax{i}.XAxis.TickLabels           = Tk;
        ax{i}.XAxis.TickLabelInterpreter = 'latex';
        clear t Tk

        t                = ax{i}.YAxis.TickValues;
        for k = 1:numel(t)
            Tk{k} = char("$\mathbf{" + string(sprintf(ax{i}.YAxis.TickLabelFormat,t(k))) + "}$");
        end
        ax{i}.YTickLabelRotation         = 0;
        ax{i}.YAxis.TickLabels           = Tk;
        ax{i}.YAxis.TickLabelInterpreter = 'latex';
        clear t Tk
        ax{i}.XTickMode  = 'manual';
        ax{i}.YTickMode  = 'manual';
        ax{i}.XLabelMode = 'manual';
        ax{i}.YLabelMode = 'manual';
    end

    saveas(gcf,char(pATh + "Results/Coeff-" + list{ceil(j/2)} + "-" + models{numel(models) - mod(j,2)} + ".png"))
    for i = 1:size(sample{1},2)
        delete(H{i});
        delete(P{i});
        for k = 1:size(G,2)
            delete(G{i,k});
        end
    end
    clear P H G
end


CO = [6 3];
for j = 1:2:numel(sample)
    for i = 1:size(sample{1},2)
        nexttile(i)
        ax{i}.XTickMode  = 'auto';
        ax{i}.YTickMode  = 'auto';
        ax{i}.XLabelMode = 'auto';
        ax{i}.YLabelMode = 'auto';

        for h = 0:1
            H{1 + h}{i} = histogram(sample{j + h}(:,i),25,'Normalization','pdf','FaceColor',coloR{CO(1 + h)},'EdgeColor',[0 0 0],'FaceAlpha',0.1,'EdgeAlpha',0.25); 
        end
        
        for h = 0:1
            [f,xi]      = ksdensity(sample{j + h}(:,i));
            P{1 + h}{i} = plot(xi,f,'LineWidth',1.25*z,'Color',coloR{CO(1 + h)});
            
            w(1 + h,:)  = prctile(sample{j + h}(:,i),[0.1 99.9]);
            xlim([min(w(:,1)) max(w(:,2))])
            fw(1 + h)   = max(f)*1.15;          
            ylim([0 max(fw)])
            
            y           = prctile(sample{j + h}(:,i),[50 2.5 97.5]);
            [f,y]       = ksdensity(sample{j + h}(:,i),y);
            for k = 1:numel(y)
                G{1 + h}{i,k} = plot(y([k k]),[0 f(k)],'LineWidth',0.75*z,'Color',coloR{CO(1 + h)},'LineStyle','-');
            end
            if isequal(h,1)
                G{1}{i,4} = plot([0 0],[0 max(fw)*1.50],'LineWidth',1.25*z,'Color','k','LineStyle','-');
            end
        end
        clear w fw
        
        if isequal(i,9)
            for k = 1:numel(models)
                leGend{k} = char("$\textbf{" + models{k} + "}$");
            end
            legend(leGend,'Interpreter','latex','FontSize',11*z,'FontAngle','oblique','Location','southoutside','NumColumns',2,'Box','off');
        end

        t                = ax{i}.XAxis.TickValues;
        for k = 1:numel(t)
            Tk{k} = char("$\mathbf{" + string(sprintf(ax{i}.XAxis.TickLabelFormat,t(k))) + "}$");
        end
        ax{i}.XTickLabelRotation         = 0;
        ax{i}.XAxis.TickLabels           = Tk;
        ax{i}.XAxis.TickLabelInterpreter = 'latex';
        clear t Tk

        t                = ax{i}.YAxis.TickValues;
        for k = 1:numel(t)
            Tk{k} = char("$\mathbf{" + string(sprintf(ax{i}.YAxis.TickLabelFormat,t(k))) + "}$");
        end
        ax{i}.YTickLabelRotation         = 0;
        ax{i}.YAxis.TickLabels           = Tk;
        ax{i}.YAxis.TickLabelInterpreter = 'latex';
        clear t Tk
        ax{i}.XTickMode  = 'manual';
        ax{i}.YTickMode  = 'manual';
        ax{i}.XLabelMode = 'manual';
        ax{i}.YLabelMode = 'manual';
    end
    
    saveas(gcf,char(pATh + "Results/Figure_A" + ((j + 1)/2 + 3) + ".png"))
    for i = 1:size(sample{1},2)
        for h = 0:1
            delete(H{1 + h}{i});
            delete(P{1 + h}{i});
            for k = 1:size(G{1 + h},2)
                delete(G{1 + h}{i,k});
            end
        end
    end
end