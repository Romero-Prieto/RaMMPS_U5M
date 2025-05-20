function tABleBAyEs(sEt,vARs,foRMaT,lABs,nOTe,leGend,dATa,R,d,W)

if min(numel(vARs{1}),numel(sEt)) == 0
    S  = [0.060 1/50];
    Sx = [0.025 0.025];
else
    S  = [0.085 1/50];
    Sx = [0.025 0.050];
end
s       = 3/200;
f       = 0.045;
r       = R/25;

n       = NaN(numel(vARs),1);
for i = 1:numel(vARs)
    n(i) = numel(vARs{i});
end
D       = R + 2*r*numel(n) + sum(n)*d;
xT      = NaN(1,numel(n));
x       = NaN(1,sum(n));
xlab    = num2cell(x);
xTab    = NaN(numel(n),2);
for i = 1:numel(n)
    for j = 1:n(i)
        x(sum(n(1:i - 1)) + j)    = (R + 2*r*(i - 1) + d*sum(n(1:i - 1)) + r + d*j)/D;
        xlab{sum(n(1:i - 1)) + j} = vARs{i}{j};
    end
    xT(i)     = (R + sum(n(1:i))*d + 2*r*i - (r + n(i)/2*d))/D;
    xTab(i,:) = (R + sum(n(1:i - 1))*d + 2*r*(i -1) + r + [0 n(i)*d])/D;
end
xTab(1) = r/D; 


n       = NaN(numel(lABs),1);
for i = 1:numel(lABs)
    n(i) = numel(lABs{i});
end
F    = sum(S) + s*(numel(n) - 1) + sum(n)*f;
y    = NaN(1,sum(n));
ylab = num2cell(x);
data = NaN(numel(y),numel(x)*3);
for i = 1:numel(n)
    for j = 1:n(i)
        y(sum(n(1:i - 1)) + j)      = 1 - (S(1) + s*(i - 1) + f*sum(n(1:i - 1)) + f*(j - 1))/F;
        ylab{sum(n(1:i - 1)) + j}   = leGend{lABs{i}{j}};
        data(sum(n(1:i - 1)) + j,:) = dATa(lABs{i}{j},:);
    end
end

if numel(W) == 0
    W      = ones(numel(y),numel(vARs));
end
for i = 1:size(W,2)
    w{1,i} = kron(W(:,i),ones(1,numel(vARs{i})));
end
W       = cell2mat(w);


pix                      = 1/37.7952755906;
fi                       = figure('Color',[1 1 1]);
fi.Position              = [0 0 27*D 15*F]/pix;
axes1                    = axes('Parent',fi,'Position',[0.00 0.00 1.00 1.00]);
hold(axes1,'on');
axis off
xlim([0 1])
ylim([0 1])

plot([0 1],.999*ones(2,1),'color','k','Linewidth',1.25);
plot([0 1],(S(2) + s + 0.005)/F*ones(2,1),'color','k','Linewidth',0.75);
plot([0 1],(S(2) + s)/F*ones(2,1),'color','k','Linewidth',1.25);
text(xTab(1,1),1 - Sx(2)/F,nOTe{1},'Interpreter','latex','HorizontalAlignment','left','FontName','Times New Roman','FontSize',10);
text(xTab(1,1),(S(2) - s/2)/F,nOTe{2},'Interpreter','latex','HorizontalAlignment','left','FontName','Times New Roman','FontSize',10);

for i = 1:numel(xT)
    if numel(sEt) > 0
        text(xT(i),1 - Sx(1)/F,sEt{i},'Interpreter','latex','HorizontalAlignment','center','FontName','Times New Roman','FontSize',13);
    end
    plot([xTab(i,1) xTab(i,2)],1 - (S(1) - 0.02)/F*ones(2,1),'color','k','Linewidth',0.75);
end

for i = 1:numel(y)
    if numel(ylab{i}) > 1
        text(r/D,y(i),ylab{i}{1},'Interpreter','latex','HorizontalAlignment','left','FontName','Times New Roman','FontSize',11);
        text(r/D + 0.005/D,y(i) - 0.015/F,ylab{i}{2},'Interpreter','latex','HorizontalAlignment','left','FontName','Times New Roman','FontSize',7.5);
    else
        text(r/D,y(i),ylab{i},'Interpreter','latex','HorizontalAlignment','left','FontName','Times New Roman','FontSize',11);
    end
    
    if numel(y) > 5
        if numel(y) < 27
            numeRaL = char(64 + i);
        else
            numeRaL = char(string(i));
        end
        text(r/4/D,y(i),numeRaL,'Interpreter','latex','HorizontalAlignment','center','FontName','Times New Roman','FontSize',5);
        text(1 - r/4/D,y(i),numeRaL,'Interpreter','latex','HorizontalAlignment','center','FontName','Times New Roman','FontSize',5);
    end
end

for i = 1:numel(x)
    text(x(i),1 - Sx(2)/F,xlab{i},'Interpreter','latex','HorizontalAlignment','right','FontName','Times New Roman','FontSize',11);
    for j = 1:numel(y)
        if isequal(W(j,i),0) 
            temp = [0.25 0.35 0.75];
        elseif sign(data(j,3*i - 1))*sign(data(j,3*i)) < 0
            temp = [0.85 0.25 0.65];
            %temp = [0.00 0.00 0.00];
        else
            temp = [0.00 0.00 0.00];
        end
        
        if ~isnan(data(j,3*i - 2))
            if isequal(floor(data(j,3*i - 2)),ceil(data(j,3*i - 2))) && data(j,3*i - 2) > 999              
                text(x(i),y(j),oBS(data(j,3*i - 2),foRMaT{1 + mod(numel(foRMaT) - 1 + i,numel(foRMaT))}),'Interpreter','latex','HorizontalAlignment','right','FontName','Times New Roman','FontSize',11,'color',temp);
            else
                text(x(i),y(j),sprintf(foRMaT{1 + mod(numel(foRMaT) - 1 + i,numel(foRMaT))},data(j,3*i - 2)),'Interpreter','latex','HorizontalAlignment','right','FontName','Times New Roman','FontSize',11,'color',temp);
            end
            
            if ~isequalwithequalnans(data(j,3*i - 1),data(j,3*i))
                text(x(i),y(j) - 0.015/F,char("$\mathit{[" + string(sprintf(foRMaT{1 + mod(numel(foRMaT) - 1 + i,numel(foRMaT))},data(j,3*i - 1))) + "," + string(sprintf(foRMaT{1 + mod(numel(foRMaT) - 1 + i,numel(foRMaT))},data(j,3*i))) + "]}$"),'Interpreter','latex','HorizontalAlignment','right','FontName','Times New Roman','FontSize',7.5,'color',temp);
            end
        end
    end
end


function N = oBS(number,format)
number = char(string(number));
n      = numel(number);
set    = cumsum([1 mod(n,3) ones(1,floor(n/3))*3]);
N      = "";
for i = 1:numel(set) - 1
    N = N + string(number(set(i):set(i + 1) - 1));
    if i < numel(set) - 1
        N = N + ",";
    end
end
s      = str2num(format(end-1)) - 1 ;
s      = string(char([32*ones(1,s) 32]));
N      = "$\mathrm{" + N + "}$";
N      = number;
