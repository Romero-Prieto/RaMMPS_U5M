function U = BraSsTruSseLl(W,B,D,T)
load('BraSsTruSseLl.mat','CD','coeff');
P = B./W;
for j = 1:numel(coeff)
    LT     = CD(j,:);
    x      = LT{1}.age;
    LT     = (table2array(LT{end - 1}) + 1.05*table2array(LT{end}))/2.05;
    [~,i]  = ismember([1 2 3 5 10 15 20],x);
    LT     = LT(i,:);
    
    k      = coeff{j}(:,2) + coeff{j}(:,3)*(P(1,:)./P(2,:)) + coeff{j}(:,4)*(P(2,:)./P(3,:));
    k      = max(k,0);
    q      = (D./B).*k;
    t      = coeff{j}(:,5) + coeff{j}(:,6)*(P(1,:)./P(2,:)) + coeff{j}(:,7)*(P(2,:)./P(3,:));    
        
    S      = 0;
    for h = 1:size(LT,2) - 1
        S      = S + (q < LT(:,h));
    end
    
    for h = 1:numel(i)
        LB      = LT(h,S(h,:) + 1);
        UB      = LT(h,S(h,:));
        s       = (q(h,:) - LB)./(UB - LB);
        q1(h,:) = LT(1,S(h,:) + 1).*s + LT(1,S(h,:)).*(1 - s);
        q5(h,:) = LT(4,S(h,:) + 1).*s + LT(4,S(h,:)).*(1 - s);
        clear LB UB s
    end
    U{j,1} = NaN(size(q1));
    U{j,2} = q1;
    U{j,3} = q5;
    U{j,4} = reshape(repelem(T,7)',7,numel(T)) - t;
    clear q1 q5 q k t S LT x i 
end