function s = CoMMaS(N,foRMaT)

s = sprintf(foRMaT,N);
k = find(ismember(s,'.'));
if numel(k) > 0
    tail = s(k:end);
    s    = s(1:k - 1);
else
    tail = '';
end

d = numel(s);
S = [mod(d,3) ones(1,floor(d/3))*3];
s = mat2cell(s,1,S);
for i = 1:numel(s) - 1
    if numel(s{i}) > 0
        s{i} = [s{i} ','];
    end
end
s = cell2mat([s tail]);