function OuTPuT = recode(InPuT,before,after)

OuTPuT = InPuT;
for j = 1:min(length(before),length(after))
    if isnan(before(j))
        set = find(isnan(InPuT));
    else
        set = find(InPuT == before(j));
    end
    OuTPuT(set) = after(j);
end