function date = eXAcTTime(date)

if isdatetime(date(1))
    y     = year(date);
    mIn   = datetime(y,1,1,0,0,0);
    mAx   = datetime(y + 1,1,1,0,0,0);
    date  = y + datenum(date - mIn)./datenum(mAx - mIn);
    date  = round(date,7);
else
    date  = round(date,7);
    y     = floor(date);
    mIn   = datetime(y,1,1,0,0,0);
    mAx   = datetime(y + 1,1,1,0,0,0);
    date  = datetime(datenum(mIn) + datenum(mAx - mIn).*(date - y),'ConvertFrom','datenum');
end