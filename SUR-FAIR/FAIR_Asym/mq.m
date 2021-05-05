function du=mq(y)
% convert monthly series to quarterly one
du=[];
if size(y,1)==1
    y=y';
end
if length(y)/3==int8(length(y)/3)
    for i=1:3:size(y,1)
        du=[du ; sum(y(i:i+2,:))/3];
    end;
else
    if ((length(y)-1)/3)==int8((length(y)-1)/3)
        for i=1:3:size(y,1)-1
            du=[du ; sum(y(i:i+2,:))/3];
        end;
    else
        for i=1:3:size(y,1)-2
            du=[du ; sum(y(i:i+2,:))/3];
        end;
    end
end

