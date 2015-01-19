function [resultQ]=insertsnvQ(Q,newid,value)

len=size(Q,1);
pos=1;
while (pos<=len)
    if (Q(pos,2)<value)
        %insert
        resultQ=[Q(1:pos-1,:);newid,value;Q(pos:len,:)];
        return;
    else
        pos=pos+1;
    end
end
resultQ=[Q(1:len,:);newid,value];