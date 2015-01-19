function [resultQ]=deletesnvQ(Q,id_list)

len=length(id_list);
resultQ=Q;
for i=1:len
    z=find(resultQ==id_list(i));
    if (~isempty(z))
        resultQ=[resultQ(1:z-1,:);resultQ(z+1:size(resultQ,1),:)];
    end
end
