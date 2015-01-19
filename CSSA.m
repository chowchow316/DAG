function [x]=CSSA(b,k,Tree)

%b:the recovered long vector d*N
%k:the sparsity
%Tree: the label hierarchy
%x:best k-tree-sparse b
%%son: Son(i,1)-num of children of i; Son(i,2:Son(i,1)+1)-id of each child of i
load('Son.mat');
N=size(b,2);
d=size(b,1);
x=sparse(d,N);
set=cell(N,1);
for i=1:N
%%initialization
    display(i);
    %set all nodes unselected
    selected=sparse(1,d);
    selected(1)=1;%select the tree root

    % set up each node as a supernode.  use an Nx4 matrix:
    %    col 1: leaf id (the orignal node that the supernode grows from)
    %    col 2: #rp  (number of root nodes)
    %    col 3: snv  (supernode value)
    %    col 4: num  (number of internal nodes)
    sn(1:d,1)=(1:d)';
    sn(1:d,2)=ones(d,1);
    sn(1:d,3)=b(1:d,i);
    sn(1:d,4)=ones(d,1);

    %set up each supernode's root node list as itself
    sl=sparse(d,1000);
    sl(1:d,1)=(1:d)';

    %%set up each supernode's internal node list as itself
    si=sparse(d,1000);
    si(1:d,1)=(1:d)';
    
    %set up supernode queue
    snvQ=sparse(d-1,2);
    snvQ(1:d-1,1)=(2:d)'; %%supernode id
    snvQ(1:d-1,2)=b(2:d,i);%%supernode value

    %initialte each node occurs in its initialized supernode.
    node=sparse(d,1000);
    node(1:d,1)=(1:d)';
    
    snvQ=sortrows(snvQ,-2);%%sort snvQ according to sn's value
%   snvQ=snvsort(snvQ);
   

    
    num=1;
    snum=d;%number of generated supernode
    count=0;
    
    while (num<k && count<=200)
        count=count+1;
%         
        if mod(count,100)==0
             display(count);
         end
    
     
    %Step1:retrieve the supernode with the largest snv and delete it from
    %snvQ;  
        max_sn=snvQ(1,:);
        Qlength=size(snvQ,1);
        snvQ=snvQ(2:Qlength,:); %delete the supernode from snvQ
        sn_id=max_sn(1);
        root_num=sn(sn_id,2);
  
    %Step2:find the rootid's parent
    %if all parent are selected,  select the supernode's all nodes
        add_list=[];
        for j=1:root_num
            root_id=sl(sn_id,j);
            parent_num=Tree(root_id,1);
            
            for l=1:parent_num
                parent_id=Tree(root_id,l+1);
                               
                if (selected(parent_id)==0)
                    if (isempty(find(node(parent_id,:)==sn_id, 1)))%% if parent is already in the supernode, do not add it                    
                       add_list=[add_list,parent_id]; %% the j-th root node's l's parent is not selected, add it to the nodes to be merged
                    end
                end
                add_list=unique(add_list);
            end
        end
        
        if (isempty(add_list)) %% all root nodes' parents are selected
            if ((sn(sn_id,4)+num)<=k)
                content=1;
                num=sn(sn_id,4)+num;
                
            else
                content=(k-num)/sn(sn_id,4);
                num=k;
            end
 
            del_queue=[];
            sn_queue=si(sn_id,:);
            sn_length=length(find(si(sn_id,:)>0));
            for j=1:sn_length
                    selected(si(sn_id,j))=content;  
                    set{i,1}=[set{i,1},si(sn_id,j)];
 
            %%delete other supernodes that share nodes with
            %%supernode_{sn_id}
                    tmp_lens=length(find(node(si(sn_id,j),:)>0));
                    del_queue=union(del_queue,node(si(sn_id,j),1:tmp_lens));
            end
           
            del_queue=setdiff(del_queue,sn_id);
            %%nodes in supernodes in del_queue
            decomp_queue=[];
            for j=1:length(del_queue)
                tmp_lens=length(find(si(del_queue(j),:)>0));
                decomp_queue=union(decomp_queue,si(del_queue(j),1:tmp_lens));
            end
            decomp_queue=unique(decomp_queue);
            z=find(selected>0);
            decomp_queue=setdiff(decomp_queue,z);
            %delete 
            snvQ=deletesnvQ(snvQ,del_queue);
            %%decompose
            lens=length(decomp_queue);
            
                    sn(decomp_queue,1)=decomp_queue';
                    sn(decomp_queue,2)=ones(lens,1);
                    sn(decomp_queue,3)=b(decomp_queue,i);
                    sn(decomp_queue,4)=ones(lens,1);
            
                    
                    for l=1:lens
                    %insert into supernode queue  
                        if isempty(find(snvQ(:,1)==decomp_queue(l), 1))
                         snvQ=insertsnvQ(snvQ,decomp_queue(l),sn(decomp_queue(l),3));
                        end
                         %initiate each node occurs in its initialized
                     %supernode and its root set as itself           
                         tmp=length(find(node(decomp_queue(l),:)>0));
                         node(decomp_queue(l),1:tmp)=0;
                         node(decomp_queue(l),1)=decomp_queue(l);
              
                         tmp=length(find(sl(decomp_queue(l),:)>0));
                         sl(decomp_queue(l),1:tmp)=0;    
                         sl(decomp_queue(l),1)=decomp_queue(l);
                         tmp=length(find(si(decomp_queue(l),:)>0));
                         si(decomp_queue(l),1:tmp)=0;
                         si(decomp_queue(l),1)=decomp_queue(l);
                    end
           
        %else form new supernodes and insert it into snvQ
        else
           %% find the parent node that has the min parent supernode value
           min_pos=0;
           for j=1:length(add_list)
               node_list=node(add_list(j),:);
               list_lens=length(find(node(add_list(j),:)>0));
               pos=0;
               for l=1:list_lens
                   cpos=find(snvQ(:,1)==node_list(l));
                   
                   if cpos>pos
                       pos=cpos;
                   end
               end
               if pos>min_pos
                   min_pos=pos;
                   min_pa=add_list(j);%%the selected parent node id
               end
           end
           %%merge sn_id with every supernode that contains min_pa
           node_list=node(min_pa,:);
           list_lens=length(find(node(min_pa,:)>0));
           for j=1:list_lens
               %%check if contain nodes in sn_id
               pa_sn=node_list(j);
  
               z1=si(pa_sn,:);
               z2=si(sn_id,:);
               z3=intersect(z1,z2);%%has no nodes in common
               z3=setdiff(z3,0);
              if isempty(z3)
               %%update rootnode num and rootnode list
                 snum=snum+1;
           
                 newid=snum;
                
                pending=[];
                rootset=sl(sn_id,:);
                rootlens=length(find(rootset>0));
                for l=1:rootlens
                    rootnum=Tree(rootset(l),1);
                    for l1=1:rootnum
                        if (Tree(rootset(l),1+l1)==min_pa)
                            pending=[pending,sl(sn_id,l)];
                            break;
                        end
                    end 
                end
                
               new_root_id_list=setdiff(union(sl(pa_sn,:),sl(sn_id,:)),pending);
       
            %%find which root node in pending has parent not in newnode
             for jj=1:length(pending)
                cur_id=pending(jj);
                flag=0;
                for l=1:Tree(cur_id,1)
                    cur_p_id=Tree(cur_id,l+1);
                    if (isempty(find(node(cur_p_id,:)==sn_id, 1)) && isempty(find(node(cur_p_id,:)==pa_sn, 1)))
                        flag=1;
                        break;
                    end
                end
                if (flag==1)
                    new_root_id_list=[new_root_id_list,pending(jj)];
                end
            end
            %%
            new_root_id_list=unique(new_root_id_list);
            new_root_id_list=setdiff(new_root_id_list,0);
            sl(newid,1:length(new_root_id_list))=new_root_id_list;
            sn(newid,1)=sn(sn_id,1);
            sn(newid,2)=length(find(sl(newid,:)));
            inter=sn(sn_id,4)+sn(pa_sn,4);
            sn_snv=(sn(sn_id,3)*sn(sn_id,4)+sn(pa_sn,3)*sn(pa_sn,4))/inter;
                
            %update SNV and inter node num
            sn(newid,3)=sn_snv;
            sn(newid,4)=inter;
        
        
            si(newid,1:sn(newid,4))=setdiff(union(si(sn_id,:),si(pa_sn,:)),0);
            
            %update node belonging to which supernode
            for l=1:sn(newid,4)
                in_id=si(newid,l);
                tmp_lens=length(find(node(in_id,:)>0));
                node(in_id,tmp_lens+1)=newid;
                
                z2=find(node(in_id,:)==sn_id);%%if node belongs to sn_id, delete sn_id from its node list
                if (~isempty(z2))
                    for jj=z2:tmp_lens
                        node(in_id,jj)=node(in_id,jj+1);
                    end
                    node(in_id,tmp_lens+1)=0;
                end
            end
            
            %insert
            snvQ=insertsnvQ(snvQ,newid,sn(newid,3));
            end
        
           end
              
        end
    end
    
    z=find(selected>0);
 %   x(:,i)=selected;
    if (length(z)>k)
        %%greedy
            root_id=1;
            pending=[];
            s=1;
            s_id=sparse(1,d);
            s_id(1)=1;
            while (s<k)
                new_id=Son(root_id,2:(Son(root_id,1)+1));
                pending=union(pending,new_id);
                pending=intersect(pending,z);
                [v,q]=sort(b1(pending,i),'descend');
                r=1;
                flag=0;
                while (flag==0 && r<=length(q))
                    root_id=pending(q(r));
                    flag=1;
                    for t=1:Tree(root_id,1)
                        if (s_id(Tree(root_id,t+1))==0)
                            flag=0;
                            break;
                        end
                    end
                    if (flag==0)
                        r=r+1;
                    else
                        s_id(root_id)=1;
                        s=s+1;
                        pending=setdiff(pending,root_id);
                    end
                end
            end
        z=find(s_id~=0);
        x(z,i)=b1(z,i);
    else
    x(z,i)=b1(z,i);
    end

end

