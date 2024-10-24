function [ allks, rsegind] = generate_k_line( num_kpts, klist )
% Reference citation: Phys. Rev. B 92, 205108 (2015).
% Ab initio tight-binding Hamiltonian for transition metal dichalcogenides
% by Shiang Fang, Rodrick Kuate Defo, Sharmila N. Shirodkar, Simon Lieu, Georgios A. Tritsaris, and Efthimios Kaxiras
% code version: July 2017

klist_size=size(klist);
num_segs=klist_size(1)-1;
allks=zeros(num_segs*num_kpts,3);
rsegind=zeros(num_segs*num_kpts,1);

comp=linspace(0,1,num_kpts);

ind=1;

seg_length(1)=0;
for inds=1:num_segs
    kstart=klist(inds,:);
    kend=klist(inds+1,:);
    deltak=kend-kstart;
    seg_length(inds+1)=sqrt(dot(deltak,deltak));
    for indk=1:num_kpts
        allks(ind,:)=kstart+comp(indk)*(kend-kstart);
        ind=ind+1;
    end
end

sum_seg_length=cumsum(seg_length)/sum(seg_length);

ind=1;
for inds=1:num_segs
    segstart=sum_seg_length(inds);
    segend=sum_seg_length(inds+1);
    
    for indk=1:num_kpts
        rsegind(ind,:)=segstart+comp(indk)*(segend-segstart);
        ind=ind+1;
    end
end

end
