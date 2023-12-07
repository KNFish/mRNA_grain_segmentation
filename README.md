# mRNA_grain_segmentation
MatLab code for segmenting mRNA grains

% replace <text> throughout with appropriate information
% mat files must be number <>_<>_<>_<> where <> are numbers

myDir = <insert directory where your mat files are lopcated>;
outDir = <insert directory where you want the mask files saved>;
outDir2 = <insert directory where you want the stats files saved>;

channels=who;
ch_label={'<channel name>'}; *add up to 4 channels

cofac=10;
cofac2=[<>]; % inserted value is multiplied by the Otsu value to determine graylevel start value
cofac_lipo=[<>]; % Ditto

cameraspecs = [<>]; % insert pixel size in on dimension followed by z-step
inc=<>; % number of gray values to jump for each segmentation interation
maxsize=<>; % insert max object size
minsize=<>; % insert min obj size


%find files in folder
[mfiles] = dir([Dir '\*.mat']);
mfiles=vertcat({mfiles(:).name})';

%for
parfor f=1:length(mfiles)
    
    subdir = [Dir mfiles{f}];
    %read in data
  
    %the x405E below is used for mask channels that are imported
    myVars= {'LipoE' 'name' 'x405E' 'x488E' 'x568E' 'x647E'};
    slide=matfile(subdir);  
    
   [finalmask,lipo,gs,puncta]=intensity_based_segmentation_rnascope...
       (slide,cofac,cofac2,cofac_lipo,ch_label, ...
       cameraspecs,inc,minsize,maxsize);
    tic
    
    writemaskout_rnascope...
        (finalmask,lipo,ch_label,slide,outDir,outDir2,gs,puncta)
    toc
end

delete(gcp);


function [finalmask,lipo,gs,puncta]=intensity_based_segmentation_rnascope...
    (slide,cofac,cofac2,cofac_lipo,ch_label,cameraspecs,...
    inc,minsize,maxsize)
% minsize is a matrix organized as channel x gate if you want to use the
% new method specify the range of gates for each channel
% eg maxsize=[0.4:0.1:2;0.4:0.1:2;0.4:0.1:2;0.4:0.1:2;];
% same goes for minsize
% cameraspecs defines the xy and z pixel sizes
% method is 'New' to activate variable gating
% pixsize=0.108; inc=50; step=0.25;

if size(maxsize,1)==1
    minsize= repmat(minsize(1,:),length(ch_label),1);
    maxsize= repmat(maxsize(1,:),length(ch_label),1);
end
tmp=fieldnames(slide);
lipoorg=slide.(tmp{~cellfun(@isempty, strfind(tmp,'Lipo'))});
threshlipo=multithresh(lipoorg);
% threshlipo=1;

%line below is where you are creating the logical data type version of
%the lipo mask
lipo = bsxfun(@lt,lipoorg,ceil(threshlipo)*cofac_lipo);
pixsize=cameraspecs(1); step=cameraspecs(2);

finalmask=cell(length(ch_label),1);
tic
for ch=1:length(ch_label)
    Tsmall = ceil(minsize(ch,1)/(pixsize*pixsize*step));
    Tbig = round(maxsize(ch,1)/(pixsize*pixsize*step));
    raw=slide.([ch_label{ch}]);
    
    %channel=bsxfun(@minus, imgaussfilt3(raw,0.7),imgaussfilt3(raw,2));
    %If you want to mask the gaussian channels uncomment above line
    %and comment below line
    channel=raw;
    channel(channel<0)=0;
    gs{ch}=channel;

    nolipo=(channel+1);
    thresh=multithresh(nolipo);
    Tstart = ceil(thresh)*cofac2(ch);

    
    thresh2=Tstart:inc:max(channel(:));
    CC=arrayfun(@(x) bwconncomp(bsxfun(@ge,channel,x)),thresh2,...
        'UniformOutput',false);
    area=arrayfun(@(x) cellfun(@numel, x{:}.PixelIdxList),CC,...
        'UniformOutput',false);
    [maskout]=grow_puncta_master(channel,area,CC,minsize(ch),...
        maxsize(ch),pixsize,step);
    fmask=zeros(size(channel));
    fmask(vertcat(maskout{:}))=1;
    finalmask{ch}=double(fmask);
    fmask=[];
    
    puncta{ch}.Connectivity=26;
    puncta{ch}.ImageSize=size(channel);
    puncta{ch}.NumObjects=length(maskout);
    puncta{ch}.PixelIdxList=maskout';    
end
lipo = bsxfun(@ge,lipoorg,ceil(threshlipo)*cofac_lipo);
toc
lipo=double(lipo);
end


function [maskout]=grow_puncta_master(channel,area,CC,min_vol,max_vol,...
    pixsize,step)

toosmall=ceil(min_vol/(pixsize^2*step));
toobig=round(max_vol/(pixsize^2*step));
npix=numel(channel);

CC=cellfun(@(x,y) x.PixelIdxList(y>=toosmall & y<=toobig),CC,area...
    ,'UniformOutput',0);
CC=CC(cellfun(@(x) ~isempty(x),CC));

[~,maxind]=cellfun(@(x) cellfun(@(y) max(channel(y)),x) ,CC,...
    'UniformOutput',0);
maxind=cellfun(@(x,y) cellfun(@(x2,y2) y2(x2),...
    num2cell(x),y),maxind,CC,'UniformOutput',0);

if npix<2^32
    maxind=cellfun(@(x) uint32(x),maxind,'UniformOutput',0);
end

[maxpuncta]=find_biggest_puncta_mask(maxind);
maskout=cellfun(@(x,y) CC{x}{y},...
    num2cell(maxpuncta(:,1)),num2cell(maxpuncta(:,2)),'UniformOutput',0);

[hungry,eatenind,lostind]=find_merges(maskout,maxpuncta,maxind,npix);

[maskout]=unmerge(maskout,maxind,hungry,CC,maxpuncta,eatenind,lostind);

end


function [maxpuncta]=find_biggest_puncta_mask(maxind)
umaxind=unique(horzcat(maxind{:}));
npuncta=numel(umaxind);

maxpuncta=zeros(npuncta,2);
tempumaxind=num2cell(umaxind);
punctaind=1:npuncta;
cnt=1;
while ~isempty(punctaind)
    found=cellfun(@(x) find(x==maxind{cnt}),tempumaxind,...
        'UniformOutput',0);
    checkfound=cellfun(@(x) ~isempty(x),found);
    found=found(checkfound);
    %save out puncta index locations
    maxpuncta(punctaind(checkfound),:)=[repmat(cnt,sum(checkfound),1),...
        vertcat(found{:})];
    %update loop variables
    checkfound=(~checkfound);
    punctaind=punctaind(checkfound);
    tempumaxind=tempumaxind(checkfound);
    cnt=cnt+1;
end
end


function [hungry,eatenind,lostind]=find_merges(maskout,maxpuncta,maxind,...
    npix)
lostind=find(maxpuncta(:,1)~=1);
lost_maxind=cellfun(@(x,y) maxind{x}(y),num2cell(maxpuncta(lostind,1)),...
    num2cell(maxpuncta(lostind,2)),'UniformOutput',0);

indmaskout=cellfun(@(x,y) repmat(x,numel(y),1),...
    num2cell(1:length(maskout))',maskout,'UniformOutput',0);
indmaskout=vertcat(indmaskout{:});
linmaskout=vertcat(maskout{:});
if npix<2^32; linmaskout=uint32(linmaskout); end

hungry=cellfun(@(x) indmaskout(bsxfun(@eq,linmaskout,x)),lost_maxind,...
    'UniformOutput',0);

hungry=cellfun(@(x,y) x(x~=y),hungry,num2cell(lostind),...
    'UniformOutput',0);
eatenind=find(cellfun(@(x) ~isempty(x),hungry));
hungry=hungry(eatenind);
end


function [maskout]=unmerge(maskout,maxind,hungry,CC,maxpuncta,eatenind,...
    lostind)
%create shortcut indexing variables
nmaxind=length(maxind);
uhungry=unique(vertcat(hungry{:}));
%find largest unmerged hungry puncta
ref=cellfun(@(x) maxpuncta(lostind(eatenind(cell2mat(...
    cellfun(@(y) any(y==x),hungry,'UniformOutput',0)))),:),...
    num2cell(uhungry),'UniformOutput',0);
start=num2cell(cellfun(@(x) max(x(:,1)),ref));
compareind=cellfun(@(x) maxind{maxpuncta(x,1)}(maxpuncta(x,2)),...
    num2cell(uhungry),'UniformOutput',0);
while ~isempty(start)
    checkmax=cellfun(@(x) x==nmaxind,start);
    if any(checkmax)
        maskout(uhungry(checkmax))={nan};
        checkmax=(~checkmax);
        start=start(checkmax);
        uhungry=uhungry(checkmax);
        compareind=compareind(checkmax);
    end
    found=cellfun(@(x,y) find(maxind{x}==y),start,compareind,...
        'UniformOutput',0);
    isfound=cellfun(@(x) ~isempty(x),found);
    %update loop vars
    maskout(uhungry(isfound))=cellfun(@(x,y) CC{x}{y},start(isfound),...
        found(isfound),'UniformOutput',0);
    isfound=(~isfound);
    start=num2cell(cellfun(@(x) x+1,start(isfound)));
    uhungry=uhungry(isfound);
    compareind=compareind(isfound);
end
maskout=maskout(cell2mat(cellfun(@(x) ~any(isnan(x)),maskout...
    ,'UniformOutput',0)));
end


function [vol,volmic,minInt,maxInt,sumInt,meanInt,c1x,c1y,c1z,PrcOverlap,...
    OverlapObjID,CenOverlap]=calc_stat_rnascope(CC1,CC2,channel)
pixsize=0.108;
obj=CC1.PixelIdxList;
cen1=regionprops(CC1,'Centroid');
cen1=vertcat({cen1.Centroid});
obj2=CC2.PixelIdxList;
cen2=regionprops(CC2,'Centroid');
cen2=vertcat(cen2.Centroid);
if ~isempty(cen1) & ~isempty(cen2)

eudist=cellfun(@(x) sqrt(sum(bsxfun(@power,bsxfun(@minus,cen2,x),2),2)),...
    cen1,'UniformOutput',false);
vol=cellfun(@numel ,obj)';

volmic=vol*((pixsize*pixsize*.25))';
minInt=cellfun(@(x) min(channel(x)),obj)';
maxInt=cellfun(@(x) max(channel(x)),obj)';
sumInt=cellfun(@(x) sum(channel(x)),obj)';
meanInt=cellfun(@(x) mean(channel(x)),obj)';
cen1=vertcat(cen1{:});
c1x=cen1(:,2);
c1y=cen1(:,1);
c1z=cen1(:,3);

mvol=max([vol; cellfun(@numel ,obj2)']);
r=(3*mvol/(4*pi))^(1/3);
PrcOverlap=num2cell(zeros(length(obj),10));
OverlapObjID=num2cell(zeros(length(obj),10));
CenOverlap=num2cell(zeros(length(obj),10));
dim=size(channel);
for i=1:length(obj)
    nind=find(eudist{i}<5*r);
    if ~isempty(nind)
    overlap=cell2mat(arrayfun(@(x) sum(ismember(x{:},obj{i})),obj2(nind),...
        'UniformOutput',false));
    if sum(overlap)>0
    lo=length(find(overlap>0));
    if lo<10
PrcOverlap(i,1:lo)=num2cell(overlap(overlap>0)/vol(i));
    OverlapObjID(i,1:lo)=num2cell(CC2.ObjectID(nind(overlap>0)));
    cor=round(cen2(nind(overlap>0),:));
    cind=sub2ind(dim,cor(:,2),cor(:,1),cor(:,3));
    CenOverlap(i,1:lo)=num2cell(double(ismember(cind,obj{i})));
 
    end
    end
    end
end
else 
    if isempty(cen1)
    vol=zeros(1);
    volmic=zeros(1);
    minInt=zeros(1);
    maxInt=zeros(1);
   sumInt=zeros(1);
   meanInt=zeros(1);
   c1x=zeros(1);
   c1y=zeros(1);
   c1z=zeros(1);
   PrcOverlap=num2cell(zeros(1,10));
   OverlapObjID=num2cell(zeros(1,10));
   CenOverlap=num2cell(zeros(1,10));
    else
       vol=cellfun(@numel ,obj)';

volmic=vol*((pixsize*pixsize*.25))';
minInt=cellfun(@(x) min(channel(x)),obj)';
maxInt=cellfun(@(x) max(channel(x)),obj)';
sumInt=cellfun(@(x) sum(channel(x)),obj)';
meanInt=cellfun(@(x) mean(channel(x)),obj)';
cen1=vertcat(cen1{:});
c1x=cen1(:,2);
c1y=cen1(:,1);
c1z=cen1(:,3);
 PrcOverlap=num2cell(zeros(length(obj),10));
OverlapObjID=num2cell(zeros(length(obj),10));
CenOverlap=num2cell(zeros(length(obj),10));
    end
end


function writemaskout_rnascope(finalmask,lipo,ch_label,slide,outDir,...
        outDir2,gs,puncta)

for i=1:length(ch_label)
    eval(['m' ,ch_label{i} ,' = finalmask{i};']);
    eval(['gs' ,ch_label{i} ,' = gs{i};']);
    switch i
        case 1
            save([outDir '\' slide.name '_results.mat'],['m' ,ch_label{i}]);
            save([outDir '\' slide.name '_results.mat'],['gs' ,...
                ch_label{i}],'-append');
        otherwise
            save([outDir '\' slide.name '_results.mat'],['m' ,...
                ch_label{i}],'-append');
            save([outDir '\' slide.name '_results.mat'],['gs' ,...
                ch_label{i}],'-append');
    end 
end

save([outDir '\' slide.name '_results.mat'],'lipo','-append');

tind=find(slide.x405E);
cellmask=slide.x405E;
for c=1:length(ch_label)
    
      tmp = puncta{c};
    ind=find((lipo+finalmask{c})==2);
    obj2del=find(cell2mat(cellfun(@(y) sum(ismember(ind,y)),...
        tmp.PixelIdxList,'UniformOutput',false)));
    
    tmp.PixelIdxList(obj2del)=[];
    obj2keep=find(cell2mat(cellfun(@(y) sum(ismember(tind,y)),...
        tmp.PixelIdxList,'UniformOutput',false)));
    tmp.PixelIdxList(setdiff(1:length(tmp.PixelIdxList),obj2keep))=[];
    tmp.NumObjects=length(tmp.PixelIdxList);
    tmp.ObjectID=(1:tmp.NumObjects);
    CC{c}=tmp;
    clearvars tmp
end


tmp=bwconncomp(cellmask);
tmp.ObjectID=1:tmp.NumObjects;
CC=[{tmp} CC];
ch_label=[{'x405E'} ch_label];

[c1 c2]=meshgrid(1:length(ch_label),1:length(ch_label));
c1=nonzeros(reshape(triu(c1,1)+tril(c1,-1),[],1));
c2=nonzeros(reshape(triu(c2,1)+tril(c2,-1),[],1));

for i=1:numel(c1)
    Stat(i).Combination=[ch_label{c1(i)} '__'  ch_label{c2(i)}];
    [Stat(i).vol, Stat(i).volmic,Stat(i).minInt, Stat(i).maxInt, ...
        Stat(i).sumInt, Stat(i).meanInt,Stat(i).c1x, Stat(i).c1y, ...
        Stat(i).c1z, Stat(i).PrcOverlap,Stat(i).OverlapObjID,...
        Stat(i).CenOverlap,]=calc_stat_rnascope(CC{c1(i)},CC{c2(i)},...
        slide.(ch_label{c1(i)}));

end

save([outDir '\' slide.name '_results.mat'],'Stat','CC','-append');
outputfile={'Species' 'Subject' 'Section' 'Layer' 'Site' 'Mask' 'Object'...
    'Volume (voxels)' 'Volume (microns)' 'MinInt' 'MaxInt' 'SumInt'...
    'MeanInt' 'COVdimx' 'COVdimy' 'COVdimz'};
for c=1:length(ch_label)
     for ii=1:10
     outputfile=cat(2,outputfile,[ch_label{c} 'OverlapObjID' num2str(ii)],...
         [ch_label{c} 'PercentOverlap' num2str(ii)],[ch_label{c} ...
         'COVoverlap(y_n)' num2str(ii)]);
     end
end

name=slide.name;
un=strfind(name,'_');

subnum=name(1:un(1)-1);
round=name(un(1)+1:un(2)-1);
layer=name(un(2)+1:un(3)-1);
site=name(un(3)+1:end);
nc=1;
tmp=[];
for c=1:length(ch_label)-1:length(Stat)
    np=length(Stat(c).vol);
  tmp=cat(1,tmp,[repmat(ch_label(nc),np,1) num2cell(1:np)']); 
  nc=nc+1;
end
start=2; stop=length(tmp)+1;
outputfile(start:stop,1) = {'Hu'};
outputfile(start:stop,2) = {subnum};
outputfile(start:stop,3) = {round};
outputfile(start:stop,4) = {layer};
outputfile(start:stop,5) = {site};
outputfile(start:stop,6:7) = tmp;
spac=1:length(ch_label)-1:length(Stat);
tmp=[];
if length(ch_label)<=2
    for c=1:length(ch_label)
        
        switch c
            case 1
                tmp=cat(2,[num2cell(Stat(c).vol), num2cell(Stat(c).volmic),...
                    num2cell(Stat(c).minInt), num2cell(Stat(c).maxInt),...
                    num2cell(Stat(c).sumInt), num2cell(Stat(c).meanInt),...
                    num2cell(Stat(c).c1x), num2cell(Stat(c).c1y),...
                    num2cell(Stat(c).c1z)],num2cell(zeros(length...
                    (Stat(c).vol),30)),Stat(c).PrcOverlap,...
                    Stat(c).OverlapObjID,Stat(c).CenOverlap);
            case 2
                tmp=cat(1,tmp,cat(2,[num2cell(Stat(c).vol), num2cell...
                    (Stat(c).volmic), num2cell(Stat(c).minInt), ...
                    num2cell(Stat(c).maxInt), num2cell(Stat(c).sumInt), ...
                    num2cell(Stat(c).meanInt), num2cell(Stat(c).c1x), ...
                    num2cell(Stat(c).c1y), num2cell(Stat(c).c1z)],...
                    Stat(c).PrcOverlap,Stat(c).OverlapObjID,...
                    Stat(c).CenOverlap,num2cell(zeros(length(...
                    Stat(c).vol),30))));
        end
    end
else
    for c=1:length(ch_label)
        switch c
            case 1
                
                tmp=cat(2,[num2cell(Stat(c).vol), num2cell(Stat(c).volmic),...
                    num2cell(Stat(c).minInt), num2cell(Stat(c).maxInt), ...
                    num2cell(Stat(c).sumInt), num2cell(Stat(c).meanInt), ...
                    num2cell(Stat(c).c1x), num2cell(Stat(c).c1y), ...
                    num2cell(Stat(c).c1z)],...
                    num2cell(zeros(length(Stat(c).vol),30)));
                for c2=1:length(ch_label)-1
                    
                    tmp=cat(2,tmp,Stat(c2).PrcOverlap,...
                        Stat(c2).OverlapObjID,Stat(c2).CenOverlap);
                end
                
            case 2
                c1=spac(c);
                tmp2=[];
                tmp2=cat(2,tmp2,[num2cell(Stat(c1).vol), ...
                    num2cell(Stat(c1).volmic), num2cell(Stat(c1).minInt), ...
                    num2cell(Stat(c1).maxInt), num2cell(Stat(c1).sumInt), ...
                    num2cell(Stat(c1).meanInt), num2cell(Stat(c1).c1x), ...
                    num2cell(Stat(c1).c1y), num2cell(Stat(c1).c1z)],...
                    Stat(c1).PrcOverlap,Stat(c1).OverlapObjID,...
                    Stat(c1).CenOverlap,num2cell...
                    (zeros(length(Stat(c1).vol),30)));
                for c2=c1+1:spac(c+1)-1
                    
                    tmp2=cat(2,tmp2,Stat(c2).PrcOverlap,...
                        Stat(c2).OverlapObjID,Stat(c2).CenOverlap);
                end
                tmp=cat(1,tmp,tmp2);
                
            case 3
                c1=spac(c);
                tmp2=[];
                tmp2=cat(2,tmp2,[num2cell(Stat(c1).vol), ...
                    num2cell(Stat(c1).volmic), num2cell(Stat(c1).minInt), ...
                    num2cell(Stat(c1).maxInt), num2cell(Stat(c1).sumInt), ...
                    num2cell(Stat(c1).meanInt), num2cell(Stat(c1).c1x), ...
                    num2cell(Stat(c1).c1y), num2cell(Stat(c1).c1z)],...
                    Stat(c1).PrcOverlap,Stat(c1).OverlapObjID,...
                    Stat(c1).CenOverlap);
                
                for c2=1:length(ch_label)-2
                    
                    tmp2=cat(2,tmp2,Stat(c2+c1).PrcOverlap,...
                        Stat(c2+c1).OverlapObjID,Stat(c2+c1).CenOverlap);
                    if c2==1
                        tmp2=cat(2,tmp2,num2cell(...
                            zeros(length(Stat(c1).vol),30)));
                    end                    
                end
                tmp=cat(1,tmp,tmp2);               
            case 4
                c1=spac(c);
               tmp2=[];

                tmp2=cat(1,tmp2,cat(2,[num2cell(Stat(c1).vol), ...
                    num2cell(Stat(c1).volmic), num2cell(Stat(c1).minInt), ...
                    num2cell(Stat(c1).maxInt), num2cell(Stat(c1).sumInt), ...
                    num2cell(Stat(c1).meanInt), num2cell(Stat(c1).c1x), ...
                    num2cell(Stat(c1).c1y), num2cell(Stat(c1).c1z)],...
                    Stat(c1).PrcOverlap,...
                    Stat(c1).OverlapObjID,Stat(c1).CenOverlap));
                for c2=1:length(ch_label)-2
                    
                    tmp2=cat(2,tmp2,Stat(c2+c1).PrcOverlap,...
                        Stat(c2+c1).OverlapObjID,Stat(c2+c1).CenOverlap);
                end
                tmp2=cat(2,tmp2,num2cell(zeros(length(Stat(c1).vol),30)));
                tmp=cat(1,tmp,tmp2);
        end
        
    end
    
end
outputfile(2:end,8:end)=tmp;


save([outDir2 '\' name '_stats2.mat'],'outputfile','name');

toc 


end
