classdef TrackingMouse_shulab
    properties
        arearadiu = 10; % cm, 如果目标中心区域边长20cm，那么填写20/2 = 10cm;
        videoduration = 600; % s， 分析的视频长度，10分钟 = 600s;
        filepath
        fileoutput
        filename
        centerpath
        headpath
        tailpath
        boxwidth
        boxheight
        boxcenter
        videoinfo
        background
        clickpointx
        clickpointy
        factor = 0.35/375;
        pathlength
        lengthincenter
        timeincenter
        incenterYN
        numframes
    end
    
    methods
        function obj = TrackingMouse_shulab(fp)
            if nargin == 0
                obj.filepath = getfilepath;
            else
                obj.filepath = getfilepath(fp);
            end
            if ~exist([obj.filepath,'.mat'],'file')
                obj.filename = getfilename(obj);
                obj.videoinfo = getvideoinfo(obj);
                obj.numframes = getnumframes(obj);
                [obj.clickpointx,obj.clickpointy,obj.boxcenter] = getclickpoint(obj);
                [obj.boxwidth,obj.boxheight] = getboxsize(obj);
                [data2,obj.background] = getvideodata(obj);
                [obj.centerpath,obj.headpath,obj.tailpath] = getpath(obj,data2);
                [obj.pathlength] = parametercalculate(obj);
                [obj.lengthincenter,obj.timeincenter,obj.incenterYN] = lengthinarea(obj);
                save([obj.filepath,'.mat'],'obj');
            else
                load([obj.filepath,'.mat']);
            end
        end
        
        
        
        function [filename] = getfilename(obj)
            temp = strsplit(obj.filepath,'\');
            filename = temp{end};
        end
        
        function [videoinfo] = getvideoinfo(obj)
            videoinfo = VideoReader(obj.filepath);
        end
        
        function [numframes] = getnumframes(obj)
            a = obj.videoinfo.NumFrames;
            b = obj.videoduration * obj.videoinfo.FrameRate;
            numframes = min([a,b]);
        end
        
        function [x,y,c] = getclickpoint(obj)
            imshow(read(obj.videoinfo,1)); hold on;
            title('click LeftUp and RightDown corner of the box');
            [x,y] = ginput(2);
            x = round(x);
            y = round(y);
            c(1) = mean(x) - x(1);
            c(2) = mean(y) - y(1);
            plot(x,y,'ro','markerfacecolor',[1 0 0]);
            pause(1);
            close all;
        end
        
        function [w,h] = getboxsize(obj)
            w = abs(obj.clickpointx(1) - obj.clickpointx(2));
            h = abs(obj.clickpointy(1) - obj.clickpointy(2));
        end
        
        function [data2,bk] = getvideodata(obj)
            tic;
            data = zeros(obj.videoinfo.Height,obj.videoinfo.Width,obj.numframes,'uint8');
            f = waitbar(0,'Loading Video...');
            x = obj.clickpointx;
            y = obj.clickpointy;
            for i = 1:obj.numframes
                tmp = read(obj.videoinfo,i);
                data(:,:,i) = tmp(:,:,1);
                if i == obj.numframes
                    data2 = data(y(1):y(2),x(1):x(2),:); % cut the videosize
                    bk = 255 - uint8(mean(data2(:,:,1:10:end),3)); % make background
                end
                waitbar(i/obj.numframes,f,'Loading Video...');
            end
            toc;
            delete(f);
        end
        
        function [ps,hs,ts] = getpath(obj,data2)
            h = imagesc(data2(:,:,1)); hold on;
            h1 = plot(0,0,'k.'); hold on;
            h2 = plot(0,0,'r.'); hold on;
            h3 = plot(0,0,'c.'); hold on;
            axis image off;
            lwsize;
            ps = zeros(obj.numframes,2);
            hs = zeros(obj.numframes,2);
            ts = zeros(obj.numframes,2);
            for i = 1:obj.numframes
                I = 255 - data2(:,:,i) - obj.background;
                I = im2bw(I,graythresh(I));
                I = bwareaopen(I,300);
                p = regionprops(I,'Centroid');
                h.CData = I;
                [hs(i,:),ts(i,:)] = tailtip(I);
                if length(p) == 0
                    ps(i,:) = ps(i-1,:);
                else
                    ps(i,:) = p(1).Centroid;
                end
                h1.XData = p(1).Centroid(1);
                h1.YData = p(1).Centroid(2);
                h2.XData = hs(i,2);
                h2.YData = hs(i,1);
                h3.XData = ts(i,2);
                h3.YData = ts(i,1);
                drawnow limitrate
                title([num2str(round(i/obj.videoinfo.FrameRate)),' s']);
            end
            close all
        end
        
        function [totalL] = parametercalculate(obj)
            d = obj.centerpath(2:end,:) - obj.centerpath(1:end-1,:);
            dL = sqrt(sum(d.^2,2));
            totalL = sum(dL)*obj.factor;
        end
        
        function [totalL2c,ttt,tidx] = lengthinarea(obj)
            d = obj.centerpath(2:end,:) - obj.centerpath(1:end-1,:);
            dL = sqrt(sum(d.^2,2));
            d2c = obj.centerpath - 0.5*[obj.boxwidth obj.boxheight];
            l2c = nthroot(sum(d2c.^10,2),10)*obj.factor;
            totalL2c = sum(dL((l2c(1:end-1)<obj.arearadiu/100)))*obj.factor;
            ttt = length(find(l2c(1:end-1)<obj.arearadiu/100))/obj.videoinfo.FrameRate;
            tidx = l2c(1:end-1)<(obj.arearadiu/100);
        end
        
        function drawpath(obj)
            figure;
            plot(obj.centerpath(:,1),obj.centerpath(:,2));hold on;
            axis equal off
            set(gca,'ydir','reverse');
            axis([0 obj.boxwidth 0 obj.boxheight]);
            xlabel('X');
            ylabel('Y');
            title(obj.filename,'Interpreter','none')
            lwsize;
        end
        
        function drawheatmap(obj)
            figure;
            oc = round(obj.centerpath/10);
            ig = zeros(max(oc(:)));
            for i = 1:length(oc)
                ig(oc(i,2),oc(i,1)) = ig(oc(i,1),oc(i,2)) + 1;
            end
            G = fspecial('gaussian', 7*[1 1], 1);
            Ig = imfilter(log2(ig+1),G,'same');
            Ig = interp2(Ig,3);
            %%lat
            imagesc(Ig);
            colormap(jet)
            axis off equal
            title(obj.filename,'Interpreter','none');
            lwsize;
        end
        
        function drawNheatmap(obj,N)
            figure;
            mc = max(obj.centerpath(:));
            cr = linspace(0,mc,N);
            [X,Y] = meshgrid(cr(1:end-1));
            X = X(:) + 0.5*cr(2);
            Y = Y(:) + 0.5*cr(2);
            for i = 1:length(X)
                temp = obj.centerpath - [X(i) Y(i)];
                d = nthroot(sum(temp.^10,2),10);
                ss(i) = length(find(d<cr(2)));
            end
            ig = reshape(ss,N-1,N-1);
            Ig = log2(ig+1);
            imagesc(Ig);
            colormap(parula)
            axis off equal
            title(obj.filename,'Interpreter','none');
            lwsize;
        end
        
        function checkvideo(obj,VMS)
            if nargin == 1
                VMS = 0;
            end
            
            figure;
            h = imagesc(read(obj.videoinfo,1));hold on;
            axis off image
            lwsize;
            set(gca,'position',[0 0 1 1]);
            plot(obj.clickpointx,obj.clickpointy,'ro','markerfacecolor',[1 0 0]);hold on;
            bc = obj.boxcenter+ [obj.clickpointx(1),obj.clickpointy(1)];
            r = obj.arearadiu/100/obj.factor;
            centeredge = [bc+[-r -r]; bc+[-r +r]; bc+[+r +r]; bc+[+r -r];bc+[-r -r]];
            l = plot(centeredge(:,1),centeredge(:,2));
            tr = find(obj.incenterYN);
            tt = cumsum(obj.incenterYN);
            time_tag = text(0,0,['CenterTime:',num2str(1)],...
                'Horizontalalignment','left',...
                'Verticalalignment','top',...
                'fontsize',6,'color',[1 1 1],...
                'fontname','Times New Roman');
            d = obj.centerpath(2:end,:) - obj.centerpath(1:end-1,:);
            dL = sqrt(sum(d.^2,2));
            dL = [dL;dL(1)];
            lt = cumsum(dL)*obj.factor;
            [strtag strtag1 strtag2 strtag3] = deal('0');
            if VMS; makevideo([obj.filepath,'_checkvideo.avi'],obj.videoinfo.FrameRate);end
            for i = 1:obj.numframes
                h.CData = read(obj.videoinfo,i);
                
                if ismember(i,tr)
                    l.Color = [1 1 1];
                    strtag = sprintf('%.2f',tt(i)/obj.videoinfo.FrameRate);
                    strtag2 = sprintf('%.2f',i/obj.videoinfo.FrameRate);
                    strtag3 = sprintf('%.2f',lt(tt(i)));
                    strtag4 = sprintf('%.2f',lt(i));
                    
                else
                    strtag2 = sprintf('%.2f',i/obj.videoinfo.FrameRate);
                    strtag4 = sprintf('%.2f',lt(i));
                    l.Color = [0 0 0];
                end
                time_tag.String = {
                    ['T_{all}: ',strtag2,' s'];
                    ['T_{center}: ',strtag,' s']
                    ['L_{all}: ',strtag4,' m']
                    ['L_{center}: ',strtag3,' m']
                    };
                drawnow;
                if VMS; makevideo;end
            end
            if VMS; makevideo(0);end
        end
        
        function [obj1] = novelobject(obj)
            objr = 0.05;
            imshow(read(obj.videoinfo,1));hold on;
            title('Click 2 objects (from L to R)');
            [x,y] = ginput(2);
            x = round(x);
            y = round(y);
            plot(x,y,'ro','markerfacecolor',[1 0 0]);
            pause(1)
            close all;
            d = obj.centerpath(2:end,:) - obj.centerpath(1:end-1,:);
            dL = sqrt(sum(d.^2,2));
            leftT = [x(1) y(1)] - [obj.clickpointx(1),obj.clickpointy(1)];
            rightT = [x(2) y(2)] - [obj.clickpointx(1),obj.clickpointy(1)];
            distanceToLeft = obj.centerpath - leftT;
            distanceToRight = obj.centerpath - rightT;
            distanceToLeft = obj.headpath - leftT;
            distanceToRight = obj.headpath - rightT;
            dTL = nthroot(sum(distanceToLeft.^10,2),10)*obj.factor;
            dTR = nthroot(sum(distanceToRight.^10,2),10)*obj.factor;
            LinL = sum(dL(dTL<objr))*obj.factor;
            LinR = sum(dL(dTR<objr))*obj.factor;
            TinL = length(find(dTL<objr))/obj.videoinfo.FrameRate;
            TinR = length(find(dTR<objr))/obj.videoinfo.FrameRate;
            Ltidx = dTL < objr;
            Rtidx = dTR < objr;
            obj1.leftT = leftT;
            obj1.rightT = rightT;
            obj1.dTL = dTL;
            obj1.dTR = dTR;
            obj1.LinL=LinL;
            obj1.LinR=LinR;
            obj1.TinL=TinL;
            obj1.TinR = TinR;
            obj1.Ltidx = Ltidx;
            obj1.Rtidx = Rtidx;
            obj2 = struct(obj);
            
            f = fieldnames(obj2);
            for i = 1:length(f)
                obj1.(f{i}) = obj2.(f{i});
                
            end
        end
        function checkvideo_NR(obj,obj2,VMS)
            %%
            if nargin == 2
                VMS = 0;
            end
            %%
            figure;
            objr = 0.05;
            h = imagesc(read(obj.videoinfo,1));hold on;
            axis off image
            lwsize;
            set(gca,'position',[0 0 1 1]);
            plot(obj2.leftT(1)+obj.clickpointx(1),obj2.leftT(2)+obj.clickpointy(1),'ro','markerfacecolor',[1 0 0]);hold on;
            plot(obj2.rightT(1)+obj.clickpointx(1),obj2.rightT(2)+obj.clickpointy(1),'ro','markerfacecolor',[1 0 0]);hold on;
            % 			obj2.leftT = obj2.leftT -
            %             obj2.rightT = obj2.rightT -
            lc = obj2.leftT + [obj.clickpointx(1),obj.clickpointy(1)];
            rc = obj2.rightT + [obj.clickpointx(1),obj.clickpointy(1)];
            r = 5/100/obj.factor;
            leftedge = lc + r*[[-1 -1]; [-1 1]; [1 1 ]; [1 -1]; [-1 -1]];
            rightedge = rc + r*[[-1 -1]; [-1 1]; [1 1 ]; [1 -1]; [-1 -1]];
            lL = plot(leftedge(:,1),leftedge(:,2)); hold on;
            lR = plot(rightedge(:,1),rightedge(:,2));hold on;
            FL = cumsum((obj2.Ltidx));
            FR = cumsum((obj2.Rtidx));
            
            ttL = cumsum(obj2.Ltidx)/obj.videoinfo.FrameRate;
            ttR = cumsum(obj2.Rtidx)/obj.videoinfo.FrameRate;
            time_tag = text(0,0,['CenterTime:',num2str(1)],...
                'Horizontalalignment','left',...
                'Verticalalignment','top',...
                'fontsize',6,'color',[1 1 1],...
                'fontname','Times New Roman');
            d = obj.centerpath(2:end,:) - obj.centerpath(1:end-1,:);
            dL = sqrt(sum(d.^2,2));
            dL = [dL;dL(1)];
            lt = cumsum(dL)*obj.factor;
            LL = cumsum(dL.*(obj2.dTL<objr))*obj.factor;
            LR = cumsum(dL.*(obj2.dTR<objr))*obj.factor;
            
            [strtag strtag1 strtag2 strtag3 strtag4 strtag5] = deal('0');
            if VMS; makevideo([obj.filepath,'_checkvideo.avi'],obj.videoinfo.FrameRate);end
            for i = 1:obj.numframes
                h.CData = read(obj.videoinfo,i);
                
                if ismember(i,find(obj2.Ltidx))
                    lL.Color = [1 1 1];
                    strtag =  sprintf('%.2f',i/obj.videoinfo.FrameRate);
                    strtag3 = sprintf('%.2f',lt(i));
                    strtag1 = sprintf('%.2f',FL(i)/obj.videoinfo.FrameRate);
                    strtag4 = sprintf('%.2f',LL(i));
                else if ismember(i,find(obj2.Rtidx))
                        lR.Color = [1 1 1];
                        strtag =  sprintf('%.2f',i/obj.videoinfo.FrameRate);
                        strtag3 = sprintf('%.2f',lt(i));
                        strtag2 = sprintf('%.2f',FR(i)/obj.videoinfo.FrameRate);
                        strtag5 = sprintf('%.2f',LR(i));
                    else
                        lL.Color = [0 0 0];
                        lR.Color = [0 0 0];
                        strtag = sprintf('%.2f',i/obj.videoinfo.FrameRate);
                        strtag3 = sprintf('%.2f',lt(i));
                    end
                end
                time_tag.String = {
                    ['T_{all}: ' ,strtag,' s'];
                    ['T_{left}: ',strtag1,' s']
                    ['T_{right}: ',strtag2,' s']
                    ['L_{all}: ' ,strtag3,' m']
                    ['L_{left}: ',strtag4,' m']
                    ['L_{right}: ',strtag5,' m']
                    };
                drawnow;
                if VMS; makevideo;end
            end
            if VMS; makevideo(0);end
            %%
        end
        
    end
end

function [filepath] = getfilepath(varargin)
if isempty(varargin{1})
    [f1,f2] = uigetfile('*.*');
else
    [filepath,name,ext] = fileparts(varargin{1});
    f2 = [filepath,'\'];
    f1 = [name,ext];
end
filepath = [f2,f1];
end

function [op,ed] = tailtip(I)
I = imclose(I,strel('disk',15));
STATS = regionprops(I);
% O = regionprops(I,'Orientation');
B = bwboundaries(I);
B{1}(:,1) = smooth(B{1}(:,1));
B{1}(:,2) = smooth(B{1}(:,2));
vLeft = circshift(B{1},7);
vMiddle = B{1};
vRight = circshift(B{1},-7);
v1 = vLeft - vMiddle;
v2 = vRight - vMiddle;
vm = v1 + v2;
dj = dot(v1,v2,2);
mv1 = sqrt(dot(v1,v1,2));
mv2 = sqrt(dot(v2,v2,2));
sc = real(asin(dj./(mv1.*mv2)))*180/pi;
[~,k] = max((sc));
A = circshift(B{1},-k);
op = A(1,:);
ed = A(round(length(B{1})/2),:);
end

function makevideo( videopath,fps )
%  for 循环前用  makevideo（视频路径，帧数);
%  每一个循环结束用  makevideo；
%  循环结束后用 makevideo（0）；
%  example:
%       makevideo('test.avi',20)
%       for i = 1:10
%           plot(i,i,'ro');
%           makevideo;
%       end
%       makevideo(0);

global Video_OBJ VMS

if VMS == 0
    return;
end

switch nargin
    case 0
        tempp = getframe(gcf);
        writeVideo(Video_OBJ,tempp.cdata);
    case 1
        close(Video_OBJ);
    case 2
        Video_OBJ = VideoWriter(videopath);
        Video_OBJ.FrameRate = fps;
        Video_OBJ.Quality = 50;
        open(Video_OBJ);
end
end

function lwsize(SizeType)
%   1 == 小尺寸（6.67*5）    2 == 中尺寸（9*6.75）   3 == 大尺寸（13.5*9）
if nargin == 0
    SizeType = 1;
end
switch SizeType
    case 1
        set(gcf,'Units','centimeters','position',[12 12 6.67 5]);
    case 2
        set(gcf,'Units','centimeters','position',[12 12 9 6.75]);
    case 3
        set(gcf,'Units','centimeters','position',[12 12 13.5 9]);
end
axn = length(findobj(gcf,'Type','Axes'));
if axn~=1
    return;
end
set(gca,'OuterPosition',[0 0 1 1]);
set(gca,'fontname','times new roman');
set(gca,'fontsize',10.5);
set(gca,'linewidth',0.5);
xlb = get(gca,'xlabel');
ylb = get(gca,'ylabel');
tlb = get(gca,'title');
set(xlb,'fontsize',10.5);
set(ylb,'fontsize',10.5);
set(tlb,'fontsize',10.5);
set(gcf,'Units','pixels');
set(gca,'Units','normalized');
set(gca,'layer','top');
end

