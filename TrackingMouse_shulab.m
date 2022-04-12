classdef TrackingMouse_shulab
    %TRACKINGMOUSE_SHULAB 此处显示有关此类的摘要
    % a tool for openfield test
    % Tracking a mouse in behavior box
    
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
        factor
        pathlength
        lengthincenter
        timeincenter
        incenterYN
        numframes
        
    end
    
    methods
        function obj = TrackingMouse_shulab(manualpath)
            if nargin ~= 0
                [filepath,name,ext] = fileparts(manualpath);
                f2 = [filepath,'\'];
                f1 = [name,ext];
            else
                [f1,f2] = uigetfile('*.*');
            end
            obj.filepath = [f2,f1];
            obj.fileoutput = [obj.filepath,'\'];
            obj.filename = f1;
            
            if ~exist([obj.filepath,'.mat'],'file')
                %%
                V = VideoReader(obj.filepath);
                factor = 0.35/375;
                tmp = read(V,1);
                imshow(tmp(:,:,1)); hold on;
                title('click LeftUp and RightDown corner of the box');
                [x,y] = ginput(2);
                x = round(x);
                y = round(y);
                W = abs(x(1) - x(2));  % Width of your cut
                H = abs(y(1) - y(2));  % Height of your cut
                plot(x,y,'ro','markerfacecolor',[1 0 0]);
                pause(1);
                close all;
                
                numframes = min([V.NumFrames,obj.videoduration*V.FrameRate]);
                obj.numframes = numframes;
                data = zeros(V.Height,V.Width,numframes,'uint8');
                data2 = zeros(V.Height,V.Width,numframes,'uint8');
                tic;
                f = waitbar(0,'Loading Video...');
                for i = 1:numframes
                    tmp = read(V,i);
                    data(:,:,i) = tmp(:,:,1);
                    if i == numframes
                        data2 = data(y(1):y(2),x(1):x(2),:); % cut the video to fit the box
                        bk = 255 - uint8(mean(data2(:,:,1:10:end),3)); % make backgound
                        centerP(1) = mean(x) - x(1);
                        centerP(2) = mean(y) - y(1);
                    end
                    waitbar(i/numframes,f,'Loading Video...');
                end
                toc;
                delete(f);
                %%
                h = imagesc(data2(:,:,1));hold on;
                h2 = plot(0,0,'k.');hold on;
                h3 = plot(1,1,'r.');hold on;
                h4 = plot(1,1,'c.');hold off;
                axis image off;
                lwsize;
                ps = zeros(size(data,3),2);
                ass = zeros(size(data,3),2);
                bss = zeros(size(data,3),2);
                for i = 1:size(data,3)
                    I = (255 - data2(:,:,i) - bk);
                    I = im2bw(I,graythresh(I));
                    I = bwareaopen(I,300);
                    p = regionprops(I,'Centroid');
                    h.CData = I;
                    [ass(i,:),bss(i,:)] = tailtip(I);
                    h3.XData = ass(i,2);
                    h3.YData = ass(i,1);
                    h4.XData = bss(i,2);
                    h4.YData = bss(i,1);
                    if length(p) == 0
                        ps(i,:) =  ps(i-1,:);
                    else
                        h2.XData = p(1).Centroid(1);
                        h2.YData = p(1).Centroid(2);
                        ps(i,:) = p(1).Centroid;
                    end
                    drawnow limitrate
                    title([num2str(round(i/V.FrameRate)),' s']);
                end
                clear I data data2 temp
                save([obj.filepath,'.mat']);
                close all
                obj.centerpath = ps;
                obj.headpath = ass;
                obj.tailpath = bss;
                obj.boxwidth = W;
                obj.boxheight = H;
                obj.boxcenter = centerP;
                obj.videoinfo = V;
                obj.background = bk;
                obj.clickpointx = x;
                obj.clickpointy = y;
                obj.factor = factor;
                %%
            else
                load([obj.filepath,'.mat']);
                obj.centerpath = ps;
                obj.headpath = ass;
                obj.tailpath = bss;
                obj.boxwidth = W;
                obj.boxheight = H;
                obj.boxcenter = centerP;
                obj.videoinfo = V;
                obj.background = bk;
                obj.clickpointx = x;
                obj.clickpointy = y;
                obj.factor = factor;
            end
            [obj.pathlength] = parametercalculate(obj);
            [obj.lengthincenter,obj.timeincenter,obj.incenterYN] = lengthinarea(obj);
            
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
        
        
        function [] = drawpath(obj)
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
            Ig = ig;
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
        
        function checkvideo(obj,VMS)
            %%
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
                'Horizontalalignment','left','Verticalalignment','top',...
                'fontsize',6,'color',[1 1 1],'fontname','Times New Roman');
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
            %%
        end
        
        
    end
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


function [ ] = makevideo( videopath,fps )
%  for 循环前用  makevideo（视频路径，帧数);每一个循环结束用  makevideo；循环结束后用 makevideo（0）；
%
%
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

function [] = lwsize(SizeType)
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
