classdef TrackingMouse_two < TrackingMouse_shulab
    
    
    properties
        
    end
    
    methods
        function obj = TrackingMouse_two(fp)
            obj = obj@TrackingMouse_shulab(fp)
            %%
            V = VideoReader(obj.filepath);
            factor = 0.35/375;
            tmp = read(V,1);
            imshow(tmp(:,:,1)); hold on;
            title('click LeftUp and RightDown corner of the box');
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
            %% rz
            h = imagesc(data2(:,:,1));hold on;
            h21 = plot(0,0,'r.');hold on;
            h22 = plot(0,0,'k.');hold on;
            h3 = plot(1,1,'r.');hold on;
            h4 = plot(1,1,'c.');hold off;
            axis image off;
            lwsize;
            ps = zeros(size(data,3),4);
            plength = [];
            xlength = [];
            last2 =1;
            for i = 1:size(data,3)
                I = (255 - data2(:,:,i) - bk);
                I = im2bw(I,graythresh(I));
                I = bwareaopen(I,300);
%                 I = imerode(I,ones(15));
                p = regionprops(I,'Centroid');
                zz = bwdist(~I);
                zz = log(zz);
                zz(zz<max(zz)*0.8) = 0;
                h.CData = zz;
                switch length(p)
                    case 0
                        ps(i,:) = [0 0 0 0];
                    case 1
                        ps(i,:) = [p.Centroid,p.Centroid];
                    case 2
                        ps(i,:) = [p.Centroid];
                    otherwise 
                        ps(i,:) = ps(i-1,:);
                end
                plength(i) = length(p);
                xlength(i) = length(find(I == 1));
                if length(p) == 2
                    last2 = i;
                end
                    
                huan = circshift(ps(i,:),2);
                if sum((huan - ps(last2,:)).^2) < sum((ps(i,:) - ps(last2,:)).^2)
                    ps(i,:) = huan;
                    1
                end
                h21.XData = ps(i,1);
                h21.YData = ps(i,2);
                h22.XData =ps(i,3);
                h22.YData =ps(i,4);
                drawnow 
                title([num2str(round(i/V.FrameRate)),' s']);
            end
        end
        
    end
