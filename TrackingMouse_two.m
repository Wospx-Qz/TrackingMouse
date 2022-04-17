classdef TrackingMouse_two < TrackingMouse_shulab

	properties
		mouse1
		mouse2
		cornerpoint
		sqpoint
		m1start
		m2start
		m1path
		m2path

	end

	methods
		function obj = TrackingMouse_two(varargin)
			obj = obj@TrackingMouse_shulab(varargin)
			obj.cornerpoint = getcornerpoint(obj);
			obj.sqpoint = getsqpoint(obj);
			data3 = getvideodata3(obj);
			data3 = getvideocutted(obj,data3);
			[obj.m1start,obj.m2start] = getstartpoint(obj);
			[obj.m1path,obj.m2path] = get2mpath(obj,data3);


		end

		function [cp] = getcornerpoint(obj)
			figure;
			imshow(read(obj.videoinfo,1)); hold on;
			title('click box 4 corner');
			[x,y] = ginput(4);
			x = round(x);
			y = round(y);
			plot(x,y,'ro','markerfacecolor',[1 0 0]);
			pause(1)
			close all
			cp = [x,y];
		end

		function [sq] = getsqpoint(obj)
			xmax = max(obj.cornerpoint(:,1));
			xmin = min(obj.cornerpoint(:,1));
			ymax = max(obj.cornerpoint(:,2));
			ymin = min(obj.cornerpoint(:,2));
			sq(1,:) = [xmin,ymin];
			sq(2,:) = [xmin,ymax];
			sq(3,:) = [xmax,ymax];
			sq(4,:) = [xmax,ymin];
		end

		function data2 = getvideodata3(obj)
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

		function data3 = getvideocutted(obj,data3)
			mb = zeros(size(data3(:,:,1)));
			x = obj.cornerpoint(:,1) - obj.clickpointx(1);
			y = obj.cornerpoint(:,2) - obj.clickpointy(1);
			f1 = fit(x(1:2),y(1:2),'poly1');
			f2 = fit(x(2:3),y(2:3),'poly1');
			f3 = fit(x(3:4),y(3:4),'poly1');
			f4 = fit(x([1,4]),y([1,4]),'poly1');
			[X,Y] = meshgrid(1:size(data3,2),1:size(data3,1))
			mb(f1(X)>Y(:)) = 1;
			mb(f2(X)>Y(:)) = 1;
			mb(f3(X)<Y(:)) = 1;
			mb(f4(X)<Y(:)) = 1;
			for i = 1:length(data3)
				data3(:,:,i) = data3(:,:,i) + uint8(255*mb);
			end
		end

		function [p1,p2] = getstartpoint(obj)
			%%
			I = read(obj.videoinfo,1);
			I = I(:,:,1);
			I = I(obj.clickpointy(1):obj.clickpointy(2),obj.clickpointx(1):obj.clickpointx(2),:);
			I = 255 - I - obj.background;
			I = im2bw(I,graythresh(I));
			I = bwareaopen(I,300);
			p = regionprops(I,'Centroid');
			switch length(p)
				case 0
					[p1,p2] = deal([1,1]);
				case 1
					pall = p(1).Centroid;
					ori = regionprops(I,'orientation');
					ori = ori.Orientation;
					edge = bwboundaries(I);
					[~,km] = max(cellfun(@(x) length(x),edge));
					edge = edge{km};
					edge = fliplr(edge);
					dis = pdist(edge);
					sf = squareform(dis);
					[~,k] = max(sf(:)); 
					[m,n] = ind2sub(size(sf),k);
					vec = [cos(ori),sin(ori)]; 
					p1 = pall + (edge(m,:)-pall)/2;
					p2 = pall + (edge(n,:)-pall)/2;
				otherwise
					p1 = p(1).Centroid;
					p2 = p(2).Centroid;
			end
			%imagesc(I); hold on;
			%plot(p1(1),p1(2),'r.'); hold on; 
			%plot(p2(1),p2(2),'g.'); hold on;
			%plot(pall(1),pall(2),'k.'); hold on;
			%plot(edge(m,1),edge(m,2),'c.'); hold on;
			%plot(edge(n,1),edge(n,2),'m.'); hold on;
			%axis image
			%%
		end

		function [p1,p2] = get2mpath(obj,data3)
			%%
			p1 = zeros(obj.numframes,2);
			p2 = zeros(obj.numframes,2);
			p1(1,:) = obj.m1start;
			p2(1,:) = obj.m2start;
			h = imagesc(data3(:,:,1));hold on;
			hc = plot(0,0,'k.');hold on;
			h1 = plot(0,0,'r.');hold on;
			h2 = plot(0,0,'b.');hold on;
			ps = zeros(obj.numframes,4);
			for i = 1:obj.numframes
				I = 255 - data3(:,:,i) - obj.background;
				I = im2bw(I,graythresh(I));
				I = bwareaopen(I,100);
				se = fspecial('disk',5);
				se(se>0.4*max(se(:))) = 1;
				se(se<0.4*max(se(:))) = 0;
				se = logical(se);
				I = imopen(I,se);
				p = regionprops(I,'Centroid');
				lengthp(i) =length(p);
				switch length(p) 
					case 1
						ps(i,:) = [p(1).Centroid,p(1).Centroid]; 
					case 2
						ps(i,:) = [p(1).Centroid,p(2).Centroid];
					otherwise 
						ps(i,:) = ps(i-1,:);
				end
				if i>1
					vec = ps(i,[1:2])/norm(ps(i,[1:2]));
					vec = cross([vec,0],[0 0 1]);
					vec = vec(1:2);
					[p1(i,:),p2(i,:)] = getpointrun(ps(i,:),p1(i-1,:),p2(i-1,:));
					[p1(i,:),p2(i,:)] = getpointsplit(vec,p1(i,:),p2(i,:));
				end
				h.CData = I;
				hc.XData = ps(i,1); 
				hc.YData = ps(i,2);
				h1.XData = p1(i,1);
				h1.YData = p1(i,2);
				h2.XData = p2(i,1); 
				h2.YData = p2(i,2);
				%drawnow
				pause(0.05);
			end
			%%
		end

		function checkvideo_two(obj,VMS)

            if nargin == 1
                VMS = 0;
            end
			figure;
			h = imagesc(read(obj.videoinfo,1));hold on;
			axis off image
			lwsize;
			set(gca,'position',[0 0 1 1]);
			h1 = plot(1,1,'r.'); hold on;
			h2 = plot(1,1,'b.'); hold on;
			fixv =  [obj.clickpointx(1),obj.clickpointy(1)];
            if VMS; makevideo([obj.filepath,'_checkvideo_two.avi'],obj.videoinfo.FrameRate);end
			for i = 1:obj.numframes
				h.CData = read(obj.videoinfo,i);
				h1.XData = obj.m1path(i,1) + fixv(1);
				h1.YData = obj.m1path(i,2) + fixv(2);
				h2.XData = obj.m2path(i,1) + fixv(1);
				h2.YData = obj.m2path(i,2) + fixv(2);
				drawnow;
                if VMS; makevideo;end
			end
            if VMS; makevideo(0);end

		end

	end

end

function [r1,r2] = getpointrun(pall,p1,p2)
	temp1 = pall(1:2);
	temp2 = pall(3:4);
	fr = 0.8;
	if all(temp1 == temp2)
		r1 = p1 + (temp1-p1)*fr;
		r2 = p2 + (temp2-p2)*fr;
	else
		dp1t1 = pdist([p1;temp1]);
		dp1t2 = pdist([p1;temp2]);
		dp2t1 = pdist([p2;temp1]);
		dp2t2 = pdist([p2;temp2]);
		if dp1t1 < dp1t2 & dp1t1 < dp2t1
			r1 = p1 + (temp1-p1)*fr;
			r2 = p2 + (temp2-p2)*fr;
		else
			r1 = p1 + (temp2-p1)*fr;
			r2 = p2 + (temp1-p2)*fr;
		end
	end
end

function [p1,p2] = getpointsplit(vec,p1,p2)
	r = 30;
	sf = 0.1;
	d = pdist([p1;p2]);
	while d < r
		t1 = p1;
		t2 = p2;
		p1 = p1 + (t1-t2)*sf;
		p2 = p2 + (t2-t1)*sf; 
		d = pdist([p1;p2]);
	end
	p1 = p1 + vec;
	p2 = p2 - vec;
end
