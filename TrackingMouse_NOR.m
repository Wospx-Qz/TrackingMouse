classdef TrackingMouse_NOR < TrackingMouse_shulab
	properties
		targetleft
		targetright
		targetradius = 0.10;
		dist2left
		dist2right
		lefttime
		righttime
		leftlength
		rightlength
		inleftidx
		inrightidx
	end

	methods
		function obj = TrackingMouse_NOR(varargin)
			obj = obj@TrackingMouse_shulab(varargin);
			[obj.targetleft,obj.targetright] = gettarget(obj);
			[obj.dist2left,obj.dist2right] = getdist2target(obj);
			[obj.inleftidx,obj.inrightidx] = getidx(obj);
			[obj.lefttime,obj.righttime] = getintime(obj);
			[obj.leftlength,obj.rightlength] = getlength(obj);
		end

		function [l,r] = gettarget(obj)
			imshow(read(obj.videoinfo,1));hold on;
			title('Click 2 objects (from L to R)');
			[x,y] = ginput(2);
			x = round(x);
			y = round(y);
			plot(x,y,'ro','markerfacecolor',[1 0 0]);
			pause(1)
			close all;
			l = [x(1),y(1)];
			r = [x(2),y(2)];
		end

		function [dl,dr] = getdist2target(obj)
			LUcorner = [obj.clickpointx(1),obj.clickpointy(1)];
			dl = nthroot(sum((obj.centerpath - obj.targetleft + LUcorner).^10,2),10);
			dr = nthroot(sum((obj.centerpath - obj.targetright + LUcorner).^10,2),10);
		end

		function [idxl,idxr] = getidx(obj)
			idxl = obj.dist2left*obj.factor < obj.targetradius;
			idxr = obj.dist2right*obj.factor < obj.targetradius;
		end

		function [tl,tr] = getintime(obj)
			tl = length(obj.inleftidx)/obj.videoinfo.FrameRate;
			tr = length(obj.inrightidx)/obj.videoinfo.FrameRate;
		end

		function [LL,RR] = getlength(obj)
			d = obj.centerpath(2:end,:) - obj.centerpath(1:end-1,:);
			dl = sqrt(sum(d.^2,2));
            dl = [dl;dl(end)];
			LL = sum(dl(obj.inleftidx));
			RR = sum(dl(obj.inrightidx));
		end

		function checkvideo_NOR(obj,VMS)

            if nargin == 1
                VMS = 0;
            end
			figure;
			h = imagesc(read(obj.videoinfo,1));hold on;
			axis off image;
			lwsize;
			set(gca,'position',[0 0 1 1]);
			plot(obj.targetleft(1),obj.targetleft(2),'ro','markerfacecolor',[1 0 0]);
			hold on;
			plot(obj.targetright(1),obj.targetright(2),'ro','markerfacecolor',[1 0 0]);
            lc = obj.targetleft;
            rc = obj.targetright;
			r = obj.targetradius/obj.factor;
            leftedge = lc + r*[[-1 -1]; [-1 1]; [1 1 ]; [1 -1]; [-1 -1]];
            rightedge = rc + r*[[-1 -1]; [-1 1]; [1 1 ]; [1 -1]; [-1 -1]];
            lL = plot(leftedge(:,1),leftedge(:,2)); hold on;
            lR = plot(rightedge(:,1),rightedge(:,2));hold on;
            FL = cumsum((obj.inleftidx));
            FR = cumsum((obj.inrightidx));
            ttL = cumsum((obj.inleftidx))/obj.videoinfo.FrameRate;
            ttR = cumsum((obj.inrightidx))/obj.videoinfo.FrameRate;
            time_tag = text(0,0,['CenterTime:',num2str(1)],...
                'Horizontalalignment','left',...
                'Verticalalignment','top',...
                'fontsize',6,'color',[1 1 1],...
                'fontname','Times New Roman');
            d = obj.centerpath(2:end,:) - obj.centerpath(1:end-1,:);
            dL = sqrt(sum(d.^2,2));
            dL = [dL;dL(1)];
            lt = cumsum(dL)*obj.factor;
            LL = cumsum(dL.*(obj.inleftidx))*obj.factor;
            LR = cumsum(dL.*(obj.inrightidx))*obj.factor;
            [strtag strtag1 strtag2 strtag3 strtag4 strtag5] = deal('0');
            if VMS; makevideo([obj.filepath,'_checkvideo.avi'],obj.videoinfo.FrameRate);end
            for i = 1:obj.numframes
                h.CData = read(obj.videoinfo,i);
                if ismember(i,find(obj.inleftidx))
                    lL.Color = [1 1 1];
                    strtag =  sprintf('%.2f',i/obj.videoinfo.FrameRate);
                    strtag3 = sprintf('%.2f',lt(i));
                    strtag1 = sprintf('%.2f',FL(i)/obj.videoinfo.FrameRate);
                    strtag4 = sprintf('%.2f',LL(i));
                else if ismember(i,find(obj.inrightidx))
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
			
			
			

		end


	end

end

