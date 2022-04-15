# TrackingMouse_shulab 

## Introduction

*TrackingMouse_shulab* is a tool for tracking in open-field test. 

Matlab (Image Processing Toolbox) is required.

## Usage 

1. Open a new blank matlab script.
2. input a variable as your wish.

   select file manually

```matlab
myvar = TrackingMouse_shulab;
```

​		or

```matlab
filepath = 'C:/MyData/micedata.avi';
myvar = TrackingMouse_shulab(filepath);
```

3. click the :arrow_upper_left: and  :arrow_lower_right:  corner of  behavior box bottom. It will return red points your click position.

![image-20220406204649197](.\img\1.png)

![image-20220406204710475](.\img\2.png)

![image-20220406204833662](.\img\3.png)

4. waiting for loading video.

![image-20220406205006591](.\img\4.png)

5. waiting for data analysis.

![image-20220406204909699](.\img\5.png)

6. after the analysis, you can use `myvar.drawpath` to draw  run path.

   for example:

```matlab
filepath = 'C:/MyData/micedata.avi';
myvar = TrackingMouse_shulab(filepath);
myvar.drawpath;
```

​	:star:**If exist `micedata.mat` in file path , `TrackingMouse_shulab` will get data from it directly rather than re-analysis video.**

`myvar` is a *TrackingMouse_shulab* class, you can use these **methods** or output data by these **properties**.

## Methods

### myvar.drawpath

![image-20220406210728686](.\img\6.png)

### myvar.drawheatmap

![image-20220406210800191](.\img\7.png)

### myvar.checkvideo(arg)

![image-20220406210826830](.\img\8.png)

if `arg` is 1, it will output a `_check.avi`  for checking the process.

```matlab
myvar.checkvideo;    %% not make a check video
myvar.checkvideo(0); %% not make a check video

myvar.checkvideo(1); %% make a check video

```
### myvar.drawNheatmap(N)

![imga](img/image_2022-04-14-01-16-31.png)
split heat map

## Properties

### myvar.arearadiu

The radiu of your ROI. For example, if your want a 22*22cm center ROI, your need set arearadiu to 11 in `TrackingMouse_shulab.m`

![image-20220406212025069](.\img\9.png)



### myvar.videoduration

The duration of your video data,  if your just want 10min, set to 600;



### myvar.filepath\filepathoutput

file path of your video data;



### myvar.centerpath\headpath\tailpath

The mouse's center\head\tail position during analysis at every frame. (pixel)



### myvar.boxwidth\boxheight

The width\height of behavior box. (pixel)



### myvar.boxcenter

The center position of behavior box. (pixel)



### myvar.videoinfo

Video information of your video data.



### myvar.background

Video background of your video data, to make a clearer mouse.



### myvar.clickpointx\clickpointy

Click positions of red points. (pixel)



### myvar.factor

convert pixel data to length data.

for example:

```matlab
a = 100; %% 100pixel 
b = a * myvar.factor; %% real length of 100 pixel, m
```

`factor` is 0.35/375 as default value.  $0.35\ m = 375\ pixel$



### myvar.pathlength

total path length. (m)



### myvar.lengthincenter

Length in Center(ROI) (m)



### myvar.timeincenter

Time in Center(ROI). (m)



### myvar.incenterYN

is mouse in center area ? Yes or No? (at every frame)



### myvar.numframes

frame number of your data, associated with `myvar.videoduration` and `myvar.videoinfo.FrameRate`



## Example

```matlab
myvar = TrackingMouse_shulab;
myvar.drawpath;
myvar.drawheatmap;
myvar.checkvideo;
```



## Others

if you want to get latency of  first time in center ROI:

```matlab
a = length(find(bwlabel(myvar.incenterYN) == 1));
b = myvar.videoinfo.FrameRate;
first_in_latency = a/b; 
```



