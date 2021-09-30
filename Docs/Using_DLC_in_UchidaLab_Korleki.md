# Applying already trained network to new videos
## 1. Starting docker: 
Login the local computer, open a new terminal window check out if the docker container is running or not
```
docker ps
```
if you don't see the name `yuxie_GPU1`, that means the docker container is not running, you should start the docker container using this command.
```
docker start yuxie_GPU1
``` 

After starting the docker, run this to get inside the docker,
```
docker exec -it yuxie_GPU1 /bin/bash
```

you should see your prompt shows `root` rather than your user name now.

## 2. Changing configuration files: 
For the new network trained for MoSeq Arena, Go to the following directory (you don't need to do this in the terminal): 
`/home/alex/Programs/DeepLabCut_new/DeepLabCut/` 
and you can find many configuration files. Open `myconfig_analysis.py`, edit the `videopath` for you new videos, and save the changes. Then save a copy of this file with name `myconfig_analysis (NETWORKNAME).py` for your future use (replace NETWORKNAME with the name of your network). 

## 3. Start extraction
Go to `Analysis-tools` subfolder inside the DeepLabCut folder in the docker terminal
```
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/Analysis-tools
``` 
start extraction by running 
```
CUDA_VISIBLE_DEVICES=0 python3 AnalyzeVideos.py
```

To do this, all your video need to be directly placed under the directory you spcified in the configeration file.
If your videos are orgnized by mice names and experiment days:
* Mouse1
    * 180925
        * video1.mp4
    * 180926
        * video2.mp4
    * ...
* Mouse2
    * 180925
        * video3.mp4
    * 180926
        * video4.mp4
    * ...

you can run
```
CUDA_VISIBLE_DEVICES=0 python3 AnalyzeVideos_yuxie_180913_recursive.py
```
to traverse the subfolders.
If you orgnize your videos in different ways, you can refer to `AnalyzeVideos_yuxie_180913_recursive.py` and see sample code to do extraction recursively.

## 4. Create labeled videos
To create a video with the DLC labels superimposed, do this:
- note: you need to have run `AnalyzeVideos.py` or `AnalyzeVideos_yuxie_180913_recursive.py` before you do this step
- edit `myconfig_analysis.py` so that videofolder matches the path to your video and DLC label files (.csv/.h5/.pickle - generated from `AnalyzeVideos.py`) then run the following
```
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/Analysis-tools
python3 MakingLabeledVideo_fast.py
```
- this will give you a file with the same name as your raw video except with 'DeepLabCutlabeled' appended
- to unlock labeled video (if needed), run this line and change videopath (see `myconfig.py` and/or `myconfig_analysis.py`)
``` 
chmod -R 777 /VIDEOPATH/
```
