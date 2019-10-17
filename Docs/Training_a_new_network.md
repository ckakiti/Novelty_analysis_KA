# Training a new neural network
#### Starting the docker environment
```
docker start yuxie_GPU1
docker exec -it yuxie_GPU1 /bin/bash
```
#### Starting the docker environment
0. Configuration of your project:

Edit `myconfig.py`
1. Selecting data to label:
```
cd Generating_a_Training_Set
python3 Step1_SelectRandomFrames_fromVideos.py
```
2. Label the frames:

 - open ImageJ or Fiji
 - File > Import > Image Sequence
 ![alt text](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Labeling_images/Screen%20Shot%202019-10-16%20at%2012.38.43.png)
 - within pop-up window navigate to folder with images to be labeled (Generating_a_Training_Set -> data-YOUR_NETWORK)
 ![alt text](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Labeling_images/Screen%20Shot%202019-10-16%20at%2012.41.15.png)
 - click first image, then click open
 - you will see window pop up named "Sequence Options"
 
   ![alt text](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Labeling_images/Screen%20Shot%202019-10-16%20at%2012.41.37.png)
 - Make sure 2 boxes are checked: "Sort names numerically" and "Use virtual stack"
 - window will pop up with all your images in stack (scroll bar at bottom)
 - in tool bar (with File, Edit, etc), click Multi-point button (to right of angle button and to left of wand button)
     - you may see this botton as "point tool" (single star). if this happens, right click and change to be multi-point
  ![alt text](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Labeling_images/Screen%20Shot%202019-10-16%20at%2012.42.48.png)
  ![alt_text](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Labeling_images/Screen%20Shot%202019-10-16%20at%2012.43.04.png)
 - click on body features in EXACT order for every image (order specified in myconfig.py - Step 2, bodyparts variable)
 
 ![alt text](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Labeling_images/Screen%20Shot%202019-10-16%20at%2012.46.30.png)
 
   (if a point can't be determined, click in the top left corner of the image, so that X and Y positions are less than 50 pixels)
   ![alt text](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Labeling_images/Screen%20Shot%202019-10-16%20at%2012.48.13.png)
 - once you get through all frames, go to Analyze -> Measure
 ![alt text](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Labeling_images/Screen%20Shot%202019-10-16%20at%2012.48.36.png)
 - window will pop up: "Results"
     - the points you labeled will appear in rows, with each column representing a different feature of that point
     - make sure that the number of rows corresponds to N x BP where N = number of frames and BP = number of body parts you label in each frame
 ![alt text](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Labeling_images/Screen%20Shot%202019-10-16%20at%2012.49.40.png)
 - save Results window as "Results.csv" in same folder as the images you're labeling
![alt text](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Labeling_images/Screen%20Shot%202019-10-16%20at%2012.50.13.png)

3. Formatting the data I:
```
python3 Step2_ConvertingLabels2DataFrame.py
```
4. Checking the formatted data:
```
python3 Step3_CheckLabels.py
```
5. Formatting the data II:
```
python3 Step4_GenerateTrainingFileFromLabelledData.py
```
6. Training the deep neural network:

If using for the first time, download a pretrained network:
```
cd pose-tensorflow/models/pretrained
./download.sh
```
Copy the two folders generated from the last step to `/pose-tensorflow/models/`
```
cp -r YOURexperimentNameTheDate-trainset95shuffle1 ../pose-tensorflow/models/
cp -r UnaugmentedDataSet_YOURexperimentNameTheDate ../pose-tensorflow/models/
```
Start training
```
cd pose-tensorflow/models/data_set_name/train
TF_CUDNN_USE_AUTOTUNE=0 CUDA_VISIBLE_DEVICES=0 python3 ../../../train.py 
```
7. Evaluate your network:
```
CUDA_VISIBLE_DEVICES=0 python3 Step1_EvaluateModelonDataset.py #to evaluate your model [needs TensorFlow]
python3 Step2_AnalysisofResults.py  #to compute test & train errors for your trained model
```

# Analyzing videos
0. Configuration of your project:

Edit: `myconfig_analysis.py`

1. AnalyzingVideos:
```
CUDA_VISIBLE_DEVICES=0 python3 AnalyzeVideos.py
```
2. Making labeled videos
```
python3 MakingLabeledVideo.py
```

# Running on the cluster
Start an interactivate session with GPU

```
srun --pty -p gpu -t 0-06:00 --mem 8000 --gres=gpu:1 /bin/bash
srun --pty --mem=8G -n 4 -N 1 -p gpu -t 60 --gres=gpu:2 bash
```
Load mudules in the interactivate session
```
module load cuda/9.0-fasrc02 cudnn/7.0_cuda9.0-fasrc01

```

Generate a singularity image from docker local image.
```
docker run -it --name for_export dlc_user/dlc_tf1.2 /bin/true
docker export for_export > dlc_tf1.2.tar
singularity build dlc_tf1.2.simg dlc_tf1.2.tar
```
