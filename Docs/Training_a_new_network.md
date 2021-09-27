# Training a new neural network

## 0. Setup/configuration of your project:

### Starting the docker environment
Open a new terminal window and type in the following:
```
docker start yuxie_GPU1
docker exec -it yuxie_GPU1 /bin/bash
```
If you're curious, you can check the status of the docker container by typing:
```
docker ps
```

### Setup
- Create a folder in /home/alex/Programs/DeepLabCut_new/DeepLabCut/videos named "Training_Videos_NETWORKNAME" (except replace NETWORKNAME with what you want your network to be named)
- In this folder, add the videos you want to train your network with (usually 3-4 videos that represent the breadth of the scenery/behavior that's in your full dataset)
- Go to /home/alex/Programs/DeepLabCut_new/DeepLabCut and edit `myconfig.py`

Important variables to modify and/or pay attention to:
 - **Task** (change this to whatever you want your network name to be)
 - **videopath** (change this to wherever you stored the video files you want to train your network with)
 - **videotype** (make sure this matches the file type of the videos in videopath folder)
 - **numframes2pick** (this will determine how many frames are extracted from each video in videopath folder; ideally you want to end up with anywhere from 50-200 frames total across all your videos)
 - **selectionalgorithm** (specifies how file will choose which frames to extract; kmeans takes much longer to extract, but extracted frames will better represent the variety of frames / breadth within each video)
 - **bodyparts** (the body parts you want to identify in each video; names do not matter so much, but must select this sequence of body parts in this exact order when labeling in later step)
 - **Scorers** (change this to whoever's doing the labeling)
 - **invisibleboundary** (when labeling, if you click a pixel with this value or lower (top left corner of image), the label will be set to NaN; use when you are not confident in or cannot identify a particular body part in the frame)
 - **date** (change this to the date that you do your labeling)
 - **scorer** (again change this to woever is doing the labeling)
 - additional variables are explained within the myconfig.py file

:exclamation: :exclamation: :exclamation: AFTER EDITING `myconfig.py`, MAKE SURE TO SAVE A COPY OF THIS FILE: i.e. `myconfig (NETWORKNAME).py` :exclamation: :exclamation: :exclamation:

## 1. Selecting data to label:
In the previously mentioned terminal window, type the following:
```
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut
cd Generating_a_Training_Set
python3 Step1_SelectRandomFrames_fromVideos.py
chmod -R 777 ./data-NETWORKNAME
```
This generates a folder in `Generating_a_Training_Set` named 'data-NETWORKNAME'. This folder should contain 1 folder for each video you added in step 0, and those folders should contain a series of images (the number of images is specified by the numframes2pick variable within `myconfig.py`
## 2. Label the frames:

 - open ImageJ or Fiji
 - File > Import > Image Sequence
 ![alt text](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Labeling_images/Screen%20Shot%202019-10-16%20at%2012.38.43.png)
 - within pop-up window navigate to folder with images to be labeled (`Generating_a_Training_Set` -> data-NETWORKNAME)
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
 - click on body features in EXACT order for every image (order specified by "bodyparts" in `myconfig.py`)
 
   ![alt text](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Labeling_images/Screen%20Shot%202019-10-16%20at%2012.46.30.png)
 
   (if a point can't be determined, click in the top left corner of the image, so that X and Y positions are less than 20 pixels -> determined by "invisibleboundary" variable in `myconfig.py`)
   ![alt text](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Labeling_images/Screen%20Shot%202019-10-16%20at%2012.48.13.png)
 - once you get through all frames, go to Analyze -> Measure
 ![alt text](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Labeling_images/Screen%20Shot%202019-10-16%20at%2012.48.36.png)
 - window will pop up: "Results"
     - the points you labeled will appear in rows, with each column representing a different feature of that point
     - make sure that the number of rows corresponds to N x BP where N = number of frames and BP = number of body parts you label in each frame
 ![alt text](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Labeling_images/Screen%20Shot%202019-10-16%20at%2012.49.40.png)
 - save Results window as "Results.csv" in same folder as the images you're labeling
![alt text](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Labeling_images/Screen%20Shot%202019-10-16%20at%2012.50.13.png)

## 3. Formatting the data I:
```
python3 Step2_ConvertingLabels2DataFrame.py
```
This step will add spreadsheet (.xls) files to each folder within 'data-NETWORKNAME' with a Results.csv file (which you've just created in step 2). The number of spreadsheet files should correspond to the number of body parts you labeled.

## 4. Checking the formatted data:
```
python3 Step3_CheckLabels.py
```
This step will add folders within 'data-NETWORKNAME' with labeled images. This simply overlays the labels from the Results.csv file in step 2 onto the corresponding image. It serves as a "sanity check" to make sure nothing was accidentally mislabeled.

## 5. Formatting the data II:
```
python3 Step4_GenerateTrainingFileFromLabelledData.py
```
This step will generate two folders under `Generating_a_Training_Set`:
 - NETWORKNAMEDATE-trainset95shuffle1
 - UnaugmentedDataSet_NETWORKNAME

Then make sure to change the file permissions so you can edit:
```
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/Generating_a_Training_Set/
chmod -R 777 ./NETWORKNAMEDATE-trainset95shuffle1
chmod -R 777 ./UnaugmentedDataSet_NETWORKNAME
```
 
## 6. Training the deep neural network:

If using for the first time, download a pretrained network:
```
cd pose-tensorflow/models/pretrained
./download.sh
```
If you want to train your network based on an already existing dataset on your machine, you need to edit the `pose_cfg.yaml` file for the new (not yet trained) network:
```
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/Generating_a_Training_Set
cd NETWORKNAMEDATE-trainset95shuffle1
cd train
```
 - Open `pose_cfg.yaml`
 - the variable you'll be changing is init_weights (but what you change it to depends on what network you want to start from)
 - to find the name of that network file, go to `pose-tensorflow/models/PREV_NETWORKNAME-trainset95shuffle1/train`
 - here you'll find a whole list of snapshot files, make note of the number of the last snapshot file in the list (usually 1030000)
 - in the end, the init_weights variable should have this format: `../../PREV_NETWORKNAMEDATE-trainset95shuffle1/train/snapshot-1030000` (or whatever number you noted from the list of snapshot files)
 
Copy the two folders generated from the last step to `/pose-tensorflow/models/`
```
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/Generating_a_Training_Set
cp -r NETWORKNAMEDATE-trainset95shuffle1 ../pose-tensorflow/models/
cp -r UnaugmentedDataSet_NETWORKNAMEDATE ../pose-tensorflow/models/
```
Start training
```
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/pose-tensorflow/models/NETWORKNAMEDATE-trainset95shuffle1/train
TF_CUDNN_USE_AUTOTUNE=0 CUDA_VISIBLE_DEVICES=0 python3 ../../../train.py 
```
## 7. Evaluate your network:
If you'd like to evaluate the performance of your newly trained network, do the following:
 - edit `myconfig.py` to match the network that you want to evaluate
 - run this code:
```
cd /home/alex/Programs/DeepLabCut_new/DeepLabCut/Evaluation-Tools
CUDA_VISIBLE_DEVICES=0 python3 Step1_EvaluateModelonDataset.py #to evaluate your model [needs TensorFlow]
python3 Step2_AnalysisofResults.py  #to compute test & train errors for your trained model
```
 - Step1 creates a .h5 file within `Evaluation-Tools/Results` that contains x/y labels and likelihood values for the images you extracted when training the network (see `Generating_a_Training_Set/data-NETWORKNAME`). Most of those are "Training" images (used to train the network in the first place) and a small subset are "Test" images (held out when training and used for testing labeling accuracy now).
 - Step2 creates a folder in `Evaluation-Tools` (LabeledImages_...) with labels superimposed on the aforementioned images
     - KEY: '+' = manual label; dot = DLC prediction w/ likelihood > p-cutoff; 'x' = DLC prediction w/ likelihood <= p-cutoff

# Analyzing videos
see [Running an existing network](https://github.com/ckakiti/Novelty_analysis_KA/blob/master/Docs/Using_DLC_in_UchidaLab_Korleki.md)

# Running on the cluster (unfinished)
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
