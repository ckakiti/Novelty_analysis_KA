# Running MoSeq on the cluster (unfinished)

1. Connect to cluster
```
ssh username@login.rc.fas.harvard.edu
```
2. start interactive session
```
srun --pty --mem=500 -n 8 -N 1 -p test,shared -t 0-1:00 bash
```
3. load necessary modules (anaconda, cuda)
```
module load Anaconda/5.0.1-fasrc01
module load cuda/8.0.61-fasrc01 cudnn/6.0_cuda8.0-fasrc01
```
to list loaded modules: ``` module list ```

4. navigate to folder with yml file containing moseq environment instructions
```
cd /net/uchdafs1/uchidafs1/share_root/Lab/kakiti
```
5. create and activate moseq environment
```
conda env create -f moseq_cluster.yml
source activate moseq_cluster
```
OR
```
conda create --name moseq_cluster2 --file moseq_cluster2.txt
source activate moseq_cluster2
```
6. navigate to folder with data
```
cd Dataset_20191007
```
