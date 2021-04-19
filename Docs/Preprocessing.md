# Preprocessing
Process raw data so it can be analysed with DLC or MoSeq

## Raw data structure
```
.
+-- Group_name
|   +-- animal1_name
|       +-- date1
|           +-- session_yyyymmddhhmmss
|               +-- depth.dat
|               +-- depth_ts.txt
|               +-- metadata.json
|               +-- rgb.mp4
|               +-- rgb_ts.txt
|       +-- date2
|       +-- date3
|   +-- animal2_name
|   +-- animal3_name
```

## Separating DLC- and MoSeq-specific files
1. copy raw data to computer that can run DLC and MoSeq (this document assumes you're using "alex" and storing larges (.dat) files on DataDrive1)
2. 
