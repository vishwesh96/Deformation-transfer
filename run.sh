#SAMPLING
# $2 target
# $1 source
# $3 num times
DIR='./run_dir'
THEA_DIR='./Thea/Code/Build/Output/bin'
src_name="${1%????}"
tar_name="${3%????}"
CORR='corr.pts'
num=${2//[!0-9]}
#SAMPLING
./Deform $3 $tar_name'.pts'
$THEA_DIR'/MeshSample' -n$4 $1 $src_name'.pts'
#Register
$THEA_DIR'/Register' --salient $tar_name'.picked' $src_name'.picked' $tar_name'.pts' $src_name'.pts' /dev/null
#DEFORM
./Deform $1 $2 $3 $CORR
python sparse_script.py -a A.txt -c c.txt -o o.txt
python arrange.py -i $3 -a o.txt -o $tar_name'_'$num'_def.off'
rm -rf ./*.txt
rm -rf ./*.pts
rm -rf ./*.mtx
$THEA_DIR'/Browse3D' $1 &
$THEA_DIR'/Browse3D' $2 &
$THEA_DIR'/Browse3D' $3 &
$THEA_DIR'/Browse3D' $tar_name'_'$num'_def.off' 
