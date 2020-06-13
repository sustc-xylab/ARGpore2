#!/bin/sh

### monitor a tab-outputed blast job by giving the approximative % done
blast=$1
query=$2
PID=$3
n=`wc -l $query | cut -f 1 -d " "`

echo "the blast out is: "$blast
echo "the fasta query is: "$query
echo ""
while true;
do
	if [ -s $blast ]; then
		curquery=`tail -1 $blast | cut -f 1`
		curline=`fgrep -n $curquery $query |  cut -f 1 -d ':'`
	else
		curline=0
	fi
	
	# report progress of metablast
	x=`echo $curline | awk '{printf("%.2f",$1)}'`  
	calc() { awk "BEGIN{print $*}"; }
	m=`calc $x/$n*100`
	echo "The blast job is about $m% finish"
	
	sleep 60
	
	# check ever one minute if the blast $PID is still runing
	if [ -s $PID ] ; then 
        for pid in `cat $PID`
        do
            kill -0 "$pid" 2>/dev/null || sed -i "/^$pid$/d" $PID  # if the monitored pid is finished then remove it from $1
        done
	else
		echo "The blast job is 100% finish"
		break
	fi
	
	# # check every 1min untill finish
	# if [ $curline -eq $n ]; then
		# break 
	# else
		# sleep 60
	# fi
done

