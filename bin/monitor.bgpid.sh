#！/usr/bin/bash
# script to monitor the progress of pid in backgroud 
# $1 is the file containing all the pid to monitor, when all the pid is finished will break the loop
while true;
  do
    if [ -s $1 ] ; then # break the loop if all the pid is finished
        for pid in `cat $1`
        do
            kill -0 "$pid" 2>/dev/null || sed -i "/^$pid$/d" $1  # if the monitored pid is finished then remove it from $1
        done
    else
        echo -e "All the backgroud jobs finished on `date +"%Y-%m-%d %T"`\n"
	break
    fi
done
