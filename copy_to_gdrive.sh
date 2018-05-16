#! /bin/sh

rsync -av --progress python /gd/MS/SCoPE/bayesian_RT --exclude-from="./rsync_exclude.txt"