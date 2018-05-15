#! /bin/sh

rsync -av --progress python /gd/MS/SCoPE/bayesian_RT/python --exclude-from="./rsync_exclude.txt"