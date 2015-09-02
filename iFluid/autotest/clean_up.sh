#! /bin/bash

echo "cleaning current directory ..."
rm -f *out
rm -f *~
rm -f FAILED
rm -f WARNING
rm -f SUCCESS
echo "cleaning autotest2d/ ..."
rm -rf autotest2d/out/*
rm -rf autotest2d/archive/*
echo "cleaning autotest3d/ ..."
rm -rf autotest3d/out/*
rm -rf autotest3d/archive/*