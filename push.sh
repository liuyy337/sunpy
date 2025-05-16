#!/bin/bash

git add .
git commit -m "Update at $(date +'%Y-%m-%d %H:%M:%S')"
git push origin main

