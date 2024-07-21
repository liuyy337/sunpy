#!/bin/bash

ffmpeg -r 10 -pattern_type glob -i 'fig_94/*.png' -c:v libx264 -pix_fmt yuv420p movie_94.mp4
# ffmpeg -r 10 -pattern_type glob -i 'fig_131/*.png' -c:v libx264 -pix_fmt yuv420p movie_131.mp4
ffmpeg -r 10 -pattern_type glob -i 'fig_171/*.png' -c:v libx264 -pix_fmt yuv420p movie_171.mp4
ffmpeg -r 10 -pattern_type glob -i 'fig_193/*.png' -c:v libx264 -pix_fmt yuv420p movie_193.mp4
ffmpeg -r 10 -pattern_type glob -i 'fig_211/*.png' -c:v libx264 -pix_fmt yuv420p movie_211.mp4
ffmpeg -r 10 -pattern_type glob -i 'fig_304/*.png' -c:v libx264 -pix_fmt yuv420p movie_304.mp4
ffmpeg -r 10 -pattern_type glob -i 'fig_335/*.png' -c:v libx264 -pix_fmt yuv420p movie_335.mp4

