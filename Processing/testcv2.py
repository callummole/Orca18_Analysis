# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import cv2
cap = cv2.VideoCapture("testeye.mp4")
print (cap.isOpened())   # True = read video successfully. False - fail to read video.

fourcc = cv2.VideoWriter_fourcc(*'XVID')
#fourcc = cv2.VideoWriter_fourcc(*'H264')
out = cv2.VideoWriter("outputeye2.mp4", fourcc, 20.0, (640, 360))
print (out.isOpened())  # True = write out video successfully. False - fail to write out video.

cap.release()
out.release()