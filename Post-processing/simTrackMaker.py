"""Module that returns course arrays for simulations"""

import numpy as np
class lineBend():

    def __init__(self, startpos, rads, size = 500,  x_dir = 1, z_dir = 1, road_width = 3.0, polywidth = None):
        """Returns a  Bend array with lines for middle and edges"""

        self.RoadStart = startpos
        self.RoadSize_Pts = size
        self.RoadWidth = road_width		
        if self.RoadWidth == 0:
            self.HalfRoadWidth = 0
        else:
            self.HalfRoadWidth = road_width/2.0		
            self.Rads = rads
            self.X_direction = x_dir

        if self.X_direction > 0:
            self.RoadArray = np.linspace(np.pi, 0.0, self.RoadSize_Pts) #right bend
        else:
            self.RoadArray = np.linspace(0.0, np.pi, self.RoadSize_Pts)  #left bend

        self.Z_direction = z_dir #[1, -1] 

        self.midline = self.LineMaker(self.Rads, self.RoadStart)
        self.OutsideLine = self.LineMaker(self.Rads + self.HalfRoadWidth, self.RoadStart)
        self.InsideLine = self.LineMaker(self.Rads - self.HalfRoadWidth, self.RoadStart)

        if polywidth is not None:
            self.OutsideLine_edge1 = self.LineMaker(self.Rads + self.HalfRoadWidth + polywidth/2, self.RoadStart)
            self.OutsideLine_edge2 = self.LineMaker(self.Rads + self.HalfRoadWidth - polywidth/2, self.RoadStart)
            self.InsideLine_edge1 = self.LineMaker(self.Rads - self.HalfRoadWidth + polywidth/2, self.RoadStart)
            self.InsideLine_edge2 = self.LineMaker(self.Rads - self.HalfRoadWidth - polywidth/2, self.RoadStart)

        translate = self.Rads * self.X_direction

        self.CurveOrigin = np.add(self.RoadStart, [translate,0])

        self.midline[:,0] = np.add(self.midline[:,0], translate)
        self.OutsideLine[:,0] = np.add(self.OutsideLine[:,0], translate)
        self.InsideLine[:,0] = np.add(self.InsideLine[:,0], translate)

        if polywidth is not None:
            self.OutsideLine_edge1[:,0] = np.add(self.OutsideLine_edge1[:,0], translate)
            self.OutsideLine_edge2[:,0] = np.add(self.OutsideLine_edge2[:,0], translate)
            self.InsideLine_edge1[:,0] = np.add(self.InsideLine_edge1[:,0], translate)
            self.InsideLine_edge2[:,0] = np.add(self.InsideLine_edge2[:,0], translate)
        
        self.RoadEnd = self.midline[-1,:]


    def LineMaker(self, Rads, startpos):
        """returns a xz array for a line"""
        #make midline        
        line = np.zeros((int(self.RoadSize_Pts),2))
        line[:,0] = Rads*np.cos(self.RoadArray) + startpos[0]
        line[:,1] = self.Z_direction*Rads*np.sin(self.RoadArray) + startpos[1]

        return line

class lineStraight():

    def __init__(self, startpos, length = 50, size = 500, z_dir = 1, road_width = 3.0, polywidth = None):

        """returns a straight, given some starting coords and length"""

        self.RoadLength = length
        self.RoadStart = startpos #2 dimensional x, z array.
        self.RoadEnd = [startpos[0],startpos[1]+(length*z_dir)] #currently only if it's north or south orientation. #2dim xz array
        self.RoadSize_Pts = size
        self.RoadWidth = road_width		
        if self.RoadWidth == 0:
            self.HalfRoadWidth = 0
        else:
            self.HalfRoadWidth = road_width/2.0		

        self.Z_direction = z_dir #[1, -1] 		
        
        #directions right, left, are from the observers point of view

        self.RightLine_Start = [self.RoadStart[0]+(self.HalfRoadWidth*z_dir), self.RoadStart[1]] 
        self.RightLine_End = [self.RoadEnd[0]+(self.HalfRoadWidth*z_dir), self.RoadEnd[1]] 
        self.LeftLine_Start = [self.RoadStart[0]-(self.HalfRoadWidth*z_dir), self.RoadStart[1]]
        self.LeftLine_End = [self.RoadEnd[0]-(self.HalfRoadWidth*z_dir), self.RoadEnd[1]]

        if polywidth is not None:
            self.RightLine_Start_1 = [self.RoadStart[0]+((self.HalfRoadWidth+polywidth/2)  *z_dir), self.RoadStart[1]] 
            self.RightLine_Start_2 = [self.RoadStart[0]+((self.HalfRoadWidth-polywidth/2)  *z_dir), self.RoadStart[1]] 
            self.RightLine_End_1 = [self.RoadEnd[0]+((self.HalfRoadWidth+ polywidth/2)*z_dir), self.RoadEnd[1]] 
            self.RightLine_End_2 = [self.RoadEnd[0]+((self.HalfRoadWidth- polywidth/2)*z_dir), self.RoadEnd[1]] 
            self.LeftLine_Start_1 = [self.RoadStart[0]-((self.HalfRoadWidth+polywidth/2)  *z_dir), self.RoadStart[1]] 
            self.LeftLine_Start_2 = [self.RoadStart[0]-((self.HalfRoadWidth-polywidth/2)  *z_dir), self.RoadStart[1]] 
            self.LeftLine_End_1 = [self.RoadEnd[0]-((self.HalfRoadWidth+ polywidth/2)*z_dir), self.RoadEnd[1]] 
            self.LeftLine_End_2 = [self.RoadEnd[0]-((self.HalfRoadWidth- polywidth/2)*z_dir), self.RoadEnd[1]] 


        self.RightLine = self.StraightLineMaker(self.RightLine_Start, self.RightLine_End)
        self.LeftLine = self.StraightLineMaker(self.LeftLine_Start, self.LeftLine_End)

        if polywidth is not None:
            self.RightLine_1 = self.StraightLineMaker(self.RightLine_Start_1, self.RightLine_End_1)
            self.RightLine_2 = self.StraightLineMaker(self.RightLine_Start_2, self.RightLine_End_2)
            self.LeftLine_1 = self.StraightLineMaker(self.LeftLine_Start_1, self.LeftLine_End_1)
            self.LeftLine_2 = self.StraightLineMaker(self.LeftLine_Start_2, self.LeftLine_End_2)
            

        
        self.midline = self.StraightLineMaker(self.RoadStart, self.RoadEnd)
    
    def StraightLineMaker(self, start, end):
        """returns midline"""
        #make midline        
            
        midline_x = np.linspace(start[0], end[0], self.RoadSize_Pts)
        midline_z = np.linspace(start[1], end[1], self.RoadSize_Pts)
        
        midline = np.column_stack((midline_x, midline_z))
            
        return midline
