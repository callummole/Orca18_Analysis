library("tidyverse")

ttlc_from_offset <- function(b, w = 1.5, r = 60, v = 8){
  b = b / 180 * pi
   
  ttlc = sqrt(w*(2*r + sign(b)*w)/(abs(b)*r*v))
}

off_tangent <- function(w = 1.5, r = 60, v = 8){
  #ttlc = sqrt(w*(2*r + sign(b)*w)/(abs(b)*r*v))
  ttlc = sqrt(w*(2*r + w)/(v^2))
}

#for lane widths:
#https://www.standardsforhighways.co.uk/dmrb/search/66f2661f-959d-4b13-8139-92fbd491cbcf

#Traffic lane widths for horizontal curvature greater than 400 metres radii shall be in accordance with
#Figures 2.1.1N1a to 2.1.1N1h. So, curvatre > 400 will have a lane width of approximately 3.65.
#Traffic lane widths for carriageways with horizontal curve radii of greater than 90 metres but below 400
#metres are given in CD 109 [Ref 14.N].
#Traffic lane widths at junctions where the horizontal curve radii are 90 metres or less are given in CD
#123 [Ref 11.N].


#for speed and radii:
#https://www.standardsforhighways.co.uk/dmrb/search/c27c55b7-2dfc-4597-923a-4d1b4bd6c9fa

#the design speed is not actually that simple to calculate. You need a Layout constraint and alignment constraint.
#But large safe rounds will have low constraints, so we can estimate a design speed of 120 kph (3.65 m lane width).
#Single carriageways have more layout constraints (from 29 - 21 - see table 2.3). This results in a range from 85B to 100A, depending on visibility)
#hard to estimate the 'typical' single carriageway, so let's pick 100 and 85.
#This means that the design radius range from 255-2040!


#for motorways. national speed limit 70 mph (112.6 kph)
mph = c(70,70,70,70,70, 70)
ms = mph /2.237
radii = c(2880, 2040, 1440, 1020, 720, 510)
yawrate = ms/radii

motorway <- data.frame("radii" = radii, "mph" = mph, "ms" = ms, "yawrate" = yawrate)
motorway <- motorway %>% 
  mutate(ttlc_tangent = off_tangent(w = 1.825, r = radii, v = ms))
print(motorway)

#for single lane carriageways. national speed limit 60 mph (96.6 kph)
mph = c(60,60,60,60,60, 60)
ms = mph /2.237
radii = c(2040, 1440, 1020, 720, 510, 360, 255) #255 is the added minimum for 85kph design.
yawrate = ms/radii

single <- data.frame("radii" = radii, "mph" = mph, "ms" = ms, "yawrate" = yawrate)
single <- singe %>% 
  mutate(ttlc_tangent_wide = off_tangent(w = 1.825, r = radii, v = ms),
         ttlc_tangent_narrow = off_tangent(w = 1.5, r = radii, v = ms),)
print(single)
