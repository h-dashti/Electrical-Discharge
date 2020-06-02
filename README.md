# Electrical-Discharge
Simulation of electrical discharge in an inhomogenous isulator uisng c++.

*compile command*: g++ Main.cpp

We develop a model of an electrical discharge in an inhomogenous isulator (e.g., lighting ). When an elctrical discharging occurs the elctrical potential phi, staistfy Laplace equation( Laplaican(phi) = 0 ).  The model specified by these steps:

1- Consider a large boundary circle of radius R and place charge source at origin. Choose the potential phi=0 at the origin(ocuiped site)
	and phi=1 for sites on circumference of the circle. The radius R should be larger than the radius of the growing pattern. 

2- We using Relaxation method to compuate phi_i for empty sites in circle.

3- Assign a random number r to each empty in circle. The random number r_i at site i represents a breakdown coefficent and the random inhomegenus nature of insulator.

4- The perimeter sites are the neigbour sites of the discharged pattern(occuiped sites). We form the product r_i*(phi_i)^a for each perimeter site i , where a is adjustable parameters.

5- The perimeter site with maximum value of r*phi^a , breaks down(means its potential equal to zero ).

6- We use the Relaxation method to recaclulate the value of potential at the remaining unoccuiped sites, and then we repeat steps (4) -(6).


Site(x,y) = 1	:	represnet occupied site.<br/>
Site(x,y) = 2	:	represnet a perimeter site.<br/>
Site(x,y) = 0	:	represnet an un tested site(empty).
