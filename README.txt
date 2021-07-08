Intelligent control with adaptive neural networks of an omnidirectional mobile robot subject to unmodeled dynamics

Gabriel S. Lima, gabriel.lima.095@ufrn.edu.br, 
Victor R. F. Moreira, victor.moreira.086@ufrn.edu.br, 
Wallace M. Bessa, wmobes@utu.fi

In this project you can find the:
(a) codes used for implementating the intelligent controller for a Robotino, an omnidirectional mobile robot made by Festo Didactic;
(b) and all the results obtained by this implementation considering the intelligent and conventional approaches.

The codes are separated in two parts:
(a) a main file with the definition of all necessary variables;
(b) a header file with all the necessary functions.

The columns of the results file are divided as follows:
[1] time
[2] x position
[3] y position
[4] angular position	
[5] x velocity	
[6] y velocity		
[7] angular velocity		
[8] x desired position
[9] y desired position
[10] angular desired position
[11] x desired velocity
[12] y desired velocity
[13] angular desired velocity
[14] x desired acceleration
[15] y desired acceleration
[16] angular desired acceleration
[17] x error position
[18] y error position
[19] angular error position
[20] x error velocity
[21] y error velocity
[22] angular error velocity
[23] s variable for x	
[24] s variable for y	
[25] s variable for the angle	
[26] torque for the wheel 1	
[27] torque for the wheel 2
[28] torque for the wheel 3
[29] angular velocity for the wheel 1	
[30] angular velocity for the wheel 2
[31] angular velocity for the wheel 3
[32] uncertainty estimation for x
[33] uncertainty estimation for y
[34] uncertainty estimation for the angle

The Armadillo library [1] is needed for the implementation.

[1] Sanderson, C. and Curtin, R., 2016. Armadillo: a template-based C++ library for linear algebra. Journal of Open Source Software, 1(2), p.26.
