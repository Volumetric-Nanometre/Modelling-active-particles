# Modelling-active-particles
MSci Project

File set up

File must be called particleInput.txt

----------------------------------

numberOfParticles \n
x y z alpha beta gamma dx dy dz dalpha dbeta dgamma\n

... \n

... \n

---------------------------------------------------

This is the order the input values are read.
The first number is the number of expected particles to allocate memory. If
negative or zero it will fail with error 0
Following rows must have 12 values to be read in per row in the order given above
with a newline (\n) at the end.  


Output file outputs as:

time x1 y1 z1 x2 y2 z2 .. xn yn zn alpha1 beta1 gamma1 ... alphan betan gamman
