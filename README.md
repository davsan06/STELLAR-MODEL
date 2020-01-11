# Numerical-Model-of-Stellar-Interior
Project developed as a Final Degree Project for Physics degree at UCM (Universidad Complutense de Madrid) related with Astrophysics.

The main target of this work is to solve numerically the differential equations that govern a star given their properties (mass and composition). First of all we try to simulate the results for a model star in which we know a priori how it is going to work. Then, once we assure our code is ready to reproduce the results we introduce our particular data.

The integration of these differential equations is made in two ways: from the surface of the star to the nuclei (up-down) and the other way round, from the nuclei to the surface (down-up). By this methos we obtain the different physical parameters in every layer of the star (temperatura, pressure, luminosity, mass). Two important things which are taken into account in this project is: 

1) The way the energy is transported, by radiation or by convection (conduction is not taken into account)
2) The way the energy is generated, we distinguish between p-p cycle and cn-cycle depending on the temperature of the layer we are studying.

Something we can do when we have the result for a certain star is to study which pair of (radius, luminosity) minimize the relative error when the integration in different ways match (and for these two values we study which temperature gets the minimum error).
